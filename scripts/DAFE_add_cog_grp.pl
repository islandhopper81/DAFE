#! /usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Carp;
use Readonly;
use Path::Class;
use Data::Dumper;
use Log::Log4perl qw(:easy);
use Log::Log4perl::CommandLine qw(:all);

# My Variables
my $help = 0;
my $man = 0;
my $list_of_dir_names;
my $cog_file;
my $db_path;

# Read in variables from the command line
GetOptions ('man'  => \$man,
            'help' => \$help,
            'dir_names|dn=s' => \$list_of_dir_names,
            'cog_data_file|cdf=s' => \$cog_file,
            'db_path|db=s' => \$db_path,
            ) || die("There was an error in the command line arguements\n");

# Use Pod usage for the Manual and Help pages
if ( $help ) { pod2usage(0) }
if ( $man )  {pod2usage(-verbose => 3) }

# Setup logging environment
my $logger = get_logger();

# Main #

#Make a hashref that contained the cog name and their group
my $cog_href = get_cog_href();

#Loop through the files and produce a new file
open my $NIN, "<", $list_of_dir_names;

foreach my $line ( <$NIN> ) {
    chomp $line;
    add_cog_group( $cog_href, $line );
}

close($NIN);

## Subroutines ##

sub get_cog_href {
    $logger->info( "Creating cog to fuction hash" );
    my %cog_hash; # Contains cog names and their functional group
    open my $FH, "<", $cog_file;
    my $skip = <$FH>;
    foreach my $line ( <$FH> ) {
        my @split_line = split /\t/,$line;
        $cog_hash{$split_line[0]} = $split_line[1];
    }
    close($FH);
    return \%cog_hash;
}

sub add_cog_group {
    my ( $href, $file_name ) = @_;
    $logger->info( "changing $file_name annote file" );
    system( "mv $db_path/$file_name/all_annote.txt $db_path/$file_name/old_all_annote.txt" );
    open my $FH, "<", "$db_path/$file_name/old_all_annote.txt";
    open my $OUT, ">", "$db_path/$file_name/all_annote.txt";
    my $first_line = <$FH>;
    my @line_a = split /\t/, $first_line;
    my $col_num = 0;
    foreach my $part ( @line_a ) {
        if ( $part =~ qr/cog/ ) {
            last;
        }
        $col_num++;
    }
    chomp $first_line;
    print $OUT $first_line . "\tcog_grp\n";
    foreach my $line ( <$FH> ) {
        chomp $line;
        my @line_array = split /\t/, $line;
        if ( ! $href->{$line_array[$col_num]} ) {
            print $OUT $line . "\tNA\n";
            next;
        }
        my $cog_group = $href->{$line_array[$col_num]};
        print $OUT $line . "\t$cog_group\n";
    }
    close($FH);
    close($OUT);
}

__END__

=head1 DAFE_add_cog_grp

This script adds a cog_grp column to all the all_annote files in the directories in the dafe_db that is passed to the script. You need to have a cog data file that links a cog to a cog_grp

=head1 Version

This documentation refers to version 0.0.1

=head1 INCLUDED MODULES

use Getopt::Long;
use Pod::Usage;
use Carp;
use Readonly;
use Path::Class;
use Data::Dumper;
use Log::Log4perl qw(:easy);
use Log::Log4perl::CommandLine qw(:all);

=head1 INHERIT

    NA
    
=head1 SYNOPSIS

    perl DAFE_add_cog_grp.pl
        -dir_names|dn genome_ids.txt
        -cog_data_file|cdf cog_info.txt
        -dafe_dir|dd dafe_db/
        
        [-help]
        [-man]
        
    -dir_names|dn       =>  Text file that contains the genomes that you want to change  the annotation files for
    
    -cog_data_file|cdf  =>  Text file that links a cogid to a cog_grp
    
    -dafe_dir|dd        =>  Directory path to the DAFE directory
    
=head1 ARGUEMENTS

    -dir_names|dn
    
    This file contains the genome names that you'd like to add the cog_grp to in the DAFE directory
    
    -cog_data_file|cdf
    
    This file must be two columns with no header. The first column needs to be a cog id and then the next column needs to be the associaed cog group
    
    -dafe_dir|dd
    
    This must be the path to the DAFE database that contains the information for every genome going to be used in downstream analysis. Every id in the dir_names file must be within this directory and have an all_annote.txt file
    
=head1 DESCRIPTION

    Given the DAFE directory and the genome ids that you would like to edit their all_annote.txt files to include cog group names. The cog group names must match to a cog_id used within the all_annote.txt file.
    
=head1 CONFIGURATION AND ENVIRONMENT

    No special configuration or environment variables needed
    
=head1 DEPENDANCIES

Getopt::Long;
Pod::Usage;
Carp;
Readonly;
Path::Class;
Data::Dumper;
Log::Log4perl qw(:easy);
Log::Log4perl::CommandLine qw(:all);

=head1 AUTHOR

Nicholas Colaianni  nick.colaianni@gmail.com

=head1 LICENCE AND COPYRIGHT

Copyright (c) 2017, Nicholas Colaianni
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met: 

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer. 
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies, 
either expressed or implied, of the FreeBSD Project.


=head1 DISCLAIMER OF WARRANTY

BECAUSE THIS SOFTWARE IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY
FOR THE SOFTWARE, TO THE EXTENT PERMITTED BY APPLICABLE LAW. EXCEPT WHEN
OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES
PROVIDE THE SOFTWARE "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER
EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE
ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE SOFTWARE IS WITH
YOU. SHOULD THE SOFTWARE PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL
NECESSARY SERVICING, REPAIR, OR CORRECTION.

IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR
REDISTRIBUTE THE SOFTWARE AS PERMITTED BY THE ABOVE LICENCE, BE
LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL,
OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE
THE SOFTWARE (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING
RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A
FAILURE OF THE SOFTWARE TO OPERATE WITH ANY OTHER SOFTWARE), EVEN IF
SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF
SUCH DAMAGES.


=cut


