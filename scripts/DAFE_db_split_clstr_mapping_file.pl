#!/usr/bin/env perl

# Splits master clstr file into individual files if DAFE db

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Carp;
use Readonly;
use version; our $VERSION = qv('0.0.1');
use Log::Log4perl qw(:easy);
use Log::Log4perl::CommandLine qw(:all);
use Scalar::Util qw(openhandle);

# Subroutines #
sub check_params;
sub _is_defined;

# Variables #
my ($dafe_db, $clstr_mapping_file, $out_file_exten, $help, $man);

my $options_okay = GetOptions (
    "dafe_db:s" => \$dafe_db,
    "clstr_mapping_file:s" => \$clstr_mapping_file,
    "out_file_exten:s" => \$out_file_exten,
    "help|h" => \$help,                  # flag
    "man" => \$man,                     # flag (print full man page)
);

# set up the logging environment
my $logger = get_logger();

# check for input errors
if ( $help ) { pod2usage(0) }
if ( $man ) { pod2usage(-verbose => 3) }
check_params();


########
# MAIN #
########
### read in the clstr_mapping_file and output based on genomes
my $current_genome = "";
open my $IN, "<", $clstr_mapping_file or
    $logger->logdie("Cannot open --clstr_mapping_file: $clstr_mapping_file");
    
my @vals = ();
my $OUT;  # output file handle
my %seen_genomes = ();
foreach my $line ( <$IN> ) {
    chomp $line;
    
    @vals = split(/\t/, $line);
    
    if ( $current_genome ne $vals[0] ) {
        $logger->debug("Changing genomes from $current_genome to $vals[0]");
        
        # close the previous file handle
        if ( openhandle($OUT) ) {
            close($OUT);
        }
        
        # check for duplicate genome entries
        if ( $seen_genomes{$vals[0]} ) {
            $logger->warn("Duplicate genome: $vals[0]");
        }
        
        # open the output file handle
        my $out_file_path = $dafe_db . "/" . $vals[0] . "/" . $vals[0] . $out_file_exten;
        open $OUT, ">", $out_file_path or
            $logger->error("Cannot open file: $out_file_path");
        
        if ( openhandle($OUT) ) {
            # reset current_genome and record it has been seen
            $current_genome = $vals[0];
            $seen_genomes{$current_genome} = 1;
            
            print $OUT $vals[2] . "\t" . $vals[1] . "\n";
        }
    }
    else {
        # print to the current file handle
        print $OUT $vals[2] . "\t" . $vals[1] . "\n";
    }
}

close($IN);


########
# Subs #
########
sub check_params {
	# check for required variables
	if ( ! defined $clstr_mapping_file) { 
		pod2usage(-message => "ERROR: required --clstr_mapping_file not defined\n\n",
					-exitval => 2); 
	}
	if ( ! defined $dafe_db ) {
		pod2usage(-message => "ERROR: required --dafe_db not defined\n\n",
					-exitval => 2);
	}

	# make sure required files are non-empty
	if ( defined $clstr_mapping_file and ! -e $clstr_mapping_file ) { 
		pod2usage(-message => "ERROR: --clstr_mapping_file $clstr_mapping_file is an empty file\n\n",
					-exitval => 2);
	}

	# make sure required directories exist
	if ( ! -d $dafe_db ) { 
		pod2usage(-message => "ERROR: --dafe_db is not a directory\n\n",
					-exitval => 2); 
	}
    
    # chekc the out_file_exten
    if ( ! defined $out_file_exten ) {
        $out_file_exten = "clstr.tab.txt";
        $logger->info("Setting --out_file_exten: $out_file_exten");
    }
    
    $logger->info("--dafe_db: $dafe_db");
    $logger->info("--clstr_mapping_file: $clstr_mapping_file");
    $logger->info("--out_file_exten: $out_file_exten");
	
	return 1;
}

sub _is_defined {
    my ($val, $default) = @_;
    
    if ( defined $val ) {
        return $val;
    }
    elsif ( defined $default ) {
        return $default;
    }
    else {
        return undef;
    }
}


__END__

# POD

=head1 NAME

DAFE_db_split_clstr_mapping_file.pl - Splits master clstr file into individual files
if DAFE db


=head1 VERSION

This documentation refers to version 0.0.1


=head1 SYNOPSIS

    DAFE_db_split_clstr_mapping_file.pl
        --dafe_db my_dafe_db/
        --clstr_mapping_file master_clstr_file.txt
        --out_file_exten ".clstr.tab.txt"
        
        [--help]
        [--man]
        [--debug]
        [--verbose]
        [--quiet]
        [--logfile logfile.log]

    --file | -f     Path to an input file
    --var | -v      Path to an input variable
    --help | -h     Prints USAGE statement
    --man           Prints the man page
    --debug	        Prints Log4perl DEBUG+ messages
    --verbose       Prints Log4perl INFO+ messages
    --quiet	        Suppress printing ERROR+ Log4perl messages
    --logfile       File to save Log4perl messages


=head1 ARGUMENTS
    
=head2 --dafe_db

Path to DAFE database root directory
    
=head2 --clstr_mapping_file

Path to master mapping cluster file.  This file should have 3 columns:
1) genomeID
2) clusterID
3) geneID

=head2 --out_file_exten

File extension for output file to be printed in each genomes dir in the dafe db.
DEFAUL: ".clstr.tab.txt"
 
=head2 [--help | -h]
    
An optional parameter to print a usage statement.

=head2 [--man]

An optional parameter to print he entire man page (i.e. all documentation)

=head2 [--debug]

Prints Log4perl DEBUG+ messages.  The plus here means it prints DEBUG
level and greater messages.

=head2 [--verbose]

Prints Log4perl INFO+ messages.  The plus here means it prints INFO level
and greater messages.

=head2 [--quiet]

Suppresses print ERROR+ Log4perl messages.  The plus here means it suppresses
ERROR level and greater messages that are automatically printed.

=head2 [--logfile]

File to save Log4perl messages.  Note that messages will also be printed to
STDERR.
    

=head1 DESCRIPTION

[FULL DESCRIPTION]

=head1 CONFIGURATION AND ENVIRONMENT
    
No special configurations or environment variables needed
    
    
=head1 DEPENDANCIES

version
Getopt::Long
Pod::Usage
Carp
Readonly
version
Log::Log4perl qw(:easy)
Log::Log4perl::CommandLine qw(:all)
Scalar::Util qw(openhandle)


=head1 AUTHOR

Scott Yourstone     scott.yourstone81@gmail.com
    
    
=head1 LICENCE AND COPYRIGHT

Copyright (c) 2015, Scott Yourstone
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
