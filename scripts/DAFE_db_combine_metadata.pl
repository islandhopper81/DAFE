#!/usr/bin/env perl

# combine the individual genomes metadata into a single large table

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Carp;
use Readonly;
use version; our $VERSION = qv('0.0.1');
use Log::Log4perl qw(:easy);
use Log::Log4perl::CommandLine qw(:all);
use UtilSY qw(:all);
use Table;

# Subroutines #
sub check_params;
sub _is_defined;

# Variables #
my ($dafe_db, $meta_file, $out_file, $help, $man);

my $options_okay = GetOptions (
    "dafe_db:s" => \$dafe_db,
	"meta_file:s" => \$meta_file,
    "out_file:s" => \$out_file,
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

# open the output file
open my $OUT, ">", $out_file or
	$logger->logdie("Cannot open --out_file: $out_file\n");

# go through each genome directory and read in the metadata table
my $tbl = Table->new();
my $file = "";
opendir my $DAFE, $dafe_db or
	$logger->logdie("Cannot open --dafe_db: $dafe_db\n");

foreach my $g ( readdir($DAFE) ) {
	# skip hidden dirs
	if ( $g =~ m/^\./g ) { next; }
	
	# skip non-dirs
	if ( ! -d "$dafe_db/$g" ) { next; }
	
	$file = "$dafe_db/$g/$meta_file";

	if ( ! -f $file ) {
		$logger->warn("Cannot find --meta_file: $file\n");
	}

	$tbl->load_from_file($file);

	# print the lines to --out_file
	foreach my $r ( @{$tbl->get_row_names()} ) {
		print $OUT "$g\t$r\t" . $tbl->get_value_at($r, "cog") . "\t" . $tbl->get_value_at($r, "clstr") . "\n";
	}
}

closedir($DAFE);

close($OUT);


########
# Subs #
########
sub check_params {
	# check for required variables
	if ( ! defined $out_file) { 
		pod2usage(-message => "ERROR: required --out_file not defined\n\n",
					-exitval => 2); 
	}
	if ( ! defined $meta_file) { 
		pod2usage(-message => "ERROR: required --meta_file not defined\n\n",
					-exitval => 2); 
	}
	if ( ! defined $dafe_db ) {
		pod2usage(-message => "ERROR: required --dafe_db not defined\n\n",
					-exitval => 2);
	}


	# make sure required directories exist
	if ( ! -d $dafe_db ) { 
		pod2usage(-message => "ERROR: --dafe_db is not a directory\n\n",
					-exitval => 2); 
	}
	
	return 1;
}


__END__

# POD

=head1 NAME

DAFE_db_combine_metadata.pl - combine the metadata from all the individual genomes


=head1 VERSION

This documentation refers to version 0.0.1


=head1 SYNOPSIS

    DAFE_db_combine_metadata.pl
        --dafe_db my_dafe_db/
        --meta_file all_annote.txt
        --out_file combined_metadata.txt
        
        [--help]
        [--man]
        [--debug]
        [--verbose]
        [--quiet]
        [--logfile logfile.log]

    --dafe_db       Path to the DAFE reference genome database
    --meta_file     Name of the metadata file for each genome
    --out_file      Path and name to output combined metadata file
    --help | -h     Prints USAGE statement
    --man           Prints the man page
    --debug	        Prints Log4perl DEBUG+ messages
    --verbose       Prints Log4perl INFO+ messages
    --quiet	        Suppress printing ERROR+ Log4perl messages
    --logfile       File to save Log4perl messages


=head1 ARGUMENTS
    
=head2 --dafe_db

Path to the DAFE reference genome database

=head2 --meta_file

Name of the metadata file for each genome.  This should be the same for
each genome.  If it is not the same then the genome that is different
will throw a warning saying the file for that genome cannot be found.
    
=head2 --out_file

Path and name of output file
 
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

I want to connect the COGs to the clusters that are encapsulated within each
COG.  To do this it will be much easier to have a file will all the genomes
metadata.  So this script goes through each genome directory, finds the
metadata file and prints the metadata to the output file.

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
UtilSY qw(:all)

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
