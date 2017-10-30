#!/usr/bin/env perl

# updates the da_vec file in dafe_out
# this speeds up the process of making a da table
# this script is specific for killdevil
# it will get called my a parent script that parallelizes it for each genome

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
use Table::Numeric;
use File::Temp qw(tempfile);
use DAFE::Utils qw(:all);

# Subroutines #
sub check_params;

# Global Variables #
Readonly::Scalar my $ENR => 1;
Readonly::Scalar my $NS => 0;
Readonly::Scalar my $DEP => -1; 
Readonly::Scalar my $LOW => -2; 
Readonly::Scalar my $UMS => -3; 
Readonly::Scalar my $ABS => -4;

# Variables #
my ($genome, $dafe_out_dir, $da_vec_file_name, $dafe_db_dir,
	$annote_file_name, $feat_col,
	$help, $man);

my $options_okay = GetOptions (
    "genome:s" => \$genome,
	"dafe_out_dir:s" => \$dafe_out_dir,
    "da_vec_file_name:s" => \$da_vec_file_name,
	"dafe_db_dir:s" => \$dafe_db_dir,
	"annote_file_name:s" => \$annote_file_name,
	"feat_col:s" => \$feat_col,
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
# Do the following opperations on the given genome

# look for the da_vec file and open it
$logger->debug("open the DA vec file");
my $da_vec_file = "$dafe_out_dir/$genome/$da_vec_file_name";
my $da_vec_tbl = Table::Numeric->new();
$da_vec_tbl->load_from_file($da_vec_file, " ", "F");

# open the all_annote
$logger->debug("open the all annote file");
my $annote_file = "$dafe_db_dir/$genome/$annote_file_name";
my $annote_tbl = Table->new();
$annote_tbl->load_from_file($annote_file);

# create a look table of all the values present in the genomes for the
# given feature column
$logger->debug("create lookup table of features in genome");
my $is_present_aref = $annote_tbl->get_col($feat_col);
my $is_present_href = aref_to_href($is_present_aref); # from UtilSY


# go through each value in the all_anote file.
my $feat;
print aref_to_str($da_vec_tbl->get_col_names()) . "\n";
foreach my $feat ( @{$da_vec_tbl->get_row_names()} ) {
	$logger->debug("Starting feature: $feat");
	# go through each feature in the da_vec to determine if it is in the
	# annotation file.  if it is that means it is an UMS feature or greater,
	# otherwise it is ABS.
	
	if ( ! defined $is_present_href->{$feat} ) {
		# If I get to this point the annote file is saying that the current
		# feature is not present in the genomes.  So I should set the value to
		# ABS
		$da_vec_tbl->set_value_at($feat, 0, $ABS);
	}
}

# output the new table to a temp file then move it to overwrite the old one
my ($fh, $filename) = tempfile();
close($fh);
$da_vec_tbl->save($filename, " ", "F");
system("mv $filename $da_vec_file");




########
# Subs #
########
sub check_params {
	# check for required variables
	if ( ! defined $genome) { 
		pod2usage(-message => "ERROR: required --genome not defined\n\n",
					-exitval => 2); 
	}
	if ( ! defined $dafe_out_dir ) {
		pod2usage(-message => "ERROR: required --dafe_out_dir not defined\n\n",
					-exitval => 2);
	}
	if ( ! defined $da_vec_file_name) { 
		pod2usage(-message => "ERROR: required --da_vec_file_name not defined\n\n",
					-exitval => 2); 
	}
	if ( ! defined $dafe_db_dir ) {
		pod2usage(-message => "ERROR: required --dafe_db_dir not defined\n\n",
					-exitval => 2);
	}
	if ( ! defined $annote_file_name) { 
		pod2usage(-message => "ERROR: required --annote_file_name not defined\n\n",
					-exitval => 2); 
	}
	if ( ! defined $feat_col ) {
		pod2usage(-message => "ERROR: required --feat_col not defined\n\n",
					-exitval => 2);
	}

	# make sure required files are non-empty
	#if ( defined $file and ! -e $file ) { 
	#	pod2usage(-message => "ERROR: --file $file is an empty file\n\n",
	#				-exitval => 2);
	#}

	# make sure required directories exist
	if ( ! -d $dafe_out_dir ) { 
		pod2usage(-message => "ERROR: --dafe_out_dir is not a directory\n\n",
					-exitval => 2); 
	}
	if ( ! -d $dafe_db_dir ) { 
		pod2usage(-message => "ERROR: --dafe_db_dir is not a directory\n\n",
					-exitval => 2); 
	}
	
	return 1;
}


__END__

# POD

=head1 NAME

update_clstr_da_vec.pl - updates _da_vec.txt file in dafe_out


=head1 VERSION

This documentation refers to version 0.0.1


=head1 SYNOPSIS

    update_clstr_da_vec.pl
        --genome 2228664006
        --dafe_out_dir ref_out_restructured/
        --da_vec_file_name gene_counts_id60_clstr_agg_da_vec.txt
        --dafe_db_dir dafe_db/
        --annote_file_name all_annote.txt
        --feat_col clstr
        
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
    
=head2 --file | -f

Path to an input file
    
=head2 --var | -v

Path to an input variable   
 
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
