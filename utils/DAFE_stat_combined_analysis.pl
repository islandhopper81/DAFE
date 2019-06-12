#!/usr/bin/env perl

# creates a combined table from the directory structure output during the DAFE mapping.
# this is an alternative to the current algorithm.  It will eventually run the statistics
# on the entire dataset and not on a single genome at a time.

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
my ($dafe_dir, $tbl_file_name, $out_dir, $out_tbl, 
	$edger_script, $edger_params, $r_libs, $help, $man);

my $options_okay = GetOptions (
    "dafe_dir:s" => \$dafe_dir,
	"tbl_file_name:s" => \$tbl_file_name,
	"out_dir:s" => \$out_dir,
    "out_tbl:s" => \$out_tbl,
	"edgeR_script:s" => \$edger_script,
	"edgeR_params:s" => \$edger_params,
	"r_libs:s" => \$r_libs,
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

###
# 0. Create the final output table
###
my $final_tbl = Table->new();

###
# 1. create the combined table
###

# go through each directory in the dafe_dir
opendir(my $DH, $dafe_dir) or
	$logger->logdie("Cannot open --dafe_dir ($dafe_dir)\n");

my $first = 1;
my $tbl = Table->new();
foreach my $d ( readdir($DH) ) {
	# skip hidden dirs
	if ( $d =~ m/^\./ ) { next; }
	
	# skip non-directories
	if ( ! -d "$dafe_dir/$d" ) { next; }

	$logger->info("Starting genome: $d\n");
	
	# look for the table file
	my $tbl_file_path = "$dafe_dir/$d/$tbl_file_name";
	if ( ! -f $tbl_file_path ) {
		$logger->warn("Cannot file tbl file $tbl_file_path\n");
	}

	# reset the table object to remove any old data from the last table
	$tbl->reset();
	$tbl->load_from_file($tbl_file_path);

	if ( $first == 1 ) {
		# this is a special case where I need to use the tbl to populate 
		# the column names on the final table
		$logger->info("Setting the col names in the final table\n");

		$final_tbl->_set_col_count(scalar @{$tbl->get_col_names()});
		$final_tbl->_set_col_names($tbl->get_col_names());
		$first = 0;
	}

	$logger->info("Adding rows\n");
	foreach my $r_name ( @{$tbl->get_row_names()} ) {
		$final_tbl->add_row("$d\_$r_name", $tbl->get_row($r_name));
	}

	$logger->info("Genome complete\n");
}

# output the final count table
$logger->info("Saving the final table\n");
$final_tbl->save("$out_dir/$out_tbl");

closedir($DH);

###
# 2. statistical tests
###
# NOTE: output files are written to the current directory
my $r_cmd = "Rscript --no-save --no-restore $edger_script ";
$r_cmd .= "params_file=\\\"$edger_params\\\" ";
$r_cmd .= "source_dir=\\\"$r_libs\\\"";
$logger->info("Running Rscript: $r_cmd\n");
system($r_cmd);
$logger->info("Finished Rscript\n");

########
# Subs #
########
sub check_params {
	# check for required variables
	if ( ! defined $dafe_dir) { 
		pod2usage(-message => "ERROR: required --dafe_dir not defined\n\n",
					-exitval => 2); 
	}
	if ( ! defined $tbl_file_name) { 
		pod2usage(-message => "ERROR: required --tbl_file_name not defined\n\n",
					-exitval => 2); 
	}
	if ( ! defined $out_dir ) {
		pod2usage(-message => "ERROR: required --out_dir not defined\n\n",
					-exitval => 2);
	}
	if ( ! defined $out_tbl) { 
		pod2usage(-message => "ERROR: required --out_tbl not defined\n\n",
					-exitval => 2); 
	}
	if ( ! defined $edger_script ) {
		pod2usage(-message => "ERROR: required --edger_script not defined\n\n",
					-exitval => 2);
	}
	if ( ! defined $edger_params ) {
		pod2usage(-message => "ERROR: required --edger_params not defined\n\n",
					-exitval => 2);
	}
	if ( ! defined $r_libs ) {
		pod2usage(-message => "ERROR: required --r_libs not defined\n\n",
					-exitval => 2);
	}

	# make sure required files are non-empty
	if ( defined $edger_script and ! -e $edger_script ) { 
		pod2usage(-message => "ERROR: --edger_script $edger_script is an empty file\n\n",
					-exitval => 2);
	}
	if ( defined $edger_params and ! -e $edger_params ) { 
		pod2usage(-message => "ERROR: --edger_params $edger_params is an empty file\n\n",
					-exitval => 2);
	}

	# make sure required directories exist
	if ( ! -d $dafe_dir ) { 
		pod2usage(-message => "ERROR: --dafe_dir is not a directory\n\n",
					-exitval => 2); 
	}
	if ( ! -d $r_libs ) { 
		pod2usage(-message => "ERROR: --r_libs is not a directory\n\n",
					-exitval => 2); 
	}
	if ( ! -d $out_dir ) {
		$logger->info("Creating dir $out_dir\n");
		`mkdir $out_dir`;
	}
	
	return 1;
}


__END__

# POD

=head1 NAME

DAFE_stat_combined_analysis.pl - does analysis on full count table


=head1 VERSION

This documentation refers to version 0.0.1


=head1 SYNOPSIS

    DAFE_stat_combined_analysis.pl
        --dafe_dir ref_out/
        --tbl_file_name abund_est.txt
        --out_dir combined_analysis_out/
        --out_tbl full_count_tbl.txt
        --edger_script scripts/edgeR_driver_combined.R
        --edger_params edger_params.yaml
        --r_libs R_lib/
        
        [--help]
        [--man]
        [--debug]
        [--verbose]
        [--quiet]
        [--logfile logfile.log]

    --dafe_dir      Path to the DAFE mapping output directory
    --tbl_file_name Name of table file to look for in the DAFE mapping dir
    --out_dir       Path to the output directory
    --out_tbl       Name of the combined count table to output
    --edger_script  Path to edgeR script that runs the tests
    --edger_params  Path to edgeR params file (yaml format)
    --r_libs        Path to R_libs that are used in edger_script
    --help | -h     Prints USAGE statement
    --man           Prints the man page
    --debug	        Prints Log4perl DEBUG+ messages
    --verbose       Prints Log4perl INFO+ messages
    --quiet	        Suppress printing ERROR+ Log4perl messages
    --logfile       File to save Log4perl messages


=head1 ARGUMENTS
    
=head2 --dafe_dir

Path to the DAFE mapping output directory.  This directory is frequently
named something like ref_out/

=head2 --tbl_file_name

Name of the table file to look for in the DAFE mapping output directires.
    
=head2 --out_dir

Path to the output directory where output files can be found.  If the
provided directory does not exist is is automatically created.

=head2 --out_tbl

Name of the combined count table to be output.  The first column in this table
will have the genomeID and geneID seperate by a dash.

=head2 --edger_script

Path to the edgeR script that runs the statistical tests on the --out_tbl.
The full command that executes this script is constructed and ran inside
the program.

=head2 --edger_params

Path to params file for running edgeR.  Thist should be in yaml file format.

=head2 --r_libs

Path to the directory that contains the R libraries required by 
--edger_script.
 
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

The original DAFE algorithm ran statistical tests on each genomes COGs 
individually. In other words, the tests were done one genome at a time.  This 
script first creates a combined table from the directory structure output 
during the DAFE mapping. It runs the statistics on the entire dataset and not
on a single genome at a time. So this is an anternative approach to the 
original algorithm.

NOTE: there are some output files that are currently being written to the 
current working directory.


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
