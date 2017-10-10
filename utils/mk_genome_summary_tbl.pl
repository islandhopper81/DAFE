#!/usr/bin/env perl

# makes a summary table with numbers about each genome
# it should have the following columns:
# genomeID, abundance estimate, # Enriched COGs, # RZ Enriched COGs, # BK Enriched COGs

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
use DAFE::AbEstimator;

# Subroutines #
sub check_params;
sub _is_defined;

# Variables #
my ($da_tbl_file, $dafe_out, $sample_meta_file, $ref_meta_file, $out_file, $help, $man);

my $options_okay = GetOptions (
    "da_tbl_file:s" => \$da_tbl_file,
	"dafe_out:s" => \$dafe_out,
	"sample_meta_file:s" => \$sample_meta_file,
	"ref_meta_file:s" => \$ref_meta_file,
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
my $out_tbl = Table->new();

# read in the DA table and calculate the enrichments for each genome
$logger->info("Enrichment summary");
my $da_tbl = Table->new();
$da_tbl->load_from_file($da_tbl_file);

foreach my $g ( @{$da_tbl->get_col_names()} ) {
	$logger->debug("Genome: $g");
	my $enr = 0;
	my $rz = 0;
	my $bk = 0;

	foreach my $val ( @{$da_tbl->get_col($g)} ) {
		if ( $val == 1 ) { $rz++; $enr++; }
		elsif ($val == -1 ) { $bk++; $enr++; }
	}

	my @names = ("Enriched", "RZ Enriched", "BK Enriched");
	my @vals = ($enr, $rz, $bk);
	$out_tbl->add_row($g, \@vals, \@names);
}

### load the sample metadata so I can split the table by fraction
my $meta_tbl = Table->new();
$meta_tbl->load_from_file($sample_meta_file);

### load the genome metadata so I can get the PA/NPA/Soil 
# designation for each genome
my $ref_meta_tbl = Table->new();
$ref_meta_tbl->load_from_file($ref_meta_file);

### calculate the abundance of each genome in each fraction
$logger->info("Calculate genome abundance");

# first create a column in the output table for this
my @rz_median_abund = (0) x $out_tbl->get_row_count();
my @bk_median_abund = (0) x $out_tbl->get_row_count();
my @rz_mean_abund = (0) x $out_tbl->get_row_count();
my @bk_mean_abund = (0) x $out_tbl->get_row_count();
$out_tbl->add_col("RZ Median Abundance Estimate", \@rz_median_abund);
$out_tbl->add_col("BK Median Abundance Estimate", \@bk_median_abund);
$out_tbl->add_col("RZ Mean Abundance Estimate", \@rz_mean_abund);
$out_tbl->add_col("BK Mean Abundance Estimate", \@bk_mean_abund);

# add another column for the PA/NPA/Soil designation
my @pa_designation = ("NA") x $out_tbl->get_row_count();
$out_tbl->add_col("PA/NPA/Soil", \@pa_designation);

opendir(my $DAFE, $dafe_out) or
	$logger->logdie("Cannot open --dafe_out: $dafe_out");

my $file;
my $abund_tbl = Table->new();
my $rz_abund_tbl = Table->new();
my $bk_abund_tbl = Table->new();
my $estimator = DAFE::AbEstimator->new();

foreach my $g ( readdir($DAFE) ) {
	if ( $g =~ m/^\./ ) { next; } #skip hidden files
	if ( ! -d "$dafe_out/$g/" ) { next; } #skip non dirs
	$logger->debug("Calculating abundance estimate for genome: $g");

	$file = "$dafe_out/$g/abund_est.txt";
	$abund_tbl->load_from_file($file);
	$rz_abund_tbl->reset();
	$bk_abund_tbl->reset();

	# split the table by fraction
	$logger->debug("Splitting table by fraction");
	foreach my $c ( @{$abund_tbl->get_col_names()} ) {
		my $col_fraction = $meta_tbl->get_value_at($c, "fraction");
		if ( $col_fraction =~ m/RZ/i ) {
			$rz_abund_tbl->add_col($c, 
				$abund_tbl->get_col($c), 
				$abund_tbl->get_row_names()
			);
		}
		elsif ( $col_fraction =~ m/BK/i ) {
			$bk_abund_tbl->add_col($c, 
				$abund_tbl->get_col($c), 
				$abund_tbl->get_row_names()
			);
		}
	}

	$logger->debug("rz dim; " . $rz_abund_tbl->get_row_count() . ", " . $rz_abund_tbl->get_col_count());
	

	# get the abundance estimate
	my $rz_mean_est = $estimator->calc_abund_est($rz_abund_tbl, "mean");
	my $rz_median_est = $estimator->calc_abund_est($rz_abund_tbl, "median");
	my $bk_mean_est = $estimator->calc_abund_est($bk_abund_tbl, "mean");
	my $bk_median_est = $estimator->calc_abund_est($bk_abund_tbl, "median");

	# add the abundance estimate to the output table
	if ( $out_tbl->has_row($g) ) {
		$logger->debug("Adding abund info for genome ($g): ($rz_median_est / $bk_median_est)");
		$out_tbl->set_value_at($g, "RZ Median Abundance Estimate", $rz_median_est);
		$out_tbl->set_value_at($g, "BK Median Abundance Estimate", $bk_median_est);
		$out_tbl->set_value_at($g, "RZ Mean Abundance Estimate", $rz_mean_est);
		$out_tbl->set_value_at($g, "BK Mean Abundance Estimate", $bk_mean_est);
		$out_tbl->set_value_at($g, "PA/NPA/Soil", $ref_meta_tbl->get_value_at($g, "Label"));
	}
	else {
		$logger->warn("Genome missing from DA table: $g");
	}
}

closedir($DAFE);

$out_tbl->save($out_file);



########
# Subs #
########
sub check_params {
	# check for required variables
	if ( ! defined $da_tbl_file) { 
		pod2usage(-message => "ERROR: required --da_tbl_file not defined\n\n",
					-exitval => 2); 
	}
	if ( ! defined $dafe_out ) {
		pod2usage(-message => "ERROR: required --dafe_out not defined\n\n",
					-exitval => 2);
	}
	if ( ! defined $sample_meta_file) { 
		pod2usage(-message => "ERROR: required --sample_meta_file not defined\n\n",
					-exitval => 2); 
	}
	if ( ! defined $ref_meta_file) { 
		pod2usage(-message => "ERROR: required --ref_meta_file not defined\n\n",
					-exitval => 2); 
	}
	if ( ! defined $out_file) { 
		pod2usage(-message => "ERROR: required --out_file not defined\n\n",
					-exitval => 2); 
	}

	# make sure required files are non-empty
	if ( defined $da_tbl_file and ! -e $da_tbl_file ) { 
		pod2usage(-message => "ERROR: --da_tbl_file $da_tbl_file is an empty file\n\n",
					-exitval => 2);
	}
	if ( defined $sample_meta_file and ! -e $sample_meta_file ) { 
		pod2usage(-message => "ERROR: --sample_meta_file $sample_meta_file is an empty file\n\n",
					-exitval => 2);
	}
	if ( defined $ref_meta_file and ! -e $ref_meta_file ) { 
		pod2usage(-message => "ERROR: --ref_meta_file $ref_meta_file is an empty file\n\n",
					-exitval => 2);
	}

	# make sure required directories exist
	if ( ! -d $dafe_out ) { 
		pod2usage(-message => "ERROR: --dafe_out is not a directory\n\n",
					-exitval => 2); 
	}
	
	return 1;
}


__END__

# POD

=head1 NAME

[NAME].pl - [DESCRIPTION]


=head1 VERSION

This documentation refers to version 0.0.1


=head1 SYNOPSIS

    [NAME].pl
        -f my_file.txt
        -v 10
        
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
