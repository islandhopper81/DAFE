#!/usr/bin/env perl

# creates a boxplot of abundance counts for a subset of genomes
# this script uses the abundance estimate

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
use DAFE::AbEstimator;
use DAFE::Utils qw(:all);

# Subroutines #
sub check_params;
sub _is_defined;

# Variables #
my ($genomes_file, $meta_file, $dafe_out, $abund_file, $out_file, $ref_meta_file, $help, $man);

my $options_okay = GetOptions (
    "genomes_file:s" => \$genomes_file,
    "meta_file:s" => \$meta_file,
	"dafe_out:s" => \$dafe_out,
	"abund_file:s" => \$abund_file,
	"out_file:s" => \$out_file,
	"ref_meta_file:s" => \$ref_meta_file,
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
#load the genomes of interest
$logger->info("Loading genomes");
my $g_aref = load_lines($genomes_file);

# load the metadata table
$logger->info("Loading metagenome metadata");
my $meta_tbl = Table->new();
$meta_tbl->load_from_file($meta_file);

# foreach gneome get the median abundace of all the marker genes
# for each sample.
$logger->info("Calculating abundance for each marker gene");
my $ab_est = DAFE::AbEstimator->new();
my $tbl = Table->new();
my $out_tbl = Table->new();

my $j = 0;
foreach my $g ( @{$g_aref} ) {
	$logger->info("Genome:$g");

	my $file = "$dafe_out/$g/$abund_file";
	check_file($file);

	$tbl->load_from_file($file);
	my $marker_avgs = $ab_est->calc_marker_abund($tbl, "median");

	for ( my $i = 0; $i < scalar @{$tbl->get_col_names()}; $i++ ) {
		my $sample = $tbl->get_col_names()->[$i];
		my $counts = $marker_avgs->[$i];
		my $fraction = $meta_tbl->get_value_at($sample, "fraction");
		my $genotype = $meta_tbl->get_value_at($sample, "genotype");
		my $age = $meta_tbl->get_value_at($sample, "age");
		my @vals = ($sample, $fraction, $genotype, $age, $g, $counts);
		my @names = ("Sample","Fraction", "Genotype", "Age", "Genome", "abund");
		$out_tbl->add_row($j, \@vals, \@names);
		$j++;
	}
}

$out_tbl->save($out_file);

# convert the genome ids to names
ids_to_names({
	in_file => $out_file,
	out_file => $out_file,
	meta_file => $ref_meta_file
});


########
# Subs #
########
sub check_params {
	# check for required variables
	if ( ! defined $genomes_file) { 
		pod2usage(-message => "ERROR: required --genomes_file not defined\n\n",
					-exitval => 2); 
	}
	if ( ! defined $meta_file) { 
		pod2usage(-message => "ERROR: required --meta_file not defined\n\n",
					-exitval => 2); 
	}
	if ( ! defined $abund_file) { 
		pod2usage(-message => "ERROR: required --abund_file not defined\n\n",
					-exitval => 2); 
	}
	if ( ! defined $dafe_out ) {
		pod2usage(-message => "ERROR: required --dafe_out not defined\n\n",
					-exitval => 2);
	}
	if ( ! defined $out_file ) {
		pod2usage(-message => "ERROR: required --out_file not defined\n\n",
					-exitval => 2);
	}
	if ( ! defined $ref_meta_file ) {
		pod2usage(-message => "ERROR: required --ref_meta_file not defined\n\n",
					-exitval => 2);
	}


	# make sure required files are non-empty
	if ( defined $genomes_file and ! -e $genomes_file ) { 
		pod2usage(-message => "ERROR: --genomes_file $genomes_file is an empty file\n\n",
					-exitval => 2);
	}
	if ( defined $meta_file and ! -e $meta_file ) { 
		pod2usage(-message => "ERROR: --meta_file $meta_file is an empty file\n\n",
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

genomes_subset_abund_boxplot.pl - colates the data for creating a boxplot of abundance for a subset of genomes


=head1 VERSION

This documentation refers to version 0.0.1


=head1 SYNOPSIS

    genome_subset_abund_boxplot.pl
        --genomes_file my_genomes.txt
        --meta_file metagenome_metadata.txt
		--abund_file abund_est.txt
		--dafe_out my_dafe_out/
		--out_file geonme_abund_boxplot_data.txt
		--ref_meta_file reference_metadata.txt
        
        [--help]
        [--man]
        [--debug]
        [--verbose]
        [--quiet]
        [--logfile logfile.log]

    --genomes_file     Path to file with list of genomes
    --meta_file      Path to sample metadata file
	--abund_file     Name of the abundance files in --dafe_out
	--dafe_out       Path to DAFE output directory
	--out_file       Path and name of output file
	--ref_meta_file  Path to reference sequence metadata
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
