#!/usr/bin/env perl

# makes a DA table when given a list of features and genomes

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
use DAFE::Utils qw(:all);

# Subroutines #
sub check_params;

# Variables #
my ($dafe_out, $dafe_db, $features_file, $feature_type, $genomes_file, 
	$tags_file, $annote_file, $out_file, $help, $man);

my $options_okay = GetOptions (
    "dafe_out:s" => \$dafe_out,
	"dafe_db:s" => \$dafe_db,
	"features_file:s" => \$features_file,
	"feature_type:s" => \$feature_type,
	"genomes_file:s" => \$genomes_file,
    "tags_file:s" => \$tags_file,
	"annote_file:s" => \$annote_file,
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
# load in the genomes
$logger->info("Loading genomes");
my $g_aref = load_lines($genomes_file);

# load in the features
my $f_aref;
eval { 
	check_file($features_file);
	$logger->info("Loading features from --features_file");
	$f_aref = load_lines($features_file);
};
if ( $@ ) {
	$logger->info("Loading all features for genomes in --genomes_file");
	$f_aref = get_genomes_features($g_aref, $dafe_db, $annote_file, $feature_type);
}

# create the output table and initialize it with -3 values
# in this table the rows are features and columns are genomes
$logger->info("Initializing output table");
my $out_tbl = Table->new();
$out_tbl->_set_col_count(scalar @{$g_aref});
$out_tbl->_set_col_names($g_aref);
foreach my $f ( @{$f_aref} ) {
	my @vals = ("-3") x scalar @{$g_aref};
	$out_tbl->add_row($f, \@vals, $g_aref);
}


# make a bunch of temporary variables that will be overwritten for
# each genome.  This saves memory and compute time.
my $tags_tbl = Table->new();
my $annote_tbl = Table->new();
my $t_file;  # full path to the tag files in dafe_out
my $a_file;  # full path to the annotation files in dafe_db
my $fdr; # fdr value from tags file for one feature
my $logFC; # logFC value from tags file for one feature
my $code;    # 0,1,-1 code based on fdr_val and logFC_val
my $feature_name; # feature name from the annote file

# go through each genome and:
# 1. look at the tags file to set the values
# 2. look at the annote file to set the -2 values
$logger->info("Building DA table");
foreach my $g ( @{$g_aref} ) {
	$logger->debug("genome: $g");

	# 1. look at the tags file to set the values (ie 1,0,-1)
	$logger->debug("Setting tag values");
	$t_file = "$dafe_out/$g/$tags_file";
	check_file($t_file);

	$tags_tbl->load_from_file($t_file, " ");
	foreach my $r ( @{$tags_tbl->get_row_names()} ) {
		# if the output table doesn't have this feature just ignore it
		if ( ! $out_tbl->has_row($r) ) {next;}

		$fdr = $tags_tbl->get_value_at($r, "FDR");
		$logFC = $tags_tbl->get_value_at($r, "logFC");
		$code = get_enr_code($logFC, $fdr);
		$out_tbl->set_value_at($r, $g, $code);
	}

	# 2. now set the -2 values (the ones that cannot be measured
	$logger->debug("Correcting -2 values");
	$a_file = "$dafe_db/$g/$annote_file";
	check_file($a_file);

	$annote_tbl->load_from_file($a_file);
	foreach my $r ( @{$annote_tbl->get_row_names()} ) {
		$feature_name = $annote_tbl->get_value_at($r, $feature_type);
		if ( $feature_name eq "NA" ) {next;}
		if ( ! $out_tbl->has_row($feature_name) ) {next;}

		if ( $out_tbl->get_value_at($feature_name, $g) == -3 ) {
			$out_tbl->set_value_at($feature_name, $g, "-2");
		}
	}
}


# print the output file
$out_tbl->save($out_file);


########
# Subs #
########
sub check_params {
	# check for required variables
	if ( ! defined $dafe_out ) { 
		pod2usage(-message => "ERROR: required --dafe_out not defined\n\n",
					-exitval => 2); 
	}
	if ( ! defined $dafe_db ) { 
		pod2usage(-message => "ERROR: required --dafe_db not defined\n\n",
					-exitval => 2); 
	}
	#if ( ! defined $features_file ) {
	#	pod2usage(-message => "ERROR: required --features_file not defined\n\n",
	#				-exitval => 2);
	#}
	if ( ! defined $feature_type ) {
		pod2usage(-message => "ERROR: required --feature_type not defined\n\n",
					-exitval => 2);
	}
	if ( ! defined $genomes_file ) {
		pod2usage(-message => "ERROR: required --genomes_file not defined\n\n",
					-exitval => 2);
	}
	if ( ! defined $tags_file ) {
		pod2usage(-message => "ERROR: required --tags_file not defined\n\n",
					-exitval => 2);
	}
	if ( ! defined $annote_file ) {
		pod2usage(-message => "ERROR: required --annote_file not defined\n\n",
					-exitval => 2);
	}
	if ( ! defined $out_file ) {
		pod2usage(-message => "ERROR: required --out_file not defined\n\n",
					-exitval => 2);
	}

	# make sure required files are non-empty
	#if ( defined $features_file and ! -e $features_file ) { 
	#	pod2usage(-message => "ERROR: --features_file $features_file is an empty file\n\n",
	#				-exitval => 2);
	#}
	# the features file is no long a required parameter
	if ( defined $genomes_file and ! -e $genomes_file ) { 
		pod2usage(-message => "ERROR: --genome_file $genomes_file is an empty file\n\n",
					-exitval => 2);
	}

	# make sure required directories exist
	if ( ! -d $dafe_db ) { 
		pod2usage(-message => "ERROR: --dafe_db is not a directory\n\n",
					-exitval => 2); 
	}
	if ( ! -d $dafe_out ) { 
		pod2usage(-message => "ERROR: --dafe_out is not a directory\n\n",
					-exitval => 2); 
	}
	
	return 1;
}


__END__

# POD

=head1 NAME

DAFE_stat_mk_da_tbl.pl - Makes a DA table given a list of genomes and features


=head1 VERSION

This documentation refers to version 0.0.1


=head1 SYNOPSIS

    DAFE_stat_mk_da_tbl.pl
        --dafe_out ref_out_restructured/
        --dafe_db dafe_db/
        --features_file my_cogs.txt
        --feature_type cog
        --genomes_file my_genomes.txt
        --tags_file gene_count_id60_tags_data.txt
        --annote_file all_annote.txt
        --out_file my_da_tbl.txt
        [--help]
        [--man]
        [--debug]
        [--verbose]
        [--quiet]
        [--logfile logfile.log]

    --dafe_out       Path to DAFE count output dir
    --dafe_db        Path to DAFE database
    --features_file  Path to file with list of features to include in DA table
    --feature_type   Name of feature (must match col header in --annote_file)
    --genomes_file   Path to file with list of genomes to include in DA table
    --tags_file      Name of tags file in --dafe_out
    --annote_file    Name of annotation file in --dafe-db
    --out_file       Path to output DA table file
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