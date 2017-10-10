#!/usr/bin/env perl

# summarizes the information in the DA table
# 1. the number of tests (ie cells > -2)
# 2. number of target enriched values (ie cells == 1)
# 3. number of reference enriched values( ie cells == -1)
# 4. number of genomes with at least one enriched value
# 5. number of genoems with at least one target enriched value
# 6. number of genomes with at least one reference enriched value
# 7. number of genomes with at least one enriched value from both groups
# 8. number of COGs with a target enrichment in at least one genome
# 9. number of COGs with a reference enrichment in at least one genome
# 10. number of COGs with at least enriched genome in either group
# 11. number of COGs with at least one enrichment from both groups in a given genome
# 12. most frequently target enriched COGs 
# 13. most frequently reference enriched COGs
# 14. genomes with the most target enrichments
# 15. genomes with the most reference enrichment 
# 16. genomes witht the most enrichments from either group

# REMEMBER: in the DA table the rows are features (ie COGs)
#			and the columns are genomes

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

# Subroutines #
sub check_params;
sub _is_defined;

# Variables #
my ($da_tbl_file, $top_x, $feat_meta_file, $out_file, $help, $man);

my $options_okay = GetOptions (
	"da_tbl_file:s" => \$da_tbl_file,
	"top_x:s" => \$top_x,
	"feat_meta_file:s" => \$feat_meta_file,
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
	$logger->logdie("Cannot open --out_file: $out_file for writing");

# read in the table
my $da_tbl = Table->new();
$da_tbl->load_from_file($da_tbl_file);

# the summary numbers that I want to gather include:
# 1. the number of tests (ie cells > -2)
my $tests = 0;

# 2. number of target enriched values (ie cells == 1)
my $tgt_cells = 0;

# 3. number of reference enriched values( ie cells == -1)
my $ref_cells = 0;

# 4. number of genomes with at least one enriched value
my $g_enr = 0;

# 5. number of genoems with at least one target enriched value
my $g_tgt_enr = 0;

# 6. number of genomes with at least one reference enriched value
my $g_ref_enr = 0;

# 7. number of genomes with at least one enriched value from both groups
my $g_bth_enr = 0;

# 8. number of COGs with a target enrichment in at least one genome
my $c_tgt_enr = 0;

# 9. number of COGs with a reference enrichment in at least one genome
my $c_ref_enr = 0;

# 10. number of COGs with at least enriched genome in either group
my $c_etr_enr = 0;

# 11. number of COGs with at least one enrichment from both groups in a given genome
my $c_bth_enr = 0;

# 12. most frequently target enriched COGs
my @tgt_cogs = ();

# 13. most frequently reference enriched COGs
my @ref_cogs = ();

# 14. genomes with the most target enrichments
my @tgt_genomes = ();

# 15. genomes with the most reference enrichment 
my $ref_genomes = ();

# 16. genomes witht the most enrichments from either group
my $etr_genomes = ();



my $val;
my $flag4 = 1;
my $flag5 = 1;
my $flag6 = 1;
my $flag7 = 1;
foreach my $c ( @{$da_tbl->get_col_names()} ) {
	foreach my $r ( @{$da_tbl->get_row_names()} ) {
		$val = $da_tbl->get_value_at($r, $c);

		# 1.
		if ( $val > -2 ) { $tests++; }
		
		#2.
		if ( $val == 1 ) { $tgt_cells++; }

		#3.
		if ( $val == -1 ) { $ref_cells++; }

		#4.
		if ( ($val == 1 or $val == -1 ) and $flag4 == 1 ) { $flag4 = 0; $g_enr++; }

		#5
		if ( $val == 1 and $flag5 == 1 ) { $flag5 = 0; $g_tgt_enr++; }
	
		#6.
		if ( $val == -1 and $flag6 == 1 ) { $flag6 = 0; $g_ref_enr++; }

		#7.
		if ( ($flag5 == 0 and $flag6 == 0) and $flag7 == 1 ) { $flag7 = 0; $g_bth_enr++ }
	}

	# reset the genome flags
	$flag4 = 1;
	$flag5 = 1;
	$flag6 = 1;
	$flag7 = 1;
}

### these are the COG based summary numbers
my $flag8 = 1;
my $flag9 = 1;
my $flag10 = 1;
my $flag11 = 1;

my $c_ref_tbl = Table::Numeric->new();
my $c_tgt_tbl = Table::Numeric->new();
my $g_ref_tbl = Table::Numeric->new();
my $g_tgt_tbl = Table::Numeric->new();

my @names = ("count");
foreach my $r ( @{$da_tbl->get_row_names()} ) { # rows are features
	my @ref_vals = (0);
	my @tgt_vals = (0);
	$c_ref_tbl->add_row($r, \@ref_vals, \@names);
	$c_tgt_tbl->add_row($r, \@tgt_vals, \@names);
}

foreach my $c ( @{$da_tbl->get_col_names()} ) { # columns are genomes
	my @ref_vals = (0);
	my @tgt_vals = (0);
	$g_ref_tbl->add_row($c, \@ref_vals, \@names);
	$g_tgt_tbl->add_row($c, \@tgt_vals, \@names);
}

foreach my $r ( @{$da_tbl->get_row_names()} ) {  # row are features
	foreach my $c ( @{$da_tbl->get_col_names()} ) {  # col are genomes
		$val = $da_tbl->get_value_at($r, $c);

		#8.
		if ( $val == 1 and $flag8 == 1 ) { $flag8 = 0; $c_tgt_enr++; }

		#9.
		if ( $val == -1 and $flag9 == 1 ) { $flag9 = 0; $c_ref_enr++; }

		#10.
		if ( ($val == 1 or $val == -1) and $flag10 == 1 ) { $flag10 = 0; $c_etr_enr++; }

		#11.
		if ( $flag8 == 0 and $flag9 == 0 and $flag11 == 1 ) { $flag11 = 0; $c_bth_enr++; }

		# for 12
		if ( $val == 1 ) { 
			$c_tgt_tbl->increment_at($r, "count");
			$g_tgt_tbl->increment_at($c, "count");
		}
		
		# for 13
		if ( $val == -1 ) {
			$c_ref_tbl->increment_at($r, "count");
			$g_ref_tbl->increment_at($c, "count");
		}
	}

	# reset the flags
	$flag8 = 1;
	$flag9 = 1;
	$flag10 = 1;
	$flag11 = 1;
}

# sort (ie order) the tables
$c_tgt_tbl->order("count", "T", "T");
$c_ref_tbl->order("count", "T", "T");
$g_tgt_tbl->order("count", "T", "T");
$g_ref_tbl->order("count", "T", "T");



# read in the feature metadata
my $feat_tbl = Table->new();
$feat_tbl->load_from_file($feat_meta_file);

# print results

print $OUT "1. Number of tests: $tests\n";
print $OUT "2. Number of target enriched values: $tgt_cells\n";
print $OUT "3. Number of reference enriched values: $ref_cells\n";
print $OUT "4. Number of genomes with at least one enriched value: $g_enr\n";
print $OUT "5. Number of genomes with at least one target (ie 1) enriched value: $g_tgt_enr\n";
print $OUT "6. Number of genomes with at least one reference (ie -1) enriched value: $g_ref_enr\n";
print $OUT "7. Number of genomes with at least one enrivhed value from both groups: $g_bth_enr\n";
print $OUT "8. Number of COGs with a target enrichment in at least one genome: $c_tgt_enr\n";
print $OUT "9. number of COGs with a reference enrichment in at least one genome: $c_ref_enr\n";
print $OUT "10. number of COGs with at least enriched genome in either the target or reference group: $c_etr_enr\n";
print $OUT "11. number of COGs with at least one enrichment from both groups in a given genome: $c_bth_enr\n";
print $OUT "12. most frequently target-enriched COGs:\n";
my $i = 0;
my $feat_str;
foreach my $r ( @{$c_tgt_tbl->get_row_names()} ) {
	if ( $feat_tbl->has_row($r) ) {
		$feat_str = join("\t", @{$feat_tbl->get_row($r)});
	}
	else { $feat_str = "NA"; }
	print $OUT "\t$r\t" . $c_tgt_tbl->get_value_at($r, "count") . "\t$feat_str\n";
	$i++;
	if ( $i >= $top_x ) { last; }
}
print $OUT "13. most frequently reference-enriched COGs: \n";
$i = 0;
foreach my $r ( @{$c_ref_tbl->get_row_names()} ) {
	if ( $feat_tbl->has_row($r) ) {
		$feat_str = join("\t", @{$feat_tbl->get_row($r)});
	}
	else { $feat_str = "NA"; }
	print $OUT "\t$r\t" . $c_ref_tbl->get_value_at($r, "count") . "\t$feat_str\n";
	$i++;
	if ( $i >= $top_x ) { last; }
}
print $OUT "14. genomes with the most target enrichments: \n";
$i = 0;
foreach my $r ( @{$g_tgt_tbl->get_row_names()} ) {
	print $OUT "\t$r\t" . $g_tgt_tbl->get_value_at($r, "count") . "\n";
	$i++;
	if ( $i >= $top_x ) { last; }
}
print $OUT "15. genomes with the most reference enrichment: \n";
$i = 0;
foreach my $r ( @{$g_ref_tbl->get_row_names()} ) {
	print $OUT "\t$r\t" . $g_ref_tbl->get_value_at($r, "count") . "\n";
	$i++;
	if ( $i >= $top_x ) { last; }
}


########
# Subs #
########
sub check_params {
	# check for required variables
	if ( ! defined $da_tbl_file) { 
		pod2usage(-message => "ERROR: required --da_tbl_file not defined\n\n",
					-exitval => 2); 
	}
	if ( ! defined $top_x) { 
		$top_x = 5;
		$logger->info("Setting --top_x to 5");
	}
	if ( ! defined $feat_meta_file) { 
		pod2usage(-message => "ERROR: required --feat_meta_file not defined\n\n",
					-exitval => 2); 
	}
	if ( ! defined $out_file ) {
		pod2usage(-message => "ERROR: required --out_file not defined\n\n",
					-exitval => 2);
	}

	# make sure required files are non-empty
	if ( defined $da_tbl_file and ! -e $da_tbl_file ) { 
		pod2usage(-message => "ERROR: --da_tbl_file $da_tbl_file is an empty file\n\n",
					-exitval => 2);
	}
	if ( defined $feat_meta_file and ! -e $feat_meta_file ) { 
		pod2usage(-message => "ERROR: --feat_meta_file $feat_meta_file is an empty file\n\n",
					-exitval => 2);
	}

	# make sure required directories exist
	#if ( ! -d $out_dir ) { 
	#	system("mkdir $out_dir");
	#	$logger->info("Making --out_dir: $out_dir");
	#}
	
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
