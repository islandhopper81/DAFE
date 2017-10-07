#!/usr/bin/env perl

# gets the metadata associated with the genes in a cluster

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
my ($clstr, $isai_clstr_file, $dafe_db, $out_file, $help, $man);

my $options_okay = GetOptions (
    "clstr:s" => \$clstr,
    "isai_clstr_file:s" => \$isai_clstr_file,
	"dafe_db:s" => \$dafe_db,
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
# read in Isai's cluster file and get all the genomes and gene 
# assigned to that cluster
open my $IN, "<", $isai_clstr_file or 
	$logger->logdie("Cannot open --isai_clstr_file: $isai_clstr_file\n");

# create a table for saving the genes in the cluster of interest
# this table can be use to figure out what genomes have cluster X
# and what genes are in cluster X
my $isai_tbl = Table->new();
$isai_tbl->_set_col_count(3);
my @names = ("genome_id", "cluster_id", "gene_id");
$isai_tbl->_set_col_names(\@names);

my @save = ();
my $row = 0;
foreach my $line ( <$IN> ) {
	chomp $line;
	my @vals = split(/\t/, $line);

	if ( $vals[1] eq $clstr ) {
		$isai_tbl->add_row($row, \@vals);
		$row++;
		#push @save, $line;
	}
}

close($IN);


# go through each line that was saved
# find the metadata in the genome of the gene
my $annote_file = "";
my $annote_tbl = Table->new();

#my $out_tbl = Table->new();
open my $OUT, ">", $out_file or 
	$logger->logdie("Cannot open $out_file");

my $tmp;
my $i = 0;
foreach my $row ( @{$isai_tbl->get_row_names()} ) {
#foreach my $saved ( @save ) {
	$logger->debug("Working on line: $row\n");
	$logger->debug("genome: " . $isai_tbl->get_value_at($row, "genome_id"));

#	@vals = split(/\t/, $saved);

#	$annote_file = "$dafe_db/" . $vals[0] . "/all_annote.txt";
	$annote_file = "$dafe_db/" . $isai_tbl->get_value_at($row, "genome_id") . "/all_annote.txt";

	if ( ! -e $annote_file ) {
		$logger->warn("all_annote.txt does not exist: $annote_file\n");
		next;
	}

	#load the annote file
	$annote_tbl->load_from_file($annote_file);

	# go through each line in the annote file and get the information about the genes
	foreach my $row2 ( @{$annote_tbl->get_row_names()} ) {
		$tmp = $annote_tbl->get_value_at($row2, "clstr");
		if ( $tmp eq $clstr ) {
			$logger->debug("Found a cluster!");
			# add this info to the output table
			my @vals = ($isai_tbl->get_value_at($row, "genome_id"), @{$annote_tbl->get_row($row2)});
			my @names = ("genome_id", @{$annote_tbl->get_col_names()});
#			$out_tbl->add_row($i, \@vals, \@names);
			$i++;

			print $OUT $isai_tbl->get_value_at($row, "genome_id") . "\t";
			print $OUT join("\t", @{$annote_tbl->get_row($row2)}) . "\n";
		}
	}

#	my $cmd = "grep -P \"^" . $vals[2] . "\" $annote_file >> $out_file";
#	$logger->debug("Runnig cmd: $cmd\n");
#	`$cmd`;
}

#$out_tbl->save($out_file);
close($OUT);


########
# Subs #
########
#my ($clstr, $isai_clstr_file, $dafe_db, $out_file, $help, $man);
sub check_params {
	# check for required variables
	if ( ! defined $clstr) { 
		pod2usage(-message => "ERROR: required --clstr not defined\n\n",
					-exitval => 2); 
	}
	if ( ! defined $isai_clstr_file ) {
		pod2usage(-message => "ERROR: required --isai_clstr_file not defined\n\n",
					-exitval => 2);
	}
	if ( ! defined $dafe_db ) {
		pod2usage(-message => "ERROR: required --dafe_db not defined\n\n",
					-exitval => 2);
	}
	if ( ! defined $out_file ) {
		pod2usage(-message => "ERROR: required --out_file not defined\n\n",
					-exitval => 2);
	}

	# make sure required files are non-empty
	if ( defined $isai_clstr_file and ! -e $isai_clstr_file ) { 
		pod2usage(-message => "ERROR: --isai_clstr_file $isai_clstr_file is an empty file\n\n",
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

[NAME].pl - [DESCRIPTION]


=head1 VERSION

This documentation refers to version 0.0.1


=head1 SYNOPSIS

    [NAME].pl
        --clstr "Cluster_37031
        --isai_clstr_file dataframe_3987_genomes_maxaccepts0_maxrejects0_to_geneid_to_clusterid_v2.tsv
        --dafe_db my_dafe_db/
        --out_file Cluster_37031_gene_metadata.txt
        
        [--help]
        [--man]
        [--debug]
        [--verbose]
        [--quiet]
        [--logfile logfile.log]

    --clstr              Name of cluster
    --isai_clstr_file    Path to Isai's cluster file with three columns
    --dafe_db            Path to DAFE database
    --out_file           Path to output file that will hold the gene metadata
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
