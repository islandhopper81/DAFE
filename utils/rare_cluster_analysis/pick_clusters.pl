#!/usr/bin/env perl

# picks X clusters from the ENR, DEP, and NS groups

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

# Variables #
my ($cluster_sum_file, $count, $out_file, $help, $man);

# GLOBALS #
Readonly::Scalar my $ENR => 1;   # Enriched
Readonly::Scalar my $NS => 0;    # not significant
Readonly::Scalar my $DEP => -1;  # depleted
Readonly::Scalar my $LOW => -2;  # too low to test
Readonly::Scalar my $UNM => -3;  # no reads map but contained in genome
Readonly::Scalar my $ABS => -4;  # not in genome

my $options_okay = GetOptions (
    "cluster_sum_file:s" => \$cluster_sum_file,
	"count:s" => \$count,
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
# load the cluster summary file
$logger->info("Loading cluster summary file");
my $clstr_tbl = Table->new();
$clstr_tbl->load_from_file($cluster_sum_file);

# get the array of all the ENR, DEP, and NS clusters
$logger->info("Seperating clusters by enrichment status");
my @enr_clusters = ();
my @dep_clusters = ();
my @ns_clusters = ();
my $val;
my $is_da = 0;
foreach my $cluster ( @{$clstr_tbl->get_row_names()} ) {
	# set a boolean marker to see if this cluster is called ENR or DEP
	$is_da = 0; # start as false
	
	$val = $clstr_tbl->get_value_at($cluster, "DA_up_count");
	if ( $val > 0 ) {
		push @enr_clusters, $cluster;
		$is_da++;
	}
	
	$val = $clstr_tbl->get_value_at($cluster, "DA_dn_count");
	if ( $val > 0 ) {
		push @dep_clusters, $cluster;
		$is_da++;
	}
	
	if ( $is_da == 0 ) {
		push @ns_clusters, $cluster;
	}
}


# create a table to save the sampled clusters
my $out_tbl = Table->new();


# sample the ENR clusters without replacement
$logger->info("Randomly sampling clusters");
my $sampled_cluster;
my @col_names = ("DA_call");
for ( 1 .. $count ) {
	my @vals = ("ENR");
	$sampled_cluster = splice @enr_clusters, rand @enr_clusters, 1;
	
	if ( ! $out_tbl->has_row($sampled_cluster) ) {
		$out_tbl->add_row($sampled_cluster, \@vals, \@col_names);
	}
}

# sample the DEP clusters without replacement
for ( 1 .. $count ) {
	my @vals = ("DEP");
	$sampled_cluster = splice @dep_clusters, rand @dep_clusters, 1;
	
	if ( ! $out_tbl->has_row($sampled_cluster) ) {
		$out_tbl->add_row($sampled_cluster, \@vals, \@col_names);
	}
}

# sample the NS clusters without replacement
for ( 1 .. $count ) {
	my @vals = ("NS");
	$sampled_cluster = splice @ns_clusters, rand @ns_clusters, 1;
	$out_tbl->add_row($sampled_cluster, \@vals, \@col_names);
}

# save the output table
$out_tbl->save($out_file, "\t", "F");  # don't save the col header

########
# Subs #
########
sub check_params {
	# check for required variables
	if ( ! defined $cluster_sum_file) { 
		pod2usage(-message => "ERROR: required --cluster_sum_file not defined\n\n",
					-exitval => 2); 
	}
	if ( ! defined $count ) {
		pod2usage(-message => "ERROR: required --count not defined\n\n",
					-exitval => 2);
	}
	if ( ! defined $out_file ) {
		pod2usage(-message => "ERROR: required --out_file not defined\n\n",
					-exitval => 2);
	}

	# make sure required files are non-empty
	if ( defined $cluster_sum_file and ! -e $cluster_sum_file ) { 
		pod2usage(-message => "ERROR: --cluster_sum_file $cluster_sum_file is an empty file\n\n",
					-exitval => 2);
	}

	# make sure required directories exist
	#if ( ! -d $dir ) { 
	#	pod2usage(-message => "ERROR: --dir is not a directory\n\n",
	#				-exitval => 2); 
	#}
	
	return 1;
}


__END__

# POD

=head1 NAME

pick_clusters.pl - picks X random clusters from the ENR, DEP and NS groups


=head1 VERSION

This documentation refers to version 0.0.1


=head1 SYNOPSIS

    pick_clusters.pl
        --cluster_sum_file cluster_summary.txt
        --count 100
        --out_file random_clusters.txt
        
        [--help]
        [--man]
        [--debug]
        [--verbose]
        [--quiet]
        [--logfile logfile.log]

    --cluster_sum_file     Path to cluster summary file
    --count                Number of random clusters to get for each group
    --out_file             Path to output file
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
