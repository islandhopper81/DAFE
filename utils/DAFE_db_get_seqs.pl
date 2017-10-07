#!/usr/bin/env perl

# given a list of genome id and gene ids get the sequences in fasta format

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
use BioUtils::FastaIO;
use BioUtils::FastaSeq;

# Subroutines #
sub check_params;
sub _is_defined;

# Variables #
my ($query_file, $out_file, $dafe_db, $help, $man);

my $options_okay = GetOptions (
    "query_file:s" => \$query_file,
    "out_file:s" => \$out_file,
	"dafe_db:s" => \$dafe_db,
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
# read in the queries.
# these should be in a table where the headers are genome_id and gene_id
my $q_tbl = Table->new();
$q_tbl->load_from_file($query_file);

if ( ! $q_tbl->has_col("genome_id") ) {
	$logger->logdie("--query_file does not have the genome_id column");
}
if ( ! $q_tbl->has_col("gene_id") ) {
	$logger->logdie("--query_file does not have the gene_id column");
}

# go through each query
my $genome;
my $gene;
my $gene_file;
my $in;
my $out = BioUtils::FastaIO->new({stream_type => '>', file => $out_file});
my $seq;
my $new_header;
foreach my $row ( @{$q_tbl->get_row_names()} ) {
	$genome = $q_tbl->get_value_at($row, "genome_id");
	$gene = $q_tbl->get_value_at($row, "gene_id");

	$logger->debug("Query genome: $genome");
	$logger->debug("Query gene: $gene");

	$gene_file = "$dafe_db/$genome/$genome\.gene.fna";
	if ( ! -e $gene_file ) {
		$logger->warn("Cannot find file for genome: $genome.  Skipping!");
		next;
	}

	$in = BioUtils::FastaIO->new({stream_type => '<', file => $gene_file});

	while ( $seq = $in->get_next_seq() ) {
		$logger->debug("Looking at seq; " . $seq->get_header());
		if ( $seq->get_header() =~ m/$gene[- ]/i ) {
			$logger->debug("Writing seq: $genome/$gene");

			# add the genome id to the header
			$new_header = $genome . "-" . $seq->get_header();
			$seq->set_header($new_header);
			$out->write_seq($seq);
			next;
		}
	}

	# If I get to this point then the sequence was not found
	$logger->warn("Sequence not found: $genome/$gene");
}



########
# Subs #
########
sub check_params {
	# check for required variables
	if ( ! defined $query_file) { 
		pod2usage(-message => "ERROR: required --query_file not defined\n\n",
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
	if ( defined $query_file and ! -e $query_file ) { 
		pod2usage(-message => "ERROR: --query_file $query_file is an empty file\n\n",
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

DAFE_db_get_seqs.pl - given a genome and gene id get the sequence(s)


=head1 VERSION

This documentation refers to version 0.0.1


=head1 SYNOPSIS

    DAFE_db_get_seqs.pl
        --query_file ids_to_get.txt
        --dafe_db my_dafe_db/
		--out_file query_seqs.fasta
        
        [--help]
        [--man]
        [--debug]
        [--verbose]
        [--quiet]
        [--logfile logfile.log]

    --query_file    Path to a file with genome ids and gene ids to get
    --dafe_db       Path to the dafe database dir
	--out_file      Path to output fasta file
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
