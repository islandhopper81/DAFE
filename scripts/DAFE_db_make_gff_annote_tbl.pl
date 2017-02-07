#!/usr/bin/env perl

# make a table like txt file with the GFF attributes

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Carp;
use Readonly;
use version; our $VERSION = qv('0.0.1');
use Log::Log4perl qw(:easy);
use Log::Log4perl::CommandLine qw(:all);

# Subroutines #
sub check_params;
sub _is_defined;

# Variables #
my ($dafe_db, $gene_file_suffix, $gff_file_suffix, $gff_feature, $genome,
	$gene_id_col, $help, $man);

my $options_okay = GetOptions (
    "dafe_db:s" => \$dafe_db,
    "gene_file_suffix:s" => \$gene_file_suffix,
	"gff_file_suffix:s" => \$gff_file_suffix,
	"gff_feature:s" => \$gff_feature,
	"genome:s" => \$genome,
	"gene_id_col:s" => \$gene_id_col,
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
# get a list of genomes to make the gff_annote.tab.txt file for
$logger->debug("Reading in genomes");
my @genome_ids = ();
if ( defined $genome ) {
    push @genome_ids, $genome;
}
else {
    opendir my $DIR, $dafe_db or
        $logger->logdie("Cannot open --dafe_db ($dafe_db)");
    
    while ( my $id = readdir($DIR) ) {
		if ( $id =~ m/^\./ ) {next;} #skip hidden dirs
		if ( ! -d "$dafe_db/$id" ) { next; } # skip non-directories
        push @genome_ids, $id;
    }
	
	closedir($DIR);
}

# foreach genome save all the ID and name data
$logger->debug("Operating on each genome");
foreach my $id ( @genome_ids ) {
	$logger->debug("Starting genome: $id");
	
	# make variables to store important info that will be printed
	my %uniq_tags = ();
	my @tbl = (); # final table stored as a list of hashes w/ tag-value pairs
	
	# parse through the GFF file to get the attribute tag-value pairs
	my $gff_file = "$dafe_db/$id/$id$gff_file_suffix";
	if ( -s $gff_file ) {
		$logger->debug("Looking at GFF file: $gff_file");
		
		# open the gff file
		open my $GFF, "<", $gff_file or
			$logger->logdie("Cannot open gff file: $gff_file");
		
		# go through each line in the gff file
		foreach my $line ( <$GFF> ) {
			chomp $line;
			
			# skip comment lines
			if ( $line =~ m/^\s*#/ ) { next; }
			
			my @vals = split(/\t/, $line);
			
			# make sure the line has 9 fields as all GFF fiels should.
			if ( @vals != 9 ) {
				$logger->warn("Line has more or less than 9 fields: $line");
			}
			
			# skip lines that are not of the correct feature
			if ( $vals[2] !~ m/$gff_feature/ ) { next; }

			
			# get the attribute field tag-value pairs
			my @pairs = split(/;/, $vals[8]);
			
			# parse each pair
			my %pair_hash = ();
			foreach my $pair ( @pairs ) {
				if ( $pair =~ m/(\S+?)=(.*)/) {
					if ( ! defined $uniq_tags{$1} ) {
						$uniq_tags{$1} = 1;
						$logger->debug("Adding unique tag: $1");
					}
					
					# save the pairs in the hash
					$pair_hash{$1} = $2;
				}
				else {
					$logger->warn("Doesn't look like a tag-value pair: $pair");
				}
			}
			
			# add all the pairs for this line to the tbl
			push @tbl, \%pair_hash;
		}
		
		close($GFF);
	}
	else {
		$logger->warn("Genome ($id) missing gff_file: $gff_file");
	}
	
	# remove the exonNumber tag because it causes problems later.  When I use
	# the output file created here in DAFE_db_make_annote_table.pl the
	# exonNumber field is wrong.
	if ( defined $uniq_tags{"exonNumber"} ) {
		delete $uniq_tags{"exonNumber"};
	}
	
	# make sure the gene_id_col is one of the keys in the uniq_tags hash
	if ( ! defined $uniq_tags{$gene_id_col} ) {
		$logger->logdie("--gene_id_col ($gene_id_col) is not a GFF attribute");
	}
	
	# set up the headers order
	my @headers = keys %uniq_tags;
	
	# put gene_id_col first.  First delete it from the current headers arr.
	# then unshift it back onto the front.
	my $index = 0;
	$index++ until $headers[$index] eq $gene_id_col;
	splice(@headers, $index, 1);
	unshift @headers, $gene_id_col;	
	
	# print out all the attributes to the output file
	my $out_file = "$dafe_db/$id/gff_annote.tab.txt";
	open my $OUT, ">", $out_file or
		$logger->logdie("Cannot open output file: $out_file");
	
	# print the headers
	my $str = join("\t", @headers);
	print $OUT $str, "\n";
	
	# print the values in the table
	foreach my $gene ( @tbl ) {
		my @vals = ();
		foreach my $tag ( @headers ) {
			if ( defined $gene->{$tag} ) {
				push @vals, $gene->{$tag};
			}
			else {
				push @vals, "NA";
			}
		}
		
		my $str = join("\t", @vals);
		print $OUT $str, "\n";
	}
	
	close($OUT);
}


########
# Subs #
########
sub check_params {
	# check for required variables
	if ( ! defined $gene_file_suffix) {
        $gene_file_suffix = ".genes.fna";
		#$logger->info("Setting --gene_file_suffix to $gene_file_suffix")
	}
	if ( ! defined $gff_file_suffix) {
        $gff_file_suffix = ".gff";
		$logger->info("Setting --gff_file_suffix to $gff_file_suffix")
	}
	
	if ( ! defined $gff_feature ) {
		$gff_feature = "CDS";
		$logger->info("Setting --gff_feature to $gff_feature");
	}
	
	if ( ! defined $gene_id_col ) {
		$gene_id_col = "ID";
		$logger->info("Setting --gene_id_col to $gene_id_col");
	}

	# make sure required directories exist
    if ( ! defined $dafe_db ) {
        pod2usage(-message => "ERROR: --dafe_db is not defined\n\n",
					-exitval => 2); 
    }
	if ( ! -d $dafe_db ) { 
		pod2usage(-message => "ERROR: --dafe_db is not a directory\n\n",
					-exitval => 2); 
	}
    
    $logger->info("--dafe_db: $dafe_db");
    #$logger->info("--gene_file_suffix: $gene_file_suffix");
	$logger->info("--gff_file_suffix: $gff_file_suffix");
	$logger->info("--gff_feature: $gff_feature");
	$logger->info("--genome $genome") if (defined $genome);
	$logger->info("--gene_id_col: $gene_id_col");
	
	return 1;
}

sub _is_defined {
    my ($val, $default) = @_;
    
    if ( defined $val ) {
        return $val;
    }
    elsif ( defined $default ) {
        return $default;
    }
    else {
        return undef;
    }
}


__END__

# POD

=head1 NAME

make_db_gff_annote_tbl.pl - make a table like txt file with the GFF attributes


=head1 VERSION

This documentation refers to version 0.0.1


=head1 SYNOPSIS

    make_db_gff_annote_tbl.pl
        --dafe_db my_dafe_db/
        [--gene_file_suffix ".gene.fna"]  DEPRECIATED!
        [--gff_file_suffix ".gff"]
        [--gff_feature "CDS"]
        [--genome 2582580847]
        [--gene_id_col "ID"]
        
        [--help]
        [--man]
        [--debug]
        [--verbose]
        [--quiet]
        [--logfile logfile.log]

    --dafe_db                   Path to DAFE database directory
    --gene_file_suffix          DEPRECIATED!
    --gff_file_suffix           Suffix of gff file in each genome dir
    --gff_feature               The feature field that you want to keep
    --genome                    Genome to make the gff_annote.tabl.txt file for
    --gene_id_col               Name of the attribute in GFF for the gene ID
    --help | -h                 Prints USAGE statement
    --man                       Prints the man page
    --debug                     Prints Log4perl DEBUG+ messages
    --verbose                   Prints Log4perl INFO+ messages
    --quiet                     Suppress printing ERROR+ Log4perl messages
    --logfile                   File to save Log4perl messages


=head1 ARGUMENTS
    
=head2 --dafe_db

Path to DAFE database directory
    
=head2 [--gene_file_suffix]

DEPRECIATED.  It is now recommended that you use the --gff_file_suffix.

In the DAFE database each genome directory should contain some type of fasta
file (faa, fna, etc) that has headers where the first value is the gene ID and
the second value is the gene name
(ie 2565986422 BR39DRAFT_0002 hypothetical protein [Caulobacter sp. UNC358MFTsu5.1]).
These files should all have a common suffix that is passed in using the
--gene_file_suffix parameter
DEFAULT: .genes.fna

=head2 [--gff_file_suffix]

In the DAFE database each genome directory should contain a GFF file.  This GFF
file will have attributes as tag-value pairs in the last column.  The output
file of this script has a column for each unique tag.  The values in the table
are the values for a given gene.
DEFAULT: .gff

=head2 [--gff_feature]

The GFF feature field that you want to keep.  All other features are ignored.
For example, if -gff_feature "CDS" is used only the CDS lines are use to
generate output (ie "gene" features are ignored).  Right now only one feature
can be used.  If you are running this as a precursor to the DAFE suite the
--gff_feature should be set to "CDS".
DEFAULT: "CDS"

=head2 [--genome]

The genome ID for which to make the gff_annote.tab.txt file.  The ID must match
the directory name in the DAFE database.  If --genome is not provided this
script a gff_annote.tab.txt file is created for each dir in the DAFE database.

=head2 [--gene_id_col]

The name of the attribute in the GFF file to use for the gene ID.  This
parameter should also match the -i parameter when running htseq-count and the
other htseq_i parameters in the DAFE suite.  It should end up being the first
column in the all_annote.txt
DEFAULT: "ID"
 
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

This script replaces the bash version (make_db_geneID_name_tbl.sh).  It is also
a little more robust than the bash version.  The bash version should not be used
from now on.

This script uses a GFF file in each directory (ie genome) of a DAFE database to
create a file with all the attribute features in a table format.  There should
be a line in this file for each gene (ie CDS) sequnece in the genome.  In the
DAFE suite this file is combined with the different gene annotation files to
create the all_annote.txt file.  The first column is the column named
--gene_id_col where --gene_id_col is a parameter and the name of an attribute in
the GFF file.  Ths column should be a uniq identifier for each gene.  It MUST be
the first column.

Previously, this script used a fasta file in each directory (ie genome) of a
DAFE database to create a file that maps the gene IDs of genes in a genome to
their gene name.  This file is required for creating the all_annote.txt file if
you are creating it using the annotation files recieved from the JGI.

If you create the all_annote.txt file using only the genbank file this script
should not be ran.  However, it is advised to use the JGI annotation files
rather than making the annotation files from the genbank files.  If you make the
all_annote.txt file using only the genbank files those annotation will not be
perfectly consistant with JGI annotations.

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
