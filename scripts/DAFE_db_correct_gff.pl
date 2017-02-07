#!/bin/usr/evn perl

# 1. ensures the strand type is either "+" or "-"
# 2. comments CRISPR lines by prepending a "#" to the line

use strict;
use warnings;

use File::Temp qw/ tempfile /;
use Log::Log4perl qw(:easy);
use Log::Log4perl::CommandLine qw(:all);

# set up the logging environment
my $logger = get_logger();


# usage and parameters
my $usage = "$0 <input gff> <htseq_i>\n";

my $input_gff = shift or die $usage;
my $htseq_i = shift;

# get the genome ID from the gff file name
my @genome = split /\//, $input_gff;
my $genome_id = $genome[scalar(@genome)-1];
$genome_id =~ s/.gff//;

# set defaults if neccessary
if ( ! defined $htseq_i ) { $htseq_i = "ID"; }

open my $IN, "<", $input_gff or die $!;

# create a temporary output file
my ($fh, $filename) = tempfile();

my @vals = ();
my $line_count = 0;
foreach my $line ( <$IN> ) {
	chomp $line;
	
	#skip empty lines
	if ( $line =~ m/^$/ ) {
		next;
	}
	
	# skip comment lines
	# this used to be in the fix strand block.  To preserve the comment line
	# print it before going on to the next line.
	if ( $line =~ m/^#/ ) {
		print $fh "$line\n";
		next;
	}
	
	# substitute out double quotes for single quotes
	# this part must go first
	if ( $line =~ m/"/ ) {
		$line =~ s/"/'/g;
	}
	
	# remove the ";" in the product tag field
	if ( $line =~ m/(.*product=(?!.*(=)).*);(.*)/ ) { 
		$line = "$1,$3";
	}

	# get all the seperate fields in the line
	@vals = split("\t", $line);
	
	# Change the scaffold name to include the genome ID
	#if ( $vals[0] !~ m/$genome_id/ ) {
	#	$vals[0] = $genome_id . "-" . $vals[0];
	#}
	
	# Fixes tag field
    if ( $vals[8] =~ m/name / ) { 
      $vals[8] =~ s/;\s/;/g;
      $vals[8] =~ s/\s/=/g;
    } 

	# fix strand
	if ($vals[6] eq "+" or $vals[6] eq "-") {
		; # skip lines that are correct
	}
	elsif ($vals[6] eq "-1") {
		$vals[6] = "-";
	}
	elsif ($vals[6] eq "1") {
		$vals[6] = "+";
	}
	else {
		warn "Bad strand type: $line";
	}

	# comment out the CRISPR lines
	if ( $line !~ m/^#/ and $vals[2] eq "CRISPR" ) {
		$vals[0] = "#" . $vals[0];
	}
	
	# comment out the "direct" lines
	# I'm not sure what those even are
	if ( $line !~ m/^#/ and $vals[2] eq "direct" ) {
		$vals[0] = "#" . $vals[0];
	}
	
	# comment out the "tandem" lines
	# I'm not sure what those even are
	if ( $line !~ m/^#/ and $vals[2] eq "tandem" ) {
		$vals[0] = "#" . $vals[0];
	}

	# comment out the "inverted" lines
	# I'm not sure what those even are
	if ( $line !~ m/^#/ and $vals[2] eq "inverted" ) {
		$vals[0] = "#" . $vals[0];
	}
	
	# check if the attribute values contain the htseq_i tag
	if ( $vals[8] !~ m/$htseq_i=/ ) {
		$vals[8] = "ID=" . $genome_id . "-" . $line_count . ";" . $vals[8];
	}
	
	# Change the ID of the genes in the gff file to include the genome ID
	#if ( $line !~ m/^#/ and $vals[8] !~ m/$htseq_i=\S+-\S+?/ ) {
	#	if ( $vals[8] =~ m/$htseq_i=(\S+?);/ ) {
	#		$vals[8] =~ s/$htseq_i=\S+?;/$htseq_i=$genome_id-$1;/;
	#	}
	#}
	
	print $fh (join("\t", @vals), "\n");
	
	$line_count++;
}
close($IN);
close($fh);


# now move the temp file to where the old gff was
# this is kind of dangerous but I know it works... :)
`mv $filename $input_gff`;

