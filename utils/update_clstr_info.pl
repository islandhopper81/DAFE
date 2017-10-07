#!/usr/bin/env perl

use strict;
use warnings;

sub find_clstr_index;

my $usage = "$0 <DAFE db> <clstr file exten> <all_annote file>\n";

my $dafe_db = shift or die $usage;
my $clstr_file_ext = shift or die $usage;
my $annote_file = shift or die $usage;


# go through each directory in the db
opendir(my $DIR, "$dafe_db") or die $!;

my $annote_path;
my $clstr_path;
foreach my $id ( readdir($DIR) ) {
	if ( $id =~ m/^\./ ) { next; } # skip hidden files

	$annote_path = "$dafe_db/$id/$annote_file";
	$clstr_path = "$dafe_db/$id/$id$clstr_file_ext";

	print "all_annote file: $annote_path\n";
	print "cluster file: $clstr_path\n";

	# read in the cluster file
	open my $CLS, "<", $clstr_path or die $!;
	my %clstr_h = ();
	my @vals = ();
	foreach my $line ( <$CLS> ) {
		chomp $line;
		
		@vals = split(/\t/, $line);
		%clstr_h{$vals[0]} = $vals[1];
	}
	close($CLS);

	# update the all_annote file
	# writes temp file to current working directory
	open my $OUT, ">", "tmp" or die $!;
	
	# open the all_annote file
	open my $ANN, "<", $annote_path or die $!;
	my $headers = <$ANN>;
	print $OUT $headers;
	my $clstr_i = find_clstr_index($headers);

	foreach my $line ( <$ANN> ) {
		chomp $line;
		
		@vals = split(/\t/, $line);
		if ( ! defined $clstr_h{$vals[0]} ) { 
			warn "Cannot find gene ID: $vals[0]";
			next;
		}

		@vals[$clstr_i] = %clstr_h{$vals[0]};

		print $OUT join(@vals), "\n";
	}

	close($ANN);
	close($OUT);
}


sub find_clstr_index {
	my ($headers) = @_;

	my @vals = split(/\t/, $headers);

	my $i = 0;
	foreach my $v ( @vals ) {
		if ( $v =~ m/clstr/ ) {
			return($i);
		}
		$i++;
	}

	warn "Cannot find clstr index in header: $headers";
	die;
}
