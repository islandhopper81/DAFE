#!/usr/bin/env perl

use strict;
use warnings;

use File::Temp qw(tempfile);
use Log::Log4perl qw(:easy);
use Log::Log4perl::CommandLine qw(:all);

# set up the logging environment
my $logger = get_logger();

my $usage = "$0 <dafe db>\n";

my $dafe_db = shift or die $usage;
my $genome = shift or undef;

my @file_ext = (
	".cog.tab.txt",
	".ko.tab.txt",
	".pfam.tab.txt",
	".tigrfam.tab.txt"
);	

my %ids = (
	"cog_id" => 1,
	"ko_id" => 1,
	"pfam_id" => 1,
	"tigrfam_id" => 1
);


# go through each dir in the dafe db directory
opendir(my $DH, $dafe_db) or die $!;

foreach my $name ( readdir($DH) ) {
	# skip hidden files
	if ( $name =~ m/^\./ ) {next;}

	# for only doing one genome
	if ( defined $genome and $genome ne $name ) {
		next;
	}
	
	my $test_file = "$dafe_db/$name/$name" . $file_ext[0];
	$logger->debug("test file: $test_file");

	open my $IN, "<", $test_file or next;

	my $first_line = <$IN>;
	close($IN);

	my @vals = split(/\t/, $first_line);

	if ( scalar @vals > 2 ) {
		$logger->info("needs correcting: $name");

		foreach my $ext ( @file_ext ) {
			# open a temp file for output
			my ($fh, $filename) = tempfile();
			#$logger->info("tmp file: $filename");

			# read in the file
			# correct only the files that I can find
			my $file = "$dafe_db/$name/$name" . $ext;
			open $IN, "<", $file or next; 
			
			# parse the header
			my $first = <$IN>;
			my @first_vals = split(/\t/, $first);
			my $i = 0;
			my $index;
			foreach my $v ( @first_vals ) {
				if ( defined $ids{$v} ) {
					$index = $i;
				}
				$i++;
			}

			if ( ! defined $index ) {
				$logger->warn("Cannot find index in header: $name, $ext");
			}

			# print the rest of the lines
			foreach my $line ( <$IN> ) {
				chomp $line;
				@vals = split(/\t/, $line);
				print $fh $vals[0], "\t", $vals[$index], "\n";
			}

			close($fh);
			system("mv $filename $file");
		}
	}
}
