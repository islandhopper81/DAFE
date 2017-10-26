package DAFE::Utils;

use warnings;
use strict;
use Readonly;
use Class::Std::Utils;
use Data::Dumper;
use List::MoreUtils qw(any);
use Log::Log4perl qw(:easy);
use Log::Log4perl::CommandLine qw(:all);
use MyX::Generic;
use version; our $VERSION = qv('0.0.2');
use UtilSY qw(:all);
use Table;
use File::Slurp;
use List::Util qw(sum);
use Exporter qw( import );
our @EXPORT_OK = qw(ids_to_names get_enr_val get_enr_code get_genomes_features is_measurable);
our %EXPORT_TAGS = (
    'all' => \@EXPORT_OK,
);

# set up the logging environment
my $logger = get_logger();

# GLOBALS #
Readonly::Scalar my $ENR => 1;   # Enriched
Readonly::Scalar my $NS => 0;    # not significant
Readonly::Scalar my $DEP => -1;  # depleted
Readonly::Scalar my $LOW => -2;  # too low to test
Readonly::Scalar my $UNM => -3;  # no reads map but contained in genome
Readonly::Scalar my $ABS => -4;  # not in genome


{
	# Usage statement
	# use each function individually

	# Attributes #
	# NA -- this object is just a set of utility functions
	
	# Functions #
	sub is_measurable;
	sub get_genomes_features;
	sub get_enr_code;
	sub get_enr_val;
	sub ids_to_names;


	###############
	# Constructor #
	###############
	# NA -- This object is just a set of utility function that can be exported
	# and used in various scripts.

	#############
	# Functions #
	#############
	sub is_measurable {
		my ($sample_meta_tbl, $test_col, $t1, $t2, $count_tbl, $feat, $min_s, $min_c) = @_;

		# sample metadata table - for spliting the samples into the test groups
		# test col - column in sample metadata use to split into test groups
		# t1 and t2 - names of the tests as specified in test_col
		# feat - the feature to test for measurability
		# min_c - min number of reads required for at least min_s samples
		# min_s - min number of samples to have at least min_c reads

		# NOTE: someday this function would be better incorperated into a DA Table object

		# get the sample names in test group 1 and 2
		my %t1_names = ();
		my %t2_names = ();
		foreach my $s ( @{$sample_meta_tbl->get_row_names()} ) {
			# each s is a sample
			if ( $t1 eq $sample_meta_tbl->get_value_at($s, $test_col) ) {
				$t1_names{$s} = 1;
			}
			elsif ( $t2 eq $sample_meta_tbl->get_value_at($s, $test_col) ) {
				$t2_names{$s} = 1;
			}
		}

		# test for measurability by looking at each col in the count table
		my $t1_pass = 0; # the number of samples with at least min_c for the t1 samples
		my $t2_pass = 0; # the number of samples with at least min_c for the t2 samples
		foreach my $c ( @{$count_tbl->get_col_names()} ) {
			if ( defined $t1_names{$c} ) {
				# this is a sample in group 1
				if ( $count_tbl->get_value_at($feat, $c) >= $min_c ) {
					$t1_pass++;
				}
			}
			elsif ( defined $t2_names{$c} ) {
				# this is a sample in group 2
				if ( $count_tbl->get_value_at($feat, $c) >= $min_c ) {
					$t2_pass++;
				}
			}
		}

		if ( $t1_pass >= $min_s and $t2_pass >= $min_s ) {
			return(1);
		}

		return(0);  # default is to return false
	}

	sub get_genomes_features {
		my ($g_aref, $dafe_db, $annote_file, $feature_type) = @_;

		my %features = ();  # hash with features
		
		my $a_file; # full path to each annotation file
		my $a_tbl = Table->new();  # Table of each genomes annotations
		my $feat;
		foreach my $g ( @{$g_aref} ) {
			$logger->debug("$g");
			$a_file = "$dafe_db/$g/$annote_file";
			$a_tbl->load_from_file($a_file);

			foreach $feat ( @{$a_tbl->get_col($feature_type)} ) {
				if ( $feat eq "NA" ) {next;}

				if ( ! defined $features{$feat} ) {
					$features{$feat} = 1;
				}
			}
		}

		my @features_arr = keys %features;

		return(\@features_arr);
	}

	sub get_enr_val {
		my ($logFC, $fdr, $up, $dn) = @_;

		# up is the factor that is up when logFC is positive
		# dn is the factor that is down when logFC is negative
		# so if I run the test BK, RZ then up == RZ and dn == BK

		# this function will return one of three values:
		# 1. $up -- enriched in $up
		# 2. $dn -- enriched in $dn
		# 3. NE -- not enriched
		
		my $code = get_enr_code($logFC, $fdr);
    
		if ( $code == 0 ) { 
			return("NE");
		}   
		elsif ( $code == 1 ) { 
			return($up);
		}   
		else {
			return($dn);
		}   
	}

	sub get_enr_code {
		my ($logFC, $fdr) = @_;

		# this function will return the code (ie 1,0,-1)
		# which indicates the enrichment value
		
		if ($fdr > 0.05 ) {
			return($NS);
		}

		if ( $logFC > 0 ) {
			return($ENR); # up
		}
		else {
			return($DEP); # down
		}
	}

	sub ids_to_names {
		my ($params_href) = @_;
		
		$logger->debug("Calling DAFE::Utils::ids_to_names");

		my $usage = "ids_to_names({[str] | [in_file, out_file], meta_file})";
		my $str = $params_href->{"str"};
		my $meta_file = $params_href->{"meta_file"};
		my $in_file = $params_href->{"in_file"};
		my $out_file = $params_href->{"out_file"};

		$logger->debug("str: $str") if defined $str;
		$logger->debug("in_file: $in_file") if defined $in_file;
		$logger->debug("out_file: $out_file") if defined $out_file;
		$logger->debug("meta_file: $meta_file");  # this is required.

		
		# load the metadata file
		check_file($meta_file);
		my $meta_tbl = Table->new();
		$meta_tbl->load_from_file($meta_file);

		if ( is_defined($str) ) {
			$logger->debug("Opperationg on string");
			return(_ids_to_names_str($str, $meta_tbl));
		}
		elsif ( is_defined($in_file) ) {
			$logger->debug("Opperating on file");
	
			check_file($in_file);

			# slurp in the in_file
			my $text = read_file($in_file);

			$text = _ids_to_names_str($text, $meta_tbl);

			open my $OUT, ">", $out_file
				or $logger->logdie("Cannot open --out_file: $out_file");

			print $OUT $text;
			
			close($OUT);		
		
			return 1;

		}
		else {
			$logger->logdie("Bad Function call: $usage");
		}

		return 0;  # on Fail
	}

	sub _ids_to_names_str {
		my ($str, $meta_tbl) = @_;
		$logger->debug("str: $str");

		foreach my $r ( @{$meta_tbl->get_row_names()} ) { 
		#	$logger->debug("Looking for ID: $r");
			my $name = $meta_tbl->get_value_at($r, "Name");
			$str =~ s/$r/$name/g;
		}
		
		return($str);	
	}

}

1; # Magic true value required at end of module
__END__

=head1 NAME

UtilSY - Scott Yourstone's utility functions


=head1 VERSION

This document describes UtilSY version 0.0.2


=head1 SYNOPSIS

use UtilSY qw(:all);
	
# Or you can load and use the functions individually. E.g.
use UtilSY qw(is_defined);
  
=head1 DESCRIPTION

These are generic functions that I tend to use in many different scripts and
modules.  It seemed safer and easier to keep them in an independent Perl object
rather than keep copying and pasting the code.


=head1 CONFIGURATION AND ENVIRONMENT
  
UtilSY requires no configuration files or environment variables.


=head1 DEPENDENCIES

	warnings
	strict
	Readonly
	Class::Std::Utils
	List::MoreUtils qw(any)
	Log::Log4perl qw(:easy)
	Log::Log4perl::CommandLine qw(:all)
	MyX::Generic
	version our $VERSION = qv('0.0.2')
	Exporter qw( import )


=head1 INCOMPATIBILITIES

None reported.


=head1 METHODS

=over

	is_defined
	check_defined
	to_bool
	check_ref
	check_file
	load_lines

=back

=head1 METHODS DESCRIPTION
	
=head2 is_defined

	Title: is_defined
	Usage: is_defined($val, $val_name)
	Function: checks if a value is defined
	Returns: boolean (ie 0 or 1)
	Args: -val => value to check
	      -val_name => name of value (for print error message)
	Throws: NA
	Comments: This function will NOT throw a warning if the value is not defined
	See Also: NA
	
=head2 check_defined

	Title: check_defined
	Usage: check_defined($val, $val_name)
	Function: checks if a value is defined
	Returns: 1 on success
	Args: -val => value to check
	      -val_name => name of value (for print error message)
	Throws: MyX::Generic::Undef::Param
	Comments: This function WILL throw a warning if the value is not defined
	See Also: NA
	
=head2 to_bool

	Title: to_bool
	Usage: to_bool($val)
	Function: converts a boolean value to either 0 or 1
	Returns: 0 or 1
	Args: -val => a boolean value
	Throws: MyX::Generic::Undef::Param
	        MyX::Generic::BadValue
	Comments: Valid boolean values include: 0, 1, T, True, Y, Yes, F, False,
			  No, and N.
	See Also: NA
	
=head2 check_ref

	Title: check_ref
	Usage: check_ref($ref, $type)
	Function: checks a value to ensure it is of the correct type
	Returns: 1 on success
	Args: -ref => an object reference
	      -type => the expected type
	Throws: MyX::Generic::Undef::Param
	        MyX::Generic::Ref::UnsupportedType
	Comments: Types are case sensitive
	See Also: perl ref function
	
=head2 check_file

	Title: check_file
	Usage: check_file($file)
	Function: checks a file to make sure it exists and is non-empty
	Returns: 1 on success
	Args: -file => file name or path
	Throws: MyX::Generic::Undef::Param
	        MyX::Generic::DoesNotExist::File
			MyX::Generic::File::Empty
	Comments: NA
	See Also: NA
	
=head2 load_lines

	Title: load_lines
	Usage: load_lines($file, $sep)
	Function: loads an array from a file
	Returns: array ref
	Args: -file => file name or path
	      [-sep => delimiter]
	Throws: MyX::Generic::Undef::Param
	        MyX::Generic::DoesNotExist::File
			MyX::Generic::File::Empty
	Comments: This function does the common operation of loading the lines of a
	          file into an array.  The $sep parameter is optional and when
			  specified it does not return anything in the line past the first
			  instance of $sep.  This can be useful when you have a table and
			  only want to load the first column (ie the row names).
	See Also: NA


=head1 BUGS AND LIMITATIONS

No bugs have been reported.

=head1 TO DO

None

=head1 AUTHOR

Scott Yourstone  C<< scott.yourstone81@gmail.com >>

=head1 LICENCE AND COPYRIGHT

Copyright (c) 2013, Scott Yourstone
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

