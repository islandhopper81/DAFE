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
use Exporter qw( import );
our @EXPORT_OK = qw(ids_to_names);
our %EXPORT_TAGS = (
    'all' => \@EXPORT_OK,
);

# set up the logging environment
my $logger = get_logger();

{
	# Usage statement
	# use each function individually

	# Attributes #
	# NA -- this object is just a set of utility functions
	
	# Functions #
	sub ids_to_names;


	###############
	# Constructor #
	###############
	# NA -- This object is just a set of utility function that can be exported
	# and used in various scripts.

	#############
	# Functions #
	#############
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

