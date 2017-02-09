package DAFE::Count::OutputRestruct;

use warnings;
use strict;
use Carp;
use Readonly;
use Class::Std::Utils;
use List::MoreUtils qw(any);
use Log::Log4perl qw(:easy);
use Log::Log4perl::CommandLine qw(:all);
use MyX::Generic;
use UtilSY qw(:all);
use version; our $VERSION = qv('0.0.1');
use Scalar::Util qw(openhandle);
use DataObj;

# set up the logging environment
my $logger = get_logger();

{
	# Usage statement
	Readonly my $NEW_USAGE => q{ new( {
		out_dir => ,
		new_out_dir => ,
		ref_aref => ,
		sample_aref => ,
		perc_ids_aref => ,
		htseq_file_prefix => ,
	})};

	# Attributes #
	my %params_of;
	
	# Getters #
	sub get_param;

	# Setters #
	sub set_param;

	# Others #
	sub restructure;
	sub _mkdir_structure;
	sub _reformat;
	sub _parse_file;
	sub _get_htseq_file;


	###############
	# Constructor #
	###############
	sub new {
		my ($class, $arg_href) = @_;

		# Croak if calling new on already blessed reference
		croak 'Constructor called on existing object instead of class'
			if ref $class;

		# Make sure the required parameters are defined
		if ( any {!defined $_}
				$arg_href->{out_dir},
				$arg_href->{new_out_dir},
				$arg_href->{ref_aref},
				$arg_href->{sample_aref},
				$arg_href->{perc_ids_aref},
				$arg_href->{htseq_file_prefix},
			) {
			MyX::Generic::Undef::Param->throw(
				error => 'Undefined parameter value',
				usage => $NEW_USAGE,
			);
		}

		# Bless a scalar to instantiate an object
		my $new_obj = bless \do{my $anon_scalar}, $class;

		# Set Attributes
		my $data_obj = DataObj->new();
		$data_obj->set_feature("out_dir", $arg_href->{out_dir});
		$data_obj->set_feature("new_out_dir", $arg_href->{new_out_dir});
		$data_obj->set_feature("ref_aref", $arg_href->{ref_aref});
		$data_obj->set_feature("sample_aref", $arg_href->{sample_aref});
		$data_obj->set_feature("perc_ids_aref", $arg_href->{perc_ids_aref});
		$data_obj->set_feature("htseq_file_prefix",
							   $arg_href->{htseq_file_prefix});
		$params_of{ident $new_obj} = $data_obj;

		return $new_obj;
	}

	###########
	# Getters #
	###########
	sub get_param {
		my ($self, $name) = @_;
		
		return $params_of{ident $self}->get_feature($name);
	}

	###########
	# Setters #
	###########
	sub set_param {
		my ($self, $name, $val) = @_;
		
		$params_of{ident $self}->set_feature($name, $val);
		
		return 1;
	}

	##########
	# Others #
	##########
	sub restructure {
		my ($self) = @_;
		
		$logger->debug("Starting Restructuring");
		
		# create the new directory stucture
		$self->_mkdir_structure();
		
		# reformat the output
		$self->_reformat();
	}
	
	sub _mkdir_structure {
		my ($self) = @_;
		
		$logger->debug("Making new dir structure");
		
		my $new_out_dir = $self->get_param("new_out_dir");
		my $ref_aref = $self->get_param("ref_aref");
		my $sample_aref = $self->get_param("sample_aref");
		
		# create the reformated output directory structure
		foreach my $ref ( @{$ref_aref} ) {
			foreach my $sample ( @{$sample_aref} ) {
				`mkdir -p $new_out_dir/$ref/$sample/`;
			}
		}
		
		return 1;
	}
	
	sub _reformat {
		my ($self) = @_;
		
		$logger->debug("Starting the reforamtting");
		
		my $sample_aref = $self->get_param("sample_aref");
		my $ref_aref = $self->get_param("ref_aref");
		my $perc_ids_aref = $self->get_param("perc_ids_aref");
		my $out_dir = $self->get_param("out_dir");
		my $new_out_dir = $self->get_param("new_out_dir");
		my $htseq_file_prefix = $self->get_param("htseq_file_prefix");
		
		foreach my $sample ( @{$sample_aref} ) {
			foreach my $perc_id ( @{$perc_ids_aref} ) {
				my $htseq_file = _get_htseq_file($out_dir, $sample,
												 $htseq_file_prefix, $perc_id);
				
				_parse_file($htseq_file, $new_out_dir, $sample,
							$htseq_file_prefix, $perc_id);
			}
		}
		
		return 1;
	}
	
	sub _parse_file {
		my ($htseq_file, $new_out_dir, $sample,
			$htseq_file_prefix, $perc_id) = @_;
		
		# open the file
		open my $HTS, "<", $htseq_file or
			$logger->logdie("Cannot open htseq file: $htseq_file");
		
		my $current_ref = "";
		my $FH;  #output file handle
		my $genome_id;
		my $gene_id;
		my $count;
		my $new_file;
		foreach my $line ( <$HTS> ) {
			if ( $line =~ m/(\S+)-(\S+)\t(\d+)/ ) {
				$genome_id = $1;
				$gene_id = $2;
				$count = $3;
				
				# check the filehandle and open a new one if needed
				if ( $genome_id ne $current_ref ) {
					if ( openhandle($FH) ) {
						close($FH);
					}
					
					$new_file = "$new_out_dir/$genome_id/$sample/";
					$new_file .= "$htseq_file_prefix\_id" . $perc_id . ".txt";
					open $FH, ">", $new_file or
						$logger->logdie("Cannot open file: $new_file");
					$current_ref = $genome_id;
				}
				
				# print the line
				if ( $genome_id =~ m/^99/ ) {
					# if the genome is one of the plate scrapes keep the
					# genome id in the name because it is in the gff file
					print $FH "$genome_id-$gene_id\t$count\n";
				}
				else {
					print $FH "$gene_id\t$count\n";
				}
			}
			else {
				$logger->warn("Bad line match: $line");
			}
		}
		
		return 1;
	}
	
	sub _get_htseq_file {
		my ($out_dir, $s, $h, $p) = @_;
		
		my $htseq_file = "$out_dir/$s/$h\_id" . $p . ".txt";
		
		# warn if the file doesn't exist
		if ( ! -s $htseq_file ) {
			$logger->warn("$htseq_file not found or empty");
		}
		
		$logger->debug("htseq_file: $htseq_file");
		
		return($htseq_file);
	}
}

1; # Magic true value required at end of module
__END__

=head1 NAME

<MODULE NAME> - [One line description of module's purpose here]


=head1 VERSION

This document describes <MODULE NAME> version 0.0.1


=head1 SYNOPSIS

    use <MODULE NAME>;

=for author to fill in:
    Brief code example(s) here showing commonest usage(s).
    This section will be as far as many users bother reading
    so make it as educational and exeplary as possible.
  
  
=head1 DESCRIPTION

=for author to fill in:
    Write a full description of the module and its features here.
    Use subsections (=head2, =head3) as appropriate.


=head1 DIAGNOSTICS

=for author to fill in:
    List every single error and warning message that the module can
    generate (even the ones that will "never happen"), with a full
    explanation of each problem, one or more likely causes, and any
    suggested remedies.

=over

=item C<< Error message here, perhaps with %s placeholders >>

[Description of error here]

=item C<< Another error message here >>

[Description of error here]

[Et cetera, et cetera]

=back


=head1 CONFIGURATION AND ENVIRONMENT

=for author to fill in:
    A full explanation of any configuration system(s) used by the
    module, including the names and locations of any configuration
    files, and the meaning of any environment variables or properties
    that can be set. These descriptions must also include details of any
    configuration language used.
  
<MODULE NAME> requires no configuration files or environment variables.


=head1 DEPENDENCIES

=for author to fill in:
    A list of all the other modules that this module relies upon,
    including any restrictions on versions, and an indication whether
    the module is part of the standard Perl distribution, part of the
    module's distribution, or must be installed separately. ]

Carp
Readonly
Class::Std::Utils
List::MoreUtils qw(any)
Log::Log4perl qw(:easy)
Log::Log4perl::CommandLine qw(:all)
MyX::Generic
version; our $VERSION = qv('0.0.1')


=head1 INCOMPATIBILITIES

=for author to fill in:
    A list of any modules that this module cannot be used in conjunction
    with. This may be due to name conflicts in the interface, or
    competition for system or program resources, or due to internal
    limitations of Perl (for example, many modules that use source code
    filters are mutually incompatible).

None reported.


=head1 METHODS

=over

=for author to fill in:
	A list of method names in the module
	
	new

=back

=head1 METHODS DESCRIPTION

=head2 new

	Title: new
	Usage: <MODULE NAME>->new({
				arg1 => $arg1,
				arg2 => $arg2
			});
	Function: Build new <MODULE NAME>
	Returns: <MODULE NAME>
	Args: -arg1 => DESCRIPTION
		  -arg2 => DESCRIPTION
	Throws: MyX::Generic::Undef::Param
	Comments: NA
	See Also: NA
	
=head2 get_arg1

	Title: get_arg1
	Usage: $obj->get_arg1()
	Function: Returns arg1
	Returns: str
	Args: NA
	Throws: NA
	Comments: NA
	See Also: NA
	
=head2 set_arg1

	Title: set_arg1
	Usage: $obj->set_arg1($arg1)
	Function: sets the arg1 value
	Returns: 1 on success
	Args: -arg1 => DESCRIPTION
	Throws: MyX::Generic::Undef::Param
	Comments: NA
	See Also: NA


=head1 BUGS AND LIMITATIONS

=for author to fill in:
    A list of known problems with the module, together with some
    indication Whether they are likely to be fixed in an upcoming
    release. Also a list of restrictions on the features the module
    does provide: data types that cannot be handled, performance issues
    and the circumstances in which they may arise, practical
    limitations on the size of data sets, special cases that are not
    (yet) handled, etc.

No bugs have been reported.

Please report any bugs or feature requests to
C<bug-<RT NAME>@rt.cpan.org>, or through the web interface at
L<http://rt.cpan.org>.

=head1 TO DO

= for author to fill in:
	Include a list of features and/or tasks that have yet to be
	implemented in this object.

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

