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
		
		return 1;
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

DAFE::Count::OutputRestruct - Restructure the output to match DAFEv1 format


=head1 VERSION

This document describes DAFE::Count::OutputRestruct version 0.0.1


=head1 SYNOPSIS

    use DAFE::Count::OutputRestruct;

	# make a new object
	my %args = ("out_dir" => $out_dir,
				"new_out_dir" => $new_out_dir,
				"ref_aref" => $ref_aref,
				"sample_ref" => $sample_aref,
				"perc_ids_aref" => $perc_ids_aref,
				"htseq_file_prefix" => $htseq_file_prefix);
	my $obj = DAFE::Count::OutputRestruct->new(\%args);
	
	# do the restructuring
	$obj->restructure();
	
	# get one of parameters
	$obj->get_param("out_dir");
	
	# set one of the parameters to something else
	$obj->set_param("out_dir", $out_dir2);
  
  
=head1 DESCRIPTION

In DAFEv1 I did the mapping to individual genomes.  In DAFEv2 I decided to
combine all the genomes into a single mapping database.  This resolves the issue
of a read mapping to multiple genomes.  In order to use the downstream scripts
(ie edgeR_driver.pl) I have to restructure the output from DAFEv2 to look like
the output from DAFEv1.  This object does that operation.


=head1 CONFIGURATION AND ENVIRONMENT
  
DAFE::Count::OutputRestruct requires no configuration files or environment
variables.


=head1 DEPENDENCIES

Carp
Readonly
Class::Std::Utils
List::MoreUtils qw(any)
Log::Log4perl qw(:easy)
Log::Log4perl::CommandLine qw(:all)
MyX::Generic
version; our $VERSION = qv('0.0.1')
UtilSY qw(:all)
Scalar::Util qw(openhandle)
DataObj


=head1 INCOMPATIBILITIES

None reported.


=head1 METHODS

=over
	
	new
	get_param
	set_param
	restructure
	_mkdir_structure
	_reformat
	_parse_file
	_get_htseq_file

=back

=head1 METHODS DESCRIPTION

=head2 new

	Title: new
	Usage: DAFE::Count::OutputRestruct->new({
				out_dir => $out_dir,
				new_out_dir => $new_out_dir,
				ref_aref => $ref_aref,
				sample_ref => $sample_aref,
				perc_ids_aref => $perc_ids_aref,
				htseq_file_prefix => $htseq_file_prefix
			});
	Function: Build new DAFE::Count::OutputRestruct
	Returns: DAFE::Count::OutputRestruct
	Args: -out_dir => output dir with DAFEv2 output
		  -new_out_dir => output dir to save DAFEv1 output
		  -ref_aref => array ref with reference (ie genome) IDs
		  -sample_aref => array ref with sample IDs
		  -perc_ids_aref => array ref with percent identities
		  -htseq_file_prefix => htseq file prefix
	Throws: MyX::Generic::Undef::Param
	Comments: NA
	See Also: NA
	
=head2 get_param

	Title: get_param
	Usage: $obj->get_param("p")
	Function: Returns value of the named parameter
	Returns: the parameter value (may be any type)
	Args: -p => the parameter name (ie "out_dir")
	Throws: NA
	Comments: I use a DataObj to store the parameters
	See Also: DataObj
	
=head2 set_param

	Title: set_param
	Usage: $obj->set_param("name", $value)
	Function: sets the value of the named parameter
	Returns: 1 on success
	Args: -name => the string name of the value (ie "out_dir")
	      -value => the value of the parameter to set
	Throws: MyX::Generic::Undef::Param
	Comments: NA
	See Also: DataObj

=head2 restructure

	Title: restructure
	Usage: $obj->restructure()
	Function: Restructure the DAFEv2 output to look like DAFEv1 output
	Returns: 1 on success
	Args: NA
	Throws: NA
	Comments: The DAFEv2 output is stored in a directory structure where each
	          sample has a directory and inside there are the bam and
			  htseq-count files.  After this function is executed that output
			  will have the same format as DAFEv1.  In DAFEv1 under the root
			  output directory each genome has a directory.  Inside each genome
			  directory is a directory for each sample.  Inside each sample
			  directory are the htseq-count files representing the number of
			  reads for the given genome-sample pair that are mapped.  Note that
			  the bam files are not included in this restucturing (only the
			  htseq-count files).  The user may elect to remove the bam files
			  anyway as they take a lot of memory.
	See Also: NA
	
=head2 _mkdir_structure

	Title: _mkdir_structure
	Usage: $obj->_mkdir_structure()
	Function: Make the directories in the DAFEv1 format
	Returns: 1 on success
	Args: NA
	Throws: NA
	Comments: This runs the system mkdir command
	See Also: NA
	
=head2 _reformat

	Title: _reformat
	Usage: $obj->_reformat()
	Function: reformats the htseq-count files
	Returns: 1 on success
	Args: NA
	Throws: NA
	Comments: NA
	See Also: NA
	
=head2 _parse_file

	Title: _parse_file
	Usage: _parse_file(...)
	Function: parses the DAFEv2 htseq-count files and prints the corresponding
	          DAFEv1 htseq-count files
	Returns: 1 on success
	Args: NA
	Throws: NA
	Comments: NA
	See Also: NA
	
=head2 _get_htseq_file

	Title: _get_htseq_file
	Usage: _get_htseq_file()
	Function: compiles the htseq-count file name and path
	Returns: 1 on success
	Args: NA
	Throws: NA
	Comments: NA
	See Also: NA


=head1 BUGS AND LIMITATIONS

No bugs have been reported.

Please report any bugs or feature requests to
C<bug-<RT NAME>@rt.cpan.org>, or through the web interface at
L<http://rt.cpan.org>.

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

