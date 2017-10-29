package DAFE::AbEstimator;

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
use Table;

# set up the logging environment
my $logger = get_logger();

{
	# Usage statement
	Readonly my $NEW_USAGE => q{ new() };

	# Attributes #
	
	# Getters #

	# Setters #

	# Others #
	sub calc_marker_abund;
	sub calc_abund_est;
	sub calc_abund_est_v2;
	sub _mean;
	sub _median;


	###############
	# Constructor #
	###############
	sub new {
		my ($class) = @_;

		# Croak if calling new on already blessed reference
		croak 'Constructor called on existing object instead of class'
			if ref $class;

		# Make sure the required parameters are defined
		if ( any {!defined $_}
			) {
			MyX::Generic::Undef::Param->throw(
				error => 'Undefined parameter value',
				usage => $NEW_USAGE,
			);
		}

		# Bless a scalar to instantiate an object
		my $new_obj = bless \do{my $anon_scalar}, $class;

		# Set Attributes

		return $new_obj;
	}

	###########
	# Getters #
	###########

	###########
	# Setters #
	###########

	##########
	# Others #
	##########
	sub calc_marker_abund {
		my ($self, $tbl, $method) = @_;

		if ( ! is_defined($method) ) {
			$method = "mean";
		}

		my @gene_avgs = ();
		foreach my $col ( @{$tbl->get_col_names()} ) {
			if ( $method eq "median" ) {
				push @gene_avgs, $self->_median($tbl->get_col($col));
			}
			else {
				push @gene_avgs, $self->_mean($tbl->get_col($col));
			}
		}

		return(\@gene_avgs);
	}

	sub calc_abund_est {
		my ($self, $tbl, $method) = @_;
		
		if ( ! is_defined($method) ) {
			$method = "mean";
		}

		my @gene_avgs = (); 
		my @col_vals = (); 

		foreach my $col ( @{$tbl->get_col_names()} ) { 
			# get the median instead of the average (ie mean)
			# if the flag is set
			if ( $method eq "median" ) { 
				push @gene_avgs, $self->_median($tbl->get_col($col));
			}   
			else {
				push @gene_avgs,  $self->_mean($tbl->get_col($col));
			}   
		}   
		
		my $est;
		if ( $method eq "median" ) {
			$est = $self->_median(\@gene_avgs);
		}
		elsif ( $method eq "mean" ) {
			$est = $self->_mean(\@gene_avgs);
		}

		return($est);
	}
	
	sub calc_abund_est_v2 {
		my ($self, $tbl, $method) = @_;
		
		if ( ! is_defined $method ) {
			$method = "mean";
		}
		
		my $est;
		my @vals = ();
		foreach my $r ( @{$tbl->get_row_names} ) {
			push @vals, @{$tbl->get_row($r)};
		}
		
		if ( $method eq "median" ) {
			$est = $self->_median(\@vals);
		}
		elsif ( $method eq "mean" ) {
			$est = $self->_mean(\@vals);
		}
		else {
			# throw an error
		}
		
		return($est);
	}

	sub _mean {
		my ($self, $aref) = @_; 

		my $len = scalar @{$aref};

		my $sum = 0;
	
		foreach my $val ( @{$aref} ) { 
			$sum = $sum + $val;
		}   

		return($sum/$len);
	}	

	sub _median {
		my ($self, $aref) = @_;

		my @vals = sort {$a <=> $b} @{$aref};
		my $len = @vals;
		if($len%2) #odd?
		{
			return $vals[int($len/2)];
		}
		else #even
		{
			return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
		}
	}
}

1; # Magic true value required at end of module
__END__

=head1 NAME

DAFE::AbEstmator - Estimates the abundance of a genome from metagenome data


=head1 VERSION

This document describes DAFE::AbEstmator version 0.0.1


=head1 SYNOPSIS

    use DAFE::AbEstmator;

	# make a new object
	my $obj = DAFE::AbEstmator->new();
	
	#
  
  
=head1 DESCRIPTION

In DAFEv1 I did the mapping to individual genomes.  In DAFEv2 I decided to
combine all the genomes into a single mapping database.  This resolves the issue
of a read mapping to multiple genomes.  In order to use the downstream scripts
(ie edgeR_driver.pl) I have to restructure the output from DAFEv2 to look like
the output from DAFEv1.  This object does that operation.


=head1 CONFIGURATION AND ENVIRONMENT
  
DAFE::AbEstmator requires no configuration files or environment variables.


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
Table


=head1 INCOMPATIBILITIES

None reported.


=head1 METHODS

=over
	
	new
	calc_marker_abund
	calc_abund_est
	calc_abund_est_v2
	_mean
	_median

=back

=head1 METHODS DESCRIPTION

=head2 new

	Title: new
	Usage: DAFE::AbEstmator->new();
	Function: Build new DAFE::AbEstmator
	Returns: DAFE::AbEstmator
	Args: NA
	Throws: MyX::Generic::Undef::Param
	Comments: NA
	See Also: NA

=head2 calc_marker_abund

	Title: calc_marker_abund
	Usage: $obj->calc_marker_abund($tbl, $method)
	Function: Calculates the mean/median abudance for each sample
	Returns: aref
	Args: -tbl => Table object of counts (rows are markers, cols are samples)
	Throws: MyX::Generic::Undef::Param
	Comments: The table should have samples represented as columns and marker
              genes as rows.  For example at cell, (1,2) there would be a number
              representing the number of reads that map to marker gene 1 in
              sample 2.

			  Insetead of reducing the full count matrix for a table (ie genome)
			  to a single value it returns a value for each sample (column).
	
=head2 calc_abund_est

	Title: calc_abund_est
	Usage: $obj->calc_abund_est($tbl, $method)
	Function: Calculates abundance estimate of the counts in the table
	Returns: int
	Args: -tbl => Table object of counts
	      -method => "mean" | "median"
	Throws: NA
	Comments: The table should have samples represented as columns and marker
	          genes as rows.  For example at cell, (1,2) there would be a number
			  representing the number of reads that map to marker gene 1 in
			  sample 2.
			  
			  This function uses the $method first on each of the columns (ie
			  samples).  Then the resulting numbers are once again submitted to
			  $method to generate the final estimate.  So when $method is "mean"
			  the estimate can be interpreted as the average number of reads
			  across all marker genes then averaged again across all samples.
	See Also: NA
	
=head2 calc_abund_est_v2

	Title: calc_abund_est_v2
	Usage: $obj->calc_abund_est_v2($tbl, $method)
	Function: Calculates abundance estimate of the counts in the table
	Returns: int
	Args: -tbl => Table object of counts
	      -method => "mean" | "median"
	Throws: NA
	Comments: The table should have samples represented as columns and marker
	          genes as rows.  For example at cell, (1,2) there would be a number
			  representing the number of reads that map to marker gene 1 in
			  sample 2.
			  
			  This function disoves the table into a single array, sorts the
			  array, and then opperates on the array by the given $method.  When
			  $method is set to "mean" it perfoms exactly the same as
			  calc_abund_est (I think).  However, when method is set to "median"
			  the outcome can be slightly different.  I would recommend using
			  calc_abund_est.
	See Also: NA
	
=head2 _mean

	Title: _mean
	Usage: $obj->_mean($aref)
	Function: Calculates the mean of the given values
	Returns: int
	Args: -aref => aref of values
	Throws: NA
	Comments: Private.  Should not be called outside this object
	See Also: NA
	
=head2 _median

	Title: _median
	Usage: $obj->_median($aref)
	Function: Calculates the median of the given values
	Returns: int
	Args: -aref => aref of values
	Throws: NA
	Comments: Private.  Should not be called outside this object
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

