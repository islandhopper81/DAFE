#! usr/bin/env perl

package DecoupleDa;

use strict;
use warnings;
use Param_handler;
use Getopt::Long;
use Pod::Usage;
use Carp;
use Readonly;
use Path::Class;
use Data::Dumper;
use XML::Simple qw(:strict);
use Class::Std::Utils;
use Scalar::Util::Numeric qw(isneg isint isfloat);
use MyX::Generic;
use Table;
use UtilSY qw(:all);

# Readonly Global Variables #
Readonly::Scalar my $ENR => 1; # enriched
Readonly::Scalar my $NS => 0;  # not significant
Readonly::Scalar my $DEP => -1; # depleted
Readonly::Scalar my $LOW => -2; #low abundance
Readonly::Scalar my $UMS => -3; # unmeasurable but present in genome
Readonly::Scalar my $ABS => -4; # absent from genome


{
    #Attributes
    my %param_handler_obj;
    my %count_file_names;
    
    ### Subroutines ###
    sub decouple; # Takes in a genome and searches its annotation file to search if the -2's from the differentially abundant statistics are b/c of no info or are not present in the genome
    sub _remake_da_with_decoupled;
    sub _check_negative_two; # DEPRECIATED
    
    sub _set_param_handler;
    sub set_count_file_name;
    sub manually_set_cnt_file_name_from_edger;
    
    sub get_param_handler;
	sub get_count_file_name;
    
    
    
    
    #################
    ## CONSTRUCTOR ##
    #################
    
    sub new {
        my ($class, $param_obj) = @_;
        
        #Bless a scalar to represent the new object
        my $new_obj = bless \do{my $anon_scalar}, $class;
        
        #Handle a passed in Param_handler
        if ( $param_obj->isa( "Param_handler" ) ) {
            $new_obj->_set_param_handler($param_obj);
            $new_obj->_set_count_file_name($param_obj);
        }
        else {
            croak "Need to pass in a Param_handler object into the constructor";
        }
        
        return $new_obj;
    }
    
    #################
    ## SUBROUTINES ##
    #################
    
    sub decouple {
        my ($self, $genome_id) = @_;
        my $param_obj = $self->get_param_handler();
        #Need to get count and dafe directories file handles
        my $count_f_name = $self->get_count_file_name();
        my $count_dir = $param_obj->get_count_dir();
        my $count_file = "$count_dir/$genome_id/$count_f_name";
        #check the count file
        if (!-e $count_file) {
            MyX::Generic::BadValue->throw(
                error => "$count_f_name is not the name of the count file from edger. Manually set it with manually_set_cnt_file_name_from_edger()",
                value => $count_f_name,
            );
        }
        
        my $dafe_dir = $param_obj->get_dafe_dir();
        my $annote_name = $param_obj->get_annote_file_name();
        my $annote_file = "$dafe_dir/$genome_id/$annote_name";
        
        #create file objects that will be read in
        my $count_file_obj = file($count_file);
        my $annote_file_obj = file($annote_file);
            
        #create a slurped array 
        my @count_file_a = $count_file_obj->slurp(chomp=>1, split=>qr/\s/);
        
		#my $slurped_annote_file = $annote_file_obj->slurp();
		my $annote_file_tbl = Table->new();
		$annote_file_tbl->load_from_file($annote_file);
		my $grp_by_col = $param_obj->get_grp_genes_by();
		my $grp_vals_href; # a lookup list of all the present features (ie COGs, COG categories, etc)
		my $grp_vals_aref;
		if ( $annote_file_tbl->has_col($grp_by_col) ) {
			$grp_vals_aref = $annote_file_tbl->get_col($grp_by_col);
			$grp_vals_href = aref_to_href($grp_vals_aref);
		}
		$annote_file_tbl->reset();
        
        #go through each feature (ie COG) and correct the DA value as needed
        foreach my $line ( @count_file_a ) {
            if ( $line->[1] eq $UMS) {
				# At this point I'm looking at a feature that is either unmeasurable
				# or completely absent from the genome.  That is what I figure out
				# below.
				if ( ! defined $grp_vals_href->{$line->[0]} ) {
					# if this feature is not defined in the grp_vals_href lookup, that means
					# it was not in the annotation file therefore it is a ABS feature
					$line->[1] = $ABS;
				}
            }
        }
        _remake_da_with_decoupled(\@count_file_a, $count_file); #print decoupled info in old file
        if ( -z "$count_file" ) {
            croak "new file that was supposed to be remade with new da values is empty";
        }
        #return a hash reference of hash references
        return \@count_file_a;
    }
    
    sub _remake_da_with_decoupled {
        my ($count_file_aref_aref, $file) = @_;
        open (my $CF_FH, ">", $file); #opens file for writing, not appending
        
        foreach my $line ( @$count_file_aref_aref ) {
            print $CF_FH join("\t", @$line), "\n";
        }
        close $CF_FH;
        
        return 1;
    }
    
    sub _check_negative_two_old {
        my ($grp_id, $slurped_annote_file) = @_;
		# DEPRECIATED

        #check to see if the group id is in the genomes annotation file
        if ( $slurped_annote_file =~ qr/$grp_id/i ) {
            return 1;
        }
        else {
            return 0;
        }
    }

	sub _check_negative_two {
		my ($grp_id, $annote_tbl, $grp_by_col) = @_;
		# DEPRECIATED

		# there was a bug in this function.  It was only looking for the $grp_id
		# in the entire annotation file.  This would cause problems when the 
		# grp_id was something like "A" as in the case when I do cog groups.
		# So I changed the parameters so that a Table object is passed in which
		# has all the data in the annotation file.  I check the correct column
		# in the annotation file for any instances of grp_id.  If there are no
		# instances then I return 0.

		# note: this function takes a ton of memory.  I don't use it anymore
		
		# get the column of values in the annoation table
		if ( $annote_tbl->has_col($grp_by_col) ) {
			my $vals = $annote_tbl->get_col($grp_by_col);

			foreach my $v ( @{$vals} ) {
				if ( $v eq $grp_id ) {
					return 1;
				}
			}
		}
		else {
			# throw some error
			MyX::Generic::BadValue->throw(
				error => "Cannot find column $grp_by_col in the annotation table"
			);
		}
		

		return 0;
	}
    
    
    ### SETTERS ###
    
    sub _set_param_handler {
        my ($self, $param_obj ) = @_;
        $param_obj->check_print_params();
        $param_obj->_check_count_file_name();
        $param_handler_obj{ident $self} = $param_obj;
        return 1;
    }
    
    sub _set_count_file_name {
        my ($self, $param_obj) = @_;

		# get the test values.  These are used so that mutliple
		# statistical tests can be ran on the same aggregated
		# count file without overwritting previous statistical
		# test outputs
		my $test_str = $param_obj->get_test();
		my $test1;
		my $test2;
		if ( $test_str =~ m/\["(.*)",\s*"(.*)"\]/ ) {
			$test1 = $1;
			$test2 = $2;
		}
		else {
			croak("Cannot find the tests in the test string: $test_str");
		}
        
        my $count_f_name = $param_obj->get_count_file_name();
        my $grp_genes_by = $param_obj->get_grp_genes_by();
        $count_f_name =~ s/\.txt//;
        $count_f_name .= "_$grp_genes_by\_agg_";
		$count_f_name .= "$test1\_v_$test2";
		$count_f_name .= "_da_vec.txt";
        $count_file_names{ident $self} = $count_f_name;
        
        return 1;
    }
    #add to documentation
    sub manually_set_cnt_file_name_from_edger {
        my ($self, $name_string ) = @_;
        $count_file_names{ident $self} = $name_string;
        return 1;
    }
    
    ### GETTERS ###
    
    sub get_param_handler {
        my ($self) = @_;
        return $param_handler_obj{ident $self};
    }
    
    sub get_count_file_name {
        my ($self) = @_;
        return $count_file_names{ident $self};
    }
}

1;

__END__

=head1  DecoupleDa Object

DecoupleDa -    takes a param handler object into it's constructor. It's main
                method, decouple, will search the annotation file of the genome of interest and check to see if the group is actually in the genome.
                
=head1  Version

This documentation refers to DecoupleDa.pm

=head1 Included Modules

use Param_handler;
use Getopt::Long;
use Pod::Usage;
use Carp;
use Readonly;
use Path::Class;
use Data::Dumper;
use XML::Simple qw(:strict);
use Class::Std::Utils;
use Scalar::Util::Numeric qw(isneg isint isfloat);
use MyX::Generic;

=head1 Inherit

    NA
    
=head1 Synopsis

    $decouple_obj = DecoupleDa->new( $param_obj );
    $grp_and_da_values_aref = $decouple_obj->decouple($genome_id);
    
=head1 Description

    This object is used to decouple the -2's that are present in the count files produced from edger analysis. The decouple() method will remake the count file with the updated DA values as well as output an array ref containing array references containing the group id in the first position and the DA value in the second position.

=head1 Methods

sub new;
sub decouple;
sub _remake_da_with_decoupled;
sub _check_negative_two;
sub _set_param_handler;
sub set_count_file_name;
sub get_param_handler;
sub get_count_file_name;

=head1 Methods Description

=head2  new()

    Title:		new
	Usage:		$decouple_obj = DecoupleDa->new( $param_handler_obj );
	Function:	Creates a new instance of a DecoupleDa object
	Returns:	DecoupleDa
	Args:		$param_handler_obj =>   Param_handler object that has the print
                                        parameters in it and checked.
	Throws:		croak
	Comments:	Need to make sure that the print params are checked in the
                Param_handler object before passing the object for the creation
                of a DecoupleDa object.
	See Also:   Param_handler

=head2  decouple()

    Title:      decouple
    Usage:      $decouple_obj->decouple( $genome_id );
    Function:   Checks each group in the genome. If there is a -2 as the Da value,
                it determines if it is there due to not enough data or because the group does not exist in the genome. After decoupling it will then update the count file with the updated decoupled information. It also returns an array reference containing array references that hold goup ids and their corresponding DA values.
    Returns:    Array Reference containing Array References
    Args:       $genome_id  =>  This is the genome id that will be used to find
                                the count file countaining all the DA values.
    Throws:     TBD
    Comments:   NA
    See Also:   NA
    
=head2  get_param_handler()

    Title:      get_param_handler
    Usage:      $decouple_obj->get_param_handler();
    Function:   Returns the Param_handler object passed into the DecoupleDa object
    Returns:    Param_handler
    Args:       $decouple_obj   =>  A blessed DecoupleDa object
    Throws:     NA
    Comments:   NA
    See Also:   Param_handler
    
=head2 get_count_file_name()
    
    Title:      get_count_file_name
    Usage:      $decouple_obj->get_count_file_name();
    Function:   Returns the name of the file found in the genomes count directory
    Returns:    String
    Args:       $decouple_obj   =>  A blessed DecoupleDa object
    Throws:     NA
    Comments:   The name will have the same extension from the data returned from
                the edger analysis.
    See Also:   NA
    
=head2 manually_set_cnt_file_name_from_edger()

    Title:      manually_set_cnt_file_name_from_edger
    Usage:      $decouple_obj->manually_set_cnt_file_name_from_edger( $name );
    Function:   Sets the file to look for in the count directory. If it is different than the standard output from edgeR
    Returns:    1
    Args:       $decouple_obj   =>  A DecoupleDa object
                $name           =>  String representing the file name from edgeR
    Throws:     NA
    Comments:   NA
    See Also:   NA


=head1 Configuration And Environment

    DecoupleDa requires no configuration files or environment variables.


=head1 Dependencies

    Requires the Param_handler module, because you must pass in a Param_handler.

=head1 Incompatibilities

    None reported.

=head1 Bugs And Limitations

    No bugs have been reported.
    
=head1 Author

    Nicholas Colaianni
    contact via C<< <ncolaian@live.unc.edu> >>
    
=head1 Licence And Copyright

Copyright (c) 2016, Nicholas Colaianni
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

=head1 Disclaimer of Warranty

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
