#! usr/bin/env perl

package DaTable;

use strict;
use warnings;

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

{
    # Attributes
    my %da_matrix; #genome arefs with arefs
    my %id_order_aref;
    my %grp_order_aref;
    my %id_order_count_href;
    my %grp_order_count_href;
    
    ### Subroutines ###
    sub _set_da_table_size;
    sub _set_grp_order_aref_and_count_hash;
    sub _set_id_order_aref_and_count_hash;
    sub set_genome;
    sub _set_grp_da_value;
    sub _check_da_file;
    sub _set_datable_from_file;
    
    sub print_full_da_table;
    sub check_for_unset_values;
    
    sub get_full_da_table;
    sub get_grp_order_aref;
    sub get_id_order_aref;
    sub _set_attributes;
    
    #################
    ## Constructor ##
    #################
    sub new {
        my ($class, $arg_href) = @_;
        
        #Bless a scalar to represent the new object
        my $new_obj = bless \do{my $anon_scalar}, $class;
        
        #check to make sure both a grp_order_aref and id_order_aref is passed
        if ( exists $arg_href->{da_file} ) {
            _check_da_file( $arg_href->{da_file} );
            $new_obj->_set_datable_from_file( $arg_href->{da_file} );
        }
        elsif (_check_arg_href($arg_href) == 1) {
            $new_obj->_set_attributes($arg_href);
        }
        else {
            croak "Must pass in a hash ref containing grp_order_aref and id_order_aref in order to create a DaTable object \n";
        }
        
        return $new_obj;
    }
    
    
    #################
    ## SUBROUTINES ##
    #################
    
    sub _check_arg_href {
        my ($arg_href) = @_;
        if ( !$arg_href->{grp_order_aref} ) {
            MyX::Generic::Undef::Param->throw(
                error => "grp_order_aref was not passed to the constructor in the hash reference",
                usage => "grp_order_aref needs to be in the href passed to the constructor",
            )
        }
        elsif ( !$arg_href->{id_order_aref} ) {
            MyX::Generic::Undef::Param->throw(
                error => "id_order_aref was not passed to the constructor in the hash reference",
                usage => "id_order_aref needs to be in the href passed to the constructor",
            )
        }
        return 1;
    }
    
    sub _set_attributes {
        my ($self, $arg_href) = @_;
        $self->_set_grp_order_aref_and_count_hash($arg_href->{grp_order_aref});
        $self->_set_id_order_aref_and_count_hash($arg_href->{id_order_aref});
        $self->_set_da_table_size(scalar( @{$self->get_grp_order_aref()} ),
                                 scalar( @{$self->get_id_order_aref()} ) );
        return 1;
    }
    
    sub _check_da_file {
        my ($file_path) = @_;
        my $file_obj = file($file_path);
        my @file_array = $file_obj->slurp(chomp=>1, split=>qw/\t/);
        
        if ( @{ $file_array[0] }[0] ne "grp_id" ) {
            MyX::Generic::BadValue->throw(
                error => "$file_path is not a DaTable file",
                value => "DaTable file",
            );
        }
        return 1;
    }
    
    sub _set_datable_from_file {
        my ($self, $file_path) = @_;
        my $file_obj = file($file_path);
        my @file_array = $file_obj->slurp(chomp=>1, split=>qw/\t/);
        #get ordered ids
        my $ordered_ids = [];
        foreach my $id ( @{ $file_array[0]} ) {
            if ( $id eq "grp_id" ) {
                next;
            }
            push @$ordered_ids, $id;
        }
        $self->_set_id_order_aref_and_count_hash($ordered_ids);
        
        #get ordered grps
        my $ordered_grps = [];
        for(my $grp_line = 1; $grp_line < scalar(@file_array); $grp_line++) {
            push @$ordered_grps, @{ $file_array[$grp_line] }[0];
        }
        $self->_set_grp_order_aref_and_count_hash($ordered_grps);
        
        #set da_table
        $self->_set_da_table_size(scalar(@$ordered_grps), scalar(@$ordered_ids));
        
        #get and set each genome
        for(my $i = 1; $i <= scalar(@$ordered_ids); $i++) {
            my $genome = $ordered_ids->[$i-1];
            my $genome_da_info = [];
            for(my $j = 1; $j <= scalar(@$ordered_grps); $j++) {
                my $grp = $ordered_grps->[$j-1];
                my $grp_da_info = $file_array[$j];
                my $grp_with_da = [ $grp, $grp_da_info->[$i] ];
                push @$genome_da_info, $grp_with_da;
            }
            $self->set_genome($genome, $genome_da_info);
        }
        return 1;
    }
    
    ##############
    ## PRINTERS ##
    ##############
    
    sub print_full_da_table {
        my ($self, $outfile) = @_;
        open (my $FDA_FH, ">", $outfile);
        
        my $ordered_ids_aref = $self->get_id_order_aref();
        my $ordered_grps_aref = $self->get_grp_order_aref();
        
        print $FDA_FH "grp_id\t", join("\t", @$ordered_ids_aref), "\n" ;
        
        for (my $i = 0; $i < scalar( @$ordered_grps_aref ); $i++ ) {
            print $FDA_FH $ordered_grps_aref->[$i], "\t";
            my $grp_da_values_aref = $self->get_group($ordered_grps_aref->[$i]);
            print $FDA_FH join("\t", @$grp_da_values_aref), "\n";
        }
        close $FDA_FH;
        return 1;
    }
    #need to add to codumentation and test
    sub check_for_unset_values {
        my ($self) = @_;
        my $da_table = $self->get_full_da_table();
        my $grp_order_aref = $self->get_grp_order_aref();
        my $id_order_aref = $self->get_id_order_aref();
        my @genome_and_grp_pairs_with_missing_values;
        for ( my $i = 0; $i < scalar @$id_order_aref; $i++) {
            for (my $j = 0; $j < scalar @$grp_order_aref; $j++) {
                if ( @{$da_table->[$i]}[$j] ne "1" &&
                     @{$da_table->[$i]}[$j] ne "-3" &&
                     @{$da_table->[$i]}[$j] ne "-2" &&
                     @{$da_table->[$i]}[$j] ne "-1" &&
                     @{$da_table->[$i]}[$j] ne "0" ) {
                    my $missing = [ $id_order_aref->[$i], $grp_order_aref->[$j] ];
                    push @genome_and_grp_pairs_with_missing_values, join("=>", @$missing);
                }
            }
        }
        if ( scalar @genome_and_grp_pairs_with_missing_values > 0 ) {
            MyX::Generic::Undef::Attribute->throw(
                error => join(", ", @genome_and_grp_pairs_with_missing_values),
                att_name => "Missing Values",
                );
        }
        
        return 1;
    }
    
    
    #############
    ## GETTERS ##
    #############
    
    sub get_grp_order_aref {
        my ($self) = @_;
        return $grp_order_aref{ident $self};
    }
    
    sub get_id_order_aref {
        my ($self) = @_;
        return $id_order_aref{ident $self};
    }
    
    sub get_genome {
        my ($self, $genome_id) = @_;
        my $da_table_aref_aref = $self->get_full_da_table();
        my $genome_counts_href = $self->get_id_order_count_href();
        if ( exists $genome_counts_href->{$genome_id} ) {
            return $da_table_aref_aref->[$genome_counts_href->{$genome_id}];
        }
        else {
            MyX::Generic::BadValue->throw(
                error => "$genome_id is not included in ordered genome ids array",
                value => "$genome_id",
            );
        }
    }
    
    sub get_group {
        my ($self, $group_id) = @_;
        my $da_table_aref_aref = $self->get_full_da_table();
        my $grp_counts_href = $self->get_grp_order_count_href();
        my @cluster_da_values;
        if (exists $grp_counts_href->{$group_id}) {
            foreach my $gene_aref (@$da_table_aref_aref) {
                push(@cluster_da_values,
                     $gene_aref->[$grp_counts_href->{$group_id}]);
            }
        }
        else {
            MyX::Generic::BadValue->throw(
                error => "$group_id was not found in the included groups",
                value => $group_id,
            );
        }
        
        return \@cluster_da_values;
    }
    
    sub get_full_da_table {
        my ($self) = @_;
        return $da_matrix{ident $self};
    }
    
    sub get_id_order_count_href {
        my ($self) = @_;
        return $id_order_count_href{ident $self};
    }
    
    sub get_grp_order_count_href {
        my ($self) = @_;
        return $grp_order_count_href{ident $self};
    }
    
    
    
    #############
    ## SETTERS ##
    #############
    
    sub _set_grp_order_aref_and_count_hash {
        my ($self, $grp_aref) = @_;
        $grp_order_aref{ident $self} = $grp_aref;
        my $count_href = {};
        for (my $i = 0; $i < scalar @$grp_aref; $i++) {
            $count_href->{$grp_aref->[$i]} = $i;
        }
        $grp_order_count_href{ident $self} = $count_href;
        return 1;
    }
    
    sub _set_id_order_aref_and_count_hash {
        my ($self, $id_aref) = @_;
        $id_order_aref{ident $self} = $id_aref;
        my $count_href = {};
        for (my $i = 0; $i < scalar @$id_aref; $i++) {
            $count_href->{$id_aref->[$i]} = $i;
        }
        $id_order_count_href{ident $self} = $count_href;
        return 1;
    }
    
    sub _set_da_table_size {
        my ($self, $grp_size, $id_size) = @_;
        my @id_array;
        #$id_array[$id_size-1] = 0; #create an array the full size
        for (my $i = 0; $i < $id_size; $i++) {
            push @id_array, [];
            for (my $j = 0; $j < $grp_size; $j++) {
                push @{ $id_array[$i] }, 4;   
            }
        }
        $da_matrix{ident $self} = \@id_array;
        return 1;
    }
    
    sub set_genome {
        my ($self, $genome_id, $da_with_grp_id_aref_aref) = @_;
        my $genome_counts_href = $self->get_id_order_count_href();
        #put an error if genome doesn't match one in the ordered genome ids
        if ( !exists $genome_counts_href->{$genome_id} ) {
            carp "$genome_id is not in the ordered genome id passed";
            return;
        }
        
        my $genome_pos = $genome_counts_href->{$genome_id};

        foreach my $grp_da_combo_aref ( @$da_with_grp_id_aref_aref ) {
            my $grp_id = $grp_da_combo_aref->[0];
            my $da_value = $grp_da_combo_aref->[1];
            
            $self->_set_grp_da_value($genome_pos, $grp_id, $da_value);
        }
        return 1;
    }
    
    sub _set_grp_da_value {
        my ($self, $genome_pos, $grp_id, $da_value) = @_;
        my $grp_counts_href = $self->get_grp_order_count_href();
        #put error for if the grp doesn't match one in the ordered grps
        if ( !exists $grp_counts_href->{$grp_id} ) {
            carp "$grp_id is in the edger file but not in the grp metadata file";
            return;
        }
        
        my $grp_pos = $grp_counts_href->{$grp_id};
        
        my $da_table_aref_aref = $self->get_full_da_table();
        #set the individual grp value equal to the da value
        @{ $da_table_aref_aref->[$genome_pos] }[$grp_pos] = $da_value;
        #set the da matrix to the updated matrix
        $da_matrix{ident $self} = $da_table_aref_aref;
        
        return 1;
    }
    
}

1;

__END__

=head1  DaTable

    DaTable -   This object is a matrix containing the da values for each group in
                each genome. Da values can be added individually or by genome.
                
=head1  Version 0.0.1

This documentation refers to DaTable.pm

=head1 Included Modules

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

my $da_table_obj = DaTable->new( { grp_order_aref => $grp_aref,
                                   id_order_aref  => $id_ref });
$da_table_obj->set_genome( $genome_id, $grp_da_value_aref_aref );
$da_table_obj->get_genome( $genome_id );
$da_table_obj->get_group( $grp_id );
$da_table_obj->print_full_da_table( $outfile );
$da_table_obj->print_filtered_da_table( $outfile, $filter_obj );


=head1 Description

This object holds the differential abundance values of all the groups in each individual genome. It takes the decoupled da data and combines it to a single table. 

=head1 Methods

sub new;
sub _check_arg_href;
sub _set_attributes;
sub print_full_da_table;
sub _get_grp_order_aref;
sub _get_id_order_aref;
sub get_genome;
sub get_group;
sub get_full_da_table;
sub _set_grp_order_aref_and_count_hash;
sub _set_id_order_aref_and_count_hash;
sub _set_da_table_size;
sub set_genome;
sub _set_grp_da_value;

=head1 Methods Description

=head2  new()
    
    Title:      new
    Usage:      $da_table_obj = DaTable->new( { grp_order_aref => $grp_order_aref
                                                id_order_aref => $id_order_aref});
    Function:   This creates a new da table object that holds a matrix      that are the dimensions of the scalars of the id_array x grp_array
    Returns:    DaTable
    Args:       HashRef containing grp_order_aref and id_order_aref
                grp_order_aref  =>  This is an array ref containing all the grps
                                    in the analysis in the correct order
                id_order_aref   =>  This is an array ref containing all the ids
                                    in the analysis, in the correct order
    Throws:     croak
                MyX::Generic::Undef::Param
                
    Comments:   The ordered arefs can, and probably should, come from Justify
                and Aggregate objects. 
    See Also:   Justify and Aggregate objects. Look at DecoupleDa program to check the output structure that is being read by this program.
                 
=head2  get_genome()
    
    Title:      get_genome
    Usage:      $da_table_obj->get_genome( $genome_id );
    Function:   Returns an array references containing all the da values from that
                genome, in the correct grp order
    Returns:    Array Reference
    Args:       $genome_id  =>  A genome id that is included in the analysis
    Throws:     MyX::General::BadValue;
    Comments:   NA
    See Also:   NA
    
=head2  get_group()

    Title:      get_group
    Usage:      $da_table_obj->get_group( $group_id );
    Function:   Returns an array reference containing the da values of that group
                for every single genome in the analysis. It will be in the correct genome order
    Returns:    Array Reference
    Args:       $group_id   =>  This is the id of the group
    Throws:     MyX::Generic::BadValue
    Comments:   Check the annotation files and grp metadata files for the exact
                spelling of the group. CAPITALIZATION MATTERS
    See Also:   NA
    
=head2  get_full_da_table()
    
    Title:      get_full_da_table
    Usage:      $da_table_obj->get_full_da_table();
    Function:   Returns an array reference of array references that represents
                the full genome x group matrix
    Returns:    Array Reference of array references
    Args:       $da_table_obj   =>  A blessed DaTable object. Table should be
                                    filled using set genome
    Throws:     NA
    Comments:   NA
    See Also:   NA
    
=head2  set_genome()

    Title:      set_genome
    Usage:      $da_table_obj->set_genome( $genome_id, $grp_da_aref_aref );
    Function:   Sets the da values of one of the genomes in the matrix contained
                in the DaTable object
    Returns:    1
    Args:       $genome_id  =>  This is the id of one of the genomes included in
                                the analysis
                $grp_da_aref_aref   =>  An array reference containing array references that pairs a da value with a grp_id. This is from a DecoupleDa object
    Throws:     MyX::Generic::BadValue
    Comments:   NA
    See Also:   DecoupleDa
    
=head2  check_for_unset_values()

    Title:      check_for_unset_values
    Usage:      $da_table_obj->check_for_unset_values()
    Function:   This will check the matrix in a DaTable object for any unset values
    Returns:    1
    Args:       $da_table_obj   => A DaTable object
    Throws:     MyX::Generic::Undef::Attribute => Containing the matching group and genome pairs that have missing values
    Comments:   NA
    See Also:   NA
    
=head2  get_grp_order_aref()

    Title:      get_grp_order_aref
    Usage:      $da_table_obj->get_grp_order_aref();
    Function:   Returns an array with the grp order that the matrix is created from
    Returns:    Array Reference
    Args:       $da_table_obj   =>  DaTable object
    Throws:     NA 
    Comments:   NA
    See Also:   NA
    
=head2  get_id_order_aref()

    Title:      get_id_order_aref
    Usage:      $da_table_obj->get_id_order_aref();
    Function:   Returns an array containing ids in the order used to create the matrix within the DaTable object.
    Returns:    Array Reference
    Args:       $da_table_obj   => DaTable object
    Throws:     NA
    Comments:   NA
    See Also:   NA
    
=head2  get_grp_order_count_href()

    Title:      get_grp_order_count_href
    Usage:      $da_table_obj->get_grp_order_count_href();
    Function:   Returns a hash reference that has the group as the keys and the order in the matrix as the value
    Returns:    Hash Reference
    Args:       $da_table_obj   =>  DaTable Object
    Throws:     NA
    Comments:   NA
    See Also:   NA
    
=head2  get_id_order_count_href()

    Title:      get_id_order_count_href
    Usage:      $da_table_obj->get_id_order_count_href();
    Function:   Returns a hash reference containing the genome id as the keys and its order in the matrix as the value
    Returns:    Hash Reference
    Args:       $da_table_obj   => DaTable object
    Throws:     NA
    Comments:   These are good for position look-ups
    See Also:   NA
    
=head2  print_full_da_table()

    Title:      print_full_da_table
    Usage:      $da_table_obj->print_full_da_table( $outfile );
    Function:   Prints the count matrix contained within the DaTable object
    Returns:    A txt file located in the area passed into the function
    Args:       $da_table_obj   => DaTable object
                $outfile        => The path to write the text file containing the count matrix
    Throws:     NA
    Comments:   NA
    See Also:   NA

=head1 Configuration And Environment

    DaTable requires no configuration files or environment variables.


=head1 Dependencies

    NA

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
