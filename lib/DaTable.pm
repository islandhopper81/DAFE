#! usr/bin/evn perl

package DaTable;

use strict;
use warnings;

#use Param_handler;
use Getopt::Long;
use Pod::Usage;
use Carp;
use Readonly;
use Path::Class;
use Data::Dumper;
use XML::Simple qw(:strict);
use Class::Std::Utils;
#use Scalar::Util::Numeric qw(isneg isint isfloat);
use lib "/proj/cdjones_lab/ncolaian/MyX-Generic-v0.0.2/lib/MyX";
#use MyX::Generic;

{
    # Attributes
    my %da_matrix; #genome arefs with arefs
    my %id_order_aref;
    my %grp_order_aref;
    my %id_order_count_href;
    my %grp_order_count_href;
    
    ### Subroutines ###
    sub set_da_table_size;
    sub set_grp_order_aref_and_count_hash;
    sub set_id_order_aref_and_count_hash;
    sub set_genome;
    sub set_grp_da_value;
    
    sub print_full_da_table;
    sub print_filtered_da_table;
    
    sub get_grp_order_aref;
    sub get_id_order_aref;
    sub get_genome;
    sub get_group;
    
    sub _check_arg_href;
    sub _set_attributes;
    
    #################
    ## Constructor ##
    #################
    sub new {
        my ($class, $arg_href) = @_;
        
        #Bless a scalar to represent the new object
        my $new_obj = bless \do{my $anon_scalar}, $class;
        
        #check to make sure both a grp_order_aref and id_order_aref is passed
        if (_check_arg_href($arg_href) == 1) {
            $new_obj->_set_attributes($arg_href);
        }
        else {
            croak "Must pass in a hash ref containing grp_order_aref and id_order_aref in order to create a DaTable object \n";
        }
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
        elsif ( !$arg_href->{grp_order_aref} ) {
            MyX::Generic::Undef::Param->throw(
                error => "id_order_aref was not passed to the constructor in the hash reference",
                usage => "id_order_aref needs to be in the href passed to the constructor",
            )
        }
        return 1;
    }
    
    sub _set_attributes {
        my ($self, $arg_href) = @_;
        $self->set_grp_order_aref_and_count_hash($arg_href->{grp_order_aref});
        $self->set_id_order_aref_and_count_hash($arg_href->{id_order_aref});
        $self->set_da_table_size(scalar($self->get_grp_order_aref()),
                                 scalar($self->get_id_order_aref()));
        return 1;
    }
    
    ##############
    ## PRINTERS ##
    ##############
    
    sub print_full_da_table {
        my ($self, $outfile) = @_;
        open my $FDA_FH, "<", "$outfile";
        
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
    
    sub print_filtered_da_table {
        my ($self, $outfile, $filter_obj) = @_;
        open my $FILTER_FH, "<", "$outfile";
        
        my $ordered_ids_aref = $self->get_id_order_aref();
        my $ordered_grps_aref = $self->get_grp_order_aref();
        
        print $FILTER_FH "grp_id\t", join("\t", @$ordered_ids_aref), "\n" ;
        #Goes through and filters each cluster and decides if it should be printed
        foreach my $grp_id ( @$ordered_grps_aref ) {
            if ( $filter_obj->filter( $self->get_group($grp_id) ) == 1 ) {
                print $FILTER_FH "$grp_id\t",
                                 join("\t",@{ $self->get_group($grp_id) });   
            }
            else {
                next;
            }
        }
        close $FILTER_FH;
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
        if ( $genome_counts_href->{$genome_id} ) {
            return $da_table_aref_aref->[$genome_counts_href->{$genome_id}];
        }
        else {
            MyX::General::BadValue->throw(
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
        if ($grp_counts_href->{$group_id}) {
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
    
    
    
    
    #############
    ## SETTERS ##
    #############
    
    sub set_grp_order_aref_and_count_hash {
        my ($self, $grp_aref);
        $grp_order_aref{ident $self} = $grp_aref;
        my $count_href = {};
        for (my $i = 0; $i < scalar @$grp_aref; $i++) {
            $count_href->{$grp_aref->[$i]} = $i;
        }
        $grp_order_count_href{ident $self} = $count_href;
        return 1;
    }
    
    sub set_id_order_aref_and_count_hash {
        my ($self, $id_aref) = @_;
        $id_order_aref{ident $self} = $id_aref;
        my $count_href = {};
        for (my $i = 0; $i < scalar @$id_aref; $i++) {
            $count_href->{$id_aref->[$i]} = $i;
        }
        $id_order_count_href{ident $self} = $count_href;
        return 1;
    }
    
    sub set_da_table_size {
        my ($self, $grp_size, $id_size) = @_;
        my @id_array;
        my $grp_da_aref = [];
        $id_array[$id_size-1] = 0; #create an array the full size
        for (my $i = 0; $i < $id_size; $i++) {
            $id_array[$i] = $grp_da_aref;
        }
        $da_matrix{ident $self} = \@id_array;
        return 1;
    }
    
    sub set_genome {
        my ($self, $genome_id, $da_with_grp_id_aref_aref) = @_;
        my $genome_counts_href = $self->get_id_order_count_href();
        #put a error if genome doesn't match one in the ordered genome ids
        
        my $genome_pos = $genome_counts_href->{$genome_id};
        
        foreach my $grp_da_combo_aref ( @$da_with_grp_id_aref_aref ) {
            my $grp_id = $grp_da_combo_aref->[0];
            my $da_value = $grp_da_combo_aref->[1];
            
            $self->set_grp_da_value($genome_pos, $grp_id, $da_value);
        }
        return 1;
    }
    
    sub set_grp_da_value {
        my ($self, $genome_pos, $grp_id, $da_value) = @_;
        my $grp_counts_href = $self->get_grp_order_count_href();
        #put error for if the grp doesn't match one in the ordered grps
        
        my $grp_pos = $grp_counts_href->{$grp_id};
        
        my $da_table_aref_aref = $self->get_full_da_table();
        #set the individual grp value equal to the da value
        ${ $da_table_aref_aref->{$genome_pos} }->{$grp_pos} = $da_value;
        #set the da matrix to the updated matrix
        $da_matrix{ident $self} = $da_table_aref_aref;
        
        return 1;
    }
    
}