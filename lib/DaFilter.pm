#! usr/bin/env perl

package DaFilter;

use strict;
use warnings;

use DaTable;
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

{
    # Attributes
    my %da_counts_per_grp_aref_href; #aref of hrefs for each groups da counts
    my %da_table_objects;
    
    #Filter Constants
    Readonly::Hash my %FILTER_NUMBERS => map { $_ => 1 } qw(
        f1
        f0
        f-1
        f-2
        f-3
    );
    
    ### Subroutines ###
    sub _set_counts;
    sub _filter_by_ones;
    sub _filter_by_negative_twos;
    sub _filter_by_negative_ones;
    sub _filter_by_negative_threes;
    sub _filter_by_zeros;
    sub filter_and_print;
    sub _find_coserved_filter_grps;
    sub _print_filtered_data;
    
    sub get_grp_count_href;
    sub get_table_obj;
    
    ##################
    ## CONSTRUCTOR ##
    ##################
    
    sub new {
        my ($class, $da_table_obj) = @_;
        
        #Bless a scalar to represent the new object
        my $new_obj = bless \do{my $anon_scalar}, $class;
        
        if ( $da_table_obj->isa("DaTable") ) {
            $da_table_obj->check_for_unset_values();
            $da_table_objects{ident $new_obj} = $da_table_obj;
            $new_obj->_set_counts();
            
        }
        else {
            croak "Must pass in a DaTable object";
        }
        
        return $new_obj;
    }
    
    #################
    ## SUBROUTINES ##
    #################
    
    sub _set_counts {
        my ($self) = @_;
        my $da_tbl_obj = $self->get_table_obj();
        my $grp_order_aref = $da_tbl_obj->get_grp_order_aref();
        my @da_count;
        
        for (my $i = 0; $i < scalar @$grp_order_aref; $i++) {
            my $group_da_aref = $da_tbl_obj->get_group($grp_order_aref->[$i]);
            #create a new href for each group
            my $href = {    1 => 0,
                                0 => 0,
                                -1 => 0,
                                -2 => 0,
                                -3 => 0,
                                'total' => 0,
                            };
            #add the da counts for the grp from every genome
            foreach my $da_value (@$group_da_aref) {
                $href->{$da_value}++;
                $href->{total}++;
            }
            push @da_count, $href;
        }
        
        $da_counts_per_grp_aref_href{ident $self} = \@da_count;
        return 1;
    }
    
    #############
    ## GETTERS ##
    #############
    
    sub get_grp_count_aref {
        my ($self) = @_;
        return $da_counts_per_grp_aref_href{ident $self};
    }
    sub get_table_obj {
        my ($self) = @_;
        return $da_table_objects{ident $self};
    }
    
    #############
    ## FILTERS ##
    #############
    
    # All filters will be filtered using percentages as whole numbers
    # Returns the groups that pass the filters
    
    sub _filter_by_ones {
        my ($self, $percentage, $filter_out_below) = @_;
        my $count_aref_hrefs = $self->get_grp_count_aref();
        my $grp_order_aref = ($self->get_table_obj())->get_grp_order_aref();
        my %passing_grps;
        
        #Get the one count and total count
        for (my $i = 0; $i < scalar @$count_aref_hrefs; $i++) {
            my $href = $count_aref_hrefs->[$i];
            my $one_count = $href->{1};
            my $total = $href->{total};
            
            my $actual_percentage = ($one_count/$total)*100;
            if ( $filter_out_below =~ qr/false/i ) {
                if ( $actual_percentage <= $percentage ) {
                    $passing_grps{$grp_order_aref->[$i]} = $i;
                }
            }
            else {
                if ( $actual_percentage >= $percentage ) {
                    $passing_grps{$grp_order_aref->[$i]} = $i;
                }
            }
        }
        return \%passing_grps;
    }
    
    sub _filter_by_negative_twos {
        my ($self, $percentage, $filter_below) = @_;
        my $count_aref_hrefs = $self->get_grp_count_aref();
        my $grp_order_aref = ($self->get_table_obj())->get_grp_order_aref();
        my %passing_grps;
        
        #Get the one count and total count
        for (my $i = 0; $i < scalar @$count_aref_hrefs; $i++) {
            my $href = $count_aref_hrefs->[$i];
            my $neg_two_count = $href->{-2};
            my $total = $href->{total};
            
            my $actual_percentage = ($neg_two_count/$total)*100;
            if ( $filter_below =~ qr/false/i ) {
                if ( $actual_percentage <= $percentage ) {
                    $passing_grps{$grp_order_aref->[$i]} = $i;
                }
            }
            else {
                if ( $actual_percentage >= $percentage ) {
                    $passing_grps{$grp_order_aref->[$i]} = $i;
                }
            }
        }
        return \%passing_grps;
    }
    
    sub _filter_by_negative_ones {
        my ($self, $percentage, $filter_out_below) = @_;
        my $count_aref_hrefs = $self->get_grp_count_aref();
        my $grp_order_aref = ($self->get_table_obj())->get_grp_order_aref();
        my %passing_grps;
        
        #Get the one count and total count
        for (my $i = 0; $i < scalar @$count_aref_hrefs; $i++) {
            my $href = $count_aref_hrefs->[$i];
            my $neg_one_count = $href->{-1};
            my $total = $href->{total};
            
            my $actual_percentage = ($neg_one_count/$total)*100;
            if ( $filter_out_below =~ qr/false/i ) {
                if ( $actual_percentage <= $percentage ) {
                    $passing_grps{$grp_order_aref->[$i]} = $i;
                }
            }
            else {
                if ( $actual_percentage >= $percentage ) {
                    $passing_grps{$grp_order_aref->[$i]} = $i;
                }
            }
        }
        return \%passing_grps;
    }
    
    sub _filter_by_negative_threes {
        my ($self, $percentage, $filter_out_below) = @_;
        my $count_aref_hrefs = $self->get_grp_count_aref();
        my $grp_order_aref = ($self->get_table_obj())->get_grp_order_aref();
        my %passing_grps;
        
        #Get the one count and total count
        for (my $i = 0; $i < scalar @$count_aref_hrefs; $i++) {
            my $href = $count_aref_hrefs->[$i];
            my $neg_three_count = $href->{-3};
            my $total = $href->{total};
            
            my $actual_percentage = ($neg_three_count/$total)*100;
            if ( $filter_out_below =~ qr/false/i ) {
                if ( $actual_percentage <= $percentage ) {
                    $passing_grps{$grp_order_aref->[$i]} = $i;
                }
            }
            else {
                if ( $actual_percentage >= $percentage ) {
                    $passing_grps{$grp_order_aref->[$i]} = $i;
                }
            }
        }
        return \%passing_grps;
    }
    
    sub _filter_by_zeros {
        my ($self, $percentage, $filter_out_below) = @_;
        my $count_aref_hrefs = $self->get_grp_count_aref();
        my $grp_order_aref = ($self->get_table_obj())->get_grp_order_aref();
        my %passing_grps;
        
        #Get the one count and total count
        for (my $i = 0; $i < scalar @$count_aref_hrefs; $i++) {
            my $href = $count_aref_hrefs->[$i];
            my $zero_count = $href->{0};
            my $total = $href->{total};
            
            my $actual_percentage = ($zero_count/$total)*100;
            if ( $filter_out_below =~ qr/false/i ) {
                if ( $actual_percentage <= $percentage ) {
                    $passing_grps{$grp_order_aref->[$i]} = $i;
                }
            }
            else {
                if ( $actual_percentage >= $percentage ) {
                    $passing_grps{$grp_order_aref->[$i]} = $i;
                }
            }
        }
        return \%passing_grps;
    }
    #Change input on documentation
    sub filter_and_print {
        my ($self, $param_handler, $outfile) = @_; #each arg ref needs filter # w/ array of percentage and true/false depending if you want to filter out below values
        my $arg_href = $param_handler->get_filter_params();
        my @numbers_not_filtered;
        my @array_passing_hrefs;
        foreach my $key (keys %FILTER_NUMBERS) {
            if ( !$arg_href->{$key} ) {
                push @numbers_not_filtered, $key; # Make sure filtering everything the user wants
            }
            else {
                if ( $key eq "f1" ) {
                    my $passing_href = $self->_filter_by_ones( @{ $arg_href->{$key} }[0], @{ $arg_href->{$key} }[1] );
                    push @array_passing_hrefs, $passing_href;
                }
                elsif ( $key eq "f0" ) {
                    my $passing_href = $self->_filter_by_zeros( @{ $arg_href->{$key} }[0], @{ $arg_href->{$key} }[1] );
                    push @array_passing_hrefs, $passing_href;
                }
                elsif ( $key eq "f-1" ) {
                    my $passing_href = $self->_filter_by_negative_ones( @{ $arg_href->{$key} }[0], @{ $arg_href->{$key} }[1] );
                    push @array_passing_hrefs, $passing_href;
                }
                elsif ( $key eq "f-2" ) {
                    my $passing_href = $self->_filter_by_negative_twos( @{ $arg_href->{$key} }[0], @{ $arg_href->{$key} }[1] );
                    push @array_passing_hrefs, $passing_href;
                }
                elsif ( $key eq "f-3" ) {
                    my $passing_href = $self->_filter_by_negative_threes( @{ $arg_href->{$key} }[0], @{ $arg_href->{$key} }[1] );
                    push @array_passing_hrefs, $passing_href;
                }
            }
        }
        my $conserved_href = _find_conserved_filter_grps( \@array_passing_hrefs );
        $self->_print_filtered_data( $conserved_href, $outfile );
        if ( scalar @numbers_not_filtered > 0 ) {
            print join(", ", @numbers_not_filtered), " - number values that were not filtered\n";
        }
        return 1;
    }
    
    sub _find_conserved_filter_grps{
        my ($passing_grps_aref_hrefs) = @_;
        my $conserved_keys = $passing_grps_aref_hrefs->[0];
        if ( scalar( @$passing_grps_aref_hrefs ) == 1) {
            return $conserved_keys;
        }
        
        #Loop through and find grp values that passed every filter
        for (my $i = 1; $i < scalar(@$passing_grps_aref_hrefs); $i++) {
             foreach my $key ( keys %$conserved_keys ) {
                my $href = $passing_grps_aref_hrefs->[$i];
                if ( !exists $href->{$key} ) {
                    delete $conserved_keys->{$key};
                }
             }
        }
        return $conserved_keys;
    }
    
    sub _print_filtered_data {
        my ($self, $conserved_href, $outfile) = @_;
        my $da_table_obj = $self->get_table_obj();
        my $ordered_grps_aref = $da_table_obj->get_grp_order_aref();
        my $ordered_ids_aref = $da_table_obj->get_id_order_aref();
        open my $FDA_FH, ">", $outfile;
        #set up the filtered file
        print $FDA_FH "grp_id\t", join("\t", @$ordered_ids_aref), "\n" ;
        
        foreach my $grp ( @$ordered_grps_aref ) {
            if ( exists $conserved_href->{$grp} ) {
                print $FDA_FH "$grp\t";
                my $grp_da_values_aref = $da_table_obj->get_group($grp);
                print $FDA_FH join("\t", @$grp_da_values_aref), "\n";
            }
        }
        close($FDA_FH);
        return 1;
    }
    
}

1;

__END__

=head1  DaFilter

    DaFilter -  This object holds a DaTable and will perform filtering on the
                DaTable's groups. You can specify the Differential Abundance values you would like to filter (1,0,-1,-2,-3). Currently the filtering is done by a percentage of DA counts for the group compared to the total.
                
=head1  Version 0.0.1

    This documentation refers to DaFilter.pm
    
=head1  Included Modules

use DaTable;
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

=head1  Inherit
    
    NA
    
=head1  Synopsis


=head1  Description

This object holds the differentially abundant counts and DaTable. You can use this object to filter out groups from the DaTable object matrix and print it out into a file. The filtering will not change the DaTable in any way so multiple different filtering operations can be done in one run. Filtering of multiple DA values at once can also be acheived.

=head1  Methods

    sub new;
    sub get_grp_count_href;
    sub get_table_obj;
    sub filter_and_print;

    sub _set_counts;
    sub _filter_by_ones;
    sub _filter_by_negative_twos;
    sub _filter_by_negative_ones;
    sub _filter_by_negative_threes;
    sub _filter_by_zeros;
    sub _find_coserved_filter_grps;
    sub _print_filtered_data;
    
=head1  Methods Description

=head2 new()
    
    Title:      new
    Usage:      DaFilter->new( $da_table_obj )
    Function:   Creates a new DaFilter object
    Returns:    DaFilter object
    Args:       $da_table_obj   => This is a DaTable object that has a matrix of count data contained within it.
    Throws:     NA
    Comments:   Make sure that the DaTable object passed in has the matrix filled correctly
    See Also:   DaTable.pm
    
=head2 get_grp_count_aref()
    
    Title:      get_grp_count_aref
    Usage:      $filter_obj->get_grp_count_aref()
    Function:   Returns an array ref containing hashes with the counts of each groups differential abundance counts
    Returns:    Array reference containing hash references
    Args:       $filter_obj => a DaFilter object
    Throws:     NA
    Comments:   NA
    See Also:   NA
    
=head2 get_table_obj()

    Title:      get_table_obj
    Usage:      $filter_obj->get_table_obj()
    Function:   returns the Da_Table object that was passed
    Returns:    DaTable object
    Args:       $filter_obj => a DaFilter object
    Throws:     NA
    Comments:   NA
    See Also:   DaTable
    
=head2 filter_and_print()

    Title:      filter_and_print
    Usage:      $filter_obj->( $param_handler, $outfile );
    Function:   Checks the ParamHandler object for the filter params and determines what to filter and how. If "false" is passed everything below the filter will be kept. Vice versa for "true". It will then print the resulting matrix to the location passed
    Returns:    1
    Args:       $filter_obj => DaFilter object
                $param_handler  => This is a Param_Handler object containing the filter parameters within it.
                $outfile    => the location where you want resulting filtered count matrix to be printed to. Will be a txt file
    Throws:     Prints out numbers that were not filtered.
    Comments:   This will produce a text file at the location passed. It will contain the filtered count data originally from the matrix passed or the Da_Table object passed.
    See Also:   NA


=head1  Configuration And Environment

    DaFilter does not require any configuration files or environment variables
    
=head1  Dependencies

    NA
    
=head1  Incompatibilities

    None Reported
    
=head1  Bugs and Limitations

    No bugs have been reported
    
=head1  Author

    Nicholas Colaianni
    contact via C<< <ncolaian@live.unc.edu> >>
    
=head1  Licence And Copyright

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
