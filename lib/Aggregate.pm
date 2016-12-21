#! usr/bin/evn perl

package Aggregate;

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

{
    #Attributes
    my %param_objects;
    my %grp_order_aref;
    my %count_header_aref;
    
    #Required Paramaters
    Readonly::Hash my %PASS_PARAMS => map { $_ => 1} qw(
        count_dir
        dafe_dir
        annote_file_name
        count_file_name
        grp_genes_by
        grp_meta_file
    );
    
    ###Subroutines###
    sub aggregate; # Main aggregation subroutine that when called will aggregate files
    sub set_attributes;
    
    #getters
    sub get_param_object; 
    sub get_count_header_aref;
    sub get_grp_order_aref;
    
    #Hidden subs
    sub _set_grp_metafile_order_array;
    sub _create_order_table;
    sub _create_connection_btwn_geneid_grpid;
    sub _set_count_header_array;
    sub _grp_count_in_count_file;
    sub _create_aggregate_file;
    
    
    #################
    ## CONSTRUCTOR ##
    #################
    sub new {
        my ($class, $param_obj) = @_;
        
        #Bless a scalar to represent the new object
        my $new_obj = bless \do{my $anon_scalar}, $class;
        
        #Handle a passed in Param_handler
        if ( $param_obj->isa( "Param_handler" ) ) {
            $new_obj->set_attributes($param_obj);
        }
        else {
            croak "Need to pass in a valid Param_handler object";
        }
        
        return $new_obj;
    }
    
    #################
    ## SUBROUTINES ##
    #################
    
    sub set_attributes {
        my ($self, $param_obj) = @_;
        
        _check_param_obj($param_obj);
        $param_objects{ident $self} = $param_obj;
        $self->_set_grp_metafile_order_array();
        $self->_set_count_header_array();
        
        return 1;
    }
    
    sub aggregate {
        my ($self, $genome_id) = @_;
        
        my $id_grp_connection_href = $self->_create_connection_btwn_geneid_grpid($genome_id);
        my $count_href = $self->_grp_count_in_count_file($id_grp_connection_href,$genome_id);
        $self->_create_aggregate_file($count_href, $genome_id);
        
        return 1;
    }
    
    ############
    ## CHECKS ##
    ############
    
    sub _check_param_obj {
        my ($param_obj) = @_;
        my $href = $param_obj->get_params_href();
        
        #check to make sure all the necessary parameters are defined
        foreach my $att_name ( keys %PASS_PARAMS ) {
            if ( !$href->{$att_name} ) {
                MyX::Generic::Undef::Param->throw(
                    error => "$att_name not defined in the Param_handler object passed",
                    usage => "Pass a Param_handler object with all attributes defined",
                );
            }
        }
        $param_obj->check_edger_params(); #Make sure the param object is good
    }
    
    #############
    ## GETTERS ##
    #############
    
    sub get_param_object {
        my ($self) = @_;
        return $param_objects{ident $self};
    }

    sub get_count_header_aref {
        my ($self) = @_;
        return $count_header_aref{ident $self};
    }
    
    sub get_grp_order_aref {
        my ($self) = @_;
        return $grp_order_aref{ident $self};
    }
    
    #############
    ## SETTERS ##
    #############
    
    #NEED TO TALK ABOUT SCOTT
        #The GRP META FILE DOES NOT HAVE SAME ID AS METAFILE
    sub _set_grp_metafile_order_array {
        my ($self) = @_;
        my $param_obj = $self->get_param_object();
        my $grp_metafile_obj = file( $param_obj->get_grp_meta_file() );
        my @slurped_grp_meta = $grp_metafile_obj->slurp( chomp=>1, split=>qr/\t/ );
        my $grp_genes_by = $param_obj->get_grp_genes_by(); #Assume the first column is the id column
        my $grp_column = 0;
        my @ordered_grp;
        
        #find the column number of the grp names
        foreach my $part ( @{ $slurped_grp_meta[0] } ) {
            if ( $part =~ qr/$grp_genes_by/i ) {
                last;
            }
            $grp_column++;
        }
        
        #Create an array that will keep track of the grp order
        foreach my $line_aref ( @slurped_grp_meta ) {
            if ( $line_aref->[$grp_column] =~ qr/$grp_genes_by/i ) {
                next;
            }
            push @ordered_grp, $line_aref->[$grp_column];
        }
        #tests to make sure the grp_meta_array_isn't empty
        if ( scalar(@ordered_grp) == 0 ) {
            croak $grp_genes_by;
        }
        
        $grp_order_aref{ident $self} = \@ordered_grp;
        return 1;
    }
    
    sub _set_count_header_array {
        my ($self) = @_;
        my $param_obj = $self->get_param_object();
        #get an id so that count headers can be found in its count file
        my $include_file_obj = file( $param_obj->get_ref_include_file() );
        my @include = $include_file_obj->slurp( chomp=>1 );
        my $id = $include[0]; #Using the first id to find headers
        #get count file
        my $count_file_name = $param_obj->get_count_file_name();
        my $count_dir = $param_obj->get_count_dir();
        my $count_file_obj = file( "$count_dir/$id/$count_file_name" );
        my @count_file = $count_file_obj->slurp( chomp=>1 );
        my @first_line = split /\t/, $count_file[0];
        
        my @count_headers; #holds the headers from the count file
        #fill in the headers
        for (my $i = 1; $i < scalar(@first_line); $i++) {
            push @count_headers, $first_line[$i];
        }
        #set count headers
        $count_header_aref{ident $self} = \@count_headers;
        return 1;
    }
    

    sub _create_connection_btwn_geneid_grpid {
        my ($self, $id) = @_;
        my $param_obj = $self->get_param_object();
        my $annote_file_name = $param_obj->get_annote_file_name();
        my $dafe_dir = $param_obj->get_dafe_dir();
        my $gene_id_col_name = $param_obj->get_gene_id_col();
        my $grp_genes_by = $param_obj->get_grp_genes_by();
        my $annote_file_obj = file( "$dafe_dir/$id/$annote_file_name" );
        my @slurped_annote_file = $annote_file_obj->slurp( chomp=>1, split=>qr/\t/ );
        
        my $gene_id_col;
        my $grp_id_col;
        my $col_num = 0;
        my %gene_id_to_grp_id;
        
        #Determine the id column for both
        foreach my $part ( @{ $slurped_annote_file[0] } ) {
            if ( $part =~ qr/$gene_id_col_name/i ) {
                $gene_id_col = $col_num;
                $col_num++;
            }
            elsif ( $part =~ qr/$grp_genes_by/i ) {
                $grp_id_col = $col_num;
                $col_num++;
            }
            else{
                $col_num++;
            }
        }
        #Go through and create a hash that has the geneid as the key and the grp id as the value
        foreach my $line_aref ( @slurped_annote_file ) {
            if ( $line_aref->[$gene_id_col] =~ qr/$gene_id_col_name/i ) {
                next;
            }
            #sets the gene id equal to a grp id in a hash
            if( !$gene_id_to_grp_id{$line_aref->[$gene_id_col]} ) {
                if ($line_aref->[$grp_id_col] =~ /\,/) { #handles genes in more than 1 group
                    my @mm_array = split /,/, $line_aref->[$grp_id_col];
                    $gene_id_to_grp_id{$line_aref->[$gene_id_col]} = $mm_array[0];
                }
                else {
                    $gene_id_to_grp_id{$line_aref->[$gene_id_col]} = $line_aref->[$grp_id_col];
                }
            }#added to make sure duplicates aren't created
        }
        return \%gene_id_to_grp_id;
    }
    
    sub _grp_count_in_count_file {
        my ($self, $connections_href, $id) = @_;
        my $param_obj = $self->get_param_object();
        my $count_dir = $param_obj->get_count_dir();
        my $count_file_name = $param_obj->get_count_file_name();
        my $count_file_obj = file( "$count_dir/$id/$count_file_name" );
        my @slurped_count_file = $count_file_obj->slurp( chomp=>1, split=>qr/\t/ );
        
        my %grp_count_hash;
        
        # Keep track of genes that show up twice or have no group
        my @genes_wo_grp;
        my %genes_used;
        my @genes_counted_multiple_times;
        
        foreach my $line_aref ( @slurped_count_file ) {
            #skip first line
            if ($line_aref->[0] =~ qr/name/i ) {
                next;
            }
            #skip the last lines
            if ( $line_aref->[0] =~ qr/feature/i ) {
                last;
            }
            
            #Go through each line and add up counts into a hash
            my $marker = 0;
            my $grp_id;
                     
            foreach my $part ( @$line_aref ) {
                if ($marker == 0) {
                    #handles genes that do not have a kog
                    if (!$connections_href->{$part}) {
                        push @genes_wo_grp, $part;
                        last;
                    }
                    #check to make sure gene hasn't been counted yet
                    if ( $genes_used{$part} ) {
                        push @genes_counted_multiple_times, $part;
                        last;
                    }
                    $genes_used{$part} = 1;
                    $grp_id = $connections_href->{$part}; #save the part
                    #create a new key with count array in hash with same size as headers
                    if ( !$grp_count_hash{$grp_id} ) {
                        $grp_count_hash{$grp_id} = [];
                        my $i = 0;
                        while ( $i < scalar(@{$self->get_count_header_aref()} ) ) {
                            push @{ $grp_count_hash{$grp_id} }, 0;
                            $i++;
                        }
                    }
                }
                #this adds the grp counts to each header
                else {
                    $grp_count_hash{$grp_id}->[($marker-1)] = ( $grp_count_hash{$grp_id}->[($marker-1)] ) + $part;
                }
                $marker++;
                    
                #Tell user something happened but it's not enough to kill prgrm
            }
        }
        if (scalar(@genes_wo_grp) > 0) {
            carp("These genes are not associated with a group: "
                . join(", ", @genes_wo_grp) . "\n");
            }
        if (scalar(@genes_counted_multiple_times) > 0) {
            carp("These genes were in the count file more than once: "
                . join(", ", @genes_counted_multiple_times) . "\nThe first "
                . "instance of the count was used and the duplicate was excluded\n");
            }
        return \%grp_count_hash;
    }

    sub _create_aggregate_file {
        my ($self, $count_href, $id) = @_;
        #get everything i'm going to need to create the file
        my $param_obj = $self->get_param_object();
        my $ordered_grp_aref = $self->get_grp_order_aref();
        my $count_headers_aref = $self->get_count_header_aref();
        #get things from param_obj
        my $grp_genes_by = $param_obj->get_grp_genes_by();
        my $count_dir = $param_obj->get_count_dir();
        my $count_file_name = $param_obj->get_count_file_name();
        
        #initializing the outfile
        $count_file_name =~ s/\.txt//;
        my $file_name .= $count_file_name . "_$grp_genes_by" . "_agg.txt";
        
        if ( !-d "$count_dir/$id" ) {
            MyX::Generic::DoesNotExist::Dir->throw(
                error => "$id is not a directory within $count_dir",
                dir_name => $id,
            );
        }
        
        open(my $OUT, ">", "$count_dir/$id/$file_name" ) ||
            croak("cannot open file: $count_dir/$id/$file_name");
        
        #print outfile header
        print $OUT $grp_genes_by, "Id\t", join( "\t", @$count_headers_aref ), "\n";
        
        #print out the grp counts in order
        foreach my $grp_id ( @$ordered_grp_aref ) {
            if ( $count_href->{$grp_id} ) {
                print $OUT "$grp_id\t", join( "\t", @{ $count_href->{$grp_id} } ), "\n";
            }
            #print out grps that did not have data -> will print 0's
            else {
                my $out_line = "$grp_id";
                foreach my $col ( @$count_headers_aref ) {
                    $out_line .= "\t0";
                }
                print $OUT "$out_line\n";
            }
        }
        close($OUT);
        
        if ( !-e "$count_dir/$id/$file_name" || -z "$count_dir/$id/$file_name" ) {
            MyX::Generic::DoesNotExist::File->throw(
                error => "$file_name could not be made in the count_dir of $id",
                file_name => "$count_dir/$id/$file_name",
            );
        }
        
        return 1;
    }
}

1;

__END__

=head1 Aggregate Object

Aggregate - Takes a Param_handler object and creates an object that can
            aggregate count data for genes to a group count of your choosing.

=head1 VERSION

This documentation refers to Aggregate 0.0.1

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

=head1 SYNOPSIS

=head2 Object Creation

    use Aggregate;
    my $aggregate_obj = Aggregate->new( $param_handler_obj );

=head2 Getters

    my $passed_param_obj = $aggregate_obj->get_param_object();
    my $count_header_aref = $aggregate_obj->get_count_header_aref();
    my $grp_order_aref = $aggregate_obj->get_grp_order_aref();

=head2 Main Aggregate Method

    $aggregate_obj->aggregate($genome_id);

=head1 DESCRIPTION

Aggregate is an object to store group information and a parameter object.
It is used to ensure that gene counts are aggregated by the same group,
and the file printed are all in the same order. These files are created by
the aggregate() function. These files are made to be compatible with edgeR for
further data processing.

The getters provided are for checking the passed parameter object and the
printing mechanics of the aggregate() function.

=head1 METHODS

=over

    new
    set_attributes
    get_param_object
    get_count_header_aref
    get_grp_order_aref
    
    _set_grp_metafile_order_array
    _create_order_table
    _create_connection_btwn_geneid_grpid
    _set_count_header_array
    _grp_count_in_count_file
    _create_aggregate_file
    
=back    

=head1 METHODS DESCRIPTION

=head2 new()

    Title:      new
    Usage:      Aggregate->new( $param_obj );
    Function:   This creates a new Aggregate object and calls set_attributes
    Returns:    Aggregate
    Args:       -$param_obj => A Param_handler object. Make sure to check the
                               edgeR params before passing
    Throws:     MyX::Generic::Undef::Param
                check_edger_param errors from Param_handler
                croak "Need to pass in a valid Param_handler object"
                croak "Empty meta-array"    
    Comments:   Make sure to check the edger params before passing the object
                in. The params are checked again but the errors could be handled
                better if they are checked before
    See Also:   Param_handler
    
=head2 set_attributes()

    Title:      set_attributes
    Usage:      $aggregate_object->set_attributes( $param_obj );
    Function:   Calls the functions _check_param_obj(), _set_grp_metafile_order_array(),
                and _set_count_header_array(). Is called with new but can be used to change
                the internals of an already instantuated Aggregate object.
    Returns:    NA
    Args:       -$aggregate_object => a blessed Aggregate object
                -$param_obj => a blessed Param_handler object
    Throws:     MyX::Generic::Undef::Param
                check_edger_param errors from Param_handler
                croak "Empty meta-array"
    Comments:   This is called automatically when new() is called
    See Also:   NA
    
=head2 aggregate()

    Title:      aggregate
    Usage:      $aggregate_object->aggregate( $genome_id );
    Function:   This method calls on _create_connection_btwn_geneid_grp_id(),
                _grp_count_in_count_file(), and _create_aggregate_file(). Finds
                the gene counts associated with that genome and aggregates them
                to the group in the saved Param_handler object.
    Returns:    Tabulated Txt File
    Args:       -$aggregate_object => a blessed Aggregate object
                -$genome_id => Genome id within the count_dir from Param_handler
                               object.
    Throws:     MyX::Generic::DoesNotExist::Dir
                MyX::Generic::DoesNotExist::File             
    Comments:   Make sure the id exists in the count_dir.
    See Also:   NA
    
=head2 get_param_object()

    Title:      get_param_object
    Usage:      $aggregate_object->get_param_object();
    Function:   This method returns the internal Param_handler object passed in
    Returns:    Param_handler
    Args:       -$aggregate_object => a blessed Aggregate object
    Throws:     NA
    Comments:   NA
    See Also:   Param_handler
    
=head2 get_count_header_aref()

    Title:      get_count_header_aref
    Usage:      $aggregate_object->get_count_header_aref();
    Function:   This method returns an array reference containing the headers of
                the count files in the count files found in the count directory
                in the Param_handler object.
    Returns:    Array Reference
    Args:       -$aggregate_object => a blessed Aggregate object
    Throws:     NA
    Comments:   NA
    See Also:   Param_handler
    
=head2 get_grp_order_aref()

    Title:      get_grp_order_aref
    Usage:      $aggregate_object->get_grp_order_aref();
    Function:   Ruturns an array reference of the ordered group ids used to
                print the files from aggregate() in the same order each time
                aggregate() is called.
    Returns:    Array Reference
    Args:       -$aggregate_object => a blessed Aggregate object
    Throws:     NA
    Comments:   NA
    See Also:   NA

=head1 CONFIGURATION AND ENVIRONMENT

    Aggregate requires no configuration files or environment variables.


=head1 DEPENDENCIES

    Requires the Param_handler module, because you must pass in a Param_handler.

=head1 INCOMPATIBILITIES

    None reported.

=head1 BUGS AND LIMITATIONS

No bugs have been reported.

Limitations =>  For any file where an Id column name is not either asked for or
                in the param handler passed, the first column is assumed to be
                the Id column

Please report any bugs or feature requests

=head1 AUTHOR

Nicholas Colaianni
contact via C<< <ncolaian@live.unc.edu> >>

=head1 LICENCE AND COPYRIGHT

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

=cut
