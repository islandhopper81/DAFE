#! usr/bin/env perl

package Justify;

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
use Bio::TreeIO;

{
    #Attributes
    my %trimmed_ordered_meta_aref;
    my %ordered_ids_aref;
    my %meta_file_headers;
    my %input_files_href;
    
    #Required Parameters or param object
    Readonly::Hash my %PASS_PARAMS => map { $_ => 1} qw(
        include_file
        exclude_file
        meta_data_file
        tree_file
        count_dir
        count_file_name
        id_col
    );
    
    #Subroutines
    sub _justify_information;
    sub spew_ordered_meta_file;
    sub spew_ordered_ids;
    
    #getters
    sub get_to_meta_aref;
    sub get_ordered_ids_aref;
    sub get_saved_href;
    
    #setters
    sub _create_trimmed_ordered_meta_a;
    sub _create_ordered_ids;
    sub _determine_included_ids;
    sub _create_ordered_tree_aref;
    
    #checks
    sub _check_count_dir; # Make sure each directory has count information
    sub _check_ids; #make sure all the included id's are in the meta_data_file
    sub _check_tree; #Need to make sure tree contains only information in the included file and nothing from excluded file
    sub _check_for_all_input;
    
    
    #############
    #CONSTRUCTOR#
    #############
    sub new {
        my ($class, $arg) = @_;

        #Bless a scalar to represet the new object
        my $new_obj = bless \do{my $anon_scalar}, $class;
        
        #Justify Information runs all checks and creates saved attributes
        if ( ref $arg eq "HASH" ) {
            $new_obj->_justify_information($arg);
        }
        elsif ( $arg->isa( "Param_handler" ) ) {
            my ($arg_href) = {
                include_file    => $arg->get_ref_include_file(),
                exclude_file    => $arg->get_ref_exclude_file(),
                meta_data_file  => $arg->get_ref_meta_file(),
                tree_file       => $arg->get_tree(),
                count_dir       => $arg->get_count_dir(),
                count_file_name => $arg->get_count_file_name(),
                id_col          => $arg->get_genome_id_col(),
            };
            $new_obj->_justify_information($arg_href);
        }
        else {
            croak "Must pass in a hash ref with all the needed information or a Param_handler object";
        }
        return $new_obj;
    }
    
    sub _justify_information {
        my ($self, $href) = @_;
        
        _check_for_all_input($href);
        _check_for_bad_tags($href);
        $href->{include_href} = _determine_included_ids($href);
        _create_ordered_tree_aref($href);
        _check_tree($href);
        _check_count_dir($href);
        _check_ids($href); # Check Id's to make sure they are in the metadata file
        $self->_set_ordered_ids($href);
        $self->_set_meta_file_save_header($href);
        $input_files_href{ident $self} = $href;
        
        $self->_check_saved_attributes;
        
        
        return 1;
    }
    
    #Getters
    sub get_passed_params_href {
        my ($self) = @_;
        
        my $href = {
            include_file => $self->get_passed_include_file(),
            exclude_file => $self->get_passed_exclude_file(),
            meta_data_file => $self->get_passed_meta_data_file(),
            tree_file => $self->get_passed_tree_file(),
            count_dir => $self->get_passed_count_dir(),
            count_file_name => $self->get_passed_count_file_name(),
            id_col => $self->get_passed_id_col(),
        };
        
        return $href;
    }
    
    sub get_trimmed_ordered_meta_aref {
        my ($self) = @_;
        return $trimmed_ordered_meta_aref{ident $self};
    }
    
    sub get_ordered_ids_aref {
        my ($self) = @_;
        return $ordered_ids_aref{ident $self};
    }
    
    sub _get_meta_file_header {
        my ($self) = @_;
        return $meta_file_headers{ident $self};
    }
    
    #Pass_params_getters
    sub get_passed_include_file {
        my ($self) = @_;
        return $input_files_href{ident $self}->{include_file};
    }
    sub get_passed_exclude_file {
        my ($self) = @_;
        return $input_files_href{ident $self}->{exclude_file};
    }
    sub get_passed_meta_data_file {
         my ($self) = @_;
        return $input_files_href{ident $self}->{meta_data_file};
    }
    sub get_passed_tree_file {
         my ($self) = @_;
        return $input_files_href{ident $self}->{tree_file};
    }
    sub get_passed_count_dir {
         my ($self) = @_;
        return $input_files_href{ident $self}->{count_dir};
    }
    sub get_passed_count_file_name {
         my ($self) = @_;
        return $input_files_href{ident $self}->{count_file_name};
    }
    sub get_passed_id_col {
         my ($self) = @_;
        return $input_files_href{ident $self}->{id_col};
    }
    
    #Printers
    sub spew_trimmed_ordered_meta_file {
        my ($self,$out_file) = @_;
        my $meta_aref = $self->get_trimmed_ordered_meta_aref();
        my $header = $self->_get_meta_file_header();
        
        open(my $OUT, ">", $out_file) ||
            croak("cannot open file: $out_file\n");
            
        print $OUT "$header\n";
        foreach my $line ( @$meta_aref ) {
            print $OUT "$line\n";
        }
        close($OUT);
        return 1;
    }
    
    sub spew_ordered_ids {
        my ($self,$out_file) = @_;
        my $id_aref = $self->get_ordered_ids_aref();
        
        open(my $OUT, ">", $out_file) || croak("cannot open file: $out_file\n");
        
        foreach my $id ( @$id_aref ) {
            print $OUT "$id\n";
        }
        close($OUT);
        return 1;
    }

    #Checks
    sub _check_for_all_input {
        my ($href) = @_;
        
        foreach my $input (keys %PASS_PARAMS) {
            if ( !$href->{$input} ) {
                #Don't throw error if the missing file is an exclude file
                if ( $input eq "exclude_file" ) {
                    next;
                }
                MyX::Generic::Undef::Param->throw(
                    error => Dumper($href) . "$input is not defined in the passed in href, and is necessary",
                    usage => $input,
                );
            }
        }
    }
    
    sub _check_for_bad_tags {
        my ($href) = @_;
        
        foreach my $tag (keys %$href) {
            if ( !$PASS_PARAMS{$tag} ) {
                MyX::Generic::Undef::Attribute->throw(
                    error => "Unknown attribute in the passed parameter hash",
                    att_name => $tag,
                );
            }
        }
    }
    
    sub _determine_included_ids {
        my ($href) = @_;
        my $include_file_obj = file( $href->{include_file} );
        my @include = $include_file_obj->slurp( chomp => 1 );
        my $include_file_href = {};
        
        if ( $href->{exclude_file} ) {
            my $exclude_amount = 0; #keep track of how many id's are excluded -> sanity check
            my $exclude_file_obj = file( $href->{exclude_file} );
            my $slurped_exclude = $exclude_file_obj->slurp();
            foreach my $id ( @include ) {
                if ( $slurped_exclude !~ qr/$id/i ) {
                    $include_file_href->{$id} = 1;
                }
                else {
                    $exclude_amount++;
                }
            }
            #Check to make sure all the exclude ids match an id in the include file
            my @exclude_ids = split( /\n/, $slurped_exclude);
            if ( scalar(@exclude_ids) != $exclude_amount ) {
                MyX::Generic::File->throw(
                    error => "Some excluded ids were not found in the include file. Check to make sure the exclude file is correct",
                );
            }
        }
        else {
            foreach my $id ( @include ) {
                $include_file_href->{$id} = 1;
            }
        }
        return $include_file_href;
    }
    
    sub _create_ordered_tree_aref {
        my ($href) = @_;
        my $include_href = $href->{include_href};
        my $tree_obj = new Bio::TreeIO( -file => $href->{tree_file},
                                        -format => "newick");
        #Parses through the newick file
        my @ordered_tree;
        while ( my $tree = $tree_obj->next_tree() ) {
            for my $leaves ($tree->get_leaf_nodes()) {
                push @ordered_tree, $leaves->{_id};
            }
        }
        
        $href->{tree_aref} = \@ordered_tree;
        return 1;
    }
    
    sub _check_count_dir {
        my ($href) = @_;
        my $file_name = $href->{count_file_name};
        my $directory = $href->{count_dir};
        my $id_href = $href->{include_href};
        
        foreach my $id (keys %$id_href) {
            if (!-e "$directory/$id/$file_name" || -z "$directory/$id/$file_name") {
                MyX::Generic::DoesNotExist::File->throw(
                    error => "$id does not have a count file within its count dir. Check count_file_name and count_file",
                    file_name => "$directory/$id/$file_name",
                );
            }
        }
        return 1;
    }
    
    sub _check_ids {
        my ($href) = @_;
        my $id_href = $href->{include_href};
        my $meta_data_file_obj = file( $href->{meta_data_file} );
        my $slurped_meta_file = $meta_data_file_obj->slurp();
        my @unmatched_ids; #Holds ids that are not found in the metadata file
        
        #Find the ids that aren't found in the meta_file
        foreach my $id ( keys %$id_href ) {
            if ( $slurped_meta_file !~ qr/$id/i ) {
                push @unmatched_ids, $id;
            }
        }
        #Throw error if any ids are unmatched
        if ( scalar(@unmatched_ids) > 0 ) {
            my $missing_ids = join(", ", @unmatched_ids);
            MyX::Generic::BadValue->throw(
                error => "$missing_ids are not found in the meta_file",
                value => $missing_ids,
            );
        }
        return 1;
    }
    
    sub _check_tree {
        my ($href) = @_;
        my $id_href = $href->{include_href};
        my $tree_aref = $href->{tree_aref};
        my @tree_ids_not_in_include; #This holds ids that are in tree but not in the included ids
        my @include_not_in_tree; #this holds ids that are supposed to be included, but are not in tree
        #check to see if id from the tree is not in the include file
        foreach my $id ( @$tree_aref ) {
            if ( !$id_href->{$id} ) {
                push(@tree_ids_not_in_include, $id);
            }
        }
        #check to see if id from include is not in the tree *****make hash with tree array reference then test if id is in hash
        foreach my $id ( keys %$id_href ) {
            my $count = 0;
            foreach my $tree_id ( @$tree_aref ) {
                if ( $tree_id eq $id ) {
                    $count = 0;
                    last;
                }
                else {
                    $count++;
                }
                if ( $count == scalar(@$tree_aref) ) {
                    push(@include_not_in_tree, $id);
                }
            }
        }
        #throw errors if ids do not match completely
        if ( scalar(@tree_ids_not_in_include) > 0 ) {
            my $ids = join(", ", @tree_ids_not_in_include);
            MyX::Generic::BadValue->throw(
                error => "$ids are found in the tree but are not in the include file",
                value => "$ids",
            );
        }
        if ( scalar(@include_not_in_tree) > 0 ) {
            my $ids = join(", ", @include_not_in_tree);
            MyX::Generic::BadValue->throw(
                error => "$ids are in the include file but are not included in the tree",
                value => $ids,
            );
        }
        return 1;
    }
    
    sub _check_saved_attributes {
        my ($self) = @_;

        #check that each attr is defined
        if ( !$trimmed_ordered_meta_aref{ident $self} ) {
            MyX::Generic::Undef::Attribute->throw(
                error => "Trimmed_ordered_meta_aref was not updated in the creation of the justify object",
                att_name => "Trimmed_ordered_meta_aref",
            );
        }
        if ( !$ordered_ids_aref{ident $self} ) {
            MyX::Generic::Undef::Attribute->throw(
                error => "ordered_ids_aref was not updated in the creation of the justify object",
                att_name => "ordered_ids_aref",
            );
        }
        if ( !$meta_file_headers{ident $self} ) {
            MyX::Generic::Undef::Attribute->throw(
                error => "meta_file_headers was not updated in the creation of the justify object",
                att_name => "meta_file_headers",
            );
        }
        if ( !$input_files_href{ident $self} ) {
            MyX::Generic::Undef::Attribute->throw(
                error => "input_files_href was not updated in the creation of the justify object",
                att_name => "input_files_href",
            );
        }
        
        return 1;
    }
   
    #Setters
    sub _set_ordered_ids {
        my ($self,$href) = @_;
        $ordered_ids_aref{ident $self} = $href->{tree_aref};
        if (!$ordered_ids_aref{ident $self}) {
            croak("setting is messed up");
        }
        return 1;
    }
    
    sub _set_meta_file_save_header {
        my ($self,$href) = @_;
        my $ordered_ids = $href->{tree_aref};
        my $meta_file_obj = file( $href->{meta_data_file} );
        my $col_id_name = $href->{id_col};
        my @slurped_meta_file = $meta_file_obj->slurp( chomp => 1, split => qr/\t/ );
        my %meta_id_hash;
        my @ordered_meta_data_array;
        my $header = join("\t", @{$slurped_meta_file[0]} );
        my $id_col_number = 0;
        
        #find id column in metadata
        foreach my $meta_col ( @{$slurped_meta_file[0]} ) {
            if ( $meta_col =~ qr/$col_id_name/i ) {
                last;
            }
            else {
                $id_col_number++;
            }
        }
        #Fill hash for easy searching
        foreach my $line_aref ( @slurped_meta_file ) {
            $meta_id_hash{ ($line_aref->[$id_col_number]) } = $line_aref;
        }
        
        #Go Through and create an ordered metadata file.
        foreach my $id ( @$ordered_ids ) {
            my $meta_line = join( "\t", @{$meta_id_hash{$id}} );
            push @ordered_meta_data_array, $meta_line;
        }
        
        #sanity_check
        if (scalar(@$ordered_ids) != scalar(@ordered_meta_data_array)) {
            croak("Error occured when creating the ordered meta file");
        }
        $trimmed_ordered_meta_aref{ident $self} = \@ordered_meta_data_array;
        $meta_file_headers{ident $self} = $header;
        
        return 1;
    }
}

1;

__END__

=head1 Justify Object

This object makes sure the include file, count_database, metafile, and tree file
all the information needed for downstream analysis.

=head1 VERSION

This documentation refers to Justify 0.0.1

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

    use Justify;
    my $justify_object = Justify->new( {include_file    => $included_ids,
                                        exclude_file    => $excluded_ids,
                                        meta_data_file  => $metadata_file,
                                        tree_file       => $tree,
                                        count_dir       => $count_directory,
                                        count_file_name => $file_name,
                                        id_col          => $name_of_id_col,
                                        }
                                     );
    OR
    my $justify_object = Justify->new( Param_handler );

=head2 Getters

    my $passed_params_href = $justify_object->get_passed_params_href(); my
    my $trimmed_ordered_metafile_aref =
            $justify_object->get_trimmed_ordered_meta_aref();
    my $ordered_ids_aref = $justify_object->get_ordered_ids_aref();
    PICK WHAT YOU WANT
    my $passed_{include_file, exclude_file, meta_data_file, tree_file,
                count_dir, count_file_name, id_col} =
                $justify_object->get_passed_{include_file, exclude_file,
                                             meta_data_file, tree_file,
                                             tree_file, count_dir,
                                             count_file_name, id_col}();

=head2 Printers

    $justify_object->spew_trimmed_ordered_meta_file( $oufile );
    $justify_object->spew_ordered_ids( $outfile );

=head1  DESCRIPTION

Passing a Param_handler or a hash reference containing 6 things can have 1
optional parameter:
1)  Include File
2)  Metadata File
3)  Tree File (In Newick Format)
4)  Count Directory
5)  Count File Name In Directory
6)  ID Column Name
Optional => Exclude File

The object will make sure the id's that are to be included in downstream
analysis are included inthe tree file, final include file, metadata file, and
have a count dir with a count file.

The object will then create and hold an ordered meta array reference and id
reference. The orderis based upon the tree order passed in. The object will also
hold the input parameters if they are needed or they would like to be looked at
as a reference or to find mistakes.

=head1 METHODS

    new
    get_passed_params_href
    get_trimmed_ordered_meta_aref
    get_ordered_ids_aref
    get_passed_include_file
    get_passed_exclude_file
    get_passed_meta_data_file
    get_passed_tree_file
    get_passed_count_dir
    get_passed_count_file_name
    get_passed_id_col
    spew_trimmed_ordered_meta_file
    spew_ordered_ids

=head1 METHODS DESCRIPTION

=head2 new()

    Title:      new
    Usage:      Justify->new( { include_file    => $included_ids,
                                exclude_file    => $excluded_ids,
                                meta_data_file  => $metadata_file,
                                tree_file       => $tree,
                                count_dir       => $count_directory,
                                count_file_name => $file_name,
                                id_col          => $name_of_id_col,
                                }
                            );
                OR
                Justify->new( $param_handler );
    Function:   Creates a new Justify object and calls on _justify_information()
    Returns:    Justify
    Args:       -   $included_ids => File containing all the ids to be included
                    in the analysis.
                -   $excluded_ids => File containing all the ids to be excluded
                    in the analysis
                -   $metadata_file => File containing all the metadata info for
                    all the ids to be included in the analysis
                -   $tree => Newick file that represents the tree of the
                    included ids
                -   $count_directory => Directory containing all the count
                    directories for each included id
                -   $file_name => The name of the count file within each id's
                    count directory
                -   $name_of_id_col => The name of the id column in the metadata
                    file
                OR
                -   $param_handler => A Param_handler object that has been
                    through check_edger_params()
    Throws:     MyX::Generic::Undef::Param
                MyX::Generic::Undef::Attribute
                MyX::Generic::File
                MyX::Generic::DoesNotExist::File
                MyX::Generic::BadValue
                croak(s)
    Comments:   This object runs through all its steps during the creation of
                the object with the internal _justify_information method.
                IMPORTANT
                If you pass in a hashref -> the hashref will contain two new
                keys with values. The keys are include_href and tree_aref. The
                values will be a hash reference and array reference respectively
    See Also:   Param_handler

=head2 get_passed_params_href()

    Title:      get_passed_params_href
    Usage:      $param_href = $justify_object->get_passed_params_href();
    Function:   returns a hash ref containing all the things used to create the
                Justify object
    Returns:    Hash Reference
    Args:       NA
    Throws:     NA
    Comments:   NA
    See Also:   NA
    
=head2 get_passed_include_file()

    Title:      get_passed_include_file
    Usage:      $include_file_handle = $justify_object->get_passed_include_file();
    Function:   Returns the file handle of the include file 
    Returns:    String
    Args:       NA
    Throws:     NA
    Comments:   NA
    See Also:   NA

=head2 get_passed_exclude_file()

    Title:      get_passed_exclude_file
    Usage:      $exclude_file_handle = $justify_object->get_passed_exclude_file();
    Function:   returns the file handle of the exclude file
    Returns:    String
    Args:       NA
    Throws:     NA
    Comments:   NA
    See Also:   NA
    
=head2 get_passed_meta_data_file()

    Title:      get_passed_meta_data_file
    Usage:      $meta_data_handle = $justify_object->get_passed_meta_data_file();
    Function:   returns the file handle of the metadata file
    Returns:    String
    Args:       NA  
    Throws:     NA
    Comments:   NA
    See Also:   NA
    
=head2 get_passed_tree_file()

    Title:      get_passed_tree_file
    Usage:      $tree_file_handle = $justify_object->get_passed_tree_file();
    Function:   Returns the file handle of the tree file
    Returns:    String
    Args:       NA
    Throws:     NA
    Comments:   NA
    See Also:   NA
    
=head2 get_passed_count_dir()

    Title:      get_passed_count_dir
    Usage:      $count_dir_path = $justify_object->get_passed_count_dir(); 
    Function:   Returns the directory path of the count directory
    Returns:    String
    Args:       NA  
    Throws:     NA
    Comments:   NA
    See Also:   NA
    
=head2 get_passed_count_file_name()

    Title:      get_passed_count_file_name
    Usage:      $name = $justify_object->get_passed_count_file_name();
    Function:   Returns the name of the count files located in the id directory
    Returns:    String
    Args:       NA
    Throws:     NA  
    Comments:   NA
    See Also:   NA
    
=head2 get_passed_id_col()

    Title:      get_passed_id_col
    Usage:      $id_column_name = $justify_object->get_passed_id_col();
    Function:   Returns the name of the column that holds the id info in the
                metadata file
    Returns:    String
    Args:       NA
    Throws:     NA
    Comments:   NA
    See Also:   NA
    
=head2 spew_trimmed_ordered_meta_file()

    Title:      spew_trimmed_ordered_meta_file
    Usage:      $justify_object->spew_trimmed_ordered_meta_file($out_file);
    Function:   Creates a file containing only the necessary metadata files and
                they are also ordered by the tree fashion
    Returns:    NA  
    Args:       -   $out_file => path of the file to be created
    Throws:     NA
    Comments:   NA
    See Also:   NA
    
=head2 spew_ordered_ids()

    Title:      spew_ordered_ids
    Usage:      $justify_object->spew_ordered_ids($out_file);
    Function:   Creates a file containing only the ids being analyzed downstream
                and they are in the order of the tree.
    Returns:    NA
    Args:       -   $out_file => path of the file to be created
    Throws:     NA
    Comments:   NA  
    See Also:   NA
    
=head2 get_ordered_ids_aref()

    Title:      get_ordered_ids_aref
    Usage:      $justify_obj->get_ordered_ids_aref();
    Function:   Returns an arra reference containing the genome ids in the analysis ordered by the tree
    Returns:    Array Reference
    Args:       $justify_obj    =>  A Justify object
    Throws:     NA
    Comments:   NA
    See Also:   NA
    
=head2 get_trimmed_ordered_meta_aref()

    Title:      get_trimmed_ordered_meta_aref
    Usage:      $justify_obj->get_trimmed_ordered_meta_aref();
    Function:   Returns an array reference of the ordered metadatafile with the genomes excluded out of the file
    Returns:    Array Reference
    Args:       $justify_obj    =>  A Justify object
    Throws:     NA
    Comments:   NA
    See Also:   NA
    
=head1 CONFIGURATION AND ENVIRONMENT

    Justify requires no configuration files or environment variables.

=head1 DEPENDENCIES

    NA
    
=head1 INCOMPATIBILITIES

    None reported.

=head1 BUGS AND LIMITATIONS

No bugs have been reported.

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
