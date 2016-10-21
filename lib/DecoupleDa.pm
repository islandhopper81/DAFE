#! usr/bin/evn perl

package DecoupleDa;

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
    #Attributes
    my %param_handler_obj;
    my %count_file_names;
    
    ### Subroutines ###
    sub decouple; # Takes in a genome and searches its annotation file to search if the -2's from the differentially abundant statistics are b/c of no info or are not present in the genome
    sub _remake_da_with_decoupled;
    sub _check_negative_two;
    
    sub _set_param_handler;
    sub set_count_file_name;
    
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
        my $grp_aref = $self->get_ordered_grp_aref();
        #Need to get count and dafe directories file handles
        my $count_f_name = $self->get_count_file_name();
        my $count_dir = $param_obj->get_count_dir();
        my $count_file = "$count_dir/$genome_id/$count_f_name";
        my $dafe_dir = $param_obj->get_dafe_dir();
        my $annote_name = $param_obj->get_annote_file_name();
        my $annote_file = "$dafe_dir/$genome_id/$annote_name";
        
        #create file objects that will be read in
        my $count_file_obj = file($count_file);
        my $annote_file_obj = file($annote_file);
        
        
        #create a slurped array 
        my @count_file_a = $count_file_obj->slurp(chomp=>1, split=>qr/\t/);
        my $slurped_annote_file = $annote_file_obj->slurp();
        
        #go through each kog and check the DA Value
        foreach my $line ( @count_file_a ) {
            if ( $line->[1] eq "-2") {
                if ( _check_negative_two($line->[0],$slurped_annote_file) == 0 ){
                    $line->[1] = -3;
                }
            }
        }
        _remake_da_with_decoupled(\@count_file_a, $count_file); #print decoupled info in old file
        #return a hash reference of hash references
        return \@count_file_a;
    }
    
    sub _remake_da_with_decoupled {
        my ($count_file_aref_aref, $count_file);
        open (my $CF_FH, ">", "$count_file"); #opens file for writing note appending
        
        foreach my $line ( @$count_file_aref_aref ) {
            print join "\t", @$line;
        }
        close $CF_FH;
        
        return 1;
    }
    
    sub _check_negative_two {
        my ($grp_id, $slurped_annote_file) = @_;
        #check to see if the group id is in the genomes annotation file
        if ( $slurped_annote_file =~ qr/$grp_id/i ) {
            return 1;
        }
        else {
            return 0;
        }
    }
    
    
    ### SETTERS ###
    
    sub _set_param_handler {
        my ($self, $param_obj ) = @_;
        $param_handler_obj{ident $self} = $param_obj;
        return 1;
    }
    
    sub _set_count_file_name {
        my ($self, $param_obj) = @_;
        
        my $count_dir = $param_obj->get_count_dir();
        my $count_f_name = $param_obj->get_annote_file_name();
        my $grp_genes_by = $param_obj->get_grp_genes_by();
        $count_f_name = $count_f_name =~ s/\.txt//;
        $count_f_name .= "_$grp_genes_by" . "_agg_da_tbl.txt";
        $count_file_names{ident $self} = $count_f_name;
        
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
