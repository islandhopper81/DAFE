#! usr/bin/evn perl

use strict;
use warnings;
use Param_handler;
use Justify;
use Aggregate;
use Getopt::Long;
use Pod::Usage;
use Carp;
use Readonly;
use Path::Class;
use Data::Dumper;
use XML::Simple qw(:strict);
use Class::Std::Utils;
use Scalar::Util::Numeric qw(isneg isint isfloat);
use lib "/proj/cdjones_lab/ncolaian/MyX-Generic-v0.0.2/lib/MyX";
use MyX::Generic;
use YAML::XS qw(LoadFile);

# My Variables
my $help = 0;
my $man = 0;
my $xml_file;
my $yml_file; # Will hold the param file passed into the driver
my $out_dir; #place to put the ordered ids and meta_file

# Read in variables from the command line
GetOptions ('man'  => \$man,
            'help' => \$help,
            'xml_file|x=s' => \$xml_file,
            'yml_file|y=s' => \$yml_file,
            ) || die("There was an error in the command line arguements\n");

#kill program if neither a xml or txt file is passed
if (!defined $xml_file && !defined $yml_file) {
    croak("Must pass in a text file with newline characters, or an xml file.
          Do this by specifying the file with either -xml_file or -txt_file");
}

# Use Pod usage for the Manual and Help pages
if ( $help ) { pod2usage(0) }
if ( $man )  {pod2usage(-verbose => 3) }

# Subroutines that will occur in the main program

# MAIN #

my $param_obj;
#Create a Param_handler object and check the edgeR Params. Put this in the param hendler
if (defined $xml_file) {
    $param_obj = Param_handler->new( { xml_file => $xml_file } );
}
elsif (defined $yml_file) {
    #my $param_href = find_param_values_from_txt_file($yml_file);
    $param_obj = Param_handler->new( { yml_file => $yml_file } );
}
$param_obj->check_edger_params();

# Need to justify the tree, metadata file, count dir, and include file and print files
my $justify_obj = Justify->new( $param_obj );
$out_dir = $param_obj->get_out_dir();
my $ids_out = $out_dir . "/ordered_ids.txt";
my $meta_out = $out_dir . "/ordered_metafile.txt";
$justify_obj->spew_trimmed_ordered_meta_file( $meta_out );
$justify_obj->spew_ordered_ids( $ids_out );

# Need to then perform the aggregation of the count data
my $aggregate_obj = Aggregate->new( $param_obj );
my $ids_aref = $justify_obj->get_ordered_ids_aref(); # Gives ids in analysis
#loop trough id's and perform the aggregation of the data
foreach my $id ( @$ids_aref ) {
    $aggregate_obj->aggregate($id);
}


__END__

=head1 edgeR_driver
    
This driver uses the Param_handler, Aggregate, and Justify object. This program will check all the passed parameters from a parameter file. It will then create
an aggregated count file by grp. This file will be printed in the count directory in each id's directory. Two files will be created in the out directory passed in to the parameter file. The files will be an ordered id and metafile text file.

=head1 VERSION

This documentation refers to edgeR_driver 0.0.1

=head1 INCLUDED MODULES

use Param_handler;
use Justify;
use Aggregate;
use Getopt::Long;
use Pod::Usage;
use Carp;
use Readonly;
use Path::Class;
use Data::Dumper;
use XML::Simple qw(:strict);
use Class::Std::Utils;
use Scalar::Util::Numeric qw(isneg isint isfloat);
use lib "/proj/cdjones_lab/ncolaian/MyX-Generic-v0.0.2/lib/MyX";
use MyX::Generic;
use YAML::XS qw(LoadFile);

=head1 INHERIT

    NA
    
=head1 SYNOPSIS

    perl edgeR_driver.pl -{y|x}
    
    -yml_file | y =>    This is a yaml parameter file that contains all the     parameters needed for the analysis
    
    -xml_file | x =>    This is a xml parameter file that contains all the parameters needed for the analysis
    
    ATTENTION: You must pass either a yaml or xml file.
    
=head1 PARAMETERS
    
    ref_meta_file       =>  File Path       ->  String
    ref_include_file    =>  File Path       ->  String
    ref_exclude_file    =>  File Path       ->  String  (Optional)
    tree                =>  File Path       ->  String
    out_dir             =>  Directory Path  ->  String 
    annote_file_name    =>  File Name       ->  String
    grp_genes_by        =>  Group Names     ->  String
    gene_id_col         =>  Column Name     ->  String
    count_dir           =>  Directory Path  ->  String
    dafe_dir            =>  Directory Path  ->  String
    genome_id_col       =>  Column Name     ->  String
    metaG_meta_file     =>  File Path       ->  String
    metaG_include_file  =>  File Path       ->  String
    metaG_exclude_file  =>  File Path       ->  String  (Optional)
    count_file_name     =>  File Name       ->  String
    min_sample_count    =>  Sample Count    ->  Number
    min_sample_cpm      =>  Counts per Million ->   Number
    test_names          =>  Test Names      ->  String
    test_col_name       =>  Column Name     ->  String
    grp_meta_file       =>  File Path       ->  String

=head1 CONFIGURATION AND ENVIRONMENT

    edgeR_driver requires no configuration files or environment variables.

=head1 DEPENDENCIES

    Param_handler
    Justify
    Aggregate
    
=head1 INCOMPATIBILITIES

    None reported.

=head1 BUGS AND LIMITATIONS

No bugs have been reported.

Please report any bugs or feature requests	
	
=head1 Author

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

