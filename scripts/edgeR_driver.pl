#! usr/bin/evn perl

use strict;
use warnings;
use Param_handler;
use Justify;
use Aggregate;
use DaTable;
use DaFilter;
use DecoupleDa;
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
use File::Temp qw(tempfile tempdir);
use Log::Log4perl qw(:easy);
use Log::Log4perl::CommandLine qw(:all);
use Cwd qw(abs_path);
use File::Basename;


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

# Use Pod usage for the Manual and Help pages
if ( $help ) { pod2usage(0) }
if ( $man )  {pod2usage(-verbose => 3) }

#kill program if neither a xml or txt file is passed
if (!defined $xml_file && !defined $yml_file) {
    croak("Must pass in a text file with newline characters, or an xml file.
          Do this by specifying the file with either -xml_file or -txt_file");
}

# Set up Logger
my $logger = get_logger();

# Subroutines that will occur in the main program

# MAIN #

my $param_obj;
#Create a Param_handler object and check the edgeR Params. Put this in the param hendler
if (defined $xml_file) {
    $logger->info( "Creating a Param object with an xml file" );
    $logger->debug( $xml_file );
    $param_obj = Param_handler->new( { xml_file => $xml_file } );
}
elsif (defined $yml_file) {
    $logger->info( "Creating a Param object with an yaml file" );
    $logger->debug( $yml_file );
    #my $param_href = find_param_values_from_txt_file($yml_file);
    $param_obj = Param_handler->new( { yml_file => $yml_file } );
}
#$param_obj->set_Rsource_dir( "../R_lib");
$logger->info( "Checking the edgeR parameters" );
$param_obj->check_edger_params();

# Need to justify the tree, metadata file, count dir, and include file and print files
my $justify_obj = Justify->new( $param_obj );
$out_dir = $param_obj->get_out_dir();
my $ids_out = $out_dir . "/ordered_ids.txt";
my $meta_out = $out_dir . "/ordered_metafile.txt";

$logger->info( "Creation of an ordered metadata file" );
$logger->debug( $meta_out );
$justify_obj->spew_trimmed_ordered_meta_file( $meta_out );

$logger->info( "Creation of an ordered id file using the tree" );
$logger->debug( $ids_out );
$justify_obj->spew_ordered_ids( $ids_out );

# Need to then perform the aggregation of the count data
$logger->info( "Creation of an aggregate object, and aggregation of all the count data" );
my $aggregate_obj = Aggregate->new( $param_obj );
my $ids_aref = $justify_obj->get_ordered_ids_aref(); # Gives ids in analysis
#loop trough id's and perform the aggregation of the data
foreach my $id ( @$ids_aref ) {
    $logger->info( "Aggregation of $id" );
    $aggregate_obj->aggregate($id);
}

$logger->info( "Getting rid of the single quotes in the temp yaml passed to edgeR" );
my $temp_yaml_file = $param_obj->get_temp_yaml_file();
$logger->debug( $temp_yaml_file );
#need to get rid of single quotes
my $temp_yaml_fo = file($temp_yaml_file);
my @temp_yaml_lines = $temp_yaml_fo->slurp( chomp=>1 );
my $temp_dir = tempdir();
my ($tfh, $filename) = tempfile();
foreach my $line (@temp_yaml_lines) {
    if ( $line =~ qr/test:/ || $line =~ qr/ref_meta_cols:/ ) {
        $line =~ s/'//g;
    }
    print $tfh $line, "\n";
}
close($tfh);

#Run edgeR
$logger->info( "Running edgeR" );
my $r_source_dir = $param_obj->get_Rsource_dir();

# get the edgeR_model path which is the same as the current script
my $path = dirname(abs_path($0));
my $cmd = "Rscript --no-save --no-restore $path/edgeR_model.R params_file=\\\"$filename\\\" source_dir=\\\"$r_source_dir\\\"";

$logger->debug( $cmd );
system($cmd);

#Make the necessary objects
$logger->info( "Creation of a Decouple and a Da_Table object" );
my $decouple_obj = DecoupleDa->new( $param_obj );
my $da_table_obj = DaTable->new( { 'grp_order_aref' => $aggregate_obj->get_grp_order_aref(),
                                   'id_order_aref' => $justify_obj->get_ordered_ids_aref(),} );

#Decouple the edgeR data and then fill in a matrix
$logger->info( "Filling in the count matrix in Da_Table");
foreach my $genome ( @{$justify_obj->get_ordered_ids_aref()} ) {
    my $ids_w_da_counts_aref_aref = $decouple_obj->decouple($genome);
    $da_table_obj->set_genome($genome, $ids_w_da_counts_aref_aref);
}

#get the outfile path
my $out_path = $param_obj->get_out_dir();

#Look to filter and print or print the full object
$logger->info( "Looking for filter params" );
my $da_filter_obj;
$da_table_obj->print_full_da_table( "$out_path/full_da_tbl_$$.txt");
$param_obj->print_yaml_file("$out_path/full_da_tbl_$$.log");

eval{ $da_filter_obj = DaFilter->new( $da_table_obj );
      $param_obj->check_filter_params() };
if ( my $err = Exception::Class->caught() ) {
    print "WARNING: See error below only if you intended to print a filtered DaTable\n$err\nThe full DA table will still be printed out. If you'd like to try and filter use the ____ program which can filter using just a full DA table file";
}
else{
    $logger->info( "Filtering the Da_Table" );
    $da_filter_obj->filter_and_print( $param_obj, "$out_path/filt_da_tbl_$$.txt");
    $param_obj->print_yaml_file("$out_path/filt_da_tbl_$$.log");
    #Append the date at the end of the log file
    open my $OUT, ">>", "$out_path/filt_da_tbl_$$.log";
    my $time_string = localtime();
    print $OUT "\nDate and time this was run: $time_string";
    close($OUT);
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
use DaFilter;
use DaTable;
use DecoupleDa;
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
use Cwd qw(abs_path)
use File::Basename

=head1 INHERIT

    NA
    
=head1 SYNOPSIS

    perl edgeR_driver.pl -{y|x}
    
    -yml_file | y =>    This is a yaml parameter file that contains all the     parameters needed for the analysis
    
    -xml_file | x =>    This is a xml parameter file that contains all the parameters needed for the analysis
    
    ATTENTION: You must pass either a yaml or xml file. For more information about passing in filter parameters check DaFilter.pm documentation
    
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
    test                =>  Test Names      ->  String
    test_col_name       =>  Column Name     ->  String
    grp_meta_file       =>  File Path       ->  String
    Rsource_dir         =>  Directory Path  ->  String
    filter_params       =>  Hash Ref        ->  Array Ref( 0-100, True/False) ^(Optional)^
    filter params example:
    filter_params:
        f1:
        -   90
        -   'false/true'
        f-1:
        -   10
        -   'false/true'
        
    That will filter both abundance values of 1 and -1 in an and fashion.
    The first param under the abundance value to filter is the percentage. True means it will keep everything above the filter percentage, and false means all the data with a percentage of abundance values less than that will be kept
    
    
=head1 CONFIGURATION AND ENVIRONMENT

    Need to load R.

=head1 DEPENDENCIES

    R -> Version 3.3.1
    Param_handler
    Justify
    Aggregate
    DaFilter
    DaTable
    DecoupleDa
    
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

