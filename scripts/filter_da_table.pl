#! usr/bin/env perl

use strict;
use warnings;
use lib "../lib";
use Param_handler;
use DaTable;
use DaFilter;
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
my $full_da_table_file; #place to put the ordered ids and meta_file

# Read in variables from the command line
GetOptions ('man'  => \$man,
            'help' => \$help,
            'xml_file|x=s' => \$xml_file,
            'yml_file|y=s' => \$yml_file,
            'da_file=s'    => \$full_da_table_file,
            ) || die("There was an error in the command line arguements\n");

# Use Pod usage for the Manual and Help pages
if ( $help ) { pod2usage(0) }
if ( $man )  {pod2usage(-verbose => 3) }

#kill program if neither a xml or txt file is passed
if (!defined $xml_file && !defined $yml_file) {
    croak("Must pass in a text file with newline characters, or an xml file.
          Do this by specifying the file with either -xml_file or -txt_file");
}

# Subroutines that will occur in the main program

# MAIN #

my $param_obj;
#Create a Param_handler object and check the Filter Param.
if (defined $xml_file) {
    $param_obj = Param_handler->new( { xml_file => $xml_file } );
}
elsif (defined $yml_file) {
    #my $param_href = find_param_values_from_txt_file($yml_file);
    $param_obj = Param_handler->new( { yml_file => $yml_file } );
}
#Just Check the filter params
$param_obj->check_filter_params();

#Create DaTable from file and then Filter obj
my $da_table_obj = DaTable->new( {'da_file' => $full_da_table_file} );
my $da_filter_obj = DaFilter->new( $da_table_obj );

#get wanted file names
my $outfile = $full_da_table_file;
$outfile =~ s/full/filtered/;
$outfile =~ s/\.txt//;
$outfile .= "_$$.txt";

#perform filter and log printing
$da_filter_obj->filter_and_print($param_obj, $outfile);

#change outfile to log file
$outfile =~ s/\.txt/\.log/;
$param_obj->print_yaml_file( $outfile );

__END__

=head1 filter_da_table
    
This program reads in a DA File and will filter the file based on the specifications in the passed in yaml file. The params must be set as...
filter_params:
    f#:
    - percentage
    - true/false
    f#:
    - percentage
    - true/false
      
You can filter multiple things at once. A log file will be created with the parametes associated with the filtered file and will be printed in the same directory as the filtered_da_table and the original da_table file. 

=head1 VERSION

This documentation refers to filter_da_table.pl 0.0.1

=head1 INCLUDED MODULES

use Param_handler;
use DaFilter;
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
use lib "/proj/cdjones_lab/ncolaian/MyX-Generic-v0.0.2/lib/MyX";
use MyX::Generic;
use YAML::XS qw(LoadFile);

=head1 INHERIT

    NA
    
=head1 SYNOPSIS

    perl edgeR_driver.pl -{y|x} -da_file
    
    -yml_file | y =>    This is a yaml parameter file that contains all the     parameters needed for the analysis
    
    -xml_file | x =>    This is a xml parameter file that contains all the parameters needed for the analysis
    
    ATTENTION: You must pass either a yaml or xml file.
    
    -da_file      =>    This is the da file that will be filtered. It will be
                        unaltered.
    
=head1 PARAMETERS
    
    

=head1 CONFIGURATION AND ENVIRONMENT

    NA

=head1 DEPENDENCIES

    Param_handler
    DaFilter
    DaTable
    
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
