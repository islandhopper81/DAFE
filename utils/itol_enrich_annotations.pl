#!/usr/bin/env perl

# makes the ITOL annotations file for enrichments

# not in this case the files have no headers lines.

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Carp;
use Readonly;
use version; our $VERSION = qv('0.0.1');
use Log::Log4perl qw(:easy);
use Log::Log4perl::CommandLine qw(:all);
use UtilSY qw(:all);
use Table;

# Subroutines #
sub check_params;

# Variables #
my ($da_file, $feat_name, $up, $dn, $out_file, $help, $man);

my $options_okay = GetOptions (
    "da_file:s" => \$da_file,
	"feat_name:s" => \$feat_name,
	"up:s" => \$up,
	"dn:s" => \$dn, 
    "out_file:s" => \$out_file,
    "help|h" => \$help,                  # flag
    "man" => \$man,                     # flag (print full man page)
);

# set up the logging environment
my $logger = get_logger();

# Global Variables #
Readonly::Scalar my $ENR => 1;
Readonly::Scalar my $NS => 0;
Readonly::Scalar my $DEP => -1; 
Readonly::Scalar my $LOW => -2; 
Readonly::Scalar my $UMS => -3; 
Readonly::Scalar my $ABS => -4; 

# check for input errors
if ( $help ) { pod2usage(0) }
if ( $man ) { pod2usage(-verbose => 3) }
check_params();


########
# MAIN #
########
# read in the reference annotation table
my $ref_tbl = Table->new();
$ref_tbl->load_from_file($da_file, "\t", "F"); # no col headers

# open the output file
open my $OUT, ">", $out_file or
	$logger->logdie("Cannot open --out_file $out_file");

# print the header info
print $OUT "DATASET_COLORSTRIP\n";
print $OUT "SEPARATOR TAB\n";
print $OUT "DATASET_LABEL\t$feat_name\n";
print $OUT "COLOR\t#ff0000\n";


# set the colors
my $ENR_col = "rgb(251,0,5)";
my $DEP_col = "rgb(0,23,252)";
my $NS_col = "rgb(255,253,29)";
my $LOW_col = "rgb(255,255,235)";
my $UMS_col = "rgb(201,201,201)";
my $ABS_col = "rgb(153,153,153)";

# print the legend information
print $OUT "LEGEND_TITLE\t$feat_name\n";
print $OUT "LEGEND_SHAPES\t1\t1\t1\t1\t1\t1\n";
print $OUT "LEGEND_COLORS\t$ENR_col\t$DEP_col\t$NS_col\t$LOW_col\t$UMS_col\t$ABS_col\n";
print $OUT "LEGEND_LABELS\t$up Enriched\t$dn Enriched\tNot Enriched\tLow Abundance\tUnmeasurable\tAbsent\n";


# print the data for each genome
print $OUT "\nDATA\n";

# each row is a genome
foreach my $g ( @{$ref_tbl->get_row_names()} ) { 

	# remove the "X" that might be on the front of the genome name
	my $name;
	if ( $g =~ m/^X(\w+)/ ) {
		$name = $1;
	}
	else { 
		$name = $g
	}
    print $OUT "$name\t";
	
	# remember that the columns names are integers starting at 0
	# because there are no explicit names in the files
	my $val = $ref_tbl->get_value_at($g, "0");
	
	if ( $val eq $ENR ) {
		print $OUT $ENR_col . "\n";
	}
	elsif ( $val eq $DEP ) {
		print $OUT $DEP_col . "\n";
	}
	elsif ( $val eq $NS ) {
		print $OUT $NS_col . "\n";
	}
	elsif ( $val eq $LOW ) {
		print $OUT $LOW_col . "\n";
	}
	elsif ( $val eq $UMS ) {
		print $OUT $UMS_col . "\n";
	}
	elsif ( $val eq $ABS ) {
		print $OUT $ABS_col . "\n";
	}
	else {
		$logger->warn("Unrecognized DA value ($val) at genome: $g");
	}
}


close($OUT);


########
# Subs #
########
sub check_params {
	# check for required variables
	if ( ! defined $da_file ) { 
		pod2usage(-message => "ERROR: required --da_file not defined\n\n",
					-exitval => 2); 
	}
	if ( ! defined $feat_name ) { 
		pod2usage(-message => "ERROR: required --feat_name not defined\n\n",
					-exitval => 2); 
	}
	if ( ! defined $up ) { 
		pod2usage(-message => "ERROR: required --up not defined\n\n",
					-exitval => 2); 
	}
	if ( ! defined $dn ) { 
		pod2usage(-message => "ERROR: required --dn not defined\n\n",
					-exitval => 2); 
	}
	if ( ! defined $out_file ) {
		pod2usage(-message => "ERROR: required --out_file not defined\n\n",
					-exitval => 2);
	}

	# make sure required files are non-empty
	if ( defined $da_file and ! -e $da_file ) { 
		pod2usage(-message => "ERROR: --da_file $da_file is an empty file\n\n",
					-exitval => 2);
	}

	# make sure required directories exist
	#if ( ! -d $dir ) { 
	#	pod2usage(-message => "ERROR: --dir is not a directory\n\n",
	#				-exitval => 2); 
	#}
	
	return 1;
}


__END__

# POD

=head1 NAME

itol_enrich_annotations.pl - creates an itol file for a single feature (ie COG)


=head1 VERSION

This documentation refers to version 0.0.1

=head1 SYNOPSIS

    itol_enrich_annotations.pl
        --da_file COG0001_da.txt
        --feat_name COG0001
        --up RZ
        --dn BK
        --out_file COG0001_itol.txt
        
        [--help]
        [--man]
        [--debug]
        [--verbose]
        [--quiet]
        [--logfile logfile.log]

    --da_file     Path to an input file of da calls for each genome for this feature
    --feat_name      Name of the feature to be use for setting the legends
    --up             Name of up enrichments for setting the legend
    --dn             Name of down enrichments for setting the legend
    --out_file       Path to output file
    --help | -h     Prints USAGE statement
    --man           Prints the man page
    --debug	        Prints Log4perl DEBUG+ messages
    --verbose       Prints Log4perl INFO+ messages
    --quiet	        Suppress printing ERROR+ Log4perl messages
    --logfile       File to save Log4perl messages


=head1 ARGUMENTS
    
=head2 --file | -f

Path to an input file
    
=head2 --var | -v

Path to an input variable   
 
=head2 [--help | -h]
    
An optional parameter to print a usage statement.

=head2 [--man]

An optional parameter to print he entire man page (i.e. all documentation)

=head2 [--debug]

Prints Log4perl DEBUG+ messages.  The plus here means it prints DEBUG
level and greater messages.

=head2 [--verbose]

Prints Log4perl INFO+ messages.  The plus here means it prints INFO level
and greater messages.

=head2 [--quiet]

Suppresses print ERROR+ Log4perl messages.  The plus here means it suppresses
ERROR level and greater messages that are automatically printed.

=head2 [--logfile]

File to save Log4perl messages.  Note that messages will also be printed to
STDERR.
    

=head1 DESCRIPTION

[FULL DESCRIPTION]

=head1 CONFIGURATION AND ENVIRONMENT
    
No special configurations or environment variables needed
    
    
=head1 DEPENDANCIES

version
Getopt::Long
Pod::Usage
Carp
Readonly
version
Log::Log4perl qw(:easy)
Log::Log4perl::CommandLine qw(:all)
UtilSY qw(:all)

=head1 AUTHOR

Scott Yourstone     scott.yourstone81@gmail.com
    
    
=head1 LICENCE AND COPYRIGHT

Copyright (c) 2015, Scott Yourstone
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
