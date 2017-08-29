#!/usr/bin/env perl

# makes the ITOL dataset disignating the PA/NPA status of the references

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
sub _is_defined;

# Variables #
my ($meta_file, $out_file, $help, $man);

my $options_okay = GetOptions (
    "meta_file:s" => \$meta_file,
	"out_file:s" => \$out_file,
    "help|h" => \$help,                  # flag
    "man" => \$man,                     # flag (print full man page)
);

# set up the logging environment
my $logger = get_logger();

# check for input errors
if ( $help ) { pod2usage(0) }
if ( $man ) { pod2usage(-verbose => 3) }
check_params();


########
# MAIN #
########
# read in the metadata
my $meta = Table->new();
$meta->load_from_file($meta_file);

my %colors = ("PA" => "#009900",
			  "NPA" => "#6600cc",
			  "Soil" => "#663300");

my @col_arr = ();

foreach my $id ( @{$meta->get_row_names()} ) {
	$logger->debug("ID: $id\n");
	my $label = $meta->get_value_at($id, "Label");
	$logger->debug("Label: $label\n");
	if ( defined $colors{$label} ) {
		push @col_arr, $colors{$label};
	}
	else {
		$logger->logdie("Cannot find label ($label) for genome ($id)\n");	
	}
}

# print the output
open my $OUT, ">", "$out_file" or
	$logger->logdie("Cannot open out_file ($out_file)\n");

print $OUT "DATASET_COLORSTRIP\n";
print $OUT "SEPARATOR TAB\n";
print $OUT "DATASET_LABEL\tPA/NPA\n";
print $OUT "COLOR\t#FFFFFF\n";
print $OUT "COLOR_BRANCHES\t0\n";
print $OUT "LEGEND_TITLE\tDataset_legend\n";
print $OUT "LEGEND_SHAPES\t1\t1\t1\n";
print $OUT "LEGEND_COLORS\t#009900\t#6600CC\ti#663300\n";
print $OUT "LEGEND_LABELS\tPA\tNPA\tSoil\n";
print $OUT "DATA\n";

my $i = 0;
foreach my $id ( @{$meta->get_row_names() }) {
	print $OUT $id, "\t";
	print $OUT $col_arr[$i], "\n";
	$i++;
}

close($OUT);

########
# Subs #
########
sub check_params {
	# check for required variables
	if ( ! defined $meta_file) { 
		pod2usage(-message => "ERROR: required --meta_file not defined\n\n",
					-exitval => 2); 
	}
	if ( ! defined $out_file ) {
		pod2usage(-message => "ERROR: required --out_file not defined\n\n",
					-exitval => 2);
	}

	# make sure required files are non-empty
	if ( defined $meta_file and ! -e $meta_file ) { 
		pod2usage(-message => "ERROR: --meta_file $meta_file is an empty file\n\n",
					-exitval => 2);
	}

	
	return 1;
}


__END__

# POD

=head1 NAME

[NAME].pl - [DESCRIPTION]


=head1 VERSION

This documentation refers to version 0.0.1


=head1 SYNOPSIS

    [NAME].pl
        -f my_file.txt
        -v 10
        
        [--help]
        [--man]
        [--debug]
        [--verbose]
        [--quiet]
        [--logfile logfile.log]

    --file | -f     Path to an input file
    --var | -v      Path to an input variable
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
