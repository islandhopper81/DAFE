#!/usr/bin/env perl

# Corrects the names in a genome fasta file.  Makes sure the name has no spaces

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Carp;
use Readonly;
use version; our $VERSION = qv('0.0.1');
use Log::Log4perl qw(:easy);
use Log::Log4perl::CommandLine qw(:all);
use File::Temp qw(tempfile);
use BioUtils::FastaIO;

# Subroutines #
sub check_params;
sub _is_defined;

# Variables #
my ($genome_fasta, $overwrite, $help, $man);

my $options_okay = GetOptions (
    "genome_fasta:s" => \$genome_fasta,
    "overwrite" => \$overwrite,
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
# make an output file
my ($fh, $filename) = tempfile();
close($fh);
my $fasta_out = BioUtils::FastaIO->new({stream_type => '>', file => $filename});

# read the genome fasta file
my $fasta_in = BioUtils::FastaIO->new({stream_type => '<', file => $genome_fasta});

while ( my $seq = $fasta_in->get_next_seq() ) {
    # set the ID as the header.  Note that the ID is defined as everything up
    # until the first space in the header string.
    $seq->set_header($seq->get_id());
    
    # print
    $fasta_out->write_seq($seq);
}

# now move the tmp file to the original genome fasta
# NOTE: this overwrite the original fasta file
$logger->info("Moving tmp file ($filename) to overwite genome fasta ($genome_fasta)");
`mv $filename $genome_fasta`;





########
# Subs #
########
sub check_params {
	# check for required variables
	if ( ! defined $genome_fasta) { 
		pod2usage(-message => "ERROR: required --genome_fasta not defined\n\n",
					-exitval => 2); 
	}

	# make sure required files are non-empty
	if ( ! -e $genome_fasta ) { 
		pod2usage(-message => "ERROR: --genome_fasta $genome_fasta is an empty file\n\n",
					-exitval => 2);
	}
    
    # check the overwrite parameter
    if ( ! _is_defined($overwrite) ) {
        pod2usage(-message => "ERROR: --overwrite flag must be given.  See docs!\n\n",
					-exitval => 2);
    }

    $logger->info("--genome_fasta: $genome_fasta");
	
	return 1;
}

sub _is_defined {
    my ($val, $default) = @_;
    
    if ( defined $val and $val eq "" ) {
        # this is basically undefined
        $val = undef;
    }
    
    if ( defined $val ) {
        return $val;
    }
    elsif ( defined $default ) {
        return $default;
    }
    else {
        return undef;
    }
}


__END__

# POD

=head1 NAME

DAFE_db_correct_genome_fasta.pl - Removes everything after the first space in
the headers of a genome fasta file.


=head1 VERSION

This documentation refers to version 0.0.1


=head1 SYNOPSIS

    DAFE_db_correct_genome_fasta.pl
        --genome_fasta my_genome.fasta
        --overwrite
        
        [--help]
        [--man]
        [--debug]
        [--verbose]
        [--quiet]
        [--logfile logfile.log]

    --genome_fasta  Path to an input genome fasta file
    --overwrite     Flag that warns user that the original genome fasta is overwritten!
    --help | -h     Prints USAGE statement
    --man           Prints the man page
    --debug	        Prints Log4perl DEBUG+ messages
    --verbose       Prints Log4perl INFO+ messages
    --quiet	        Suppress printing ERROR+ Log4perl messages
    --logfile       File to save Log4perl messages


=head1 ARGUMENTS
    
=head2 --genome_fasta

Path to an input genome fasta file
    
=head2 --overwrite

Flag that warns user that the original genome fasta file is overwritten with the
new fasta file with corrected names.  The orignal fasta cannot be recovered!  So
only run this script if you can afford to lose the header inforamtion after the
first space!
 
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

Note that the original genome fasta file is overwritten with this script!  It
cannot be recovered.

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
File::Temp qw(tempfile)
BioUtils::FastaIO


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
