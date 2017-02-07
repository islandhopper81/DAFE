#!/usr/bin/env perl

# Make a bbmap database with all the combined genomes
# cats all genome.fna files
# mergeds all gff files

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
use Config::Std;
use Data::Dumper;
use IPC::Cmd qw(can_run run);
use File::Basename;
use File::Temp qw(tempfile);
use BioUtils::FastaIO;

# Subroutines #
sub combine_genomes;
sub correct_genome;
sub correct_gff;
sub make_bbmap_db;
sub check_params;
sub _is_defined;

# Variables #
my ($config_file, $dafe_db_dir, $ref_names_file, $sample_names_file,
	$genome_file_exten, $gff_file_exten, $combined_db_name, $htseq_i,
    $lsf_threads, $lsf_mem, $lsf_queue, $lsf_out_file, $lsf_err_file,
    $lsf_job_name, 
    $help, $man);

my $options_okay = GetOptions (
    "config_file:s" => \$config_file,
    "dafe_db_dir:s" => \$dafe_db_dir,
    "ref_names_file:s" => \$ref_names_file,
    "sample_names_file:s" => \$sample_names_file,
    "genome_file_exten:s" => \$genome_file_exten,
    "gff_file_exten:s" => \$gff_file_exten,
	"combined_db_name:s" => \$combined_db_name,
	"htseq_i:s" => \$htseq_i,
    "lsf_threads:i" => \$lsf_threads,
    "lsf_mem:i" => \$lsf_mem,
    "lsf_queue:s" => \$lsf_queue,
    "lsf_out_file:s" => \$lsf_out_file,
    "lsf_err_file:s" => \$lsf_err_file,
    "lsf_job_name:s" => \$lsf_job_name,
    "help" => \$help,                   # flag
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
# combine the genome fasta files for all the genomes in ref_names_file
my $combined_fasta = combine_genomes($dafe_db_dir, $ref_names_file,
									 $combined_db_name);

# make the bbmap database
make_bbmap_db($combined_fasta, $combined_db_name, $dafe_db_dir);

########
# Subs #
########
sub combine_genomes() {
	my ($dafe_db_dir, $ref_names_file, $combined_db_name) = @_;
	
	# read in and store the genomes in the ref_names_file
	my @genomes = ();
	open my $IN, "<", $ref_names_file or
		$logger->logdie("Cannot open --ref_names_file: $ref_names_file");
	
	foreach my $line ( <$IN> ) {
		chomp $line;
		push @genomes, $line;
	}
	
	# for each genome in the @genome arr cat them into a single file
	my $genome_out_file = dirname($dafe_db_dir) . "/" . $combined_db_name . ".fna";
	my $gff_out_file = dirname($dafe_db_dir) . "/" . $combined_db_name . ".gff";
	#open my $OUT, "$out_file" or
	#	$logger->logdie("Cannot open combined genome fasta file: $out_file");
	
	my $genome_path;
	my $gff_path;
	foreach my $g ( @genomes ) {
		$genome_path = "$dafe_db_dir/$g/$g" . $genome_file_exten;
		$gff_path = "$dafe_db_dir/$g/$g" . $gff_file_exten;
		
		# get the corrected files
		$genome_path = correct_genome($genome_path, $g);
		$gff_path = correct_gff($gff_path, $g);
		
		# concatenate the files
		`cat $genome_path >> $genome_out_file`;
		`cat $gff_path >> $gff_out_file`;
	}
	
	return($genome_out_file);
}

sub correct_genome {
	my ($file, $genome_id) = @_;
	
	# correct the mistakes -- you know it is necessary at this point.
	$logger->info("Genome fasta file needs correction: $file");				
 
	# make an output file		
	my ($fh, $filename) = tempfile();		
	close($fh);		
	my $fasta_out = BioUtils::FastaIO->new({stream_type => '>', file => $filename});		

	# read the genome fasta file		
	my $fasta_in = BioUtils::FastaIO->new({stream_type => '<', file => $file});		
	
	while ( my $seq = $fasta_in->get_next_seq() ) {
		# set the ID as the header.  Note that the ID is defined as everything up
		# until the first space in the header string.
		$seq->set_header($seq->get_id());
				
		if ( $seq->get_id() !~ m/$genome_id-\S+/ ) {		
			$seq->set_header($genome_id . "-" . $seq->get_id());		
		}		
		
		# print		
		$fasta_out->write_seq($seq);		
	}
	
  	return($filename);
}

sub correct_gff {
	my ($file, $genome_id) = @_;
	
	# create a temp file for outputing the corrected gff file
	my ($fh, $filename) = tempfile();
	
	# open the gff file
	open my $GFF, "<", $file or
		$logger->warn("Cannot open gff: $file");
	
	my @vals = ();
	foreach my $line ( <$GFF> ) {
		chomp $line;
		
		#skip empty lines
		if ( $line =~ m/^$/ ) {
			next;
		}
		
		# skip comment lines
		# this used to be in the fix strand block.  To preserve the comment line
		# print it before going on to the next line.
		if ( $line =~ m/^#/ ) {
			print $fh "$line\n";
			next;
		}
		
		# get all the seperate fields in the line
		@vals = split("\t", $line);
		
		# Change the scaffold name to include the genome ID
		if ( $vals[0] !~ m/$genome_id/ ) {
			$vals[0] = $genome_id . "-" . $vals[0];
		}
		
		# Change the ID of the genes in the gff file to include the genome ID
		if ( $vals[8] !~ m/$htseq_i=$genome_id-\S+?/ ) {
			if ( $vals[8] =~ m/$htseq_i=(\S+?);/ ) {
				my $new_id = $htseq_i . "=" . $genome_id . "-" . $1;
				$vals[8] =~ s/$htseq_i=\S+?;/$new_id;/;
			}
		}
		
		print $fh (join("\t", @vals), "\n");
	}
	close($fh);
	close($GFF);
	
	return($filename);
}

sub make_bbmap_db {
	my ($combined_fasta, $combined_db_name, $dafe_db_dir) = @_;
	
	my $out_dir = dirname($dafe_db_dir);
	
	`bbmap.sh path="$out_dir/$combined_db_name/" ref=$combined_fasta`;
	
	return 1;
}

sub check_params {
    # check the config file
    if ( defined $config_file ) {
        # load the params in the config file
        load_config($config_file);
    }
    
    # check the dafe db dir
	if ( ! defined $dafe_db_dir ) {
		pod2usage(-message => "ERROR: required --dafe_db_dir not defined\n\n",
					-exitval => 2);
	}
    if ( ! -d $dafe_db_dir ) {
        pod2usage(-message => "ERROR: --dafe_db_dir ($dafe_db_dir) is not a directory\n\n",
                    -exitval => 2); 
    }
	if ( ! -r dirname($dafe_db_dir) ) {
		my $msg = "--dafe_db_dir must have write permissions in parent dir\n\n";
		pod2usage(-message => "ERROR: $msg",
					-exitval => 2);
	}
    
    # check ref_names_file
    if ( defined $ref_names_file and ! -e $ref_names_file ) {
        pod2usage(-message => "ERROR: --ref_names_file is an empty file\n\n",
					-exitval => 2); 
    }
    
    # check sample_names_file
    if ( defined $sample_names_file and ! -e $sample_names_file ) {
        pod2usage(-message => "ERROR: --sample_names_file is an empty file\n\n",
					-exitval => 2); 
    }
    
    # check genome_file_exten
    if ( ! defined $genome_file_exten ) {
        $genome_file_exten = ".fna";
    }
    
    # check gff_file_exten
    if ( ! defined $gff_file_exten ) {
        $gff_file_exten = ".gff";
    }
	
	# check combined_db_name
	if ( ! defined $combined_db_name ) {
		$combined_db_name = "all_genomes";
	}
	
	# check htseq_i
	if ( ! defined $htseq_i ) {
		$htseq_i = "ID";
	}
    
    # check lsf_threads
    if ( ! defined $lsf_threads ) {
        $lsf_threads = 2;
    }
    if ( $lsf_threads < 2 ) {
        pod2usage(-message => "ERROR: --lsf_threads must be >= 2\n\n",
					-exitval => 2); 
    }
    
    # check lsf_mem
    # there are many more checks I should do here
    if ( ! defined $lsf_mem ) {
        $lsf_mem = 4;
    }
    
    # check lsf_queue
    # there are many more checks I should do here
    if ( ! defined $lsf_queue ) {
        $lsf_queue = "week";
    }
    
    # check lsf_out_file
    if ( ! defined $lsf_out_file ) {
        $lsf_out_file = "lsf.out";
    }
    
    # check lsf_err_file
    if ( ! defined $lsf_err_file ) {
        $lsf_err_file = "lsf.err";
    }
    
    # check lsf_job_name
    if ( ! defined $lsf_job_name ) {
        $lsf_job_name = "DAFE_count";
    }
    
    # log the parameters
    $logger->info("Parameters");
    $logger->info("--config_file: $config_file") if (defined $config_file);
    $logger->info("--dafe_db_dir: $dafe_db_dir");
    $logger->info("--ref_names_file: $ref_names_file");
    $logger->info("--sample_names_file: $sample_names_file");
    $logger->info("--genome_file_exten: $genome_file_exten");
    $logger->info("--gff_file_exten: $gff_file_exten");
	$logger->info("--combined_db_name: $combined_db_name");
    $logger->info("--lsf_threads: $lsf_threads");
    $logger->info("--lsf_mem: $lsf_mem");
    $logger->info("--lsf_queue: $lsf_queue");
    $logger->info("--lsf_out_file: $lsf_out_file");
    $logger->info("--lsf_err_file: $lsf_err_file");
    $logger->info("--lsf_job_name: $lsf_job_name");
	
	return 1;
}

sub load_config {
    my ($file) = @_;
    
    my %params = ();
    read_config($file, %params);
    $logger->debug("Config Data:");
    $logger->debug(Dumper(%params));
    
    # NOTE: the parameters in the config file overide the parameters given on
    #       the command line except in the case of the skip_* parameters and
    #       the lsf* parameters
    
    $dafe_db_dir = _is_defined($params{''}{dafe_db_dir}, "dafe_db_dir");
    $ref_names_file = _is_defined($params{''}{ref_names_file}, "ref_names_file");
    $sample_names_file = _is_defined($params{''}{sample_names_file}, "sample_names_file");
    $genome_file_exten = _is_defined($params{''}{genome_file_exten}, "genome_file_exten", ".fna");
    $gff_file_exten = _is_defined($params{''}{gff_file_exten}, "gff_file_exten", ".gff");
	$combined_db_name = _is_defined($params{''}{combined_db_name}, "combined_db_name", "all_genomes");
    $lsf_threads = _is_defined($params{''}{lsf_threads}, "lsf_threads", $lsf_threads);
    $lsf_mem = _is_defined($params{''}{lsf_mem}, "lsf_mem", $lsf_mem);
    $lsf_queue = _is_defined($params{''}{lsf_queue}, "lsf_queue", $lsf_queue);
    $lsf_out_file = _is_defined($params{''}{lsf_out_file}, "lsf_out_file", $lsf_out_file);
    $lsf_err_file = _is_defined($params{''}{lsf_err_file}, "lsf_err_file", $lsf_err_file);
    $lsf_job_name = _is_defined($params{''}{lsf_job_name}, "lsf_job_name", $lsf_job_name);
    
    return 1;
}

sub check_env {
    # this subroutines checks for required executables
    # you may need to load certain modules to pass this check
    
    if ( ! can_run('bbmap.sh') ) {
        $logger->logdie("Cannot find required bbmap program!");
    }
    
    return 1;
}

sub _is_defined {
	my ($val, $name, $default) = @_;

	# set the default if needed
	if ( ! defined $val and defined $default ) {
		$val = $default;
	}

	# make sure it is defined
	if ( UtilSY::check_defined($val, $name) ) {
		return($val);
	}
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
