#!/usr/bin/env perl

# Runs the making of DA tables across multiple processors to speed up the process
# currently

use strict;
use warnings;
use Cwd;  # used in getting the current directory
use Getopt::Long;
use Pod::Usage;
use Carp;
use Readonly;
use version; our $VERSION = qv('0.0.1');
use Log::Log4perl qw(:easy);
use Log::Log4perl::CommandLine qw(:all);
use File::Temp qw(tempfile);
use POSIX qw(ceil);
use UtilSY qw(:all);
use Table;


# Subroutines #
sub split_files;
sub run_commands;
sub stall;
sub merge_results;
sub clean_up;
sub build_fastxIO;
sub get_fastx_seq_count;
sub check_params;

# set up the logging environment
my $logger = get_logger();

# Parallelization Variables #
my ($exe, $queue, $jobs, $out_dir, $keep_tmp, $job_name, $help, $man);

# Command Variable #
my ($dafe_out, $dafe_db, $features_file, $feature_type, $genomes_file, 
	$tags_file, $annote_file, $count_file, $sample_meta_file, $test_col, $t1, $t2,
	$out_file, $min_s, $min_c,
	$use_da_vec, $da_vec_file_name);

my $options_okay = GetOptions (
    "dafe_out:s" => \$dafe_out,
	"dafe_db:s" => \$dafe_db,
	"features_file:s" => \$features_file,
	"feature_type:s" => \$feature_type,
	"genomes_file:s" => \$genomes_file,
    "tags_file:s" => \$tags_file,
	"annote_file:s" => \$annote_file,
	"count_file:s" => \$count_file,
	"sample_meta_file:s" => \$sample_meta_file,
	"test_col:s" => \$test_col,
	"t1:s" => \$t1,
	"t2:s" => \$t2,
	"out_file:s" => \$out_file,
	"min_s:s" => \$min_s,
	"min_c:s" => \$min_c,
	"use_da_vec:s" => \$use_da_vec,
	"da_vec_file_name:s" => \$da_vec_file_name,
    "exe:s" => \$exe,
    "queue:s" => \$queue,
    "jobs:i" => \$jobs,
    "out_dir:s" => \$out_dir,
    "keep_tmp" => \$keep_tmp,          # flag
    "job_name:s" => \$job_name,
    "help" => \$help,                  # flag
    "man" => \$man,                     # flag (print full man page)
);

# check for input errors
if ( $help ) { pod2usage(0) }
if ( $man ) { pod2usage(-verbose => 3) }
check_params();


### MAIN ###
my @time = localtime time;
$logger->info("Time: " . $time[2] . ":" . $time[1] . ":" . $time[0]);
$logger->info("Day: " . $time[3] . ":" . $time[4] . ":" . ($time[5]+1900));

# split the $genome_file into multiple files
$logger->info("Spliting genome files");
my $genomes_aref = split_files($genomes_file);

system("ls -alh /tmp/");

# run all the commands
$logger->info("Running Commands");
my $out_files_aref = run_commands($genomes_aref);
my $file = $genomes_aref->[0];
print("print file: $file\n");
system("cat $file");

# stall until all the jobs are finished
$logger->info("Stalling");
stall();

# merge the results
$logger->info("Merging Results");
merge_results($out_files_aref);

# clean up the tmp files
clean_up($genomes_aref, $out_files_aref);

@time = localtime time;
$logger->info("Time: " . $time[2] . ":" . $time[1] . ":" . $time[0]);
$logger->info("Day: " . $time[3] . ":" . $time[4] . ":" . ($time[5]+1900));






# Sub routines #
sub split_files {
    my ($genomes_file) = @_;
    
    #### splits the genomes file into $splits different files
    
    # read in the genomes file
    my $g_aref = load_lines($genomes_file);
    
    # calculate the number of lines that should be in each subset file
    my $subset_size = ceil( (scalar @{$g_aref}) / $jobs );
    
    # create a genomes file for each job
    my @genome_files = ();
    for ( my $i = 0; $i < $jobs; $i++ ) {
        # create a temporary file
        #my ($fh, $filename) = tempfile();
        my $filename = ".tmp" . $i;
        open my $fh, ">", $filename or
            $logger->logdie("Cannot print to tmp file $filename");
        $logger->debug("Tmp genome file: $filename");
        
        # print the subset
        my @subset = splice(@{$g_aref}, 0, $subset_size);
        print "print to file\n";
        print $fh aref_to_str(\@subset);
        print "closing file\n";
        close($fh);
        
        # save the file name
        push @genome_files, $filename;
    }
    
    # split the files based on the number of jobs needed
    
    return \@genome_files;
}

sub run_commands {
    my ($genomes_aref) = @_;
    
    # make an output files for each genome input file
    my @out_files = ();
    my $i = 0;
    foreach my $g_file ( @{$genomes_aref} ) {
        my $filename = ".tmp.merge" . $i;
        open my $fh, ">", $filename or
            $logger->logdie("Cannot open tmp file: $filename");
        push @out_files, $filename;
        $i++;
    }
    
    for ( my $i = 0; $i < scalar @out_files; $i++ ) {
        # suppress warnings of uninitialized values because not all
        # parameters might be used
        no warnings 'uninitialized';
        
        my $com = "bsub -q $queue -o mk_da_tbl.out -e mk_da_tbl.err -J $job_name ";
        $com .= "$exe ";
        $com .= "--genomes_file " . $genomes_aref->[$i] . " ";
        $com .= "--out_file " . $out_files[$i] . " ";
        $com .= "--dafe_out $dafe_out " if defined $dafe_out;
        $com .= "--dafe_db $dafe_db " if defined $dafe_db;
        $com .= "--features_file $features_file " if defined $features_file;
        $com .= "--feature_type $feature_type " if defined $feature_type; 
        $com .= "--tags_file $tags_file " if defined $tags_file;
        $com .= "--annote_file $annote_file " if defined $annote_file;
        $com .= "--count_file $count_file " if defined $count_file;
        $com .= "--sample_meta_file $sample_meta_file " if defined $sample_meta_file;
        $com .= "--test_col $test_col " if defined $test_col;
        $com .= "--t1 $t1 " if defined $t1;
        $com .= "--t2 $t2 " if defined $t2;
        $com .= "--min_s $min_s " if defined $min_s;
        $com .= "--min_c $min_c " if defined $min_c;
        $com .= "--use_da_vec $use_da_vec " if defined $use_da_vec;
        $com .= "--da_vec_file_name $da_vec_file_name " if defined $da_vec_file_name;
        $com .= "--debug";
        
        $logger->info("Running Command: $com");
        system($com);
    }
    
    return(\@out_files);
}

sub stall {
    my $com = 'bjobs > ' . "$out_dir/ACTIVE_JOBS" . '
                while [ -s ' . "$out_dir/ACTIVE_JOBS" . ' ]
                do
                    sleep 10
                    bjobs -J ' . $job_name . '\* > ' . "$out_dir/ACTIVE_JOBS" . '
                    done
                    rm ' . "$out_dir/ACTIVE_JOBS";
    
    `$com`;
}

sub merge_results {
    my ($out_files_aref) = @_;
    
    my $final_tbl = Table->new();
    my $first = 1;
    my $tmp_tbl = Table->new();
    
    foreach my $tmp_file ( @{$out_files_aref} ) {
        $logger->debug("cbind tmp file: $tmp_file");
        $tmp_tbl->load_from_file($tmp_file);
        
        if ( $first == 1 ) {
            $final_tbl = $tmp_tbl->copy();
            $first = 0;
        }
        else {
            $final_tbl->cbind($tmp_tbl);
        }
    }
    
    # output file final talbe
    $final_tbl->save($out_file);
}

sub clean_up {
    my ($genomes_aref, $out_files_aref) = @_;
    
    foreach my $file ( @{$genomes_aref} ) {
        system("rm $file");
    }
    
    foreach my $file ( @{$out_files_aref} ) {
        system("rm $file");
    }
}

sub check_params {
	# check for required variables
	if ( ! defined $exe) { 
		pod2usage(-message => "ERROR: required --exe not defined\n\n",
					-exitval => 2); 
	}
	if ( ! defined $out_file ) {
		pod2usage(-message => "ERROR: required --out_file not defined\n\n",
					-exitval => 2);
	}
	if ( ! defined $dafe_out ) { 
		pod2usage(-message => "ERROR: required --dafe_out not defined\n\n",
					-exitval => 2); 
	}
	if ( ! defined $genomes_file ) {
		pod2usage(-message => "ERROR: required --genomes_file not defined\n\n",
					-exitval => 2);
	}
    
	## make sure required files are non-empty
   if ( defined $genomes_file and ! -e $genomes_file ) { 
		pod2usage(-message => "ERROR: --genome_file $genomes_file is an empty file\n\n",
					-exitval => 2);
	}
    
    # set some defualts
    if ( ! defined $queue ) { $queue = 'week'; }
    if ( ! defined $jobs ) { $jobs = 1; }  # serial
    if ( ! defined $out_dir ) { $out_dir = getcwd(); }
    if ( ! defined $keep_tmp ) { $keep_tmp = 0; }
    if ( ! defined $job_name ) { $job_name = "JOB"; }

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

<NAME> - <DESC>


=head1 VERSION

This documentation refers to <NAME> version 0.0.1


=head1 SYNOPSIS

    <NAME>  --queue $q
            --jobs $j
            --out_dir my_output
            --keep_tmp
            --command "echo helloworld"
            [--help]

    --queue  = The LSF queue to submit the jobs to
    --jobs   = The number of jobs to submit
    --out_dir = the output directory
    --keep_tmp = A flag to set when you want to keep the tmp files
    --command   = The command to run for each job
    --help  = Prints USAGE statement


=head1 ARGUMENTS
    
=head2 --queue

The LSF queue to submit the jobs to
    
=head2 --jobs

The number of jobs to submit

=head2 --out_dir

The output directory

=head2 --keep_tmp

A flag to set when you want to keep the tmp files

=head2 --command

The command to run for each job

=head2 [--help]
    
An optional parameter to print a usage statement.
    

=head1 DESCRIPTION

Runs process_mate_pairs.pl in parallel.


=head1 CONFIGURATION AND ENVIRONMENT
    
No special configurations or environment variables needed
    
    
=head1 DEPENDANCIES

Getopt::Long
Pod::Usage
Carp
Cwd
BioUtils::FastqIO
BioUtils::FastaIO



=head1 AUTHOR

Scott Yourstone     scott.yourstone81@gmail.com
    
    
=head1 LICENCE AND COPYRIGHT

Copyright (c) 2013, Scott Yourstone
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
