#!/usr/bin/env perl

# Runs read mapping, htseq-count, and bam filtering on the UNC cluster

use strict;
use warnings;

# record the start time
BEGIN { our $start_time = time(); }

use Getopt::Long;
use Pod::Usage;
use Carp;
use Readonly;
use version; our $VERSION = qv('0.0.1');
use Log::Log4perl qw(:easy);
use Log::Log4perl::CommandLine qw(:all);
use Config::Std;
use IPC::Cmd qw(can_run run);
use Data::Dumper;
use List::Util qw(min);
use FindBin qw($Bin);
use DAFE::Count::OutputRestruct;

# Subroutines #
sub check_params;

# Variables #
my ($config_file, $dafe_db_dir, $combined_db_name, $reads_dir, $out_dir,
    $restruct_out_dir, $ref_names_file, $sample_names_file, $perc_ids,
    $perc_ids_aref, $min_perc_id, $read_file_exten, $genome_file_exten,
    $gff_file_exten, $bam_file_prefix, $cov_stat_file_name, $htseq_file_prefix,
    $htseq_i, $make_count_tbls_exe, $pairs_file, $skip_mapping, $skip_filtering,
    $skip_htseq, $keep_bam_files, $lsf_threads, $lsf_mem, $lsf_queue,
    $lsf_out_file, $lsf_err_file, $lsf_job_name, $queue_batch_count,
    $node_batch_count, $queue_max, $runtime_max, $jobs_file,
    $help, $man);

# NOTE: perc_ids_aref and min_perc_id are not a parameters passed in via the
#       command line.  They are populated when the parameters are checked using
#       $perc_ids variable.

# NOTE: Wehn $perc_ids_aref is populated I remove the smallest value and assign
#       that value to $min_perc_id.  However, I add $min_perc_id back to the
#       $perc_ids_aref in later in the get_commands() subroutines.

# NOTE: $jobs_file is a hidden parameter.  The users doesn't see it when the
#       docs are printed.  And the user should be using it.  It is only used
#       when this script is required to restart because it reached --runtime_max

my $options_okay = GetOptions (
    "config_file:s" => \$config_file,
    "dafe_db_dir:s" => \$dafe_db_dir,
    "combined_db_name:s" => \$combined_db_name,
    "read_dir:s" => \$reads_dir,
    "out_dir:s" => \$out_dir,
    "restruct_out_dir:s" => \$restruct_out_dir,
    "ref_names_file:s" => \$ref_names_file,
    "sample_names_file:s" => \$sample_names_file,
    "perc_ids:s" => \$perc_ids,
    "read_file_exten:s" => \$read_file_exten,
    "genome_file_exten:s" => \$genome_file_exten,
    "gff_file_exten:s" => \$gff_file_exten,
    "bam_file_prefix:s" => \$bam_file_prefix,
    "cov_stat_file_name:s" => \$cov_stat_file_name,
    "htseq_file_prefix:s" => \$htseq_file_prefix,
    "htseq_i:s" => \$htseq_i,
    "make_count_tbls_exe:s" => \$make_count_tbls_exe,
    "pairs_file:s" => \$pairs_file,
    "skip_mapping" => \$skip_mapping,
    "skip_filtering" => \$skip_filtering,
    "skip_htseq" => \$skip_htseq,
    "keep_bam_files" => \$keep_bam_files,
    "lsf_threads:i" => \$lsf_threads,
    "lsf_mem:i" => \$lsf_mem,
    "lsf_queue:s" => \$lsf_queue,
    "lsf_out_file:s" => \$lsf_out_file,
    "lsf_err_file:s" => \$lsf_err_file,
    "lsf_job_name:s" => \$lsf_job_name,
    "queue_batch_count:i" => \$queue_batch_count,
    "node_batch_count:i" => \$node_batch_count,
    "queue_max:i" => \$queue_max,
    "runtime_max:i" => \$runtime_max,
    "jobs_file:s" => \$jobs_file,       # a hidden parameter
    "help" => \$help,                   # flag
    "man" => \$man,                     # flag (print full man page)
);

# set up the logging environment
my $logger = get_logger();

# check for input errors
if ( $help ) { pod2usage(2) }
if ( $man ) { pod2usage(-verbose => 3) }
check_params();
check_env();

########
# MAIN #
########
# get all the jobs
my $jobs_aref;
if ( defined $jobs_file ) {
    $jobs_aref = load_jobs($jobs_file);
}
else {
    $jobs_aref = get_jobs();
}

# submit the commands
while ( @{$jobs_aref} ) {
    # check the queue to see if it is full
    if ( full_queue($queue_max) ) {
        $logger->debug("Full queue");
        sleep 100;
        
        # check the current runtime to see if this script needs to restart
        check_runtime(time() - our $start_time);
    }
    else {
        submit_batch($jobs_aref, $queue_batch_count);
        
        $logger->info("Number of remaining jobs: " . scalar(@{$jobs_aref}));
    }
}

# stall until all $lsf_job_name jobs are complete
while ( stall($lsf_job_name) ) {
    sleep 100;
    
    # check the current runtime to see if this script needs to restart
    check_runtime(time() - our $start_time);
}

# check the output for error.  correct if possible
check_output_files();

# reformat the output to look like the old output structure with a directory
# for each genome and inside that directory and directory for each sample with
# count files
reformat_output_structure();

# make the gene count tables for each genome
make_gene_count_tbls();


$logger->info("COMPLETE!!");







########
# Subs #
########
sub reformat_output_structure {
    # reformat the output to look like the old output structure with a directory
    # for each genome and inside that directory and directory for each sample with
    # count files
    
    # get the list of samples and references (ie genomes)
    my ($ref_aref, $sample_aref);
    if ( defined $pairs_file and -s $pairs_file ) {
        ($ref_aref, $sample_aref) = get_pairs($pairs_file);
    }
    else {
        # this is actually the standard way where everything in the
        # ref_names_file and sample_names_file is ran
        $ref_aref = get_names($ref_names_file);
        $sample_aref = get_names($sample_names_file);
    }
    
    # get a list of all the percent identies
    # NOTE: $perc_ids_aref does not include the smallest value
    #       (ie $min_perc_id) so I make a temporary aref and add back
    #       $min_perc_id.
    my $tmp_perc_ids_aref = [($min_perc_id), @{$perc_ids_aref}];
    
    # make the OutputRestruct object
    my %args = ("out_dir" => $out_dir,
                "new_out_dir" => $restruct_out_dir,
                "ref_aref" => $ref_aref,
                "sample_aref" => $sample_aref,
                "perc_ids_aref" => $tmp_perc_ids_aref,
                "htseq_file_prefix" => $htseq_file_prefix);
    my $output_restruct = DAFE::Count::OutputRestruct->new(\%args);
    $output_restruct->restructure();
    
    return 1;
}

sub make_gene_count_tbls {
    # this subroutine makes a table for each reference (ie genome) of the
    # counts for each genomic feature (ie gene) across each of the samples
    # NOTE: this operation is done on all the files in the out_dir not just
    #       the ones that are specified in the ref_names_file and
    #       sample_names_file.
    
    $logger->info("Making gene count tables");
    
    # NOTE: $perc_ids_aref does not include the smallest value
    #       (ie $min_perc_id) so I make a temporary aref and add back
    #       $min_perc_id.
    my $tmp_perc_ids_aref = [($min_perc_id), @{$perc_ids_aref}];

    
    # build a loop for all the different percent identities
    foreach my $perc_id ( @{$tmp_perc_ids_aref} ) {
        my $command = "bash $make_count_tbls_exe ";
        $command .= "-b $restruct_out_dir ";
        $command .= "-i $ref_names_file ";
        $command .= "-s $sample_names_file ";
        $command .= "-c $htseq_file_prefix\_id" . $perc_id . ".txt ";
        $command .= "-o gene_counts_id" . $perc_id . ".txt";
        
        $logger->info("Running make_gene_count_tbls command: $command");
        `$command`;
    }
    
    return 1;
}

sub check_runtime {
    my ($time) = @_;
    
    if ( $time > $runtime_max ) {
        restart();  # the program will exit and restart at this point
    }
    
    return 1;
}

sub restart {
    # no parameters
    
    $logger->info("Restarting");
    
    my $config_file = save_config();
    $logger->debug("Restart config file: $config_file");
    
    my $jobs_file = save_jobs();
    $logger->debug("Restart jobs file: $jobs_file");
    
    # resubmit the master job (ie this script)
    my $command = "bsub -q week -o lsf.out -e lsf.err -J restart ";
    $command .= "perl ~/scripts/comparitive_metaG/DAFE_count.pl ";
    $command .= "--config_file $config_file ";
    $command .= "--jobs_file $jobs_file ";
    $command .= "--debug ";
    $logger->debug("Restart command: $command");
    
    `$command`;
    
    exit 0;  # die but with success.  :)
}

sub save_config {
    # no parameters
    
    my $file = "$out_dir/.restart.conf";
    open my $CON, ">", $file or
        $logger->fatal("Cannot open file: $file");
    
    print $CON "dafe_db_dir=$dafe_db_dir\n";
    print $CON "reads_dir=$reads_dir\n";
    print $CON "out_dir=$out_dir\n";
    print $CON "ref_names_file=$ref_names_file\n";
    print $CON "sample_names_file=$sample_names_file\n";
    print $CON "perc_ids=$perc_ids\n";
    print $CON "read_file_exten=$read_file_exten\n" if is_defined($read_file_exten);
    print $CON "bam_file_prefix=$bam_file_prefix\n" if is_defined($bam_file_prefix);
    print $CON "cov_stat_file_name=$cov_stat_file_name\n" if is_defined($cov_stat_file_name);
    print $CON "htseq_file_prefix=$htseq_file_prefix\n" if is_defined($htseq_file_prefix);
    print $CON "htseq_i=$htseq_i\n" if is_defined($htseq_i);
    print $CON "pairs_file=$pairs_file\n" if is_defined($pairs_file);
    print $CON "skip_mapping=$skip_mapping\n" if is_defined($skip_mapping);
    print $CON "skip_filtering=$skip_filtering\n" if is_defined($skip_filtering);
    print $CON "skip_htseq=$skip_htseq\n" if is_defined($skip_htseq);
    print $CON "lsf_threads=$lsf_threads\n" if is_defined($lsf_threads);
    print $CON "lsf_mem=$lsf_mem\n" if is_defined($lsf_mem);
    print $CON "lsf_queue=$lsf_queue\n" if is_defined($lsf_queue);
    print $CON "lsf_out_file=$lsf_out_file\n" if is_defined($lsf_out_file);
    print $CON "lsf_err_file=$lsf_err_file\n" if is_defined($lsf_err_file);
    print $CON "lsf_job_name=$lsf_job_name\n" if is_defined($lsf_job_name);
    print $CON "queue_batch_count=$queue_batch_count\n" if is_defined($queue_batch_count);
    print $CON "node_batch_count=$node_batch_count\n" if is_defined($node_batch_count);
    print $CON "queue_max=$queue_max\n" if is_defined($queue_max);
    print $CON "runtime_max=$runtime_max\n" if is_defined($runtime_max);
    
    close($CON);
    
    return($file);
}

sub save_jobs {
    # no parameters
    
    my $file = "$out_dir/.restart.jobs";
    open my $JOBS, ">", $file or
        $logger->fatal("Cannot open file: $file");
    
    # NOTE: here the $jobs_aref variable is a global variable
    my $job_str = join("\n", @{$jobs_aref});
    print $JOBS $job_str;
    
    close($JOBS);
    
    return($file);
}

sub load_jobs {
    my ($file) = @_;
    
    my @jobs_arr = ();
    
    open my $JOBS, "<", $file or
        $logger->fatal("Cannot open jobs_file ($file) after restart");
    
    foreach my $line ( <$JOBS> ) {
        chomp $line;
        push @jobs_arr, $line;
    }
    
    close($JOBS);
    
    return(\@jobs_arr);
}

sub check_output_files {
    $logger->info("Checking output files for correctness");
    
    my ($ref_aref, $sample_aref);
    
    # NOTE: $perc_ids_aref does not include the smallest value
    #       (ie $min_perc_id) so I make a temporary aref and add back
    #       $min_perc_id.
    my $tmp_perc_ids_aref = [($min_perc_id), @{$perc_ids_aref}];

    
    # if the pairs file is provided use that to check the output
    if ( defined $pairs_file and -s $pairs_file ) {
        ($ref_aref, $sample_aref) = get_pairs($pairs_file);
        my $len = scalar @{$ref_aref};
        
        for ( my $i = 0; $i < $len; $i++ ) {
            check_bam_file($ref_aref->[$i], $sample_aref->[$i], $min_perc_id);
            foreach my $perc_id ( @{$tmp_perc_ids_aref} ) {
                check_htseq_file($ref_aref->[$i], $sample_aref->[$i], $perc_id,
                                 $htseq_i);
            }
        }
    }
    else {
        # this is actually the standard way where everything in the
        # ref_names_file and sample_names_file is ran
        
        $ref_aref = get_names($ref_names_file);
        $sample_aref = get_names($sample_names_file);
        
        foreach my $sample ( @{$sample_aref} ) {
            check_bam_file($sample, $min_perc_id);
            foreach my $perc_id ( @{$tmp_perc_ids_aref} ) {
                check_htseq_file($sample, $perc_id, $htseq_i);
            }
        }
    }
    
    $logger->info("Finished checking output files");
    
    return 1;
}

sub check_bam_file {
    my ($sample, $perc_id) = @_;
    
    my $file_name = "$out_dir/$sample/$bam_file_prefix\_id" . $perc_id . ".bam";
    
    if ( is_false($keep_bam_files) ) {
        # if keep_bam_files is false the bam files are deleted before I get to
        # this function which checks if they were created.
        return 1;
    }
    
    if ( ! -s $file_name ) {
        $logger->warn("Bam file is missing or empty: $file_name");
        # At this point it might be useful to check if the bam file does not
        # exist because the genome fna file is missing in the db or empty.
    }
    
    return 1;
}

sub check_htseq_file {
    my ($sample, $perc_id, $htseq_i) = @_;
    
    my $file_name = "$out_dir/$sample/$htseq_file_prefix\_id" . $perc_id . ".txt";
    
    if ( ! -e $file_name ) {
        # if the file is missing it's likely the bam file is also missing
        $logger->warn("htseq-count file is missing: $file_name");
    }
    
    if ( ! -s $file_name ) {
        # this means the file exists but is empty.  It likely means there were
        # not reads that mapped to this reference sequence
        $logger->warn("htseq-count file is empty: $file_name");
        $logger->warn("Trying to make an htseq-count file with all 0 values");
        make_htseq_file($file_name, $htseq_i);
    }
    
    return 1;
}

sub make_htseq_file {
    my ($file_name, $htseq_i) = @_;
    
    # first make sure the gff file exists and is non-empty in the dafe_db
    my $gff_file = "$dafe_db_dir/../$combined_db_name" . ".gff";
    if ( ! -s $gff_file ) {
        $logger->warn("Cannot make htseq-count file because missing gff: $gff_file");
    }
    else {
        # get all the CDS feature names
        my $cds_names_aref = parse_gff($gff_file, $htseq_i);
        
        # open the output file
        open my $HTSEQ, ">", $file_name or
            $logger->fatal("Cannot make htseq-count file because " .
                           "cannot open output file $file_name");
        foreach my $cds_name ( @{$cds_names_aref} ) {
            print $HTSEQ $cds_name . "\t0\n";
        }
        print $HTSEQ "__no_feature\t0\n";
        print $HTSEQ "__ambiguous\t0\n";
        print $HTSEQ "__too_lo_aQual\t0\n";
        print $HTSEQ "__not_aligned\t0\n";
        print $HTSEQ "__alignment_not_unique\t0";
        
        close($HTSEQ);
    }
    
    return 1;
}

sub parse_gff {
    my ($gff_file, $htseq_i) = @_;
    
    # this funciton simply returns a list of gene ID tags found in the CDS
    # features in the given GFF file.  The gene ID tag is specified in the
    # htseq_i variable.  Usually it will be "ID", but it could also be
    # "locus_tag".
    
    my @htseq_i_arr = ();
    
    open my $GFF, "<", $gff_file or
        $logger->fatal("Cannot make htseq-count file because cannot open " .
                       "GFF file: $gff_file");
    
    foreach my $line ( <$GFF> ) {
        chomp $line;
        
        my @vals = split(/\t/, $line);
        
        # skip features that are not CDS
        if ( $vals[2] !~ m/CDS/i ) {
            next;
        }
        
        my @tags = split(/;/, $vals[8]);
        foreach my $tag ( @tags ) {
            if ( $tag =~ m/$htseq_i=(.*)/i ) {
                push @htseq_i_arr, $1;
            }
        }
    }
    
    # sort the features
    my @sorted_htseq_i_arr = sort @htseq_i_arr;
    
    return(\@sorted_htseq_i_arr);
}

sub full_queue {
    my ($max) = @_;
    
    # checks to see if the number of jobs in the queue is less than $max
    
    # NOTE: I chose NOT to check only the number of jobs associated with this
    #       script in the queue when determining if the queue is full.  The
    #       reason for this is that I don't want to overload the queue if there
    #       are already hundreds of other jobs running.
    
    $logger->debug("Checking if the queue is full");
    
    my $command = "bjobs | wc -l";
    
    $logger->debug("Running command: $command");
    
    my $job_count = `$command`;
    $job_count--;  # the header line is not counted.
    
    $logger->debug("Current number of jobs: $job_count");
    
    if ( $job_count < $max ) {
        return 0;
    }
    
    return 1;
}

sub stall {
    my ($job_name) = @_;
    
    $logger->info("Stalling until <$job_name> jobs finish");
    
    my $command = 'bjobs -J ' . $job_name . " > $out_dir/ACTIVE_JOBS";
    `$command`; # run the command
    
    if ( -s "$out_dir/ACTIVE_JOBS" ) {
        return 1;  # send stall singnal
    }
    
    # remove the ACTIVE_JOBS file before sending the continue signal
    `rm "$out_dir/ACTIVE_JOBS"`;
    
    return 0;  # send continue signal
}

sub bash_stall {
    # DEPRECIATED
    my ($job_name) = @_;
    
    $logger->info("Stalling until <$job_name> jobs finish");
    
    # stalls this script while waiting for jobs in the cluster to complete
    
    my $command = 'bjobs -J ' . $job_name . ' > ' .  "$out_dir/ACTIVE_JOBS" . '
                          while [ -s ' . "$out_dir/ACTIVE_JOBS" . ' ]
                          do
                              sleep 10
                              bjobs -J ' . $job_name . ' > ' . "$out_dir/ACTIVE_JOBS" . '
                          done
                          rm ' . "$out_dir/ACTIVE_JOBS";
    
    $logger->debug("Running stall command: $command");
    
    `$command`;
    
    $logger->info("Stall finished.  All $job_name jobs are complete");
    
    return 1;
}

sub submit_batch {
    my ($jobs_aref, $count) = @_;
    
    # submits $count number of jobs to the LSF queue
    
    for (my $i = 0; $i < $count; $i++ ) {
        # exit loop if there are no jobs left to submit
        if ( scalar @{$jobs_aref} < 1 ) {
            last;
        }
        
        my $command = shift @{$jobs_aref};
        $logger->info("Submitting command: $command");
        `$command`;
    }
    
    return 1;
}

sub get_jobs {
    my @jobs = (); 
    my ($ref_aref, $sample_aref);
    my $command = get_bsub_command() . "\"";
    my $batch_count = 1;
    
    $sample_aref = get_names($sample_names_file);
    
    foreach my $sample ( @{$sample_aref} ) {
        my $tmp_command = get_command($sample);
        
        if ( $batch_count == $node_batch_count ) {
            # end the command and add it to the jobs array
            $command .= $tmp_command . "\";";
            push @jobs, $command;
            
            # start a new command the bsub part and a quote
            $command = get_bsub_command() . "\"";
            $batch_count = 1;
        }
        else {
            # keep adding to the current command
            $command .= $tmp_command;
            $batch_count++;
        }
        
        $logger->debug("Building command for: $sample");
        $logger->debug("Built command: $tmp_command");
    }

    return(\@jobs);
}

sub get_command {
    my ($sample) = @_;
    
    my $command = "";
    
    if ( is_false($skip_mapping) ) {
        $command .= get_mapping_command($sample, $min_perc_id);
    }
    
    if ( is_false($skip_filtering) ) {
        $command .= get_filter_command($sample, $perc_ids_aref);
    }
    
    if ( is_false($skip_htseq) ) {
        # NOTE: $perc_ids_aref does not include the smallest value
        #       (ie $min_perc_id) so I make a temporary aref and add back
        #       $min_perc_id.
        my $tmp_perc_ids_aref = [($min_perc_id), @{$perc_ids_aref}];
        
        $command .= get_htseq_command($sample, $tmp_perc_ids_aref, $htseq_i);
    }
    
    if ( is_false($keep_bam_files) ) {
        # NOTE: $perc_ids_aref does not include the smallest value
        #       (ie $min_perc_id) so I make a temporary aref and add back
        #       $min_perc_id.
        my $tmp_perc_ids_aref = [($min_perc_id), @{$perc_ids_aref}];
        
        $command .= get_rm_bam_command($sample, $tmp_perc_ids_aref);
    }
    
    return($command);
}

sub get_rm_bam_command {
    my ($sample, $perc_id_aref) = @_;
    
    my $command = " ";
    my $bam_dir = "$out_dir/$sample/";
    
    foreach my $perc_id ( @{$perc_id_aref} ) {
        my $bam_file = $bam_dir . $bam_file_prefix . "_id" . $perc_id . ".bam";
        
        # generate the command
        $command .= "rm $bam_file;";
    }
    
    return($command);
}

sub get_filter_command {
    my ($sample, $perc_id_aref) = @_;
    
    my $command = " ";
    my $bam_dir = "$out_dir/$sample/";
    
    foreach my $perc_id ( @{$perc_id_aref} ) {
    
        my ($in_bam_file, $out_bam_file, $json_file);
        
        # set the input bam file to be filtered   
        $in_bam_file = $bam_dir . $bam_file_prefix . "_id" . $min_perc_id . ".bam";
        
        # set the output bam file
        $out_bam_file = $bam_dir . $bam_file_prefix . "_id" . $perc_id . ".bam";
        
        # make the json file
        $json_file = make_json_file($sample, $perc_id, $bam_dir);
        
        # generate the command
        $command .= "bamtools filter ";
        $command .= "-in $in_bam_file ";
        $command .= "-out $out_bam_file ";
        $command .= "-script $json_file; ";
    }
    
    return($command);
}

sub make_json_file {
    my ($sample, $perc_id, $bam_dir) = @_;
    
    my $json_file = $bam_dir . "id" . $perc_id . ".json";
    open my $JSON, ">", $json_file or
        $logger->fatal("Cannot open json file ($json_file)");
    
    print $JSON "{\n";
    print $JSON "\t\"tag\":\"YI:>" . $perc_id . "\"\n";
    print $JSON "}";
    
    close($JSON);
    
    return($json_file);
}

sub get_htseq_command {
    my ($sample, $perc_id_aref, $htseq_i) = @_;
    
    my $command = " ";
    my $bam_dir = "$out_dir/$sample/";
    
    foreach my $perc_id ( @{$perc_id_aref} ) {
    
        my ($bam_file, $gff_file, $out_file);
        
        # set the bam file
        $bam_file = $bam_dir . $bam_file_prefix . "_id" . $perc_id . ".bam";
        
        # set the gff file
        $gff_file = "$dafe_db_dir/../$combined_db_name" . ".gff";
        
        # set the output file
        $out_file = "$bam_dir/$htseq_file_prefix" . "_id" . $perc_id . ".txt";
        
        # generate the command
        {
            $command .= "if [ `samtools view -c $bam_file` != 0 ]; then ";
            $command .= "htseq-count ";
            $command .= "-s no ";
            $command .= "-t CDS ";
            $command .= "-i $htseq_i ";
            $command .= "-f bam ";
            $command .= "-a 0 ";
            $command .= "$bam_file ";
            $command .= "$gff_file ";
            $command .= "> $out_file; ";
            $command .= "fi; ";
        }
    }
    
    return($command);
}

sub get_mapping_command {
    my ($sample, $perc_id) = @_;
    
    my ($reads_file, $reads_file_regx, $db_path, $out_bam_dir,
        $out_bam_file, $command);
            
    # find the sample reads file
    {
        $reads_file_regx = "$reads_dir/";
        $reads_file_regx .= $sample . "/";
        $reads_file_regx .= $sample . "*" . $read_file_exten;
        
        $reads_file = `find $reads_file_regx`;
        chomp $reads_file;
    }
    
    # generate the reference file name
    {
        $db_path = "$dafe_db_dir/../$combined_db_name/";
    }
    
    # generate the output bam file
    # this also creates the parent directories if they do not exist
    {
        $out_bam_dir = "$out_dir/";
        $out_bam_dir .= $sample . "/";
        
        if ( ! -d $out_bam_dir ) {
            `mkdir -p $out_bam_dir`;
        }
        
        $out_bam_file = $out_bam_dir . $bam_file_prefix . "_id" . $perc_id . ".bam";
    }
    
    # generate the command
    {
        my $t = $lsf_threads - 1;
        $perc_id = $perc_id / 100;
        $logger->debug("Converting mapping percent ID to decimal: $perc_id");
        
        $command = "bbmap.sh ";
        $command .= "path=$db_path ";
        $command .= "in=$reads_file ";
        $command .= "interleaved=false ";
        $command .= "ambiguous=random ";
        $command .= "minid=$perc_id ";
        $command .= "idtag=t ";
        $command .= "out=$out_bam_file ";
        $command .= "outputunmapped=f ";
        $command .= "sam=1.3 ";
        $command .= "nodisk=t ";
        $command .= "threads=$t; ";
        
        if ( defined $cov_stat_file_name and length $cov_stat_file_name > 0 ) {
            # an option to make a coverage stats output file in addition
            # to the output bam file
            my $cov_stat_file = $out_bam_dir . $cov_stat_file_name;
            $command .= "covstats=$cov_stat_file ";
        }
    }
    
    return($command);
}

sub get_bsub_command {
    
    # prepend the bsub stuff to the command
    
    my $bsub = "bsub ";
    $bsub .= "-M $lsf_mem ";
    $bsub .= "-q $lsf_queue ";
    $bsub .= "-n $lsf_threads ";
    $bsub .= "-R \"span[hosts=1]\" ";
    $bsub .= "-o $out_dir/$lsf_out_file ";
    $bsub .= "-e $out_dir/$lsf_err_file ";
    $bsub .= "-J $lsf_job_name ";
    
    return($bsub);
}

sub get_names {
    my ($file) = @_;
    
    open my $IN, "<", $file or
        $logger->fatal("Cannot open ref_names_file ($ref_names_file)");
    
    my @names = ();
    
    # these files should NOT have headers
    # I assume the names are in the first column
    my @vals = ();
    foreach my $line ( <$IN> ) {
        chomp $line;
        
        @vals = split(/\t/, $line);
        push @names, $vals[0];
        $logger->debug("Adding name ($vals[0]) from $file")
    }
    
    return(\@names);
}

sub get_pairs {
    my ($file) = @_;
    
    # gets the reference, sample pairs and retuns as two seperate arrays
    open my $IN, "<", $file or
        $logger->logdie("Cannot open pairs_file ($file)");
    
    my @ref_ids = ();
    my @sample_names = ();
    my @vals = ();
    foreach my $line ( <$IN> ) {
        chomp $line;
        
        @vals = split(/\t/, $line);
        push @ref_ids, $vals[0];
        push @sample_names, $vals[1];
    }
    
    return(\@ref_ids, \@sample_names);
}

sub load_modules {
    # NOTE: I couldn't get this code to work.  Module must be loaded manually.
    # this runs a few systems commands to load the modules I will need
    
    #$logger->info("Start load_modules()");
    #
    #`bash /nas02/apps/Modules/default/init/bash`;
    #`module add bbmap`;
    #`module add samtools`;
    #
    #$logger->info("End load_modules()");
    
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
        pod2usage(-message => "ERROR: --dafe_db_dir is not a directory\n\n",
                    -exitval => 2); 
    }
    
    # check the combined_db_name
    if ( ! defined $combined_db_name ) {
        pod2usage(-message => "ERROR: required --combined_db_name not defined\n\n",
					-exitval => 2);
    }
    if ( ! -d "$dafe_db_dir/$combined_db_name" ) {
        my $msg = "--combined_db_name is not a directory in $dafe_db_dir\n\n";
        pod2usage(-message => "ERROR: $msg",
                    -exitval => 2); 
    }
    
    # check the reads dir
    if ( ! defined $reads_dir ) {
        pod2usage(-message => "ERROR: required --reads_dir not defined\n\n",
					-exitval => 2);
	}
    if ( ! -d $dafe_db_dir ) {
        pod2usage(-message => "ERROR: --read_dir is not a directory\n\n",
                    -exitval => 2); 
    }
    
    # check the output dir
    if ( ! defined $out_dir ) {
        pod2usage(-message => "ERROR: required --out_dir not defined\n\n",
					-exitval => 2);
    }
    if ( ! -d $out_dir ) {
        `mkdir $out_dir`;
        $logger->info("Creating out_dir ($out_dir)");
    }
    
    # check restruct_out_dir
    if ( ! defined $restruct_out_dir ) {
        my $tmp = $out_dir;
        $tmp =~ s/\/$//;
        $restruct_out_dir = $tmp . "_restructured/";
    }
    if ( ! -d $restruct_out_dir ) {
        `mkdir $restruct_out_dir`;
        $logger->info("Creating restructured out dir ($restruct_out_dir)");
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
    
    # check perc_ids
    if ( ! defined $perc_ids ) {
        pod2usage(-message => "ERROR: required --perc_ids not defined\n\n",
					-exitval => 2);
    }
    else {
        my @ids = split(/,/, $perc_ids);
        my @sorted = sort { $a <=> $b } @ids;
        $min_perc_id = shift @sorted;
        $perc_ids_aref = \@sorted;
        
        $logger->debug("Setting min_perc_id to: $min_perc_id");
        
        # NOTE: I take off the minimum value and assign that to $min_perc_id.
        #       Later in the code I want to do operations on all the files
        #       created at different percent identites (ie htseq-count).  So
        #       in the get_command subroutine I add min_perc_id back to the
        #       $perc_ids_aref.  Sorry this is kind of poor design.
    }
    
    # check for at least some type of sample,ref input file
    if ( ! defined $pairs_file and
         ! defined $sample_names_file and
         ! defined $ref_names_file ) {
        pod2usage(-message => "ERROR: must provide [--pairs_file] OR [--sample_name_file and --ref_names_file]\n\n",
					-exitval => 2); 
    }
    
    # check the read_file_exten
    if ( ! defined $read_file_exten ) {
        $read_file_exten = "non_contaminants.fastq.gz";
    }
    
    # check genome_file_exten
    if ( ! defined $genome_file_exten ) {
        $genome_file_exten = ".fna";
    }
    
    # check gff_file_exten
    if ( ! defined $gff_file_exten ) {
        $gff_file_exten = ".gff";
    }
    
    # check the bam file name
    if ( ! defined $bam_file_prefix ) {
        $bam_file_prefix = "mapped";
    }
    
    # check cov_stat_file_name
    # --no checks required
    
    # check the htseq_file_prefix
    if ( ! defined $htseq_file_prefix ) {
        $htseq_file_prefix = "htseq_count";
    }
    
    # check the htseq_i value
    if ( ! defined $htseq_i ) {
        $htseq_i = "ID";
    }
    
    # check make_count_tbls_exe
    if ( ! defined $make_count_tbls_exe ) {
        $make_count_tbls_exe = "$Bin/DAFE_make_gene_count_tbls.sh";
    }
    
    # check the pairs file
	if ( defined $pairs_file and ! -e $pairs_file) { 
		pod2usage(-message => "ERROR: --pairs_file is an empty file\n\n",
					-exitval => 2); 
	}
    
    # check skip_mapping -- no checks required  It's a flag
    
    # check skip_filtering -- no checks required.  It's a flag
    
    # check skip_htseq -- no checks required  It's a flag
    
    # check keep_bam_files
    if ( ! defined $keep_bam_files ) {
        $keep_bam_files = 0;  # False
    }
    else {
        $keep_bam_files = 1;  # True
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
    
    # check queue_batch_count
    # I should do more checks here to make sure this is a digit
    if ( ! defined $queue_batch_count ) {
        $queue_batch_count = 100;
    }
    
    # check node_batch_count
    # I should do more checks here to make sure this is a digit
    if ( ! defined $node_batch_count ) {
        # when set to "1" each job is submitted individually
        $node_batch_count = 1;
    }
    
    # check queue_max
    # I should do more checks here to make sure this is a digit
    if ( ! defined $queue_max ) {
        $queue_max = 1000;
    }
    
    # check runtime_max
    if ( ! defined $runtime_max ) {
        $runtime_max = 561600;  # about 6.5 days
    }
    
    # log the parameters
    $logger->info("Parameters");
    $logger->info("--config_file: $config_file") if (defined $config_file);
    $logger->info("--dafe_db_dir: $dafe_db_dir");
    $logger->info("--combined_db_name: $combined_db_name");
    $logger->info("--reads_dir: $reads_dir");
    $logger->info("--out_dir: $out_dir");
    $logger->info("--restruct_out_dir: $restruct_out_dir");
    $logger->info("--ref_names_file: $ref_names_file");
    $logger->info("--sample_names_file: $sample_names_file");
    $logger->info("--perc_ids: $perc_ids");
    $logger->info("--read_file_exten: $read_file_exten");
    $logger->info("--genome_file_exten: $genome_file_exten");
    $logger->info("--gff_file_exten: $gff_file_exten");
    $logger->info("--bam_file_prefix: $bam_file_prefix");
    $logger->info("--cov_stat_file_name: $cov_stat_file_name") if (defined $cov_stat_file_name);
    $logger->info("--htseq_file_prefix: $htseq_file_prefix");
    $logger->info("--htseq_i: $htseq_i");
    $logger->info("--make_count_tbls_exe $make_count_tbls_exe");
    $logger->info("--pairs_file: $pairs_file") if (defined $pairs_file);
    $logger->info("--skip_mapping: $skip_mapping") if (defined $skip_mapping);
    $logger->info("--skip_filtering: $skip_filtering") if (defined $skip_filtering);
    $logger->info("--skip_htseq: $skip_htseq") if (defined $skip_htseq);
    $logger->info("--keep_bam_files") if ( is_true($keep_bam_files) );
    $logger->info("--lsf_threads: $lsf_threads");
    $logger->info("--lsf_mem: $lsf_mem");
    $logger->info("--lsf_queue: $lsf_queue");
    $logger->info("--lsf_out_file: $lsf_out_file");
    $logger->info("--lsf_err_file: $lsf_err_file");
    $logger->info("--lsf_job_name: $lsf_job_name");
    $logger->info("--queue_batch_count: $queue_batch_count");
    $logger->info("--node_batch_count: $node_batch_count");
    $logger->info("--queue_max: $queue_max");
    $logger->info("--runtime_max: $runtime_max");
	
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
    
    $dafe_db_dir = is_defined($params{''}{dafe_db_dir});
    $combined_db_name = is_defined($params{''}{combined_db_name});
    $reads_dir = is_defined($params{''}{reads_dir});
    $out_dir = is_defined($params{''}{out_dir});
    $restruct_out_dir = is_defined($params{''}{restruct_out_dir});
    $ref_names_file = is_defined($params{''}{ref_names_file});
    $sample_names_file = is_defined($params{''}{sample_names_file});
    $perc_ids = is_defined($params{''}{perc_ids});
    $read_file_exten = is_defined($params{''}{read_file_exten});
    $genome_file_exten = is_defined($params{''}{genome_file_exten});
    $gff_file_exten = is_defined($params{''}{gff_file_exten});
    $bam_file_prefix = is_defined($params{''}{bam_file_prefix});
    $cov_stat_file_name = is_defined($params{''}{cov_stat_file_name});
    $htseq_file_prefix = is_defined($params{''}{htseq_file_prefix});
    $htseq_i = is_defined($params{''}{htseq_i});
    $pairs_file = is_defined($params{''}{pairs_file});
    $skip_mapping = is_defined($params{''}{skip_mapping}, $skip_mapping);
    $skip_filtering = is_defined($params{''}{skip_filtering}, $skip_filtering);
    $skip_htseq = is_defined($params{''}{skip_htseq}, $skip_htseq);
    $keep_bam_files = is_defined($params{''}{keep_bam_files}, $keep_bam_files);
    $lsf_threads = is_defined($params{''}{lsf_threads}, $lsf_threads);
    $lsf_mem = is_defined($params{''}{lsf_mem}, $lsf_mem);
    $lsf_queue = is_defined($params{''}{lsf_queue}, $lsf_queue);
    $lsf_out_file = is_defined($params{''}{lsf_out_file}, $lsf_out_file);
    $lsf_err_file = is_defined($params{''}{lsf_err_file}, $lsf_err_file);
    $lsf_job_name = is_defined($params{''}{lsf_job_name}, $lsf_job_name);
    $queue_batch_count = is_defined($params{''}{queue_batch_count});
    $node_batch_count = is_defined($params{''}{node_batch_count});
    $queue_max = is_defined($params{''}{queue_max});
    $runtime_max = is_defined($params{''}{runtime_max});
    
    return 1;
}

sub check_env {
    # this subroutines checks for required executables
    # you may need to load certain modules to pass this check
    
    if ( ! can_run('bbmap.sh') and is_false($skip_mapping) ) {
        $logger->logdie("Cannot find required bbmap program!");
    }
    
    if ( ! can_run('samtools') and is_false($skip_mapping) ) {
        $logger->logdie("Cannot find required samtools program!");
    }
    
    if ( ! can_run('bamtools') and is_false($skip_filtering) ) {
        $logger->logdie("Cannot find required bamtools program!");
    }
    
    if ( ! can_run('htseq-count' ) and is_false($skip_htseq) ) {
        $logger->logdie("Cannot find required htseq-count program!");
    }
    
    # can_run does not work here unless the file is set with executable
    # permission.  For some reason can_run doesn't recognize the
    # #!usr/bin/env bash line but it does recognize the perl one.  So this
    # simply checks to make sure there is a non-zero sized file in the variable
    if ( ! -s $make_count_tbls_exe ) {
        $logger->logdie("Cannot find DAFE_make_gene_count_tbls.sh: $make_count_tbls_exe");
    }
    
    return 1;
}

sub is_true {
    my ($val) = @_;
    
    my $is_true = 1;  # default is true
    
    if ( ! defined $val ) { $is_true = 0; }
    elsif ( $val eq "" ) { $is_true = 0; }
    elsif ( $val == 0 ) { $is_true = 0; }
    elsif ( $val =~ m/f|false/i ) { $is_true = 0; }
    elsif ( $val =~ m/t|true/i ) { $is_true = 1; }
    elsif ( $val == 1 ) { $is_true = 1; }
    
    return($is_true);
}

sub is_false {
    my ($val) = @_;
    
    return( ! is_true($val) );
}

sub is_defined {
    my ($val, $default) = @_;
    
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

DAFE_count.pl - Runs read mapping, htseq-count, and bam filtering on the UNC cluster


=head1 VERSION

This documentation refers to version 0.0.1


=head1 SYNOPSIS

    DAFE_count.pl
        [--config_file DAFE_count.conf]
        --dafe_db_dir my_dafe_db/
        --combined_db_name all_genomes
        --read_dir my_metagenome_samples/
        --out_dir my_out_dir/
        [--restruct_out_dir restructured_out_dir/]
        --ref_names_file names.txt
        --sample_names_file samples.names.txt
        --perc_ids "60,70,80,90,95"
        [--read_file_exten "non_contaminants.fastq.gz"]
        [--genome_file_exten ".fna"]
        [--gff_file_exten ".gff"]
        [--bam_file_prefix "mapped"]
        [--cov_stat_file "cov_stat.txt"]
        [--htseq_file_prefix "htseq.count.txt"]
        [--htseq_i "ID"]
        [--make_count_tbls_exe "DAFE_make_gene_count_tbls.sh"]
        [--pairs_file run_these_pairs.txt]
        [--skip_mapping]
        [--skip_filtering]
        [--skip_htseq]
        [--keep_bam_files]
        [--lsf_threads 2]
        [--lsf_mem 4]
        [--lsf_queue week]
        [--lsf_out_file lsf.out]
        [--lsf_err_file lsf.err]
        [--lsf_job_name DAFE_count]
        [--queue_batch_count 100]
        [--node_batch_count 1]
        [--queue_max 1000]
        [--runtime_max 561600]
        
        [--help]
        [--man]
        [--debug]
        [--verbose]
        [--quiet]
        [--logfile logfile.log]

    --config_file = Path to config file with any of the described parameters
    --dafe_db_dir = Path to DAFE database directory
    --combined_db_name = Name of the directory with the bbmap reference db
    --reads_dir = Path to directory with metagenome sample reads
    --out_dir = Path to where output files will be stored
    --restruct_out_dir = Path to where restructured output files are stored
    --ref_names_file = Path to file with reference names
    --sample_names_file = Path to file with sample names
    --perc_ids = Comma seperated list of percent IDs to using when filtering
    --read_file_exten = File extension of sample read files
    --genome_file_exten = File extension for the genome fasta in the DAFE db
    --gff_file_exten = File extension for the gff file in the DAFE db
    --bam_file_prefix = Prefix for output bam file
    --cov_stat_file_name = File name for cov_stat_file created by bbmap
    --htseq_file_prefix = Prefix for output htseq-count file
    --htseq_i = GFF attribute tag name for htseq-count -i parameter
    --make_count_tbls_exe = Path to the bash script for making gene count tables
    --pairs_file  = Path to an input file with two columns: 1)ref 2)sample
    --skip_mapping = Flag that skips the mapping command when it is set
    --skip_filtering = Flag that skips the filtering step when it isset
    --skip_htseq = Flag that skips the htseq-count step when it is set
    --keep_bam_files = Flag to prevent bam files from being deleted
    --lsf_threads = Number of threads to use for read mapping
    --lsf_mem = Number of GB of memory to allocate
    --lsf_queue = UNC compute cluster queue to submit jobs to
    --lsf_out_file = Name of lsf output file
    --lsf_err_file = Name of lsf error file
    --lsf_job_name = Job name for lsf jobs
    --queue_batch_count = Number of jobs to batch submit to --lsf_queue
    --node_batch_count = Number of jobs to submit to a single node
    --queue_max = Number of jobs allowed in queue (plus --queue_batch_count)
    --runtime_max = Number of seconds the master job will run before it restarts
    --help  = Prints USAGE statement
    --man   = Prints the man page
    --debug	= Prints Log4perl DEBUG+ messages
    --verbose	= Prints Log4perl INFO+ messages
    --quiet	= Suppress printing ERROR+ Log4perl messages
    --logfile	= File to save Log4perl messages
    
    Remember to load the following module if running on Kure/Killdevil
    - bbmap
    - htseq-count
    - bamtools
    - samtools


=head1 ARGUMENTS

=head2 [--config_file]

The config file is option, but I would STRONGLY recommend using it as opposed
to passing individual parameter values.

=head2 --dafe_db_dir

Path to DAFE database directory.  In the DAFE database directory root there
should exist a directory for each reference sequence (eg genome).  Each of these
reference directories are required to have a fasta genome file (*.fna), GFF file
(*.gff), ......  The script ....... can be used to check the DAFE database for
missing data or potential errors.

=head2 --combined_db_name

Name of the directory with the bbmap reference database.  This directory must
be found in --dafe_db_dir.  It can be created by running
DAFE_db_make_combined_db.pl.

=head2 --reads_dir

Path to directory with metagneome sample reads.  Within --reads_dir there should
exist a directory with a name matching that in the --sample_names_file.  Within
each of those sample directories should be a gzipped fastq file of metagenome
reads associated with that sample.  The fastq.gz file should be named something
like "sample_name*.fastq.gz" where sample_name is the same name as found in the
--sample_names_file.

=head2 --out_dir

Path to where output files will be stored.  If --out_dir does not exist as a
directory it will be created.

=head2 [--restruct_out_dir]

Path to where restrucutred output files are stored.  When I changed the
algorithm to use a combined database I restrucutre the output files to have the
same format as when I use the multi-genome database scheme.  This way downstream
scripts will continue to work.
DEFAULT: --out_dir . "_restructured"

=head2 --ref_names_file

Path to file with reference names.  The names (ie JGI Genome OIDs) must be in
the first column.  They must also correspond to the names/IDs used in the DAFE
database.  If the --pairs_file parameter is given this file is ignored.  This
file should NOT have headers.

=head2 --sample_names_file

Path to file with sample names.  The names must be in the first column.  If the
--pairs_file paramter is given this file is ignored.  This file should NOT have
headers.

=head2 --perc_ids

Comma seperated list of percent IDs to using when filtering.  The smallest value
is used for the initial mapping step.  Other bam files are created using the
remaining percent IDs by running "bamtools filter".  Also, htseq-count is ran on
each bam file.  So there will be a bam file and htseq-count file for each of the
given percent IDs.  This parameter can be a single number if the user wants to
only analyize mappings at a single percent identity (ie "60").

=head2 [--read_file_exten]

File extension of sample read files.
DEFAULT: "non_contaminants.fastq.gz"

=head2 [--genome_file_exten]

File extension for the genome fasta file in the DAFE database.
DEFAULT: ".fna"

=head2 [--gff_file_exten]

File extension for the gff file in the DAFE database
DEFAULT: ".gff"

=head2 [--bam_file_prefix]

Prefix for output bam file.  The actual output file name will append "_id60.bam"
where "60" is the percent identity used to create that bam file.
DEFAULT: "mapped"

=head2 [--cov_stat_file_name]

A file name for the cog_stat_file create by bbmap.  If this parameter is present
bbmap will output the cov_stat_file which gives coverage statistics.  See bbmap
documetation for details.
DEFAUL: NULL

=head2 [--htseq_file_prefix]

Prefix for htseq-count output file.  The actual output file name will append
"_id60.txt" where "60" is the percetn identity used to create that count txt
file
DEFAULT: "htseq_count"

=head2 [--htseq_i]

GFF tag name for the -i parameter when running htseq-count.  This parameter
tells htseq-count what identifier to use when outputing the counts of each gene.
This tag name MUST be in each gene (ie CDS entry), and it MUST be unique across
all genes in the genome.  I was previously using "locus_tag" for --htseq_i but
found that not every gene has a locus_tag attribute in the GFF file.
DEFAULT: "ID"

=head2 [--make_count_tbls_exe]

Path to the bash script for making gene count tables.  Eventually I want to
incorperate the commands in this bash script into DAFE_count.pl.  By default it
looks in the same directory as this perl script for a bash script named
DAFE_make_gene_count_tbls.sh
DEFAULT: DAFE_make_gene_count_tbls.sh
    
=head2 [--pairs_file]

Path to an input file with two columns.  The first column must be reference IDs
and the second must be sample names.  When using this file only those pairs are
ran.  This is a good option for only running very specific sample,reference
pairs.  There should be no headers on this file.

=head2 [--skip_mapping]

If this flag is used the bbmap mapping command is skipped.  This flag should
only be used if the mapping command was correctly executed in a previous run.
This option become useful if a step after the mapping command fails and the
program needs to be restarted, or if the user wants to rerun the steps after
mapping with different parameter values.

=head2 [--skip_filtering]

If this flag is used the bamtools filtering command is skipped.  This flag
should only be used if the filtering command was correctly executed in a
previous run. This option become useful if a step after the filtering command
fails and the program needs to be restarted, or if the user wants to rerun the
steps after filtering with different parameter values.

=head2 [--skip_htseq]

If this flag is used the htseq-count command is skipped.  This flag should
only be used if the htseq-count command was correctly executed in a previous run.
This option become useful if a step after the htseq-count command fails and the
program needs to be restarted, or if the user wants to rerun the steps after
htseq-count with different parameter values.

=head2 [--keep_bam_files]

If this flag is used all the bam files (ie for each percent identity) are kept.
For large runs it is recommended to NOT use this flag as the bam files can take
up a lot of space.  In an older version I kept only the lowest percent identity
bam file.  All other bam files could easily be recreated by running the
filtering commands.  However, even with only one bam file it was taking up way
to much space (more than my quota would allow--8TB).

=head2 [--lsf_threads]

Number of threads to use for read mapping step.  Must be greater than or equal
to 2.
DEFAUL: 2

=head2 [--lsf_mem]

Number of GB of memory to allocate.  Must be a nunmber > 0.
DEFAULT: 4

=head2 [--lsf_queue]

Queue to submit LSF jobs to
DEFAULT: week

=head2 [--lsf_out_file]

Name of LSF output file
DEFAULT: lsf.out

=head2 [--lsf_err_file]

Name of LSF error file
DEFAULT: lsf.err

=head2 [--lsf_job_name]

Name of jobs submitted to queue
DEFAULT: DAFE_count

=head2 [--queue_batch_count]

The number of jobs submitted to the queue at one time.  When there are fewer
than --queue_max jobs in the queue --queue_batch_count jobs are submitted until
there are now more jobs to submit.
DEFAULT: 100

=head2 [--node_batch_count]

The number of jobs to submit to a single node.  The default is "1" meaning that
each jobs gets it's own node.  Remember a job consists of mapping reads from one
sample to one reference sequence (eg genome).  If there are very few reads in
the sample it may be best to batch sevearl of these jobs on a single node.  You
just need to be sure that the batched group of jobs will finish before the queue
runtime limit.
DEFAULT: 1

=head2 [--queue_max]

The max number of jobs allowed in the queue (plus --queue_batch_count).  When the
number of jobs in the queue falls bellow --queue_max a new batch of
--queue_batch_count jobs may be submitted.
DEFAUL: 1000

=head2 [--runtime_max]

Number of seconds the master job will run before it restarts.
DEFAULT: 561600 (about 6.5 days)

=head2 [--help]
    
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

Each section below describes various aspects of this program.

=head General Notes

There are three main commands executed by this script on each reference, sample
pair.  These commands include "bbmap" which maps reads to individual reference
sequences, "bamtools filter" which creates mapping files at the given percent
identities greater than the percent identity used in the original mapping, and
"htseq-count" which counts the number of reads mapping to each [gene] feature in
each reference sequecne for each percent identity used in the filter step.

This scirpt is meant to be ran on the UNC clusters Kure or KillDevil.  It
incorperates its own job scheduling system.  The job scheduling is controled by
several parameters including --queue_batch_count, --queue_max, etc.

There is a parameter purposfully not included in the documentation
(--jobs_file).  This parameter is only used when this script kills itself and
restarts.  This restart is required if the jobs cannot be completed within the
runtime limit of the cluster (usually a week).

=head2 Config File

=head2 Mapping

The lowest value in the commma seperated list of --perc_ids is used in the
original mapping.

=head2 Filtering

=head2 htseq

Right now the name of the filtered bam files are hard coded as id70_mapped.bam
where "70" represents the percent identity of the filter.

If you try filtering below the percent identity of the original mapping you will
get the same bam file as the original mapping.  So it only makes sense to filter
at higher percent idenities than the original mapping percent identity.

=head2 To Do

- Split this scirpt into a collection of perl objects (params_obj,
lsf_control_obj, command_master_obj, check_output_obj, etc)
- Add features to resubmit the restart with the original lsf parameters (ie
queue, mem, out files, etc.)

=head1 CONFIGURATION AND ENVIRONMENT
    
The commands "bbmap", "samtools", "bamtools", "htseq-count", and
"DAFE_make_gene_count_tbls.sh" must be executable.  They are called as system
commands in this program.
    
    
=head1 DEPENDANCIES

version
Getopt::Long
Pod::Usage
Carp
Readonly
version
Log::Log4perl qw(:easy)
Log::Log4perl::CommandLine qw(:all)
Config::Std
IPC::Cmd qw(can_run run)
Data::Dumper
List::Util qw(min)
FindBin qw($Bin)


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
