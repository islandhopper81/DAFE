#!/usr/bin/env perl

# makes the full DA table by parsing the tag_data files output by edgeR analysis

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Carp;
use Readonly;
use version; our $VERSION = qv('0.0.1');
use Log::Log4perl qw(:easy);
use Log::Log4perl::CommandLine qw(:all);
use Array::Transpose;
use Scalar::Util qw(looks_like_number);

# Subroutines #
sub guess_pa_index;
sub get_clstr_col;
sub check_params;
sub _is_defined;

# Variables #
my ($stat_dir, $dafe_db, $clstr_metadata, $ref_meta_file, $clstr_mapping_file,
    $tag_file_name, $da_col_file, $annote_file_name, $clstr_summary_file,
    $da_filt_col_file, $final_tbl_file, $MIN_PERC, $MIN_FDR,
    $help, $man);

my $options_okay = GetOptions (
    "stat_dir:s" => \$stat_dir,
    "dafe_db:s" => \$dafe_db,
    "clstr_metadata:s" => \$clstr_metadata,
    "ref_meta_file:s" => \$ref_meta_file,
    "clstr_mapping_file:s" => \$clstr_mapping_file,
    "tag_file_name:s" => \$tag_file_name,
    "da_col_file:s" => \$da_col_file,
    "annote_file_name:s" => \$annote_file_name,
    "clstr_summary_file:s" => \$clstr_summary_file,
    "da_filt_col_file:s" => \$da_filt_col_file,
    "final_tbl_file:s" => \$final_tbl_file,
    "min_perc:i" => \$MIN_PERC,
    "min_fdr:i" => \$MIN_FDR,
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


# read in the group annotation file to get a list of all possible groups.  Each
# line in the annotation file is a group.  Groups can be defined like COGs, KOs,
# or even clusters generated using Isai and Asaf's clustering approaches.
# These names are used to get a compresensive list of all possible groups and
# determine their order

open my $META, "<", $clstr_metadata or
    $logger->logdie("Cannot open --clstr_metadata: $clstr_metadata");

my @clstr_names = ();  # names of all the clusters
my %index = ();  # index of each cluster key => cluster; value => index
my @vals = ();
my $i = 0;
foreach my $line ( <$META> ) {
    chomp $line;
    
    @vals = split(/\t/, $line);
    
    push @clstr_names, $vals[0];
    $index{$vals[0]} = $i;
    $i++;
}

close($META);


# go through each of the genomes in the stat_dir
opendir my $DIR, "$stat_dir" or
    $logger->logdie("Cannot open --stat_dir: $stat_dir");

my @genome_names = ();  # an array of genomes names used while printing.
my $annote_file = "";
my $tag_file = "";
my $first_line = 1;  # boolean for testing if we are on the first line
my $clstr_col_index;  # clstr column number in annotation file
my $clstr_id;  # not this gets set twice in this loop. 1)annote file, 2)tag file
my $logfc; # set when looking at tag file
my $fdr; # set when looking at tag file
my $out_path;
my $DA_val = -4; # note this is used in a couple loops
my $num_genomes_in_clstr = init_num_genomes_in_clstr(\@clstr_names);
my %seen_clstr = ();  # records when we see a clstr in a genome
my %DA_count = ();
my %DA_up_count = ();
my %DA_dn_count = ();


# loop through all the genomes
# I will have to do this a second time after I figure out what clusters
# I want to keep
foreach my $dir ( readdir($DIR) ) {
    # skip hidden dirs
    if ( $dir =~ m/^\./ ) {next;}
    push @genome_names, $dir;
    
    $logger->debug("Starting dir: $dir");
    
    # create an array to hold the DA values
    # the array has position for each possible cluster
    # the -3 is the DA value indicating the genome does NOT contain the cluster.
    # these values will be change appropriately later in the algorithm.
    my @full_DA_col = (-3) x (scalar @clstr_names);
    
    
    ### Determine the clusters that are present in the genome (will get a -2 DA
    # value) by looking in the genome annotation file.  Also, record the
    # number of genomes in each cluster.
    $annote_file = "$dafe_db/$dir/$annote_file_name";
    open my $AN, "<", "$annote_file" or
        $logger->logdie("Cannot open genome annotation file: $annote_file");
    
    %seen_clstr = ();
    $first_line = 1;
    foreach my $line ( <$AN> ){
        chomp $line;
        if ( $first_line ) {
            $clstr_col_index = get_clstr_col($line);
            $first_line=0;
            next; # skip the header line.
        }
        
        @vals = split(/\t/, $line);
        $clstr_id = $vals[$clstr_col_index];
        
        # the clstr value might be NA if it was not included in a cluster
        if ( $clstr_id eq "NA" ) { next; }
        
        # make sure the cluster is defined.  If it is not defined it must not
        # be in the master cluster metadata file
        if ( ! defined $index{$clstr_id} ) {
            $logger->warn("Cluster not defined in annotation file: $clstr_id");
            next;
        }
        
        # set all the known clusters in the genome to -2.  In the next step
        # they will be set to -1, 0, 1 if they are in the DA tag output
        $full_DA_col[$index{$clstr_id}] = -2;
        
        # record this cluster as being found in this genome
        # but only count clusters that we have not yet seen (ie don't double count)
        if ( ! $seen_clstr{$clstr_id} ) {
            $num_genomes_in_clstr->{$clstr_id}++;
            $seen_clstr{$clstr_id} = 1;  # now we have seen this cluster
        }
    }
    close($AN);
    
    
    ### Determine the clusters that are DA (-1, 0, 1) by looking in the tag file
    $tag_file = "$stat_dir/$dir/$tag_file_name";
    open my $TAG, "<", "$tag_file" or
        $logger->warn("Cannot open tag file: $tag_file");
    
    # if the file handle is not open skip this one
    if ( tell($TAG) == -1 ) {next;}
    
    $first_line = 1;
    foreach my $line ( <$TAG> ) {
        chomp $line;
        if ( $first_line ) { $first_line=0; next; } # skip the header line.
        
        @vals = split(/ /, $line);
        $clstr_id = $vals[0];
        $logfc = $vals[1];
        $fdr = $vals[4];
        
        # make sure the cluster is defined.  If it is not defined it must not
        # be in the master cluster metadata file
        if ( ! defined $index{$clstr_id} ) {
            $logger->warn("Cluster not defined in tag file: $clstr_id");
            next;
        }
        
        # determine if the cluster is DA
        $DA_val = get_DA_call($fdr);
        $full_DA_col[$index{$clstr_id}] = $DA_val;
        count_DA_val($DA_val, $clstr_id);
       
    }
    close($TAG);
    
    
    ### print the outputs
    
    # print the DA column
    print_full_DA_col($dir, \@clstr_names, \@full_DA_col)

}
closedir($DIR);

# build a lookup has that can translate a cluster to a list of genomes
# into the number of PA genomes
# I need the genome-to-cluster-to-gene table file and the reference (ie genome)
# metadata file.

# read in the reference metadata
my %genome_to_pa_call = ();
open my $REF, "<", $ref_meta_file or
    $logger->logdie("Cannot open --ref_meta_file: $ref_meta_file");
    
# I only consder the genomes for the genomes that are in the stat_dir
my %genome_names = ();
foreach my $g_name ( @genome_names ) {
    $genome_names{$g_name} = 1;
}

# go through the reference metadata
$first_line = 1;
my $pa_index;
foreach my $line ( <$REF> ) {
    chomp $line;
    
    if ( $first_line ) {
        $pa_index = guess_pa_index($line);
        $first_line = 0;
        next;
    }
    
    @vals = split(/\t/, $line);
    if ( defined $genome_names{$vals[0]} ) {
        $genome_to_pa_call{$vals[0]} = $vals[$pa_index];
    }
}
close($REF);

# read in the genome-to-cluster-to-gene table
my %clstr_pa_count = ();
open my $GTCTG, "<", $clstr_mapping_file or
    $logger->logide("Cannot open --clstr_mapping_file: $clstr_mapping_file");

my $current_genome = "";
my %seen_clusters = ();
# no headers in this file
foreach my $line ( <$GTCTG> ) {
    chomp $line;

    @vals = split(/\t/, $line);
    
    # make sure the genome is in the genome_to_pa_call hash.  Otherwise I can
    # just skip it
    if ( ! defined $genome_to_pa_call{$vals[0]} ) { next; }
    
    # reset the seen clusters if on a different genome
    if ( $current_genome ne $vals[0] ) {
        $current_genome = $vals[0];
        %seen_clusters = ();
    }
    
    # don't count the cluster if I've already seen it for this genomes
    if ( defined $seen_clusters{$vals[1]} ) { next; }
    
    # Count it as PA if it is
    if ( defined $clstr_pa_count{$vals[1]} ) {
        if ( $genome_to_pa_call{$vals[0]} eq "PA" ) {
            $clstr_pa_count{$vals[1]}++;
        }
    }
    else {
        if ( $genome_to_pa_call{$vals[0]} eq "PA" ) {
            $clstr_pa_count{$vals[1]} = 1;
        }
    }
    
    # remember that we have seen this cluster
    $seen_clusters{$vals[1]} = 1;
}

close($GTCTG);

### Print the file with all the cluster information (ie the clstr summary file)
open my $FULL, ">", $clstr_summary_file or
    $logger->logdie("Cannot open full DA file: $clstr_summary_file");

print $FULL "clstr_id\tnum_genomes_in_clstr\tDA_count\tDA_up_count\tDA_dn_count\tPA_count\n";

foreach my $clstr_id ( keys %{$num_genomes_in_clstr} ) {
    print $FULL $clstr_id, "\t";
    print $FULL $num_genomes_in_clstr->{$clstr_id}, "\t";
    if ( $DA_count{$clstr_id} ) {
        print $FULL $DA_count{$clstr_id};
    }
    else {
        print $FULL "0";
    }
    print $FULL "\t";
    
    if ( $DA_up_count{$clstr_id} ) {
        print $FULL $DA_up_count{$clstr_id};
    }
    else {
        print $FULL "0";
    }
    print $FULL "\t";
    
    if ( $DA_dn_count{$clstr_id} ) {
        print $FULL $DA_dn_count{$clstr_id};
    }
    else {
        print $FULL "0";
    }
    print $FULL "\t";
    
    if ( $clstr_pa_count{$clstr_id} ) {
        print $FULL $clstr_pa_count{$clstr_id};
    }
    else {
        print $FULL "0";
    }
    print $FULL "\n";
}

close($FULL);


### At this point I have all the information I need to filter out clusters that are
# uninformative, but I have to go through each dir again.  But first I need to get
# a list of clusters that I want to keep based on the filtering parameters

# get a list of clusters I want to keep
# I want to filter based on the percent of genomes in which a cluster is called DA
my $perc_DA;
my %clstrs_to_keep = ();
my $genome_count = scalar(@genome_names);

foreach my $clstr_id ( keys %DA_count ) {
    $perc_DA = ($DA_count{$clstr_id} / $genome_count) * 100;
    if ( $perc_DA >= $MIN_PERC ) {
        $clstrs_to_keep{$clstr_id} = 1;
    }
}

# go through each dir again to print the filtered table
# also save each genome's column to print the filtered table
opendir $DIR, "$stat_dir" or
    $logger->logdie("Cannot open --stat_dir: $stat_dir");

my $full_DA_col_path;
my $COL;
$first_line = 1;
my %filtered_tbl = ();
my $filt_DA_col_path;
my $FILT;

foreach my $dir ( readdir($DIR) ) {
    if ( $dir =~ m/^\./ ) { next; } # skip hidden files
    
    # go through the full_DA_col_file
    $full_DA_col_path = "$stat_dir/$dir/$da_col_file";
    $filt_DA_col_path = "$stat_dir/$dir/$da_filt_col_file";
    
    open $COL, "<", $full_DA_col_path or
        $logger->logdie("Cannot open recently created full_DA_col: $full_DA_col_path");
        
    open $FILT, ">", $filt_DA_col_path or
        $logger->logdie("Cannot open filtered DA column path: $filt_DA_col_path");
    
    foreach my $line ( <$COL> ) {
        chomp $line;
        #if ( $first_line ) { $first_line = 0; next; } #skip header line
        
        @vals = split(/\t/, $line);
        
        if ( defined $clstrs_to_keep{$vals[0]} ) {
            # this is a good cluster
            # print it and save it in the final table
            
            # remember $vals[0] is the cluster ID
            print $FILT $vals[0], "\t", $vals[1], "\n";
            $filtered_tbl{$vals[0]}{$dir} = $vals[1];
        }
    }
    close($COL);
    close($FILT);
}
close($DIR);


# print the final filtered table with information from all final clusters
# for each genome.
open my $FINAL, ">", $final_tbl_file or
    $logger->logdie("Cannot open final table file: $final_tbl_file");

# print the headers (ie genome names)
print $FINAL join("\t", keys %filtered_tbl), "\n";

foreach my $genome ( @genome_names ) {
    @vals = ($genome);
    foreach my $clstr_id ( keys %filtered_tbl ) {
        push @vals, $filtered_tbl{$clstr_id}{$genome};
    }
    print $FINAL join("\t", @vals), "\n";
}
close($FINAL);



########
# Subs #
########
sub guess_pa_index {
    my ($headers) = @_;
    
    my @vals = split(/\t/, $headers);
    
    my $i = 0;
    foreach my $val ( @vals ) {
        if ( $val =~ m/Label/i ) {
            return($i);
        }
        $i++;
    }
    
    return(-1);
}

sub count_DA_val {
    my ($DA_val, $clstr_id) = @_;
    
    if ( defined $DA_count{$clstr_id} ) {
        if ( $DA_val == 1 ) {
            $DA_count{$clstr_id}++;
            $DA_up_count{$clstr_id}++;
        }
        elsif ( $DA_val == -1 ) {
            $DA_count{$clstr_id}++;
            $DA_dn_count{$clstr_id}++;
        }
    }
    else {
        if ( $DA_val == 1 ) {
            $DA_count{$clstr_id}++;
            $DA_up_count{$clstr_id} = 1;
        }
        elsif ( $DA_val == -1 ) {
            $DA_count{$clstr_id}++;
            $DA_dn_count{$clstr_id} = 1;
        }
    }
    
    return(1);
}

sub print_full_DA_col {
    my ($dir, $clstr_names_aref, $full_DA_col_aref) = @_;
    
    $out_path = "$stat_dir/$dir/$da_col_file";
    open my $OUT, ">", $out_path or
        $logger->logdie("Cannot open out file: $out_path");
    
    foreach my $clstr_name ( @$clstr_names_aref ) {
        $DA_val = $full_DA_col_aref->[$index{$clstr_name}];
        print $OUT $clstr_name, "\t", $DA_val, "\n";
    }
    
    close($OUT);
    
    return(1);
}

# initializes the num_genomes_in_clstr hash.  This hash counts the number of
# genomes that have at least one sequence in a given cluster.
sub init_num_genomes_in_clstr {
    my ($clstr_names_aref) = @_;
    
    my %num_genomes_in_clstr = ();
    
    foreach my $clstr_name ( @$clstr_names_aref ) {
        $num_genomes_in_clstr{$clstr_name} = 0;
    }
    
    return(\%num_genomes_in_clstr);
}

# For getting the DA call (-1, 0, 1) given a genes FDR value
sub get_DA_call {
    my ($fdr, $clstr) = @_;
    
    my $DA_val;
    
    if ( $fdr < $MIN_FDR ) {
        if ( $logfc > 0 ) {
            $DA_val = 1;
        }
        else {
            $DA_val = -1;
        }
    }
    else {
        $DA_val = 0;
    }
    
    return($DA_val);
}

# For guessing the column in the annotation file that corresponds to the
# cluster ID for each gene.
sub get_clstr_col {
    my ($headers) = @_;
    
    my @vals = split(/\t/, $headers);
    my $i = 0;
    foreach my $v ( @vals ) {
        if ( $v =~ m/clstr/i ) {
            return($i);
        }
        $i++;
    }
    
    $logger->warn("Clstr column not found in annotation file headers");
    return(9)
}

sub check_params {
	# check for required variables
    if ( ! defined $stat_dir ) {
		pod2usage(-message => "ERROR: required --stat_dir not defined\n\n",
					-exitval => 2);
	}
    if ( ! defined $dafe_db ) {
        pod2usage(-message => "ERROR: required --dafe_db not defined\n\n",
					-exitval => 2);
    }
    if ( ! defined $clstr_metadata ) {
        pod2usage(-message => "ERROR: required --clstr_metadata not defined\n\n",
					-exitval => 2);
    }
    if ( ! defined $ref_meta_file ) {
        $ref_meta_file = "/netscr/yourston/compMetaG/metadata/reference_metadata_20160519.txt";
        $logger->info("Setting --ref_meta_file: $ref_meta_file");
    }
    if ( ! defined $clstr_mapping_file ) {
        $clstr_mapping_file = "/netscr/yourston/compMetaG/metadata/gene_to_cdhit_clstr_mapping_file.txt";
        $logger->info("Setting -- clstr_mapping_file: $clstr_mapping_file");
    }
    if ( ! defined $tag_file_name ) {
        $tag_file_name = "tags_data.txt";
        $logger->info("Setting --tag_file_name: $tag_file_name");
    }
	if ( ! defined $da_col_file) {
        $da_col_file = "dafe_full_DA_col.txt";
        $logger->info("Setting --da_col_file: $da_col_file");
	}
    if ( ! defined $annote_file_name ) {
        $annote_file_name = "all_annote.txt";
        $logger->info("Setting --annote_file_name: $annote_file_name");
    }
    if ( ! defined $clstr_summary_file ) {
        $clstr_summary_file = "clstr_summary.txt";
        $logger->info("Setting --clstr_summary_file: $clstr_summary_file");
    }
    if ( ! defined $da_filt_col_file ) {
        $da_filt_col_file = "dafe_filt_DA_col.txt";
        $logger->info("Setting --da_filt_col_file: $da_filt_col_file");
    }
    if ( ! defined $final_tbl_file ) {
        $final_tbl_file = "dafe_filt_tbl.txt";
        $logger->info("Setting --final_tbl_file: $final_tbl_file");
    }
    if ( ! defined $MIN_PERC ) {
        $MIN_PERC = 2;
        $logger->info("Setting --min_perc: $MIN_PERC");
    }
    if ( ! defined $MIN_FDR ) {
        $MIN_FDR = 0.05;
        $logger->info("Setting --min_fd: $MIN_FDR");
    }

	# make sure required directories exist
	if ( ! -d $stat_dir ) { 
		pod2usage(-message => "ERROR: --stat_dir is not a directory\n\n",
					-exitval => 2); 
	}
    if ( ! -d $dafe_db ) { 
		pod2usage(-message => "ERROR: --dafe_db is not a directory\n\n",
					-exitval => 2); 
	}
    
    # check the numbers
    if ( ! looks_like_number($MIN_PERC) ) {
        pod2usage(-message => "ERROR: --min_perc must be a digit\n\n",
                  -exitval => 2);
    }
    if ( ! looks_like_number($MIN_FDR) ) {
        pod2usage(-message => "ERROR: --min_fdr must be digit\n\n",
                  -exitval=> 2);
    }
    
    
    $logger->info("Parameters");
    $logger->info("--stat_dir: $stat_dir");
    $logger->info("--dafe_db $dafe_db");
    $logger->info("--clstr_metadata: $clstr_metadata");
    $logger->info("--ref_meta_file: $ref_meta_file");
    $logger->info("--clstr_mapping_file: $clstr_mapping_file");
    $logger->info("--tag_file_name: $tag_file_name");
    $logger->info("--da_col_file: $da_col_file");
    $logger->info("--annote_file_name: $annote_file_name");
    $logger->info("--clstr_summary_file: $clstr_summary_file");
    $logger->info("--da_filt_col_file: $da_filt_col_file");
    $logger->info("--final_tbl_file: $final_tbl_file");
    $logger->info("--min_perc: $MIN_PERC");
    $logger->info("--min_fdr: $MIN_FDR");
	
	return 1;
}

sub _is_defined {
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

DAFE_stat_mk_DA_tbl.pl - makes a DA tbl from cluster output


=head1 VERSION

This documentation refers to version 0.0.1


=head1 SYNOPSIS

    DAFE_stat_mk_DA_tbl.pl
        --stat_dir clstr_ref_out
        --dafe_db my_dafe_db/
        --clstr_metadata clstr_metadata_file.txt
        [--ref_meta_flie eference_metadata_20160519.txt]
        [--clstr_mapping_file gene_to_cdhit_clstr_mapping_file.txt]
        [--annote_file_name all_annote.txt]
        [--tag_file_name tags_data.txt]
        [--da_col_file dafe_full_DA_col.txt]
        [--clstr_summary_file clstr_summary.txt]
        [--da_filt_col_file dafe_filt_DA_col.txt]
        [--final_tbl_file dafe_filt_tbl.txt]
        [--min_perc 2]
        [--min_fdr 0.05]
        [--help]
        [--man]
        [--debug]
        [--verbose]
        [--quiet]
        [--logfile logfile.log]

    --stat_dir             Dir with DAFE_stat output
    --dafe_db              Dir with dafe database
    --clstr_metadata       File with each cluster's metedata
    --ref_meta_file        File with reference genome metadata
    --clstr_mapping_file   File mapping cluster IDs to genes and genomes
    --annote_file          Name of annotation file in dafe db
    --tag_file_name        Name of file with tag info for each genome
    --da_col_file          Name of file to print full DA column info
    --clstr_summary_file   Name of file to print the cluster summary info
    --da_filt_col_file     Name of file to print filtered DA column info
    --final_tbl_file       Name of file to print the final filtered table
    --min_perc             Min percentage of genomes a cluster must be DA to keep it
    --min_fdr              Min FDR to call a cluster DA
    --help | -h     Prints USAGE statement
    --man           Prints the man page
    --debug	        Prints Log4perl DEBUG+ messages
    --verbose       Prints Log4perl INFO+ messages
    --quiet	        Suppress printing ERROR+ Log4perl messages
    --logfile       File to save Log4perl messages


=head1 ARGUMENTS
    
=head2 --stat_dir

Direcotory where outputs from running the DAFE stat edgeR_model.R script.
This directory should have a subdirectory for each genome.  Inside each
genome's directory should be a file (--tag_file_name) that has information
about each cluster that had enough counts to make a statistical inference
about its differential abundance (DA).
    
=head2 --dafe_db

Path to DAFE database.  I use the DAFE database to look at each genomes
annotation file to determine all the clusters that are present in a genome.

=head2 --clstr_metadata

Path to the cluster metadata file.  This file should have a line for each
cluster.  The first value in the line should be the cluster ID.  The
cluster ID is the only thing used in this script.

=head2 [--ref_meta_file]

File with reference genome metadata.  The only piece of metadata currently being
used in the "Label" column.  This column designates each genome as either PA,
NPA, or Soil.
DEFAULT: "/netscr/yourston/compMetaG/metadata/reference_metadata_20160519.txt"

=head2 [--clstr_mapping_file]

File mapping cluster IDs to genome and gene IDs.  It has three columns and no
headers.  The first column is the genomeID.  The second is the cluster ID.  And
the third is the gene ID.  This file is currently only used to map cluster IDs
to genomes ID.  This allows me to count the number of PA genomes in each cluster
by cross referencing with the data in --ref_meta_file.
DEFAULT: "/netscr/yourston/compMetaG/metadata/gene_to_cdhit_clstr_mapping_file.txt"

=head2 [--annote_file]

Name of annotation file in DAFE database.  This file is used to determine
which clusters are present in a genome.  I can't simply use the tag file
because not all clusters have enough reads to make it into the tag file.
So I would have no way of identifying clusters that are present in the
genome, but simply don't have mapped reads.
DEFAULT: all_annote.txt

=head2 [--tag_file_name]

The name of the tag file that should be found in the --stat_dir within
each genome's subdirectory.  This file should have a line for each
cluster that had enough reads to infer its differential abundance (DA)
status.  This file should be a tab delemited file with headers.  The
first column should be the cluster ID.  The second column should be the
log fold change which is used to determine the DA direction (ie -1 or 1).
The fifth column should be an FDR value which is used to determine if
the cluster id DA.
DEFAULT: tags_data.txt

=head2 [--da_col_file]

The name of the file where the full column of clusters is printed and
their corresponding differential abuance value (-3, -2, -1, 0, 1).
This file has each possible cluster in the analysis.  This file is
printed for each genome.
DEFAULT: dafe_full_DA_col.txt

=head2 [--clstr_summary_file]

The name of the file where the summary information for all the clusters
can be printed.  This will be a tab delemited file with 5 columns: 1)
cluster ID, 2) number of genomes contained in that cluster (ie a measure
of conservation), 3) the number of DA sequences in that cluster, 4) the
number of DA up sequences in that cluster, and 5) the number of DA down
sequences in that cluster.
DEFAULT: clstr_summary.txt

=head2 [--da_filt_col_file]

The name of the file where the DA information for only the clusters that
are interested is printed for each genome.  Interesting clusters are
determined after all the cluster information for each genome is read.
Because there are so many clusters only the most interested in kept in
this file.
DEFAULT: dafe_filt_DA_col.txt

=head2 [--final_tbl_file]

The name of the output file that will contain a table where rows are
only the clusters of interest and columns are each genome.  The values
are the DA values (-3,-2,-1,0,1).
DEFAULT: dafe_filt_tbl.txt

=head2 [--min_perc]

The min percentage of genomes where a cluster is called DA to be kept
during the filtering process.  I want to keep only the most interesting
cluster and I define these as the clusters that are DA in the most
genomes.
DEFAULT: 2

=head2 [--min_fdr]

The min false discovery rate to call a cluster DA
DEFAULT: 0.05
 
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

When I run the DAFE stat script edgeR_model.R using clusters as the grouping
factor (remember the clusters are made using uclust using all ~3800 genomes)
there are over 2 million clusters.  Making a table of DA values for each
genome for all 2 million clusters would take some very creative memory tricks
(ie use sparse table).  Alternatively, clusters that are not of interest can
be removed.  The primary purpose of this script is to remove clusters that
are not of interest and print a genome X cluster table for the remaining
clusters and their DA values in each genome.

The DA value code is as follows:
-3 -- not present in the genome
-2 -- present in the genome but has very few or no reads to map to it
-1 -- statistically DA down
0 -- statistically equal
1 -- statistically DA up

For the final cluster_summary table to make sense this scrtip should be ran on
the full database.  

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
