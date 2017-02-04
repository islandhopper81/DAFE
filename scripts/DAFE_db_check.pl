#!/usr/bin/env perl

# Checks a DAFE database for potential errors and tries to correct when possible

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Carp;
use Readonly;
use version; our $VERSION = qv('0.0.2');
use Log::Log4perl qw(:easy);
use Log::Log4perl::CommandLine qw(:all);
use FindBin qw($Bin);
use BioUtils::FastaIO;
use File::Temp qw(tempfile);

# Subroutines #
sub check_dir_exists;
sub check_required_files;
sub check_genome_fasta;
sub correct_genome_fasta;		
sub check_gene_faa;
sub correct_gene_faa;		
sub check_gene_fna;
sub correct_gene_fna;		
sub correct_fasta_IDs;		
sub check_gff;
sub check_for_bad_lines;
sub check_annote;

# Variables #
my ($dafe_db, $genome, $meta_data_file,
    $genome_fasta_ext, $gene_faa_ext, $gene_fna_ext, $gff_ext, $annote_file,
	$annote_method, $annote_jgi_exe, $annote_gbk_exe, $annote_gbk_params,
	$correct_gff_exe, $correct_gff_htseq_i, $correct_genome_exe,
    $help, $man);

my $options_okay = GetOptions (
    "dafe_db:s" => \$dafe_db,
	"genome:s" => \$genome,
    "meta_data_file:s" => \$meta_data_file,
    "genome_fasta_ext:s" => \$genome_fasta_ext,
    "gene_faa_ext:s" => \$gene_faa_ext,
    "gene_fna_ext:s" => \$gene_fna_ext,
    "gff_ext:s" => \$gff_ext,
	"annote_file:s" => \$annote_file,
	"annote_method:s" => \$annote_method,
	"annote_jgi_exe:s" => \$annote_jgi_exe,
	"annote_gbk_exe:s" => \$annote_gbk_exe,
	"annote_gbk_params:s" => \$annote_gbk_params,
	"correct_gff_exe:s" => \$correct_gff_exe,
	"correct_gff_htseq_i:s" => \$correct_gff_htseq_i,
	"correct_genome_exe:s" => \$correct_genome_exe,
    "help" => \$help,                  # flag
    "man" => \$man,                     # flag (print full man page)
);

# set up the logging environment
my $logger = get_logger();

# check for input errors
if ( $help ) { pod2usage(2) }
if ( $man ) { pod2usage(-verbose => 3) }
check_params();


########
# MAIN #
########
my $ids_aref;
if ( defined $genome ) {
	$ids_aref = [$genome];
}
else {
	# get the ref IDs in the meta_data_file
	$ids_aref = get_IDs($meta_data_file);
}

foreach my $id ( @{$ids_aref} ) {
    # check if a directory exists for this id
    check_dir_exists($id);
    
    # checks for all the required files
    check_required_files($id);
}


########
# Subs #
########
sub check_dir_exists {
    my ($id) = @_;
    
    if ( ! -d "$dafe_db/$id" ) {
        $logger->warn("$id -- Missing dir");
    }
    
    return 1;
}

sub check_required_files {
    my ($id) = @_;
	
	# NOTE that the order of these calls is important.  For example, check_gff
	# must come before check_annote because check_annote assumes a valid GFF
	# file.
    
    check_genome_fasta($id);
    check_gene_faa($id);
	check_gene_fna($id);
    check_gff($id);  # this must come before check_annote
	check_annote($id);
    
    # do I need these for anything?
    #check_gene_fna($id);
    
}

sub check_genome_fasta {
    my ($id) = @_;
    
    my $file = "$dafe_db/$id/$id" . "$genome_fasta_ext";
    if ( ! -e $file ) {
        $logger->warn("$id -- Missing genome fasta");
		return 0;
    }
    elsif ( ! -s $file ) {
        $logger->warn("$id -- Zero size genome fasta");
		return 0;
    }
	
	# check if I need to correct the names in the genome file
	open my $FNA, "<", $file
		or $logger->warn("Cannot open genome fasta: $file");
	
	# NOTE: this loop looks like it goes through each line, but notice the
	#		"last" statement at the end.  This effectively looks only at the
	#		first non-comment line.
	foreach my $line ( <$FNA> ) {
		chomp $line;
		
		if ( $line !~ m/^>/ ) {
			$logger->warn("This file( $file) may not be a fasta at line: $line");
		}
		
		# this looks for space in the name.  There should be none
		if ( $line =~ m/ / ) {
			correct_genome_fasta($file, $id);
		}
		
		# looks for the genome id in the name
		if ( $line !~ m/$id-\s+/ ) {
			correct_genome_fasta($file, $id);
		}
		
		last;
	}
	
	close($FNA);
    
    return($file);
}

sub correct_genome_fasta {
	my ($file, $id) = @_;
	
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
				
		if ( $seq->get_id() !~ m/$id-\S+/ ) {		
			$seq->set_header($id . "-" . $seq->get_id());		
		}		
				
		# print		
		$fasta_out->write_seq($seq);		
	}		
	
	# now move the tmp file to the original genome fasta		
	# NOTE: this overwrite the original fasta file		
	$logger->info("Moving tmp file ($filename) to overwite genome fasta ($file)");		
	`mv $filename $file`;		
	
  	return 1;
}

sub check_gene_faa {
    my ($id) = @_;
	
	# NOTE this file is no longer needed
	# but I'm keeping it because it might be useful in the future
    
    my $file = "$dafe_db/$id/$id" . "$gene_faa_ext";
    if ( ! -e $file ) {
        $logger->warn("$id -- Missing gene faa");
		return 0;
    }
    elsif ( ! -s $file ) {
        $logger->warn("$id -- Zero size gene faa");
		return 0;
    }
	
	# check that the gene ID has the patter: genomeID-geneID
	my $faa_in = BioUtils::FastaIO->new({stream_type => '<', file => $file});
	my $first_seq = $faa_in->get_next_seq();
	if ( $first_seq->get_id() !~ m/$id-\S+/ ) {
		$logger->info("Gene faa file needs correction: $file");
		correct_gene_faa($id, $file, $faa_in, $first_seq);
	}
    
    return($file);
}

sub correct_gene_faa {
	my ($id, $file, $faa_in, $first_seq) = @_;
	
	$logger->info("Correcting gene faa file");
	
	correct_fasta_IDs($id, $file, $faa_in, $first_seq);
	
	return 1;
}

sub check_gene_fna {
    my ($id) = @_;
    
    my $file = "$dafe_db/$id/$id" . "$gene_fna_ext";
    if ( ! -e $file ) {
        $logger->warn("$id -- Missing gene fna");
		return 0;
    }
    elsif ( ! -s $file ) {
        $logger->warn("$id -- Zero size gene fna");
		return 0;
    }
	
	# check that the gene ID has the patter: genomeID-geneID
	my $fna_in = BioUtils::FastaIO->new({stream_type => '<', file => $file});
	my $first_seq = $fna_in->get_next_seq();
	if ( $first_seq->get_id() !~ m/$id-\S+/ ) {
		$logger->info("Gene fna file needs correction: $file");
		correct_gene_faa($id, $file, $fna_in, $first_seq);
	}
    
    return($file);
}

sub correct_gene_fna {
	my ($id, $file, $fna_in, $first_seq) = @_;
	
	$logger->info("Correcting gene fna file");
	
	correct_fasta_IDs($id, $file, $fna_in, $first_seq);
	
	return 1;
}

sub correct_fasta_IDs {
	my ($id, $file, $in, $first_seq) = @_;
	
	# correct the fasta IDs by adding te genomeID to the beginning of the
	# geneID in the header
	
	# make a temp file to put the correct seqs
	my ($tmp_fh, $tmp_file) = tempfile();
	close($tmp_fh);
	
	# make a BioUtilsIO object that points to the temp file
	my $out = BioUtils::FastaIO->new({stream_type => '>', file => $tmp_file});
	
	# update and add the first sequence
	# here I assume that the first value in the header is the gene ID
	my $header = $first_seq->get_header();
	$first_seq->set_header("$id-$header");
	$out->write_seq($first_seq);
	
	# update the remaining sequences in the file
	while ( my $seq = $in->get_next_seq() ) {
		$header = $seq->get_header();
		$seq->set_header("$id-$header");
		$out->write_seq($seq);
	}
	
	# replace the original file with the tmp file
	`mv $tmp_file $file`;
	
	return 1;
}

sub check_gff {
    my ($id) = @_;
    
    my $file = "$dafe_db/$id/$id" . "$gff_ext";
    if ( ! -e $file ) {
        $logger->warn("$id -- Missing gff");
		return 0;
    }
    elsif ( ! -s $file ) {
        $logger->warn("$id -- Zero size gff");
		return 0;
    }
	
	# check if I need to correct the GFF file
	open my $GFF, "<", $file
		or $logger->warn("Cannot open gff: $file");
	
	# NOTE: this loop looks like it goes through each line, but notice the
	#		"last" statement at the end.  This essentially only looks only at the
	#		first non-comment line.
	foreach my $line ( <$GFF> ) {
		chomp $line;
		
		if ( $line =~ m/^#/ ) {next;}  # ignore comment lines
		
		my @vals = split(/\t/, $line);
		
		# look for the scaffold name to include the genome ID
		if ( $vals[0] !~ m/$id/ ) {
			correct_gff($file);
		}
		
		# look for bad strand field
		if ($vals[6] eq "+" or $vals[6] eq "-") {
			; # doesn't need correction
		}
		elsif ($vals[6] eq "-1") {
			# needs correction
			correct_gff($file);
		}
		elsif ($vals[6] eq "1") {
			# needs correction
			correct_gff($file);
		}
		
		# look for bad attribute tag fields
		if ( $vals[8] =~ m/name / ) { 
			correct_gff($file);
		}
		
		# look for the htseq_i field in the tags
		if ( $vals[8] !~ m/$correct_gff_htseq_i=/ ) {
			correct_gff($file);
		}
        
        # look for genome ID in the gene id's
        if ( $vals[scalar(@vals) - 1] !~ qr/$id/ ) {
            correct_gff($file);
        }
		
		last;
	}
	close($GFF);
	
	# check for the bad lines that are at the end of the file
	check_for_bad_lines($file, qr/\tCRISPR\t/);
	check_for_bad_lines($file, qr/\tdirect\t/);
	check_for_bad_lines($file, qr/\ttandem\t/);
	check_for_bad_lines($file, qr/\tinverted\t/);
	
	# check if there are double quotes (which) are an illegal character
	my $has_double_quotes = `grep '"' $file`;
	if ( $has_double_quotes ) {
		correct_gff($file);
	}
	
	# check if there are semi colons in the product tag information.
	# these are illegal characters because each tag is seperated by ";"
	my $has_semi_colon = `grep -P 'product=(?!.*(=)).*;' $file`;
	if ( $has_semi_colon ) {
		correct_gff($file);
	}
    
    return($file);
}

sub check_for_bad_lines {
	my ($file, $regex) = @_;
	
	# figure out how to pass in the regex so it can be included in the
	# grep statement and match statement
	
	my $direct_lines = `grep -P "$regex" $file`;
	my @lines = split(/\n/, $direct_lines);
	foreach my $line ( @lines ) {
		chomp $line;
		
		if ( $line !~ m/^#/i and $line =~ m/$regex/i ) {
			correct_gff($file);
			last;
		}
	}
}

sub correct_gff {
	my ($file) = @_;
	
	$logger->info("GFF file needs correction: $file");
	my $command = "perl $correct_gff_exe $file";
	$logger->info("Running correct_gff command: $command");
	`$command`;
	
	return 1;
}

sub check_annote {
	my ($id) = @_;
	
	my $file = "$dafe_db/$id/$annote_file";
	if ( ! -e $file ) {
		$logger->warn("$id -- Missing annotation file");
		make_annote($id, $file);
	}
	elsif ( ! -s $file ) {
		$logger->warn("$id -- Zero size annotation file");
		make_annote($id, $file);
	}
	
	return($file);
}

sub make_annote {
	my ($id, $file) = @_;
	
	if ( $annote_method =~ m/none/i ) {
		$logger->warn("Not trying to make an annotaiton file");
		return 0;  # Early return here so other logger messages are ignored
	}
	
	$logger->info("Trying to make annotaiton file");
	
	# There are two ways to make an annotation file.  Option number 1 is the
	# default and recommended method.
	# 1) combine the IMG downloaded files (cog.tab.txt, ko.tab.txt, etc.) by
	#	 running the make_isolate_annote_tables.pl script
	# 2) use the genbank file and blast to CDD database by running the
	#	 genbank_to_annotation_tbl.pl script
	
	# try either method starting with jgi and then trying gbk if needed
	if ( $annote_method =~ m/either/i ) {
		# first try using the jgi method to make the all_annote.txt file
		jgi_annotation($id);
		
		# if that didn't work (ie there is still no alL_annote.txt file
		# then try the gbk method
		if ( ! -s $file ) {
			# now try the gbk method
			$logger->info("JGI annotation method FAILED");
			$logger->info("Trying to annotate using gbk method");
			gbk_annotation($id);
		}
	}
	elsif ( $annote_method =~ m/jgi/i ) {
		# try only the jgi method
		jgi_annotation($id);
	}
	elsif ( $annote_method =~ m/gbk/i ) {
		# try only the gbk method
		gbk_annotation($id);
	}
	else {
		$logger->warn("Unknown --make_annote_method ($annote_method)");
	}
	
	# check to see if it worked
	if ( ! -s $file ) {
		$logger->warn("Cannot make $annote_file file for genome ($id)");
		return 0;
	}
	else {
		$logger->info("Successfully made $annote_file for genome ($id)");
	}
	
	return 1;
}

sub jgi_annotation {
	my ($id) = @_;
	
	# DAFE_db_make_annote_tables.pl
	my $jgi_command = "perl $annote_jgi_exe ";
	$jgi_command .= "--dafe_db $dafe_db ";
	$jgi_command .= "--genome $id ";
	$jgi_command .= "--annote_file_name $annote_file ";
	$jgi_command .= "--gene_id_col $correct_gff_htseq_i ";
	$jgi_command .= "--verbose";
	
	$logger->info("Trying to make annote file using jgi command: $jgi_command");
	`$jgi_command`;
	
	return 1;
}

sub gbk_annotation {
	my ($id) = @_;
	
	# look for the gbk file
	my $gbk_file = "$dafe_db/$id/$id\.gbk";
	if ( ! -s $gbk_file ) {
		$logger->warn("No gkb file for genome($id)");
		$logger->warn("Cannot annotate using gbk method");
	}
	else {
		# the gbk file does exist
		my $gbk_command = "perl $annote_gbk_exe ";
		$gbk_command .= " --gbk_file $gbk_file";
		$gbk_command .= " --out_dir $dafe_db/$id/";
		$gbk_command .= " --prefix $id";
		$gbk_command .= " --faa_suffix $gene_faa_ext";
		$gbk_command .= " --tbl_file $annote_file";
		
		if ( _is_defined($annote_gbk_params) ) {
			$gbk_command .= " $annote_gbk_params ";
		}
		
		# look for the faa file.  if it is there skip making a new faa
		if ( check_gene_faa($id) != 0 ) {
			$gbk_command .= " --skip_faa ";
		}
		
		$gbk_command .= " --verbose ";
		
		$logger->info("Trying to make annote file using gkb command: $gbk_command");
		`$gbk_command`;
	}
	
	return 1;
}

sub get_IDs {
    my ($file) = @_;
    
    my @ids = ();
    
    # the ID should be the first column in the file
    open my $IN, "<", $file or
        $logger->fatal("Cannot open meta_data_file ($file)");
    
    my $first = 1;
    foreach my $line ( <$IN> ) {
        chomp $line;
        
        # skip the header line.
        if ( $first == 1 ) {
            $first = 0;
            next;
        }
        
        my @vals = split(/\t/, $line);
        push @ids, $vals[0];
    }
    
    close($IN);
    
    return \@ids;
}

sub check_params {
	# check for required variables
	if ( ! defined $dafe_db) { 
		pod2usage(-message => "ERROR: required --dafe_db not defined\n\n",
					-exitval => 2); 
	}
	if ( ! defined $meta_data_file ) {
		pod2usage(-message => "ERROR: required --meta_data_file not defined\n\n",
					-exitval => 2);
	}
	
	# check the genome parameter
	# no checks required as it is completely optional

	# make sure required files are non-empty
	if ( ! -e $meta_data_file ) { 
		pod2usage(-message => "ERROR: --meta_data_file $meta_data_file is an empty file\n\n",
					-exitval => 2);
	}

	# make sure required directories exist
	if ( ! -d $dafe_db ) { 
		pod2usage(-message => "ERROR: --dafe_db is not a directory\n\n",
					-exitval => 2); 
	}
    
    # set defaults for file extensions
    if ( ! defined $genome_fasta_ext ) {
        $genome_fasta_ext = ".fna";
    }
    if ( ! defined $gene_faa_ext ) {
        $gene_faa_ext = ".gene.faa";
    }
    if ( ! defined $gene_fna_ext ) {
        $gene_fna_ext = ".gene.fna";
    }
    if ( ! defined $gff_ext ) {
        $gff_ext = ".gff";
    }
	
	# set annotation file default
	if ( ! defined $annote_file ) {
		$annote_file = "all_annote.txt";
	}
	
	# check annote_method
	if ( ! defined $annote_method ) {
		$annote_method = "jgi";
		$logger->info("Setting --annote_method to $annote_method");
	}
	if ( $annote_method !~ m/either|none|jgi|gbk/i ) {
		pod2usage(-message => "ERROR: --annote_method is not an accepted value\n\n",
					-exitval => 2); 
	}
	
	# check annote_jgi_exe
	if ( ! defined $annote_jgi_exe ) {
		$annote_jgi_exe = "$Bin/DAFE_db_make_annote_tables.pl";
		$logger->info("Setting --annote_jgi_exe to $annote_jgi_exe");
	}
	
	# check annote_gkb_exe
	if ( ! defined $annote_gbk_exe ) {
		$annote_gbk_exe = "genbank_to_annotation_tbl.pl";
		$logger->info("Setting --annote_gbk_exe to $annote_gbk_exe");
	}
	
	# check annote_gbk_params
	# no checks required because this is completely optional
	
	# check correct_gff_exe
	if ( ! defined $correct_gff_exe ) {
		$correct_gff_exe = "$Bin/DAFE_db_correct_gff.pl";
		$logger->info("Setting --correct_gff_exe $correct_gff_exe");
	}
	
	# check correct_gff_htseq_i
	if ( ! defined $correct_gff_htseq_i ) {
		$correct_gff_htseq_i = "ID";
		$logger->info("Setting --correct_gff_htseq_i $correct_gff_htseq_i");
	}
	
	# check correct_genome_exe
	if ( ! defined $correct_genome_exe ) {
		$correct_genome_exe = "$Bin/DAFE_db_correct_genome_fasta.pl";
		$logger->info("Setting --correct_genome_exe: $correct_genome_exe");
	}
    
    $logger->info("Parameters");
    $logger->info("--dafe_db: $dafe_db");
	$logger->info("--genome $genome") if (defined $genome);
    $logger->info("--meta_data_file: $meta_data_file");
    $logger->info("--genome_fasta_ext: $genome_fasta_ext");
    $logger->info("--gene_faa_ext: $gene_faa_ext");
    $logger->info("--gene_fna_ext: $gene_fna_ext");
    $logger->info("--gff_ext: $gff_ext");
	$logger->info("--annote_file: $annote_file");
	$logger->info("--annote_method: $annote_method");
	$logger->info("--annote_jgi_exe: $annote_jgi_exe");
	$logger->info("--annote_gbk_exe: $annote_gbk_exe");
	if (defined $annote_gbk_params) {
		$logger->info("--annote_gbk_params: $annote_gbk_params")
	}
	$logger->info("--correct_gff_exe: $correct_gff_exe");
	$logger->info("--correct_gff_htseq_i: $correct_gff_htseq_i");
	$logger->info("--correct_genome_exe: $correct_genome_exe");
	
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

DAFE_db_check.pl.pl - Checks for required files in DAFE database


=head1 VERSION

This documentation refers to version 0.0.2


=head1 SYNOPSIS

    DAFE_db_check.pl
        --dafe_db my_dafe_db/
        [--genome 646564591]
        --meta_data_file reference_metadata_20160412.txt
        [--genome_fasta_ext ".fna"]
        [--gene_faa_ext ".gene.faa"]
        [--gene_fna_ext ".gene.fna"] -- DEPRECIATED!
        [--gff_ext ".gff"]
        [--annote_file "all_annote.txt"]
        [--annote_method either | none | jgi | gbk]
        [--annote_jgi_exe DAFE_db_make_annote_tables.pl]
        [--annote_gbk_exe genbank_to_annotation_tbl.pl]
        [--annote_gbk_params "--cog_pfam_tbl cog_pfam_correspondance_tbl.txt"]
        [--correct_gff_exe DAFE_db_correct_gff.pl]
        [--correct_gff_htseq_i "ID"]
        [--correct_genome_exe DAFE_db_correct_genome_fasta.pl]
        
        [--help]
        [--man]
        [--debug]
        [--verbose]
        [--quiet]
        [--logfile logfile.log]

    --dafe_db               Path to DAFE database directory
    [--genome]              Genome ID to run checks on
    --meta_data_file        Path to txt file with all reference metadata
    [--genome_fasta_ext]    File extension for genome fasta in DAFE db
    [--gene_faa_ext]        File extension for gene faa file in DAFE db
    [--gene_fna_ext]        DEPRECIATED!
    [--gff_ext]             File extension for gff file in DAFE db
    [--annote_file]         Name of annotation txt file in DAFE db
    [--annote_method]       Name of method to use for making annotation file
    [--annote_jgi_exe]      Path to Perl script for the JGI annotation method
    [--annote_gbk_exe]      Path to Perl script for genbank annotation method
    [--annote_gbk_params]   Extra parameters to pass to --annote_gbk-exe
    [--correct_gff_exe]     Path to Perl script for correcting GFF files
    [--correct_gff_htseq_i] Tag name specified in the htseq-count -i parameter
    [--correct_genome_exe]  Path to Perl script for correcting the genome fastas
    --help                  Prints USAGE statement
    --man                   Prints the man page
    --debug                 Prints Log4perl DEBUG+ messages
    --verbose               Prints Log4perl INFO+ messages
    --quiet                 Suppress printing ERROR+ Log4perl messages
    --logfile               File to save Log4perl messages


=head1 ARGUMENTS
    
=head2 --dafe_db

Path to DAFE database
    
=head2 [--genome]

Genome ID to run checks on
DEFAULT: ""

=head2 --meta_data_file

Path to metadata table with information about each reference in the DAFE
database.  This script uses the table to get all the names of the reference.  So
if there is a reference in the DAFE database structure that is not in the
metadata file it will not be checked.  The other values in this script are not
currently being used.

=head2 [--genome_fasta_ext]

File extension for the genome fasta file.  This file is required in the DAFE
database because it is what the metagenome reads are mapped to.  The headers of
this file should stop after the first space.  Otherwise the names of the
scaffolds as they appear in the BAM file created by BBMap will not match the
sequence headers in the file.  This will cause htseq-count to classify all the
mapped reads to unknown features.
DEFAULT: ".fna"

=head2 [--gene_faa_ext]

File extension for the gene faa file.  This file is only used if there is no
annotation file associated with a reference.  When using the gbk method this
file is used to BLAST against a database of annotated proteins.  This file was
previously used to create the geneID_name_tbl.txt file.  This functionality was
replace by using the GFF file to create the gff_annote.tab.txt file.
DEFAULT: ".gene.faa"

=head2 [--gene_fna_ext]

DEPRECIATED!

File extension for the gene fna file.  This file was previously used to create
the geneID_name_tbl.txt file.  This functionality was replace by using the GFF
file to create the gff_annote.tab.txt file.
DEFAULT: ".gene.fna"

=head2 [--gff_ext]

File extention for the gff file.  The gff file is used by htseq-count to count
the number of reads assigned to each genome feature (ie gene).  
DEFAULT: ".gff"

=head2 [--annote_file]

Name of the annotation file.  This file is required for downstream analysis
after the mapping has all be completed.  It is used to group genes into COGs or
other functional groups contained in the annotation file.  If this file does not
exist the script attempts to create it using one of two methods.  The first
method combines files downloaded from IMG that contain all the annotations for
each gene (ie cog.tab.txt).  The second method uses a genbank file to create
a gene faa file and then BLAST those genes against a database of annotated
proteins.
DEFAULT: "all_annote.txt"

=head2 [--annote_method]

Name of method used for making annotation file.  Allowed methods include:
either, none, jgi, gbk.  The "either" method first tries "jgi" (which is the
recommended method) and then tries "gbk" if the "jgi" method fails because the
required files are not present.  The "jgi" method requires all the annotation
files to be inseperate files (ie cog.tab.txt, pfam.tab.txt, etc).  It also
requires a mapping file of gene IDs to gene names.  If the mapping file is not
present it is created using either the gene.faa or gene.fna file.  In those
files the headers must have the format ">geneID geneName ...." where the first
value is the gene ID and the next value is the gene name.  The "gbk" method uses
the genbank file to create a gene.faa file (if the gene.faa does not already
exist).  It then uses that file to BLAST against an annotation database
(ie CDS).  The database can be specified in the --annote_gbk_params string.  If
the attempt to make the annotation file succeeds a new file will appear in each
reference directory called --annote_file.
DEFAULT: jgi

=head2 [--annote_jgi_exe]

Path to Perl script for JGI annotation method.  I wrote a script to do the JGI
annotation method.  It is called DAFE_db_make_annote_tables.pl and is likely
located in the same location as this script (DAFE_db_check.pl).  The default
automatically looks for the DAFE_db_make_annote_tables.pl script in the same
directory as where the current script is ran from.
DEFAULT: $Bin/DAFE_db_make_annote_tables.pl

=head2 [--annote_gbk_exe]

Path to Perl script for genbank annotation method.  I have a script that I include
in my SeqTools called genbank_to_annotation_tbl.pl.  This should be the script
passed in this parameter.  The default will attempt to run
genbank_to_annotation_tbl.pl assuming it is located somewhere defined in the $PATH
environment variable.
DEFAULT: genbank_to_annotation_tbl.pl

=head2 [--annote_gbk_params]

A string of parameters that are passed to genbank_to_annotation_tbl.pl.
Parameters that are automatically set in DAFE_db_check.pl include --gbk_file,
--out_dir, --prefix, --faa_suffix, and --tbl_file.  Other common parameters that
can be passed in this situation include --cog_pfam_tbl, --gbk_id_tag --blast_db,
--blast_eval, --min_cog_pid, and --min_cog_aln_frac. See the documentation for
genbank_to_annotation_tbl.pl for details about each of these parameters and other
parameter that can be passed to the genbnak_to_annotation_tbl.pl script.
DEFAULT: ""

=head2 [--correct_gff_exe]

Path to Perl script for correcting GFF files.  These corrections include
changing strand designations from "-1" and "1" to "-" and "+".  It also comments
the CRISPR line.  The default automatically looks for the DAFE_db_correct_gff.pl
script in the same directory as where the current script is ran from.
DEFAULT: $Bin/DAFE_db_correct_gff.pl

=head2 [--correct_gff_htseq_i]

Tag name specified when running htseq-count as the -i parameter.  This parameter
tells htseq-count what to name the gene in it's output file.  Each gene in the
GFF file should have one of these tags and a corresponding unique attribute
specify the gene indentifier.
DEFAULT: "ID"

=head2 [--correct_genome_exe]

Path to Perl script for correcting genome fasta files.  The headers in the fasta
cannot have any spaces.  Mostly they just need to match what is in the BAM file
which cuts off everything after the first space (for BBMap).  Note that the
original fasta headers are lost (ie unrecoverable).  Make sure there is nothing
in the fasta headers that is important!  The default automatically looks for the
DAFE_db_correct_genome_fasta.pl script in the same directory as where the
current script is ran from.
DEFAULT: $Bin/DAFE_db_correct_genome_fasta.pl
 
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

The purpose of this script is to check the DAFE database by looking for all the
required files.  For some files, if they do not exist, DAFE_db_check.pl will
create them.  Also, some files may be formatted incorrectly (ie gff) in which
case DAFE_db_check.pl attempts to convert them to the correct format.

Here are the files that are required:

[ID].gff -- required for running htseq-count
[ID].fna -- required for running bbmap
all_annote.txt -- required for running the DAFE_stat* scripts

These files are optional in some cases but required in others:

[ID].gene.fna -- required for creating the all_annote.txt file using blast method
[ID].gene.faa -- required for creating the all_annote.txt file using blast method

If these files are missing this script attempts to create them:

all_annote.txt

When looking for gff formatting mistakes DAFE_db_check.pl does the following:

1. checks strand characters are anything other than "+" and "-".
2. checks for lines added by some JGI program that tend to be incorrecty
formatted (CRISPR, direct, tantdem, inverted).  
3. checks for double quotes from anywhere in the file
4. checks for semicolon characters in the "product" tag value

If any of the above are true DAFE_db_check.pl runs the --correct_gff_exe program
(default is DAFE_db_correct_gff.pl).  Note that the GFF file must contain a gene
indentifier as one of the tag values (ie. ID or locus_tag).  Each gene in all
the gff files must have this tag and a unique value used to indentify the gene.

Also, when this script checks the GFF and genome fasta files it looks at the
first line of the files to see if they need to be corrected.  If there are
mistakes in later lines these will be missed.  This is a limitation, but I don't
want this script to have to look at every line looking for mistakes.

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
FindBin qw($Bin)
BioUtils::FastaIO
File::Temp qw(tempfile)


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
