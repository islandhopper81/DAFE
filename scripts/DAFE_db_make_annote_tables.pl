#!/usr/bin/env perl

# Given each of the annotation files for cog, ko, etc make all_annote.txt

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Carp;
use Readonly;
use version; our $VERSION = qv('0.0.1');
use Log::Log4perl qw(:easy);
use Log::Log4perl::CommandLine qw(:all);
use FindBin qw($Bin);

# Subroutines #
sub check_params;

# Variables #
my ($dafe_db, $genome_id, $gene_file_ext, $gff_file_ext, $mk_gff_annote_tbl_exe,
	$annote_file_name, $gene_id_col, $help, $man);
my @gff_file_ext_arr = (".gff", ".gff2", ".gff3", ".GFF");
my @gene_file_ext_arr = (".genes.fna", ".genes.faa", ".gene.faa", ".gene.fna");
my %annotes_hash = (cog => ".cog.tab.txt",
					ko => ".ko.tab.txt",
					pfam => ".pfam.tab.txt",
					tigrfam => ".tigrfam.tab.txt",
					clstr => ".clstr.tab.txt",
					go => ".go.tab.txt",
					kegg => ".kegg.tab.txt",
					kog => ".kog.tab.txt"
					);

my $options_okay = GetOptions (
    "dafe_db|d:s" => \$dafe_db,
    "genome|g:s" => \$genome_id,
    "mk_gff_annote_tbl_exe:s" => \$mk_gff_annote_tbl_exe,
	"gff_file_ext:s" => \$gff_file_ext,
    "gene_file_ext:s" => \$gene_file_ext,
	"annote_file_name:s" => \$annote_file_name,
	"gene_id_col:s" => \$gene_id_col,
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
# get a list of genomes to make the all_annote.txt file for
$logger->debug("Reading in genomes");
my @genome_ids = ();
if ( defined $genome_id ) {
    push @genome_ids, $genome_id;
}
else {
    opendir my $DIR, $dafe_db or
        $logger->logdie("Cannot open --dafe_db ($dafe_db)");
    
    while ( my $id = readdir($DIR) ) {
        if ( $id =~ m/^\./ ) {next;} #skip hidden dirs
		if ( ! -d "$dafe_db/$id" ) { next; } # skip non-directories
        push @genome_ids, $id;
    }
    
    closedir($DIR);
}

# foreach genome save all the annotation data
$logger->debug("Operating on each genome");
foreach my $genome_id ( @genome_ids ) {
    $logger->debug("Starting genome: $genome_id");
    
    # create a table for storing all the info from each file
    my %tbl = ();
    
    # check for the gff_annote.tab.txt file
    my $gff_annote_file = "$dafe_db/$genome_id/gff_annote.tab.txt";
    if (! -s $gff_annote_file) {
        my $msg = "No gff_annote.tab.txt file for genome ($genome_id).";
        $logger->warn($msg);
        
        my $unix_exit_status = make_gff_annote_tbl($genome_id, $gene_id_col);
		
        if ( $unix_exit_status > 0 ) {  # remember in unix 0 means success
            $logger->warn("Cannot make gff_annote.tab.txt.  Skipping $genome_id");
            next;
        }
        else {
            $logger->info("gff_annote.tab.txt file successfully made for genome ($genome_id)");
        }
    }
    
    # save the gff_annote.tab.txt file info
	# this must be outside the loop below because I make the gff_annote_file
	# above.
    add_to_tbl(\%tbl, $gff_annote_file, "gff_annote", $genome_id);
    
	# save all the annotation info in the table
	foreach my $ann ( keys %annotes_hash ) {
		my $exten = $annotes_hash{$ann};
		my $file = "$dafe_db/$genome_id/$genome_id$exten";
		add_to_tbl(\%tbl, $file, $ann, $genome_id);
	}
    
    # print the all_annote.txt table file.
	# I had to add the gff_annote_file as a sperate variable because I made
	# that variable in the above code.  I need the headers from that file
	# when I print the table.
    print_tbl(\%tbl, $genome_id, $gff_annote_file);
}

#2582580847

########
# Subs #
########
sub make_gff_annote_tbl {
    my ($genome_id, $gene_id_col) = @_;
    
    $logger->warn("Trying to make the gff_annote.tab.txt for genome ($genome_id)");
    
    # do some checks for the $mk_gff_annote_tbl_exe parameter
    if ( ! defined $mk_gff_annote_tbl_exe ) {
        $mk_gff_annote_tbl_exe = "$Bin/DAFE_db_make_gff_annote_tbl.pl";
        $logger->warn("setting --mk_ID_gff_annote_exe to $mk_gff_annote_tbl_exe");
    }
    if ( ! -e $mk_gff_annote_tbl_exe ) {
        $logger->warn("--mk_gff_annote_tbl_exe ($mk_gff_annote_tbl_exe) does not exist");
        return(1);  # return fail signal (unix exist status)
    }
    
    my $suffix = find_gff_file_suffix($genome_id);
    if ( $suffix eq 0 ) { return(1); }  # exit status "fail"
    
    my $command  = "perl $mk_gff_annote_tbl_exe ";
    $command .= "--dafe_db $dafe_db ";
    $command .= "--gff_file_suffix \"$suffix\" ";
	$command .= "--gene_id_col $gene_id_col ";
	$command .= "--genome $genome_id ";
    $command .= "--verbose";
    
    $logger->info("Running command: $command");
    my $exit_status = system($command);
    
    return($exit_status);
}

sub find_gff_file_suffix {
	my ($id) = @_;
	
	if ( defined $gff_file_ext ) {
		return($gff_file_ext);
	}
	
	$logger->info("Looking for gff file");
	foreach my $suffix ( @gff_file_ext_arr ) {
		my $file_name = "$dafe_db/$id/$id" . $suffix;
		$logger->debug("Looking for potential gff file: $file_name");
		if ( -s $file_name ) {
			$logger->info("Found gff file ($file_name) for genome ($id)");
			return($suffix);
		}
	}
	
	$logger->warn("GFF file for genome ($id) NOT found!");
    
    return(0);  # fail!
}

# DEPRECIATED
sub find_gene_file_suffix {
    my ($id) = @_;
    
    if ( defined $gene_file_ext ) {
        return($gene_file_ext);
    }
    
    $logger->info("Looking for gene file");
    foreach my $suffix ( @gene_file_ext_arr ) {
        my $file_name = "$dafe_db/$id/$id" . $suffix;
        $logger->debug("Looing for potential gene file: $file_name");
        if ( -s $file_name ) {
            $logger->info("Found gene file ($file_name) for genome ($id)");
            return($suffix);
        }
    }
    
    $logger->warn("Gene file for genome ($id) NOT found!");
    
    return(0);  # fail!
}

sub add_to_tbl {
    my ($tbl_href, $file, $type, $genome_id) = @_;
    
    if ( ! -s $file) {
        my $msg = "No $type annotation file ($file) for genome ($genome_id).  ";
        $msg .= "$type set to NA for genome ($genome_id)";
        $logger->warn($msg);
        return 0; 
    }
    
    # save the file info
    open my $IN, "<", $file or
        $logger->fatal("Cannot open annotation file: $file");
    
	my $first = 1;
    my @vals = ();
    foreach my $line ( <$IN> ) {
        chomp $line;
		
		# skip the header line of only the gff_annote file
		# it is the only one that is allowed to have headers.
		if ( $first and $type eq "gff_annote" ) { $first = 0; next; }
		
		# skip the header line if the line starts with the gene_id_col varaible
		if ( $first and $line =~ m/$gene_id_col/i ) { $first = 0; next; }
        
		# this code assumes that the first column is a unique value for each
		# line (ie gene) in the file.  This assumption is violated for the
		# fungal gff_annote.tab.txt files.
        
		# get the values in the line.
		@vals = split(/\t/, $line);
		
		# if there are more than 2 columns (ie 2 values in @vals) then there
		# must be a header line which will need to be skipped
		# UPDATE: I'm only using two column files now so the line below is
		# irrelevant
		#if ( $first and scalar @vals > 2 ) { $first = 0; next; }
		
		# get the gene_id value which MUST be in the first column
        my $gene_id = shift @vals;
        if ($type ne "gff_annote") {
            $gene_id = $gene_id;
        }
		#my $annote_id = shift @vals;
		
        $tbl_href->{$gene_id}{$type} = join("\t", @vals);
    }
    
    close($IN);
    
    return 1;
}

sub print_tbl {
    my ($tbl_href, $genome_id, $gff_annote_file) = @_;
    
    # open the output file
    my $out_file = "$dafe_db/$genome_id/$annote_file_name";
    open my $OUT, ">", $out_file or
        $logger->fatal("Cannot print to otuput file: $out_file");
    
	$logger->debug("Printing $annote_file_name file: $out_file");
	
	# get the headers from the gff_annote.tab.txt file
	my $gff_annote_file_headers = `head -n 1 $gff_annote_file`;
	chomp $gff_annote_file_headers;
	my @headers = split(/\t/, $gff_annote_file_headers);
	
	# add the headers from all the annotation columns in the tbl
	push @headers, keys %annotes_hash;
	
    # print the headers to the output file
    print $OUT join("\t", @headers), "\n";
    
	# print the table values
    foreach my $gene_id ( keys %{$tbl_href} ) {
        $logger->debug("Getting info for gene_id: $gene_id");
		print $OUT $gene_id, "\t";
        
		# print all the values from the gff_annote file
        if ( defined $tbl_href->{$gene_id}{"gff_annote"} ) {
            print $OUT $tbl_href->{$gene_id}{"gff_annote"};
        }
        else {
            $logger->warn("Missing gene name for gene ($gene_id) in genome ($genome_id)");
			print $OUT "NA";
        }
		print $OUT "\t";
		
		# print the values from the other annotation files
		my @ann_vals = ();
		foreach my $ann ( keys %annotes_hash ) {
			if ( defined $tbl_href->{$gene_id}{$ann} ) {
				push @ann_vals, remove_white_space($tbl_href->{$gene_id}{$ann});
			}
			else { push @ann_vals, "NA"; }
		}
        my $ann_str = join("\t", @ann_vals);
		print $OUT $ann_str, "\n";
    }
    
    close($OUT);
}

sub remove_white_space {
    my ($val) = @_;
    
    if ( $val =~ m/^(\S*)\s*/ ) {
        return $1
    }
    
    return($val);
}

sub check_params {
	# make sure required directories exist
    if ( ! defined $dafe_db ) {
        pod2usage(-message => "ERROR: --dafe_db is not defined\n\n",
					-exitval => 2); 
    }
	if ( ! -d $dafe_db ) { 
		pod2usage(-message => "ERROR: --dafe_db is not a directory\n\n",
					-exitval => 2); 
	}
    
    # check mk_gff_annote_tbl_exe
    # I actually do these checks in the make_gff_annote_tbl subroutine because
    # I only care about this parameter if the gff_annote_tbl needs to be made.
    
    
	# check gff_file_ext
    # this is a little different test.  if it is defined add it to the
    # gff_file_ext_arr and remove the default values
    if ( defined $gff_file_ext ) {
        @gff_file_ext_arr = ($gff_file_ext);
    }
	
    # check gene_file_ext
	# DEPRECIATED
    # this is a little different test.  if it is defined add it to the
    # gene_file_ext_arr and remove the default values
    if ( defined $gene_file_ext ) {
        @gene_file_ext_arr = ($gene_file_ext);
    }
	
	# check the annote_file_name
	if ( ! defined $annote_file_name ) {
		$annote_file_name = "all_annote.txt";
		$logger->info("Setting --annote_file_name to $annote_file_name");
	}
    
    $logger->info("--dafe_db $dafe_db");
    $logger->info("--genome $genome_id") if (defined $genome_id);
    $logger->info("--mk_gff_annote_tbl_exe $mk_gff_annote_tbl_exe")
		if (defined $mk_gff_annote_tbl_exe);
	$logger->info("--gff_file_ext $gff_file_ext") if (defined $gff_file_ext);
    $logger->info("--gene_file_ext $gene_file_ext") if (defined $gene_file_ext);
	$logger->info("--annote_file_name $annote_file_name");
	
	return 1;
}


__END__

# POD

=head1 NAME

make_isolate_annote_tables.pl - make all_annote.txt file(s)


=head1 VERSION

This documentation refers to version 0.0.1


=head1 SYNOPSIS

    make_isolate_annote_tables.pl
        --dafe_db dafe_db/
        [--genome 2623620901]
        [--gff_file_ext ".gff"]
        [--gene_file_ext ".gene.faa"]  DEPRECIATED!
        [--annote_file_name "all_annote.txt"]
        [--mk_gff_annote_tbl_exe DAFE_db_make_gff_annote_tbl.pl]
        [--gene_id_col "ID"]
        
        [--help]
        [--man]
        [--debug]
        [--verbose]
        [--quiet]
        [--logfile logfile.log]

    --dafe_db            	  Path to DAFE database directory
    [--genome]     	          Genome ID for only doing one genome
    [--gff_file_ext]          GFF file extension
    [--gene_file_ext]   	  DEPRECIATED
    [--annote_file_name]	  The name to use for the annotation output file
    [--mk_gff_annote_tbl_exe] Perl executable for making gff_annote_tab.txt
    [--gene_id_col]           Name of GFF attribute to use as the uniq gene ID
    --help | -h         	  Prints USAGE statement
    --man               	  Prints the man page
    --debug	            	  Prints Log4perl DEBUG+ messages
    --verbose           	  Prints Log4perl INFO+ messages
    --quiet	            	  Suppress printing ERROR+ Log4perl messages
    --logfile           	  File to save Log4perl messages


=head1 ARGUMENTS
    
=head2 --dafe_db

Path to an DAFE database directory.  The default is to make the all_annote.txt
file for all the genomes in the DAFE database.  If --genome is use only that
genome is done.

=head2 [--genome | -g]

The genome ID for which to make the all_annote.txt file.  The ID must match the
directory name in the DAFE database.  If --genome is not provided this script a
gff_annote.tab.txt file is created for each dir in the DAFE database.

=head2 [--gff_file_ext]

The file extension for the GFF file.  The GFF file is used to create the
gff_annote.tab.txt file.  If not value is given several extensions will be
searched for (see default values)
DEFAULT: (".gff", ".gff2", ".gff3", ".GFF")

=head2 [--gene_file_ext]

DEPRECIATED!

The file extension for the gene file.  The gene file is used to create the
geneID_namt_tbl file if it does not already exist.  If no value is given
several extensions will be searched for (see default values)
DEFAULT: (".genes.fna", "genes.faa", "gene.faa", "gene.fna")

=head2 [--annote_file_name]

The name to use for the output annotation file.
DEFAULT: "all_annote.txt"

=head2 [--mk_gff_annote_tbl_exe]

Perl executable for makeing the gff_annote.tab.txt file.
DEFAULT: $Bin/DAFE_db_make_gff_annote_tbl.pl

=head2 [--gene_id_col]

Name of GFF attribute to use as the unique gene identifier.  This attribute
will be the first column in the output all_annote.txt file.  It should also
match the htseq-count -i parameter and htseq_i paramters in the other DAFE suite
scripts.
 
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

Given each of the annotation files for cog, ko, etc make all_annote.txt.

These are the files that MUST be in the DAFE database at before running this
script:

gff_annote.tab.txt

These are the optional files that will be used to create the all_annote.txt
file:

.ko.tab.txt
.cog.tab.txt
.ipr.tab.txt
.tigrfam
.pfam

This script is probably incomplete because it doesn't include the IDs and names
of the annotations.  It only includes the IDs.  I think this all that I
currently use, but that might change.  


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
