#!/usr/bin/env perl

# adds the cog group to the all_annote.txt file

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
use Scalar::Util qw(looks_like_number);

# Subroutines #
sub check_params;
sub _check_seed;
sub _check_method;
sub _check_over;

# Variables #
my ($db, $all_annote_file, $cog_grp_file, $names_file, $seed, $method, $over_write, $help, $man);

my $options_okay = GetOptions (
    "db|d:s" => \$db,
    "all_annote_file|f:s" => \$all_annote_file,
	"cog_grp_file|c:s" => \$cog_grp_file,
	"names_file|n:s" => \$names_file,
	"seed|s:s" => \$seed,
	"method|m:s" => \$method,
	"over:s" => \$over_write,
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
# read in the cog_grp_file
# this file should match each COG to a cog group.  It will be the
# lookup table.
my $cog_to_grp_tbl = Table->new();
$cog_to_grp_tbl->load_from_file($cog_grp_file);


# load in all the genomes
open my $NAME, "<", $names_file or
	$logger->logdie("Cannot open --names_file ($names_file)\n");

my @names = ();
foreach my $line ( <$NAME> ) {
	chomp $line;
	
	push @names, $line;
}

close($NAME);


# foreach name find it's annotation file, read in the file, and
# edit it by adding the cog grp
my $cog_id;
my $cog_grp;
foreach my $n ( @names ) {
	$logger->debug("genome: $n\n");
	my $ann_file = "$db/$n/$all_annote_file";

	if ( ! -f $ann_file ) {
		$logger->warn("Cannot find annotation file ($ann_file)\n");
		next;
	}

	my $ann_tbl = Table->new();
	$ann_tbl->load_from_file($ann_file);

	# if the annotation table already has the column cog_grp
	# the use must pass the flag "--over".  
	if ( $ann_tbl->has_col("cog_grp") ) {
		if ( $over_write == 1 ) {
			$ann_tbl->drop_col("cog_grp");
		}
		else {
			my $msg = "Annotation file ($ann_file) already has cog_grp column\n";
			$msg .= "You can set \"--over T\" to overwrite the existing column.";
			$msg .= "  However, the existing column will be PERMANENTLY deleted\n";
			$logger->logdie($msg);
		}
	}

	my @cog_grps = ();
	foreach my $gene_id ( @{$ann_tbl->get_row_names()} ) {
		$cog_id = $ann_tbl->get_value_at($gene_id, "cog");

		if ( $cog_id eq "NA" ) {
			push @cog_grps, "NA";
			next;
		}
		
		eval {
			$cog_grp = $cog_to_grp_tbl->get_value_at($cog_id, "func");
		};
		if ( my $e = MyX::Table::Row::UndefName->caught() ) {
       		# Do something to handle the exception like print an error message
			$logger->warn("COG id not found: $cog_id\n");
			$cog_grp = "NA";
    	}

		$cog_grp = get_cog_grp($cog_grp, $method, $seed);
		
		push @cog_grps, $cog_grp;
	}

	$ann_tbl->add_col("cog_grp", \@cog_grps);

	# save the table and overwrite the old one
	$ann_tbl->save($ann_file);
}

########
# Subs #
########
sub get_cog_grp {
	my ($grps, $method, $seed) = @_;

	if ( $grps eq "NA" ) {
		return("NA");
	}

	if ( $method eq "all" ) {
		return($grps);
	}

	my @vals = split(//, $grps);

	if ( $method eq "random" ) {
		my $len = scalar @vals;
		my $i = int(rand($len - 1));

		return($vals[$i]);
	}
	elsif ( $method eq "first" ) {
		return($vals[0]);
	}
	else { 
		$logger->warn("WARN: Bad method value: $method\n");
	}
}

sub _check_method {
	my ($method) = @_;

	# set the default
	if ( ! is_defined($method, "method") ) {
		return("all");
	}

	# check for a legal value: random|all|first
	if ( $method =~ m/random/i ) {
		return("random");
	}
	elsif ( $method =~ m/all/i ) {
		return("all");
	}
	elsif ( $method =~ m/first/i ) {
		return("first");
	}
	else {
		pod2usage(-message => "ERROR: --method must be random|all|first\n\n",
					-exitval => 2);
	}

	return("all");
}

sub _check_seed {
	my ($seed) = @_;

	# set the default to 10
	if ( ! is_defined($seed, "seed") ) {
		return(10);
	}

	# makes sure it looks like a number
	if ( looks_like_number($seed) ) {
		return($seed);
	}	
	else {
		pod2usage(-message => "ERROR: --seed must be a number\n\n",
					-exitval => 2);
	}
}

sub _check_over {
	my ($over) = @_;

	if ( ! is_defined($over, "over") ) {
		$over = 0;  # set to FALSE
	}

	$over = to_bool($over);

	return($over);
}

sub check_params {
	# check for required variables
	if ( ! defined $db) { 
		pod2usage(-message => "ERROR: required --db not defined\n\n",
					-exitval => 2); 
	}
	if ( ! defined $cog_grp_file ) {
		pod2usage(-message => "ERROR: required --cog_grp_file not defined\n\n",
					-exitval => 2);
	}
	if ( ! defined $names_file ) {
		pod2usage(-message => "ERROR: required --names_file not defined\n\n",
					-exitval => 2);
	}

	# make sure required files are non-empty
	if ( defined $cog_grp_file and ! -e $cog_grp_file ) { 
		pod2usage(-message => "ERROR: --cog_grp_file $cog_grp_file is an empty file\n\n",
					-exitval => 2);
	}
	if ( defined $names_file and ! -e $names_file ) { 
		pod2usage(-message => "ERROR: --names_file $names_file is an empty file\n\n",
					-exitval => 2);
	}

	# make sure required directories exist
	if ( ! -d $db ) { 
		pod2usage(-message => "ERROR: --db is not a directory\n\n",
					-exitval => 2); 
	}

	# set the default all_annote_file if it is not defined by the user
	if ( ! defined $all_annote_file ) {
		$all_annote_file = "all_annote.txt";
		$logger->info("Setting --all_annote_file to all_annote.txt\n");
	}

	# check seed
	$seed = _check_seed($seed);
	
	# check method
	$method = _check_method($method);

	# check over
	$over_write = _check_over($over_write);

	$logger->info("--db: $db\n");
	$logger->info("--cog_grp_file: $cog_grp_file\n");
	$logger->info("--names_file: $names_file\n");
	$logger->info("--method: $method\n");
	$logger->info("--seed: $seed\n");
	$logger->info("--over: $over_write\n");
	
	return 1;
}


__END__

# POD

=head1 NAME

DAFE_db_add_cog_grp.pl - add the cog grp info to each all_annote.txt file in the database


=head1 VERSION

This documentation refers to version 0.0.1


=head1 SYNOPSIS

    DAFE_db_add_cog_grp.pl
        --db dafe_db_path/
        --all_annote_file all_annote.txt
        --cog_grp_file cognames2003-2014.tab
        --names_file genomes.txt
        --seed 10
        --method all
        --over F
        
        [--help]
        [--man]
        [--debug]
        [--verbose]
        [--quiet]
        [--logfile logfile.log]

    --db | -d          Path to the DAFE database
    --all_annote_file  Name of the annotation files in the DAFE database
    --cog_grp_file     Path to file with each COG and the matching grp
    --name_file        Path to file with each genome to add to in the database
    --seed             A random seed (ie integer)
    --method           Method for choosing between multiple groups (all|random|first)
    --over             Bool indiciating it's safe to overwrite an existing cog_grp column        
    --help | -h        Prints USAGE statement
    --man              Prints the man page
    --debug	           Prints Log4perl DEBUG+ messages
    --verbose          Prints Log4perl INFO+ messages
    --quiet	           Suppress printing ERROR+ Log4perl messages
    --logfile          File to save Log4perl messages


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
