#! usr/bin/evn perl
sdf
use strict;
use warnings;
use Bio::TreeIO;
use Bio::Tree::TreeFunctionsI;
use Getopt::Long;
use Pod::Usage;
use Carp;
use Readonly;
use Path::Class;

#My Variables
my $help = 0;
my $man = 0;
my $orig_tree;
my $remove_file;
my $outfile;

# Read in variables from the command line
GetOptions ('man'  => \$man,
            'help' => \$help,
            'tree=s' => \$orig_tree,
            'rm_file=s' => \$remove_file,
            'out_tree|o=s' => \$outfile,
            ) || die("There was an error in the command line arguements\n");

# Use Pod usage for the Manual and Help pages
if ( $help ) { pod2usage(0) }
if ( $man )  {pod2usage(-verbose => 3) }

# MAIN #

my $rmv_fo = file($remove_file);
my @remove_ids = $rmv_fo->slurp( chomp=>1 );

my $input = Bio::TreeIO->new(-file => $orig_tree, -format => "newick");
my $tree = $input->next_tree;

foreach my $id ( @remove_ids ) {
    my $bool = $tree->remove_Node($id);
    if (!$bool) {
        print "$id not not removed from tree\n";
    }
}

my $tree_out = Bio::TreeIO->new(-file => ">$outfile", -format => "newick");
$tree_out->write_tree($tree);

