use strict;
use warnings;

use Cwd qw/ abs_path /;
use File::Basename;
use Param_handler;
use Justify;
use Aggregate;
use DecoupleDa;
use Data::Dumper;
use XML::Simple qw(:strict);
use File::Temp qw/ tempfile tempdir /;
use Test::Exception;
use Test::Warn;
use Test::More tests => 6; # Need to place the number of tests you run here.

BEGIN { use_ok( 'DecoupleDa' ); }

#Global Variables

my $decda_to;
my $dummy_param;
my $param_handler_obj;
my $abs_path = dirname(abs_path($0));

###########
## TESTS ##
###########

#test new
{
    lives_ok( sub{ $decda_to = DecoupleDa->new( get_fake_param_obj_with_testable_data()) }, "new() with Param_handler" );
    
    dies_ok( sub{ $decda_to = DecoupleDa->new() }, "new() - no params dies" );
    
    $param_handler_obj = get_fake_param_obj_with_testable_data();
    my $href = $param_handler_obj->get_params_href();
    $href->{ref_meta_cols} = "Bad";
    $dummy_param = Param_handler->new( {href => $href} );
    #throw the bad print parameters error
    throws_ok( sub{ $decda_to = DecoupleDa->new( $dummy_param ) }, "MyX::Generic::BadValue", "Bad print params" );
}
#test decouple
{
    $decda_to = DecoupleDa->new( get_fake_param_obj_with_testable_data() );
    lives_ok( sub{ $decda_to->decouple(10001) }, "decouple lives" );
    #has bad count file name and it needs to manually be set
    $param_handler_obj = get_fake_param_obj_with_testable_data();
    my $href = $param_handler_obj->get_params_href();
    $href->{count_file_name} = "Bad";
    throws_ok( sub{ $decda_to = DecoupleDa->new ( Param_handler->new({href => $href}) ) }, "MyX::Generic::BadValue", "bad count file name" );
}



# Easily testable data in Param_handler object
sub get_fake_param_obj_with_testable_data {
    my $xml_href = {
       annote_file_name => "all_annote.txt", #Each genome should have one
       grp_genes_by => "kog", # Each annote file should have this column
       gene_id_col => "proteinId", #Check for a column that matches this name
       count_dir => "$abs_path/../t/test_dir/test_aggregate/count_dir",
       dafe_dir => "$abs_path/../t/test_dir/test_aggregate/dafe_dir",
       # ^The database directory^
       out_dir => "$abs_path/../t/test_dir/full_test_dir/test_R_stats",
       ref_meta_file => "$abs_path/../t/test_dir/full_test_dir/fungi_test.txt", #Need to make sure this annotation file exists
       ref_meta_cols => '["Fraction", "Source", "Label"]', #Columns used in heatmap creation and each should be a column in the metadata file
       ref_include_file => "$abs_path/../t/test_dir/test_aggregate/test_include_file.txt", #The ids of genomes that should be included
       ref_exclude_file => "$abs_path/../t/test_dir/full_test_dir/empty_file.txt",
       # ^The exclude file should be optional ^
       genome_id_col => "ID", 
       metaG_meta_file => "$abs_path/../t/test_dir/full_test_dir/metagenome_metadata.txt",
       metaG_include_file => "$abs_path/../t/test_dir/full_test_dir/sample.names.txt",
       metaG_exclude_file => "$abs_path/../t/test_dir/full_test_dir/metaG_exclude.txt", #optional parameter
       tree => "$abs_path/../t/test_dir/full_test_dir/Rooted_test_newick.nwk",
       grp_meta_file => "$abs_path/../t/test_dir/full_test_dir/kog_metadata.txt", #file containing all the metadata for the features being used
       count_file_name => "gene_counts_id60.txt", #this is the default value
       min_sample_count => 3, #Must be a positive integer
       min_sample_cpm => 0.03, #Must be positive and greater than 0
       test => '["BK", "RZ"]', #defines the two sample groups to compare
       test_col_name => "fraction", #where to look at the MetaG meta file to identify which group each experiment is from
       Rsource_dir => "$abs_path/../R_lib/",
    };

    my $param_obj = Param_handler->new( {href=>$xml_href} );
    return $param_obj;
}