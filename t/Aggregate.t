use strict;
use warnings;

use Cwd qw/ abs_path /;
use File::Basename;
use MyX::Generic;
use Param_handler;
use Aggregate;
use XML::Simple;
use File::Temp qw/ tempfile tempdir /;
use Test::Exception;
use Test::Warn;
use Test::More tests => 16; # Need to place the number of tests you run here.

BEGIN { use_ok( 'Param_handler' );
        use_ok( 'Aggregate' ); }
        
#Global Variables to be used in tests
my $agg_to;
my $dummy_param;
my $abs_path = dirname(abs_path($0));

###########
## TESTS ##
###########

#test new
{
    lives_ok( sub{ $agg_to = Aggregate->new( get_test_param_obj() ) },
             "new() with aggregate object");
    dies_ok( sub{ $agg_to = Aggregate->new() }, "new() no params - dies");
}

#test a bad Param_handler object
{
    $dummy_param = get_test_param_obj();
    my $href = $dummy_param->get_params_href();
    delete($href->{count_dir});
    $dummy_param = Param_handler->new( {href => $href} );
    throws_ok( sub{ $agg_to = Aggregate->new( $dummy_param ) }, "MyX::Generic::Undef::Param", "bad Param_handler object" );
}

#test attribute getters
{
    $agg_to = Aggregate->new( get_test_param_obj() );
    $dummy_param = get_test_param_obj();
    is_deeply( ($agg_to->get_param_object())->get_params_href(), $dummy_param->get_params_href(), "get_param_object()");
    
    is( ref($agg_to->get_grp_order_aref()), "ARRAY", "get_grp_order_aref() -returns an aref" );
    my $grp_meta_file = ($agg_to->get_param_object())->get_grp_meta_file();
    is( scalar(@{$agg_to->get_grp_order_aref()}),
        ( (`cat $grp_meta_file | wc -l`) - 1 ),
        "get_grp_order_aref() -correct size");
    
    is( ref($agg_to->get_count_header_aref()), "ARRAY", "get_count_header_aref() -returns an aref" );
    is ( scalar(@{$agg_to->get_count_header_aref()}), 6, "get_count_header_aref() -correct size" );
}

#test aggregate componants and then aggregate
{
    # test _create_connection_btwn_geneid_grpid
    $agg_to = Aggregate->new( get_test_param_obj());
    my $connection_href = $agg_to->_create_connection_btwn_geneid_grpid("10001");
    ok( %$connection_href, "creates a non empty connection hash" );
    
    #test _grp_count_in_count_file
    my $count_href = $agg_to->_grp_count_in_count_file($connection_href,"10001");
    ok( %$count_href, "_grp_count_in_count_file - non empty hash created" );
    
    #test the creation of the aggregate file
    lives_ok( sub{ $agg_to->_create_aggregate_file($count_href, "10001") }, "_create_aggregate_file - creates a non-empty file");
    
}

#test aggregate componants on a fake param object.
{
    $agg_to = Aggregate->new( get_fake_param_obj_with_testable_data() );
    my $test_connection_href =   {
                            '101' => 'KOG01',
                            '102' => 'KOG04', #Tests the actual ordering
                            '103' => 'KOG02',
                            '104' => 'KOG03',
                            '105' => 'KOG05',
                            '106' => 'KOG06',
                            '107' => 'KOG01', #WILL TEST ADDITIONS
                            };
    my $connection_href = $agg_to->_create_connection_btwn_geneid_grpid("10001");
    is_deeply($connection_href, $test_connection_href, "_connection_btwn_grp_gene");
    
    my $test_count_href =   {
                                'KOG01' => [2,2],
                                'KOG02' => [1,1],
                                'KOG03' => [1,1],
                                'KOG04' => [1,1],
                                'KOG05' => [1,1],
                                'KOG06' => [1,1],
                            };
    my $count_href = $agg_to->_grp_count_in_count_file($connection_href,"10001");
    is_deeply($count_href, $test_count_href,"_grp_count_in_count_file()");
    
    $agg_to->_create_aggregate_file( $count_href, "10001");
    
    #test full aggregate() method. look in test directory for the output file
    lives_ok( sub{ $agg_to->_create_aggregate_file($count_href, "10001") },
             "aggregate() on test data");
}






sub get_test_param_obj {
    #put all test data into the hash ref
    my $xml_href = {
       Rsource_dir => "$abs_path/../R_lib",
       annote_file_name => "all_annote.txt", #Each genome should have one
       grp_genes_by => "kog", # Each annote file should have this column
       gene_id_col => "proteinId", #Check for a column that matches this name
       count_dir => "$abs_path/../t/test_dir/full_test_dir/DAFE_cnt_results",
       dafe_dir => "$abs_path/../t/test_dir/full_test_dir/JGI_fungal_db",
       # ^The database directory^
       out_dir => "$abs_path/../t/test_dir/full_test_dir/test_R_stats",
       ref_meta_file => "$abs_path/../t/test_dir/full_test_dir/fungi_test.txt", #Need to make sure this annotation file exists
       ref_meta_cols => '["Fraction", "Source", "Label"]', #Columns used in heatmap creation and each should be a column in the metadata file
       ref_include_file => "$abs_path/../t/test_dir/full_test_dir/names.txt", #The ids of genomes that should be included
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
    };

    my $param_obj = Param_handler->new( {href=>$xml_href} );
    return $param_obj;
}

#this contains easily testable data specifically for the aggregate function.
#The directory can be found /proj/cdjones_lab/ncolaian/tests/test_aggregate/
sub get_fake_param_obj_with_testable_data {
    my $xml_href = {
       Rsource_dir => "$abs_path/../R_lib",
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
       grp_meta_file => "$abs_path/../t/test_dir/test_aggregate/fake_kog_metadata.txt", #file containing all the metadata for the features being used
       count_file_name => "gene_counts_id60.txt", #this is the default value
       min_sample_count => 3, #Must be a positive integer
       min_sample_cpm => 0.03, #Must be positive and greater than 0
       test => '["BK", "RZ"]', #defines the two sample groups to compare
       test_col_name => "fraction", #where to look at the MetaG meta file to identify which group each experiment is from
    };

    my $param_obj = Param_handler->new( {href=>$xml_href} );
    return $param_obj;
}
