use strict;
use warnings;

use lib("../lib");
use DaFilter;
use DaTable;
use Param_handler;
use XML::Simple;
use File::Temp qw/ tempfile tempdir /;
use Test::Exception;
use Test::Warn;
use Test::More tests => 15; # Need to place the number of tests you run here.

BEGIN { use_ok( 'DaFilter' );
        use_ok( 'DaTable' );
        use_ok( 'Param_handler' ); }

#Global variable to use in tests
my $da_filter_to;
my $dummy_param;
my $test_count_info = [ { 1 => 0, 0 => 2, -1 => 0, -2 => 0, -3 => 0,'total' => 2,},
                        { 1 => 2, 0 => 0, -1 => 0, -2 => 0, -3 => 0,'total' => 2,},
                        { 1 => 0, 0 => 0, -1 => 2, -2 => 0, -3 => 0,'total' => 2,},
                        { 1 => 0, 0 => 0, -1 => 0, -2 => 2, -3 => 0,'total' => 2,},
                        { 1 => 0, 0 => 0, -1 => 0, -2 => 0, -3 => 2,'total' => 2,},
                        { 1 => 2, 0 => 0, -1 => 0, -2 => 0, -3 => 0,'total' => 2,},
                        { 1 => 0, 0 => 2, -1 => 0, -2 => 0, -3 => 0,'total' => 2,}, ];

###TESTS###

#new()
{
    $dummy_param = get_testable_da_table();
    set_test_da_table($dummy_param);
    lives_ok( sub{ $da_filter_to = DaFilter->new( $dummy_param )},
             "new()");
    dies_ok( sub{ $da_filter_to = DaFilter->new()}, "new w/o da_table_obj");
    $dummy_param->set_genome(10001, [['KOG01',4]]);
    throws_ok( sub{$da_filter_to = DaFilter->new( get_testable_da_table())}, "MyX::Generic::Undef::Attribute", "unset/bad value");
}

#getters
{
    $dummy_param = get_testable_da_table();
    set_test_da_table($dummy_param);
    $da_filter_to = DaFilter->new( $dummy_param );
    #get_grp_count_aref
    is_deeply( $da_filter_to->get_grp_count_aref(), $test_count_info, "correct count arref of hrefs");
    #get_table_obj
    is_deeply( $da_filter_to->get_table_obj(), $dummy_param, "get_table_obj");
}

#filter_and_print
{
    $dummy_param = get_testable_da_table();
    set_test_da_table($dummy_param);
    $da_filter_to = DaFilter->new( $dummy_param );
    is_deeply( $da_filter_to->_filter_by_ones(90, 'true'), {'KOG02'=>1, 'KOG06'=>5}, "f1");
    is_deeply( $da_filter_to->_filter_by_zeros(90, 'true'), {'KOG01'=>0, 'KOG07'=>6}, "f0");
    is_deeply( $da_filter_to->_filter_by_negative_ones(90, 'true'), {'KOG03'=>2}, "f-1");
    is_deeply( $da_filter_to->_filter_by_negative_twos(90, 'true'), {'KOG04'=>3}, "f-2");
    is_deeply( $da_filter_to->_filter_by_negative_threes(90, 'true'), {'KOG05'=>4}, "f-3");
    my $conserved_href;
    push @$conserved_href, $da_filter_to->_filter_by_ones(90, 'false');
    push @$conserved_href, $da_filter_to->_filter_by_zeros(90, 'false');
    push @$conserved_href, $da_filter_to->_filter_by_negative_ones(90, 'false');
    is_deeply(DaFilter::_find_conserved_filter_grps($conserved_href),
              {'KOG04'=>3, 'KOG05'=>4}, "find conserved");
    
    my $param_handl_obj = get_fake_param_obj_with_testable_data();
    lives_ok( sub{$da_filter_to->filter_and_print($param_handl_obj, "../t/test_dir/test_filter.txt")}, "filter_and_print" );
}


sub get_fake_param_obj_with_testable_data {
    my $xml_href = {
       annote_file_name => "all_annote.txt", #Each genome should have one
       grp_genes_by => "kog", # Each annote file should have this column
       gene_id_col => "proteinId", #Check for a column that matches this name
       count_dir => "../t/test_dir/test_aggregate/count_dir",
       dafe_dir => "../t/test_dir/test_aggregate/dafe_dir",
       # ^The database directory^
       out_dir => "../t/test_dir/full_test_dir/test_R_stats",
       ref_meta_file => "../t/test_dir/full_test_dir/fungi_test.txt", #Need to make sure this annotation file exists
       ref_meta_cols => '["Fraction", "Source", "Label"]', #Columns used in heatmap creation and each should be a column in the metadata file
       ref_include_file => "../t/test_dir/test_aggregate/test_include_file.txt", #The ids of genomes that should be included
       ref_exclude_file => "../t/test_dir/full_test_dir/empty_file.txt",
       # ^The exclude file should be optional ^
       genome_id_col => "ID", 
       metaG_meta_file => "../t/test_dir/full_test_dir/metagenome_metadata.txt",
       metaG_include_file => "../t/test_dir/full_test_dir/sample.names.txt",
       metaG_exclude_file => "../t/test_dir/full_test_dir/metaG_exclude.txt", #optional parameter
       tree => "../t/test_dir/full_test_dir/Rooted_test_newick.nwk",
       grp_meta_file => "../t/test_dir/test_aggregate/fake_kog_metadata.txt", #file containing all the metadata for the features being used
       count_file_name => "gene_counts_id60.txt", #this is the default value
       min_sample_count => 3, #Must be a positive integer
       min_sample_cpm => 0.03, #Must be positive and greater than 0
       test => '["BK", "RZ"]', #defines the two sample groups to compare
       test_col_name => "fraction", #where to look at the MetaG meta file to identify which group each experiment is from
       heat_filter => "FALSE", #Do the columns of the heatmap need to be filtered
       p3_height => 8, #Must be a number, but determines the plot height
       Rsource_dir => "/netscr/yourston/compMetaG_R_dev/R/",
       filter_params => { "f1" => [90, "false"],
						"f-1" => [90,"false"], }
    };

    my $param_obj = Param_handler->new( {href=>$xml_href} );
    return $param_obj;
}

sub get_testable_da_table {
    my $href = {};
    my $grp_aref = [ 'KOG01', 'KOG02', 'KOG03', 'KOG04', 'KOG05', 'KOG06', 'KOG07' ];
    my $id_aref = [ '10001', '10002' ];
    $href->{grp_order_aref} = $grp_aref;
    $href->{id_order_aref} = $id_aref;
    
    my $test_da_table = DaTable->new( $href );
    
    return $test_da_table;
}

sub set_test_da_table {
    my ($da_test_table_obj) = @_;
    my $da_values = [ ['KOG01', 0], ['KOG02', 1], ['KOG03', -1], ['KOG04', -2], ['KOG05', -3], ['KOG06', 1], ['KOG07', 0] ];
    $da_test_table_obj->set_genome(10001, $da_values);
    $da_test_table_obj->set_genome(10002, $da_values);
    return;
}

