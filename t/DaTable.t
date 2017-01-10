use strict;
use warnings;

use Cwd qw/ abs_path /;
use File::Basename;
use DaTable;
use Data::Dumper;
use XML::Simple qw(:strict);
use File::Temp qw/ tempfile tempdir /;
use Test::Exception;
use Test::Warn;
use Test::More tests => 20; # Need to place the number of tests you run here.

BEGIN { use_ok( 'DaTable' )}

#Global Variables
my $da_tbl_to;
my $dummy_param;
my $abs_path = dirname(abs_path($0));

###########
## TESTS ##
###########

#test new
{
    lives_ok( sub{ $da_tbl_to = DaTable->new( get_grp_and_id_order_arrays_in_href() ) }, "new with correct info");
    dies_ok( sub{ $da_tbl_to = DaTable->new() }, "new w/ no parameters dies" );
    
    $dummy_param = get_grp_and_id_order_arrays_in_href();
    delete $dummy_param->{grp_order_aref};
    throws_ok( sub{ $da_tbl_to = DaTable->new( $dummy_param ) }, "MyX::Generic::Undef::Param", "new w/o one of the parameters" );
}

#check the setters
{
    #check the correct creation of the count hashes
    $da_tbl_to = DaTable->new( get_grp_and_id_order_arrays_in_href() );
    $dummy_param = $da_tbl_to->get_grp_order_count_href();
    is( $dummy_param->{KOG03}, 2, "grp order href made correctly");
    $dummy_param = $da_tbl_to->get_id_order_count_href();
    is( $dummy_param->{10001}, 0, "id order href made correctly");
    
    #check the da table size
    $dummy_param = $da_tbl_to->get_full_da_table();
    is( scalar @$dummy_param, scalar( @{$da_tbl_to->get_id_order_aref()}), "the table has right # of genomes" );
    is( scalar( @{$dummy_param->[1]} ), scalar( @{$da_tbl_to->get_grp_order_aref()} ), "the table has the right # of groups" );
}

#test setting the first genome
#test::warn to test for carped warnings
{
    $da_tbl_to = DaTable->new( get_grp_and_id_order_arrays_in_href() );
    $dummy_param = [ ['KOG01', 0], ['KOG02', 1], ['KOG03', -1], ['KOG04', -2], ['KOG05', -3], ['KOG06', 1], ['KOG07', 0] ];
    lives_ok( sub{ $da_tbl_to->set_genome(10001, $dummy_param) }, "set_genome lives" );
    is_deeply( $da_tbl_to->get_genome(10001), [0,1,-1,-2,-3,1,0], "correct table addition via get genome" );
    my $data_table = $da_tbl_to->get_full_da_table();
    is_deeply( $data_table->[0], [0,1,-1,-2,-3,1,0], "Genome set in right spot on da table" );
    $da_tbl_to->set_genome(10002, $dummy_param);
    is_deeply( $da_tbl_to->get_group('KOG01'), [0,0], "Get group works correctly" );
    
    warnings_are { $da_tbl_to->set_genome(10003, $dummy_param) } [{carped => '10003 is not in the ordered genome id passed'}], "try to add bad genome";
    warnings_are { $da_tbl_to->set_genome(10002, [ ['KOG09', 1] ]) } [{carped => "KOG09 is in the edger file but not in the grp metadata file" }], "try to add bad group";
    
    #Test printer
    my $tempdir = tempdir();
	my ($fh, $filename) = tempfile();
    lives_ok( sub{ $da_tbl_to->print_full_da_table($filename) }, "full print lives" );
    close($fh);
    
    #test the passing of a da_table file
    lives_ok( sub{ $dummy_param = DaTable->new( {'da_file' => "$abs_path/../t/test_dir/test_full_da_tbl.txt",} ) }, "create DaTable with file");
    throws_ok( sub{ my $fail = DaTable->new( {'da_file' => "$abs_path/../t/test_dir/test_aggregate/test_include_file.txt"} ) }, "MyX::Generic::BadValue", "pass bad file to constructor");
    
    is_deeply( $dummy_param->get_full_da_table(), $da_tbl_to->get_full_da_table(),"creation of DaTable from file is correct");
    
    
    lives_ok( sub{ $da_tbl_to->check_for_unset_values()}, "check good da_table");
    $da_tbl_to->set_genome(10002, [['KOG01', 4]]);
    throws_ok( sub{ $da_tbl_to->check_for_unset_values()}, "MyX::Generic::Undef::Attribute", "throwsfor bad set value");
}



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
       grp_meta_file => "$abs_path/../t/test_dir/test_aggregate/fake_kog_metadata.txt", #file containing all the metadata for the features being used
       count_file_name => "gene_counts_id60.txt", #this is the default value
       min_sample_count => 3, #Must be a positive integer
       min_sample_cpm => 0.03, #Must be positive and greater than 0
       test => '["BK", "RZ"]', #defines the two sample groups to compare
       test_col_name => "fraction", #where to look at the MetaG meta file to identify which group each experiment is from
       heat_filter => "FALSE", #Do the columns of the heatmap need to be filtered
       Rsource_dir => "$abs_path/../R_lib/",
    };

    my $param_obj = Param_handler->new( {href=>$xml_href} );
    return $param_obj;
}

sub get_grp_and_id_order_arrays_in_href {
    my $href = {};
    my $grp_aref = [ 'KOG01', 'KOG02', 'KOG03', 'KOG04', 'KOG05', 'KOG06', 'KOG07' ];
    my $id_aref = [ '10001', '10002' ];
    $href->{grp_order_aref} = $grp_aref;
    $href->{id_order_aref} = $id_aref;
    
    return $href;
}