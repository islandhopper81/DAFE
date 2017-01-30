use strict;
use warnings;

use Cwd qw/ abs_path /;
use File::Basename;
use Param_handler;
use Justify;
use XML::Simple;
use File::Temp qw/ tempfile tempdir /;
use Test::Exception;
use Test::Warn;
use Test::More tests => 17; # Need to place the number of tests you run here.

BEGIN { use_ok( 'Param_handler' );
        use_ok( 'Justify' ); }

#Global variable to use in tests
my $justify_to;
my $justify_thref;
my $dummy_param;
my $abs_path = dirname(abs_path($0));
my $empty_file = "$abs_path/../t/test_dir/full_test_dir/empty_file.txt";

###TESTS###

#Test new
{
    lives_ok( sub{ $justify_to = Justify->new( get_test_href() ) }, "new() with href" );
    lives_ok( sub{ $justify_to = Justify->new( get_test_param_obj() ) }, "new() with Param_handler object" );
    dies_ok( sub{ $justify_to = Justify->new() }, "new() - die when no params are passed" );
}

#test _justify_information
{
    #test if input must all be there
    $justify_thref = get_test_href();
    delete $justify_thref->{count_dir};
    throws_ok( sub{ $justify_to = Justify->new( $justify_thref ) }, "MyX::Generic::Undef::Param", "_check_for_all_input" );
    
    #test if an unknown tag is passed
    $justify_thref = get_test_href();
    $justify_thref->{bad_tag} = 1;
    throws_ok( sub{ $justify_to = Justify->new( $justify_thref ) }, "MyX::Generic::Undef::Attribute", "_check_for_unknown_tags" );
    
    #test a bad exclude file that has an id that does not match an id in the exclude file
    $justify_thref = get_test_href();
    $justify_thref->{exclude_file} = "$abs_path/../t/test_dir/test_bad/sample.names.txt";
    throws_ok( sub{ $justify_to = Justify->new( $justify_thref ) }, "MyX::Generic::File", "_determine_included_ids" );
    
    #test check tree and create ordered tree href
    #Need to look at the possibility of id's being strings and not numbers. This will currently cause problems
    $justify_thref = get_test_href();
    $justify_thref->{include_file} = "$abs_path/../t/test_dir/test_bad/full_db_names.txt";
    throws_ok( sub{ $justify_to = Justify->new( $justify_thref ) }, "MyX::Generic::BadValue", "bad include file passed" );
    $justify_thref = get_test_href();
    $justify_thref->{tree_file} = "$abs_path/../t/test_dir/test_bad/trial_raxml.nwk";
    throws_ok( sub{ $justify_to = Justify->new( $justify_thref ) }, "MyX::Generic::BadValue", "bad tree file" );
    
    #test check_count_dir
    $justify_thref = get_test_href();
    $justify_thref->{count_dir} = "bad/directory";
    throws_ok( sub{ $justify_to = Justify->new( $justify_thref ) }, "MyX::Generic::DoesNotExist::File", "bad count_dir" );
    $justify_thref = get_test_href();
    $justify_thref->{count_file_name} = "bad_file_name";
    throws_ok( sub{ $justify_to = Justify->new( $justify_thref ) }, "MyX::Generic::DoesNotExist::File", "bad count_file_name" );
    
    #test check ids by passing in a bad metadata file
    $justify_thref = get_test_href();
    $justify_thref->{meta_data_file} = "$abs_path/../t/test_dir/full_test_dir/metagenome_metadata.txt";
    throws_ok( sub{ $justify_to = Justify->new( $justify_thref ) }, "MyX::Generic::BadValue", "bad metadata file" );
}

###GETTERS###

#test get_passed_params_href plus all the get_passed_single_params
{
    $justify_thref = get_test_href();
    $justify_to = Justify->new( $justify_thref );
    $justify_thref = get_test_href(); #The creation of the object changes the passed in hash ref - Need to ask scott if this is a problem
    is_deeply( $justify_to->get_passed_params_href(), $justify_thref, "get_passed_params_href()")
}

{
    #test get_ordered_ids_aref
    $justify_thref = get_test_href();
    $justify_to = Justify->new( $justify_thref );
    my $test_id_aref = [10003,10013,10011,10008,10001,10006,10002,10012,10010,10007,10004,10015,10005,10009,10184];
    is_deeply( $justify_to->get_ordered_ids_aref(), $test_id_aref, "creation of ordered id aref");
}

#test get_trimmed_ordered_aref (Eyeball the spewed meta file in the tester below)

###PRINTERS###

#test spew_trimmed_ordered_meta_file and spew_ordered_ids ( Tests the gets too since they are used in these methods)
{
    $justify_to = Justify->new( get_test_href() );
    unlink("$abs_path/../t/test_dir/full_test_dir/test_R_stats/test_meta.xml");
	$justify_to->spew_trimmed_ordered_meta_file("$abs_path/../t/test_dir/full_test_dir/test_R_stats/test_meta.xml");
    cmp_ok(-s "$abs_path/../t/test_dir/full_test_dir/test_R_stats/test_meta.xml", '>', 0);
    unlink("$abs_path/../t/test_dir/full_test_dir/test_R_stats/test.xml");
    $justify_to->spew_ordered_ids("$abs_path/../t/test_dir/full_test_dir/test_R_stats/test.xml");
    cmp_ok(-s "$abs_path/../t/test_dir/full_test_dir/test_R_stats/test.xml", '>', 0);
}

### SUBROUTINES ###

#Subroutine to get a test hashref that has test variables inside of it
sub get_test_param_obj {
    #put all test data into the hash ref
    my $xml_href = {
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
       Rsource_dir => "$abs_path/../R_lib/",
    };

    my $param_obj = Param_handler->new( {href=>$xml_href} );
    return $param_obj;
}

sub get_test_href {
    
    my $href = {
        include_file    => "$abs_path/../t/test_dir/full_test_dir/names.txt",
        exclude_file    => "$abs_path/../t/test_dir/full_test_dir/empty_file.txt",
        meta_data_file  => "$abs_path/../t/test_dir/full_test_dir/fungi_test.txt",
        tree_file       => "$abs_path/../t/test_dir/full_test_dir/Rooted_test_newick.nwk",
        count_dir       => "$abs_path/../t/test_dir/full_test_dir/DAFE_cnt_results",
        count_file_name => "gene_counts_id60.txt",
        id_col          => "ID",
    };
    
    return $href;
}