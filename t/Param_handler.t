use strict;
use warnings;

use lib("../lib");
use Param_handler;
use XML::Simple;
use File::Temp qw/ tempfile tempdir /;
use Test::Exception;
use Test::Warn;
use Test::More tests => 159; # Need to place the number of tests you run here.


BEGIN { use_ok( 'Param_handler' ); }

sub get_test_href; #Simple hashref to pass into tests
sub get_test_file; #Simple file to pass into tests

#Global variable to use in tests
my $params_to;
my $params_thref;
my $dummy_param;
my $empty_file = "../test_dir/full_test_dir/empty_file.txt";
my $yml_file = "../test_dir/test_edgeR_driver.yaml";

### TESTS ###


# Test new()
{
	lives_ok( sub{$params_to = Param_handler->new()},
		  "New() w/ no params lives" );
	lives_ok( sub{$params_to = Param_handler->new({href => get_test_href()})},
		  "new(href) - lives" );
	lives_ok( sub{$params_to = Param_handler->new({xml_file => get_test_file()})},
		  "new(xml_file) - lives" );
	$params_to = new_ok( "Param_handler" );
	lives_ok( sub{$params_to = Param_handler->new({yml_file => $yml_file})}, "new(yml_file)");
}

# Test set_params()
{
	#ask scott if we need to have the set params create a new object
	$params_to = Param_handler->new();
	lives_ok( sub{$params_to->set_params({href => get_test_href()})},
		  "set_params(href) - lives" );
	lives_ok( sub{$params_to->set_params({xml_file => get_test_file()})},
		  "set_params(xml_file) - lives" );
}

# Test get_params_href()
{
	$params_thref = get_test_href();
	$params_to = Param_handler->new({href => $params_thref});
	is_deeply( $params_thref, $params_to->get_params_href(), "get_params_href()" );
	$dummy_param = $params_to->get_test_names();
	$params_to = Param_handler->new({xml_file => get_test_file()});
	is_deeply( $params_thref, $params_to->get_params_href(), "get_params_href() with xml_file" );
	is($params_to->get_test_names(), $dummy_param, "test quoted variable test");
}

#test check_edger_params()
{
	#Need to test if it works without print params
	$params_thref = get_test_href();
	$params_thref->{p3_height} = "";
	$params_thref->{ref_meta_cols} = "";
	$params_thref->{heat_filter} = "";
	$params_to = Param_handler->new( {href => $params_thref});
	warnings_are { $params_to -> check_edger_params() } [], "just global and edger params";
}

#test check_print_params()
{
	#Check the print parameters without the edger parameters
	$params_thref->{p3_height} = 3;
	$params_thref->{ref_meta_cols} = '["Fraction", "Source", "Label"]';
	$params_thref->{heat_filter} = "FALSE";
	$params_thref->{ref_meta_file} = "../test_dir/full_test_dir/fungi_test.txt";
	$params_thref->{ref_include_file} = "../test_dir/full_test_dir/names.txt";
	$params_thref->{ref_exclude_file} = "../test_dir/full_test_dir/empty_file.txt";
	$params_thref->{tree} = "../test_dir/full_test_dir/Rooted_test_newick.nwk";
	$params_thref->{out_dir} = "../test_dir/full_test_dir/test_R_stats";
	
	$params_to = Param_handler->new( {href => $params_thref});
	warnings_are { $params_to -> check_print_params() } [], "just global and print params";
}

#test _check_global_params()
{
	$params_thref->{ref_meta_file} = "../test_dir/full_test_dir/fungi_test.txt";
	$params_thref->{ref_include_file} = "../test_dir/full_test_dir/names.txt";
	$params_thref->{ref_exclude_file} = "../test_dir/full_test_dir/empty_file.txt";
	$params_thref->{tree} = "../test_dir/full_test_dir/Rooted_test_newick.nwk";
	$params_thref->{out_dir} = "../test_dir/full_test_dir/test_R_stats";
	
	$params_to = Param_handler->new( {href => $params_thref});
	warnings_are { $params_to -> _check_global_params() } [], "just global params";

}

#test spews (spew_xml_file() & spew_out())
{
	$params_to = Param_handler->new( {href => get_test_href()} );
	$params_to->spew_xml_file("../test_dir/full_test_dir/test_R_stats/test.xml");
	$params_to->spew_out();
	cmp_ok(-s "../test_dir/full_test_dir/test_R_stats/test.xml", '>', 0);
}

#test _find_unknown_tags
{
	$params_thref = get_test_href();
	$params_thref->{bad_tag} = "hi";
	$params_to = Param_handler->new( {href => $params_thref});
	throws_ok( sub{ $params_to->check_print_params() }, 'MyX::Generic::Undef::Attribute', "_find_unknown_tags");
}

#TEST CHECKS

#_check_ref_meta_file
{
	# Check if it works if correct param is passed
	$params_to = Param_handler->new( {href => get_test_href} );
	warnings_are { $params_to->_check_ref_meta_file() } [], "_check_ref_meta_file -- no warnings";

	# Check if MyX::Generic::Undef::Param is thrown if param isn't passed
	$params_to = Param_handler->new();
	throws_ok( sub{ $params_to->_check_ref_meta_file }, 'MyX::Generic::Undef::Param', "_check_ref_meta_file -- check when ref_meta_file isn't passed" );

	# Check if a file that doesn't exist or an empty file is passed
	$dummy_param = "file";
	$params_thref = { ref_meta_file => $dummy_param };
	$params_to = Param_handler->new( {href => $params_thref} );
	throws_ok( sub{ $params_to->_check_ref_meta_file }, 'MyX::Generic::DoesNotExist::File', "_check_ref_meta_file -- file doesn't exist" );

	$params_thref = { ref_meta_file => $empty_file };
	$params_to = Param_handler->new( {href => $params_thref} );
	throws_ok( sub{ $params_to->_check_ref_meta_file }, 'MyX::Generic::DoesNotExist::File', "_check_ref_meta_file -- given an empty file" );
}

#_check_ref_include_file
{
	$params_to = Param_handler->new( {href => get_test_href} );
	warnings_are { $params_to->_check_ref_include_file() } [], "_check_ref_include_file -- no warnings";

	#checks if error is thrown when param isn't passed
	$params_to = Param_handler->new();
	throws_ok( sub{ $params_to->_check_ref_include_file }, 'MyX::Generic::Undef::Param', "_check_ref_include_file -- check when ref_include_file isn't passed" );

	# Check if a file that doesn't exist or an empty file is passed
	$dummy_param = "file";
	$params_thref = { ref_include_file => $dummy_param };
	$params_to = Param_handler->new( {href => $params_thref} );
	throws_ok( sub{ $params_to->_check_ref_include_file }, 'MyX::Generic::DoesNotExist::File', "_check_ref_include_file -- file doesn't exist" );

	$params_thref = { ref_include_file => $empty_file };
	$params_to = Param_handler->new( {href => $params_thref} );
	throws_ok( sub{ $params_to->_check_ref_include_file }, 'MyX::Generic::DoesNotExist::File', "_check_ref_include_file -- given an empty file" );
	
	# Check to make sure the refs to include are in the metadata file
	$params_thref = get_test_href();
	$params_thref->{ref_include_file} = "../test_dir/full_test_dir/sample.names.txt";
	$params_to = Param_handler->new( {href => $params_thref} );
	throws_ok( sub{ $params_to->_check_ref_include_file }, 'MyX::Generic::BadValue', "_check_ref_include_file -- given wrong inclusion file" );
}

#_check_ref_exclude_file
{
	$params_to = Param_handler->new( {href => get_test_href} );
	warnings_are { $params_to->_check_ref_exclude_file() } [], "_check_ref_exclude_file -- no warnings";
	
	$params_thref = get_test_href();
	$params_thref->{ref_exclude_file} = "umm";
	$params_to = Param_handler->new( {href => $params_thref} );
	throws_ok( sub{ $params_to->_check_ref_exclude_file }, 'MyX::Generic::DoesNotExist::File', "_check_ref_exclude_file -- file doesn't exist" );

	# Check to make sure exclude values are in include file
	$params_thref->{ref_exclude_file} = "../test_dir/full_test_dir/sample.names.txt";
	$params_to = Param_handler->new( {href => $params_thref} );
	throws_ok( sub{ $params_to->_check_ref_exclude_file }, 'MyX::Generic::BadValue', "_check_ref_exclude_file -- given wrong exclusion file" );
}

#_check_tree
{
	# Check if it works if correct param is passed
	$params_to = Param_handler->new( {href => get_test_href} );
	warnings_are { $params_to->_check_tree() } [], "_check_tree -- no warnings";

	# Check if MyX::Generic::Undef::Param is thrown if param isn't passed
	$params_to = Param_handler->new();
	throws_ok( sub{ $params_to->_check_tree }, 'MyX::Generic::Undef::Param', "_check_tree -- check when file isn't passed" );

	# Check if a file that doesn't exist or an empty file is passed
	$dummy_param = "file";
	$params_thref = { tree => $dummy_param };
	$params_to = Param_handler->new( {href => $params_thref} );
	throws_ok( sub{ $params_to->_check_tree }, 'MyX::Generic::DoesNotExist::File', "_check_tree -- file doesn't exist" );

	$params_thref = { tree => $empty_file };
	$params_to = Param_handler->new( {href => $params_thref} );
	throws_ok( sub{ $params_to->_check_tree }, 'MyX::Generic::DoesNotExist::File', "_check_tree -- given an empty file" );
}

#_check_out_dir
{
	# Check if it works if correct param is passed
	$params_to = Param_handler->new( {href => get_test_href} );
	warnings_are { $params_to->_check_out_dir() } [], "_check_out_dir -- no warnings";
	
	# Check if MyX::Generic::Undef::Param is thrown if param isn't passed
	$params_to = Param_handler->new();
	throws_ok( sub{ $params_to->_check_out_dir }, 'MyX::Generic::Undef::Param', "_check_out_dir -- check when the parameter isn't passed" );
}

#_check_annote_file_name
{
	# Check if it works if correct param is passed
	$params_to = Param_handler->new( {href => get_test_href} );
	warnings_are { $params_to->_check_annote_file_name() } [], "_check_annote_file_name -- no warnings";

	# Check if MyX::Generic::Undef::Param is thrown if param isn't passed
	$params_to = Param_handler->new();
	throws_ok( sub{ $params_to->_check_annote_file_name() }, 'MyX::Generic::Undef::Param', "_check_annote_file_name -- check when name isn't passed" );

	# Check for bad annote_file_name
	$params_thref = get_test_href();
	$params_thref->{annote_file_name} = "not_good";
	$params_to = Param_handler->new( {href => $params_thref} );
	throws_ok( sub{ $params_to->_check_annote_file_name() }, 'MyX::Generic::BadValue', "_check_annote_file_name -- wrong name is passed" );
}

#_check_grp_genes_by
{
	# Check if it works if correct param is passed
	$params_to = Param_handler->new( {href => get_test_href} );
	warnings_are { $params_to->_check_grp_genes_by() } [], "_check_grp_genes_by -- no warnings";

	#Check if MyX::Generic::Undef::Param is thrown if param isn't passed
	$params_to = Param_handler->new();
	throws_ok( sub{ $params_to->_check_grp_genes_by() }, 'MyX::Generic::Undef::Param', "_check_grp_genes_by -- check when name isn't passed" );

	# Check annotation file for grp_genes_by column
	$params_thref = get_test_href();
	$params_thref->{grp_genes_by} = "not_good";
	$params_to = Param_handler->new( {href => $params_thref} );
	throws_ok( sub{ $params_to->_check_grp_genes_by() }, 'MyX::Generic::BadValue', "_check_grp_genes_by -- wrong name is passed" );
}

#_check_gene_id_col
{
	# Check if it works if correct param is passed
	$params_to = Param_handler->new( {href => get_test_href} );
	warnings_are { $params_to->_check_gene_id_col() } [], "_check_gene_id_col -- no warnings";

	#Check if MyX::Generic::Undef::Param is thrown if param isn't passed
	$params_to = Param_handler->new();
	throws_ok( sub{ $params_to->_check_gene_id_col() }, 'MyX::Generic::Undef::Param', "_check_gene_id_col() -- check when name isn't passed" );

	# Check annotation file for wrong column name
	 $params_thref = get_test_href();
	 $params_thref->{gene_id_col} = "not_good";
	 $params_to = Param_handler->new( {href => $params_thref} );
	 throws_ok( sub{ $params_to->_check_gene_id_col() }, 'MyX::Generic::BadValue', "_check_gene_id_col() -- wrong name is passed" );

}

#_check_count_dir
{
	# Check if it works if correct param is passed
	$params_to = Param_handler->new( {href => get_test_href} );
	warnings_are { $params_to->_check_count_dir() } [], "_check_count_dir() -- no warnings";

	#Check if MyX::Generic::Undef::Param is thrown if param isn't passed
	$params_to = Param_handler->new();
	throws_ok( sub{ $params_to->_check_count_dir() }, 'MyX::Generic::Undef::Param', "_check_count_dir() -- check when name isn't passed" );

	#Check if count dir doesnt exist
	$params_thref = get_test_href();
	$params_thref->{count_dir} = "not_good";
	$params_to = Param_handler->new( {href => $params_thref} );
	throws_ok( sub{ $params_to->_check_count_dir() }, 'MyX::Generic::DoesNotExist::Dir', "_check_count_dir() -- dir does not exist" );
}

#_check_dafe_dir
{
	# Check if it works if correct param is passed
	$params_to = Param_handler->new( {href => get_test_href} );
	warnings_are { $params_to->_check_dafe_dir() } [], "_check_dafe_dir() -- no warnings";

	#Check if MyX::Generic::Undef::Param is thrown if param isn't passed
	$params_to = Param_handler->new();
	throws_ok( sub{ $params_to->_check_dafe_dir() }, 'MyX::Generic::Undef::Param', "_check_dafe_dir() -- check when name isn't passed" );

	#Check if dafe dir doesnt exist
	$params_thref = get_test_href();
	$params_thref->{dafe_dir} = "not_good";
	$params_to = Param_handler->new( {href => $params_thref} );
	throws_ok( sub{ $params_to->_check_dafe_dir() }, 'MyX::Generic::DoesNotExist::Dir', "_check_dafe_dir() -- dir does not exist" );
}

#_check_genome_id_col
{
	# Check if it works if correct param is passed
	$params_to = Param_handler->new( {href => get_test_href} );
	warnings_are { $params_to->_check_genome_id_col() } [], "_check_genome_id_col() -- no warnings";

	#Check if MyX::Generic::Undef::Param is thrown if param isn't passed
	$params_to = Param_handler->new();
	throws_ok( sub{ $params_to->_check_genome_id_col() }, 'MyX::Generic::Undef::Param', "_check_genome_id_col() -- check when name isn't passed" );
	
	#Check if bad col name is passed
	$params_thref = get_test_href();
	$params_thref->{genome_id_col} = "not_good";
	$params_to = Param_handler->new( {href => $params_thref} );
	throws_ok( sub{ $params_to->_check_genome_id_col() }, 'MyX::Generic::BadValue', "_check_genome_id_col() -- dir does not exist" );
}

#_check_metaG_meta_file
{
	# Check if it works if correct param is passed
	$params_to = Param_handler->new( {href => get_test_href} );
	warnings_are { $params_to->_check_metaG_meta_file() } [], "_check_metaG_meta_file() -- no warnings";

	# Check if MyX::Generic::Undef::Param is thrown if param isn't passed
	$params_to = Param_handler->new();
	throws_ok( sub{ $params_to->_check_metaG_meta_file }, 'MyX::Generic::Undef::Param', "_check_metaG_meta_file() -- check when file isn't passed" );

	# Check if a file that doesn't exist or an empty file is passed
	$dummy_param = "file";
	$params_thref = { metaG_meta_file => $dummy_param };
	$params_to = Param_handler->new( {href => $params_thref} );
	throws_ok( sub{ $params_to->_check_metaG_meta_file }, 'MyX::Generic::DoesNotExist::File', "_check_metaG_meta_file() -- file doesn't exist" );

	$params_thref = { metaG_meta_file => $empty_file };
	$params_to = Param_handler->new( {href => $params_thref} );
	throws_ok( sub{ $params_to->_check_metaG_meta_file }, 'MyX::Generic::DoesNotExist::File', "_check_metaG_meta_file() -- given an empty file" );
}

#_check_metaG_include_file
{
	# Check if it works if correct param is passed
	$params_to = Param_handler->new( {href => get_test_href} );
	warnings_are { $params_to->_check_metaG_include_file() } [], "_check_metaG_include_file() -- no warnings";

	# Check if MyX::Generic::Undef::Param is thrown if param isn't passed
	$params_to = Param_handler->new();
	throws_ok( sub{ $params_to->_check_metaG_include_file }, 'MyX::Generic::Undef::Param', "_check_metaG_include_file() -- check when file isn't passed" );

	# Check if a file that doesn't exist or an empty file is passed
	$dummy_param = "file";
	$params_thref = { metaG_include_file => $dummy_param };
	$params_to = Param_handler->new( {href => $params_thref} );
	throws_ok( sub{ $params_to->_check_metaG_include_file() }, 'MyX::Generic::DoesNotExist::File', "_check_metaG_include_file() -- file doesn't exist" );

	$params_thref = { metaG_include_file => $empty_file };
	$params_to = Param_handler->new( {href => $params_thref} );
	throws_ok( sub{ $params_to->_check_metaG_include_file() }, 'MyX::Generic::DoesNotExist::File', "_check_metaG_include_file() -- given an empty file" );

	#Check to see error is given when bad include file is given
	$params_thref = get_test_href();
	$params_thref->{metaG_include_file} = "../test_dir/test_bad/full_db_names.txt";
	$params_to = Param_handler->new( {href => $params_thref} );
	throws_ok( sub{ $params_to->_check_metaG_include_file() }, 'MyX::Generic::BadValue', "_check_metaG_include_file -- give bad inclusion file" );
}

#_check_metaG_exclude_file
{
	# Check if it works if correct param is passed
	$params_to = Param_handler->new( {href => get_test_href} );
	warnings_are { $params_to->_check_metaG_exclude_file() } [], "_check_metaG_exclude_file() -- no warnings";

	# Check if a file that doesn't exist or an empty file is passed
	$dummy_param = "DNE";
	$params_thref = { metaG_exclude_file => $dummy_param };
	$params_to = Param_handler->new( {href => $params_thref} );
	throws_ok( sub{ $params_to->_check_metaG_exclude_file() }, 'MyX::Generic::DoesNotExist::File', "_check_metaG_exclude_file() -- file doesn't exist" );

	# Check to make sure the exclude info is in the include file
	$params_thref = get_test_href();
	$params_thref->{metaG_exclude_file} = "../test_dir/full_test_dir/kog_metadata.txt";
	$params_to = Param_handler->new( {href => $params_thref} );
	throws_ok( sub{ $params_to->_check_metaG_exclude_file() }, 'MyX::Generic::BadValue', "_check_metaG_exclude_file() -- wrong exclude file" );
}

#_check_count_file_name
{
	# Check if it works if correct param is passed
	$params_to = Param_handler->new( {href => get_test_href} );
	warnings_are { $params_to->_check_count_file_name() } [], "_check_count_file_name() -- no warnings";

	# Check if MyX::Generic::Undef::Param is thrown if param isn't passed
	$params_to = Param_handler->new();
	throws_ok( sub{ $params_to->_check_count_file_name() }, 'MyX::Generic::Undef::Param', "_check_count_file_name() -- check when file isn't passed" );

	#Check for error when bad file name is passed
	$params_thref = get_test_href();
	$params_thref->{count_file_name} = "not_good";
	$params_to = Param_handler->new( {href => $params_thref} );
	throws_ok( sub{ $params_to->_check_count_file_name() }, 'MyX::Generic::BadValue', "_check_count_file_name() -- dir does not exist" );
}

#_check_min_sample_count
{
	# Check if it works if correct param is passed
	$params_to = Param_handler->new( {href => get_test_href} );
	warnings_are { $params_to->_check_min_sample_count() } [], "_check_min_sample_count() -- no warnings";

	# Check if MyX::Generic::Undef::Param is thrown if param isn't passed
	$params_to = Param_handler->new();
	throws_ok( sub{ $params_to->_check_min_sample_count() }, 'MyX::Generic::Undef::Param', "_check_min_sample_count() -- check when file isn't passed" );

	#Check if it will throw right error when a negative or non int is passed
	$params_thref->{min_sample_count} = -2;
	$params_to = Param_handler->new( {href => $params_thref} );
	throws_ok( sub{ $params_to->_check_min_sample_count() }, 'MyX::Generic::Digit', "_check_min_sample_count() -- negative number" );

	$params_thref->{min_sample_count} = 0.02;
	$params_to = Param_handler->new( {href => $params_thref} );
	throws_ok( sub{ $params_to->_check_min_sample_count() }, 'MyX::Generic::Digit', "_check_min_sample_count() -- non-integer passed" );
}

#_check_min_sample_cpm
{
	# Check if it works if correct param is passed
	$params_to = Param_handler->new( {href => get_test_href} );
	warnings_are { $params_to->_check_min_sample_cpm() } [], "_check_min_sample_cpm() -- no warnings";

	# Check if MyX::Generic::Undef::Param is thrown if param isn't passed
	$params_to = Param_handler->new();
	throws_ok( sub{ $params_to->_check_min_sample_cpm() }, 'MyX::Generic::Undef::Param', "_check_min_sample_cpm() -- check when file isn't passed" );

	# Throw error when a negative value is given NEED TO COME BACK TO THIS
	$params_thref->{min_sample_cpm} = -2;
	$params_to = Param_handler->new( {href => $params_thref} );
	throws_ok( sub{ $params_to->_check_min_sample_cpm() }, 'MyX::Generic::Digit::TooSmall', "_check_min_sample_cpm() -- negative number" );
}

#_check_test_names
{
	# Check if it works if correct param is passed
	$params_to = Param_handler->new( {href => get_test_href} );
	warnings_are { $params_to->_check_test_names() } [], "_check_test_names() -- no warnings";

	# Check if MyX::Generic::Undef::Param is thrown if param isn't passed
	$params_to = Param_handler->new();
	throws_ok( sub{ $params_to->_check_test_names() }, 'MyX::Generic::Undef::Param', "_check_test_names() -- check when file isn't passed" );

	#check if wrong test_names are given
	$params_thref = get_test_href();
	$params_thref->{test_names} = "not_good";
	$params_to = Param_handler->new( {href => $params_thref} );
	throws_ok( sub{ $params_to->_check_test_names() }, 'MyX::Generic::BadValue', "_check_test_names() -- wrong names are passed" );	
}

#_check_test_col_name
{
	# Check if it works if correct param is passed
	$params_to = Param_handler->new( {href => get_test_href} );
	warnings_are { $params_to->_check_test_col_name() } [], "_check_test_col_name() -- no warnings";

	# Check if MyX::Generic::Undef::Param is thrown if param isn't passed
	$params_to = Param_handler->new();
	throws_ok( sub{ $params_to->_check_test_col_name() }, 'MyX::Generic::Undef::Param', "_check_test_col_name() -- check when file isn't passed" );
	
	#check if wrong test_col_name is given
	$params_thref = get_test_href();
	$params_thref->{test_col_name} = "not_good";
	$params_to = Param_handler->new( {href => $params_thref} );
	throws_ok( sub{ $params_to->_check_test_col_name() }, 'MyX::Generic::BadValue', "_check_test_col_name() -- wrong names are passed" );
}

#_check_grp_meta_file
{
	# Check if it works if correct param is passed
	$params_to = Param_handler->new( {href => get_test_href} );
	warnings_are { $params_to->_check_grp_meta_file() } [], "_check_grp_meta_file -- no warnings";

	# Check if MyX::Generic::Undef::Param is thrown if param isn't passed
	$params_to = Param_handler->new();
	throws_ok( sub{ $params_to->_check_grp_meta_file }, 'MyX::Generic::Undef::Param', "_check_grp_meta_file -- check when file isn't passed" );

	# Check if a file that doesn't exist or an empty file is passed
	$dummy_param = "file";
	$params_thref = { grp_meta_file => $dummy_param };
	$params_to = Param_handler->new( {href => $params_thref} );
	throws_ok( sub{ $params_to->_check_grp_meta_file }, 'MyX::Generic::DoesNotExist::File', "_check_grp_meta_file -- file doesn't exist" );

	$params_thref = { grp_meta_file => $empty_file };
	$params_to = Param_handler->new( {href => $params_thref} );
	throws_ok( sub{ $params_to->_check_grp_meta_file }, 'MyX::Generic::DoesNotExist::File', "_check_grp_meta_file -- given an empty file" );
}

#_check_p3_height
{
	# Check if it works if correct param is passed
	$params_to = Param_handler->new( {href => get_test_href} );
	warnings_are { $params_to->_check_p3_height() } [], "_check_p3_height -- no warnings";

	# Check if MyX::Generic::Undef::Param is thrown if param isn't passed
	$params_to = Param_handler->new();
	throws_ok( sub{ $params_to->_check_p3_height }, 'MyX::Generic::Undef::Param', "_check_p3_height -- check when file isn't passed" );
}

#_check_ref_meta_cols
{
	# Check if it works if correct param is passed
	$params_to = Param_handler->new( {href => get_test_href} );
	warnings_are { $params_to->_check_ref_meta_cols() } [], "_check_ref_meta_cols -- no warnings";

	# Check if MyX::Generic::Undef::Param is thrown if param isn't passed
	$params_to = Param_handler->new();
	throws_ok( sub{ $params_to->_check_ref_meta_cols }, 'MyX::Generic::Undef::Param', "_check_ref_meta_cols -- check when file isn't passed" );

	#check if wrong ref_meta_cols are given
	$params_thref = get_test_href();
	$params_thref->{ref_meta_cols} = '["not_good","bad"]';
	$params_to = Param_handler->new( {href => $params_thref} );
	throws_ok( sub{ $params_to->_check_ref_meta_cols() }, 'MyX::Generic::BadValue', "_check_ref_meta_cols() -- wrong names are passed" );
}

#_check_heat_filter
{
	# Check if it works if correct param is passed
	$params_to = Param_handler->new( {href => get_test_href()} );
	warnings_are { $params_to->_check_heat_filter() } [], "_check_heat_filter -- no warnings";

	# Check if MyX::Generic::Undef::Param is thrown if param isn't passed
	$params_to = Param_handler->new();
	throws_ok( sub{ $params_to->_check_heat_filter() }, 'MyX::Generic::Undef::Param', "_check_heat_filter -- check when filter isn't passed" );

	#check to see if lowercase single letter can be passed
	$params_thref = get_test_href();
	$params_thref->{heat_filter} = "f";
	$params_to = Param_handler->new( {href => $params_thref} );
	warnings_are { $params_to->_check_heat_filter() } [], "_check_heat_filter -- no warnings when pass either a t or f";
	
	#check to see if wrong value throws error
	$params_thref = get_test_href();
	$params_thref->{heat_filter} = "no";
	$params_to = Param_handler->new( {href => $params_thref} );
	throws_ok( sub { $params_to->_check_heat_filter() }, 'MyX::Generic::BadValue', "_check_heat_filter -- no warnings when pass either a 0 or 1");
}
	

#TEST GETTERS#

#get_ref_meta_file
{
	$dummy_param = "file";
	$params_thref = { ref_meta_file => $dummy_param };
	$params_to = Param_handler->new({href => $params_thref});
	is($params_to->get_ref_meta_file(),$dummy_param,"get_ref_meta_file()");
}

#get_ref_include_file
{
	$dummy_param = "file";
	$params_thref = { ref_include_file => $dummy_param };
	$params_to = Param_handler->new({href => $params_thref});
	is($params_to->get_ref_include_file(),$dummy_param,"get_ref_include_file()");
}

#get_ref_exclude_file
{
	$dummy_param = "file";
	$params_thref = { ref_exclude_file => $dummy_param };
	$params_to = Param_handler->new({href => $params_thref});
	is($params_to->get_ref_exclude_file(),$dummy_param,"get_ref_exclude_file()");
}

#get_tree
{
	$dummy_param = "file";
	$params_thref = { tree => $dummy_param };
	$params_to = Param_handler->new({href => $params_thref});
	is($params_to->get_tree(),$dummy_param,"get_tree()");
}

#get_out_dir
{
	$dummy_param = "file";
	$params_thref = { out_dir => $dummy_param };
	$params_to = Param_handler->new({href => $params_thref});
	is($params_to->get_out_dir(),$dummy_param,"get_out_dir()");
}

#get_annote_file
{
	$dummy_param = "file";
	$params_thref = { annote_file_name => $dummy_param };
	$params_to = Param_handler->new({href => $params_thref});
	is($params_to->get_annote_file_name(),$dummy_param,"get_annote_file_name()");
}

#get_grp_genes_by
{
	$dummy_param = "file";
	$params_thref = { grp_genes_by => $dummy_param };
	$params_to = Param_handler->new({href => $params_thref});
	is($params_to->get_grp_genes_by(),$dummy_param,"get_grp_genes_by()");
}

#get_gene_id_col
{
	$dummy_param = "file";
	$params_thref = { gene_id_col => $dummy_param };
	$params_to = Param_handler->new({href => $params_thref});
	is($params_to->get_gene_id_col(),$dummy_param,"get_gene_id_col()");
}

#get_count_dir
{
	$dummy_param = "file";
	$params_thref = { count_dir => $dummy_param };
	$params_to = Param_handler->new({href => $params_thref});
	is($params_to->get_count_dir(),$dummy_param,"get_count_dir()");
}

#get_dafe_dir
{
	$dummy_param = "file";
	$params_thref = { dafe_dir => $dummy_param };
	$params_to = Param_handler->new({href => $params_thref});
	is($params_to->get_dafe_dir(),$dummy_param,"get_dafe_dir()");
}

#get_genome_id_col
{
	$dummy_param = "file";
	$params_thref = { genome_id_col => $dummy_param };
	$params_to = Param_handler->new({href => $params_thref});
	is($params_to->get_genome_id_col(),$dummy_param,"get_genome_id_col()");
}

#get_metaG_meta_file
{
	$dummy_param = "file";
	$params_thref = { metaG_meta_file => $dummy_param };
	$params_to = Param_handler->new({href => $params_thref});
	is($params_to->get_metaG_meta_file(),$dummy_param,"get_metaG_meta_file()");
}

#get_metaG_include_file
{
	$dummy_param = "file";
	$params_thref = { metaG_include_file => $dummy_param };
	$params_to = Param_handler->new({href => $params_thref});
	is($params_to->get_metaG_include_file(),$dummy_param,"get_metaG_include_file()");
}

#get_metaG_exclude_file
{
	$dummy_param = "files";
	$params_thref = { metaG_exclude_file => $dummy_param };
	$params_to = Param_handler->new({href => $params_thref});
	is($params_to->get_metaG_exclude_file(),$dummy_param,"get_metaG_exclude_file()");
}

#get_count_file_name
{
	$dummy_param = "file";
	$params_thref = { count_file_name => $dummy_param };
	$params_to = Param_handler->new({href => $params_thref});
	is($params_to->get_count_file_name(),$dummy_param,"get_count_file_name()");
}

#get_min_sample_count
{
	$dummy_param = "file";
	$params_thref = { min_sample_count => $dummy_param };
	$params_to = Param_handler->new({href => $params_thref});
	is($params_to->get_min_sample_count(),$dummy_param,"get_min_sample_count()");
}

#get_min_sample_cpm
{
	$dummy_param = "file";
	$params_thref = { min_sample_cpm => $dummy_param };
	$params_to = Param_handler->new({href => $params_thref});
	is($params_to->get_min_sample_cpm(),$dummy_param,"get_min_sample_cpm()");
}

#get_test_names
{
	$dummy_param = "file";
	$params_thref = { test_names => $dummy_param };
	$params_to = Param_handler->new({href => $params_thref});
	is($params_to->get_test_names(),$dummy_param,"get_test_names()");
}

#get_test_col_name
{
	$dummy_param = "file";
	$params_thref = { test_col_name => $dummy_param };
	$params_to = Param_handler->new({href => $params_thref});
	is($params_to->get_test_col_name(),$dummy_param,"get_test_col_name()");
}

#get_grp_meta_file
{
	$dummy_param = "file";
	$params_thref = { grp_meta_file => $dummy_param };
	$params_to = Param_handler->new({href => $params_thref});
	is($params_to->get_grp_meta_file(),$dummy_param,"get_grp_meta_file()");
}

#get_p3_height
{
	$dummy_param = "file";
	$params_thref = { p3_height => $dummy_param };
	$params_to = Param_handler->new({href => $params_thref});
	is($params_to->get_p3_height(),$dummy_param,"get_p3_height()");
}

#get_ref_meta_cols
{
	$dummy_param = "file";
	$params_thref = { ref_meta_cols => $dummy_param };
	$params_to = Param_handler->new({href => $params_thref});
	is($params_to->get_ref_meta_cols(),$dummy_param,"get_ref_meta_cols()");
}

#get_heat_filter
{
	$dummy_param = "file";
	$params_thref = { heat_filter => $dummy_param };
	$params_to = Param_handler->new({href => $params_thref});
	is($params_to->get_heat_filter(),$dummy_param,"get_heat_filter()");
}

#TEST SETTERS

#set_ref_meta_file
{
	#check to see that the set will work with a good value
	$params_to = Param_handler->new();
	$dummy_param = "../test_dir/full_test_dir/fungi_test.txt";
	$params_to->set_ref_meta_file($dummy_param);
	is($params_to->get_ref_meta_file(),$dummy_param, "set_ref_meta_file() -- good param passed");

	#if a bad value is given to the setter it will revert to the old value
	$dummy_param = "bad_value";
	$params_thref = get_test_href();
	$params_to = Param_handler->new( {href => $params_thref} );
	my $old_var = $params_to->get_ref_meta_file();
	$params_to->set_ref_meta_file($dummy_param);
	is($params_to->get_ref_meta_file(),$old_var, "set_ref_meta_file() -- bad parameter passed");
}

#set_ref_include_file
{
	$params_thref = get_test_href();
	my $old_value = $params_thref->{ref_include_file};
	$params_thref->{ref_include_file} = "filler";
	$params_to = Param_handler->new( {href => $params_thref} );
	$params_to->set_ref_include_file($old_value);
	is($params_to->get_ref_include_file(),$old_value, "set_ref_include_file() -- good param passed");
	
	#if a bad value is given to the setter it will revert to the old value
	$dummy_param = "bad_value";
	$params_to->set_ref_include_file($dummy_param);
	is($params_to->get_ref_include_file(),$old_value, "set_ref_include_file() -- bad parameter passed");
}

#set_ref_exclude_file
{
	$params_thref = get_test_href();
	my $old_value = $params_thref->{ref_exclude_file};
	$params_thref->{ref_exclude_file} = "filler";
	$params_to = Param_handler->new( {href => $params_thref} );
	$params_to->set_ref_exclude_file($old_value);
	is($params_to->get_ref_exclude_file(),$old_value, "set_ref_exclude_file() -- good param passed");
	
	#if a bad value is given to the setter it will revert to the old value
	$dummy_param = "bad_value";
	$params_to->set_ref_exclude_file($dummy_param);
	is($params_to->get_ref_exclude_file(),$old_value, "set_ref_exclude_file() -- bad parameter passed");
}

#set_tree
{
	$params_thref = get_test_href();
	my $old_value = $params_thref->{tree};
	$params_thref->{tree} = "filler";
	$params_to = Param_handler->new( {href => $params_thref} );
	$params_to->set_tree($old_value);
	is($params_to->get_tree(),$old_value, "set_tree() -- good param passed");
	
	#if a bad value is given to the setter it will revert to the old value
	$dummy_param = "bad_value";
	$params_to->set_tree($dummy_param);
	is($params_to->get_tree(),$old_value, "set_tree() -- bad parameter passed");
}

#set_annote_file_name
{
	$params_thref = get_test_href();
	my $old_value = $params_thref->{annote_file_name};
	$params_thref->{annote_file_name} = "filler";
	$params_to = Param_handler->new( {href => $params_thref} );
	$params_to->set_annote_file_name($old_value);
	is($params_to->get_annote_file_name(),$old_value, "set_annote_file_name() -- good param passed");
	
	#if a bad value is given to the setter it will revert to the old value
	$dummy_param = "bad_value";
	$params_to->set_annote_file_name($dummy_param);
	is($params_to->get_annote_file_name(),$old_value, "set_annote_file_name() -- bad parameter passed");
}

#set_grp_genes_by
{
	$params_thref = get_test_href();
	my $old_value = $params_thref->{grp_genes_by};
	$params_thref->{grp_genes_by} = "filler";
	$params_to = Param_handler->new( {href => $params_thref} );
	$params_to->set_grp_genes_by($old_value);
	is($params_to->get_grp_genes_by(),$old_value, "set_grp_genes_by() -- good param passed");
	
	#if a bad value is given to the setter it will revert to the old value
	$dummy_param = "bad_value";
	$params_to->set_grp_genes_by($dummy_param);
	is($params_to->get_grp_genes_by(),$old_value, "set_grp_genes_by() -- bad parameter passed");
}

#set_gene_id_col
{
	$params_thref = get_test_href();
	my $old_value = $params_thref->{gene_id_col};
	$params_thref->{gene_id_col} = "filler";
	$params_to = Param_handler->new( {href => $params_thref} );
	$params_to->set_gene_id_col($old_value);
	is($params_to->get_gene_id_col(),$old_value, "set_gene_id_col() -- good param passed");
	
	#if a bad value is given to the setter it will revert to the old value
	$dummy_param = "bad_value";
	$params_to->set_gene_id_col($dummy_param);
	is($params_to->get_gene_id_col(),$old_value, "set_gene_id_col() -- bad parameter passed");
}

#set_count_dir
{
	$params_thref = get_test_href();
	my $old_value = $params_thref->{count_dir};
	$params_thref->{count_dir} = "filler";
	$params_to = Param_handler->new( {href => $params_thref} );
	$params_to->set_count_dir($old_value);
	is($params_to->get_count_dir(),$old_value, "set_count_dir() -- good param passed");
	
	#if a bad value is given to the setter it will revert to the old value
	$dummy_param = "bad_value";
	$params_to->set_count_dir($dummy_param);
	is($params_to->get_count_dir(),$old_value, "set_count_dir() -- bad parameter passed");
}

#set_dafe_dir
{
	$params_thref = get_test_href();
	my $old_value = $params_thref->{dafe_dir};
	$params_thref->{dafe_dir} = "filler";
	$params_to = Param_handler->new( {href => $params_thref} );
	$params_to->set_dafe_dir($old_value);
	is($params_to->get_dafe_dir(),$old_value, "set_dafe_dir() -- good param passed");
	
	#if a bad value is given to the setter it will revert to the old value
	$dummy_param = "bad_value";
	$params_to->set_dafe_dir($dummy_param);
	is($params_to->get_dafe_dir(),$old_value, "set_dafe_dir() -- bad parameter passed");
}

#set_genome_id_col
{
	$params_thref = get_test_href();
	my $old_value = $params_thref->{genome_id_col};
	$params_thref->{genome_id_col} = "filler";
	$params_to = Param_handler->new( {href => $params_thref} );
	$params_to->set_genome_id_col($old_value);
	is($params_to->get_genome_id_col(),$old_value, "set_genome_id_col() -- good param passed");
	
	#if a bad value is given to the setter it will revert to the old value
	$dummy_param = "bad_value";
	$params_to->set_genome_id_col($dummy_param);
	is($params_to->get_genome_id_col(),$old_value, "set_genome_id_col() -- bad parameter passed");
}

#set_metaG_meta_file
{
	$params_thref = get_test_href();
	my $old_value = $params_thref->{tree};
	$params_thref->{metaG_meta_file} = "filler";
	$params_to = Param_handler->new( {href => $params_thref} );
	$params_to->set_metaG_meta_file($old_value);
	is($params_to->get_metaG_meta_file(),$old_value, "set_metaG_meta_file() -- good param passed");
	
	#if a bad value is given to the setter it will revert to the old value
	$dummy_param = "bad_value";
	$params_to->set_metaG_meta_file($dummy_param);
	is($params_to->get_metaG_meta_file(),$old_value, "set_metaG_meta_file() -- bad parameter passed");
}

#set_metaG_include_file
{
	$params_thref = get_test_href();
	my $old_value = $params_thref->{metaG_include_file};
	$params_thref->{metaG_include_file} = "filler";
	$params_to = Param_handler->new( {href => $params_thref} );
	$params_to->set_metaG_include_file($old_value);
	is($params_to->get_metaG_include_file(),$old_value, "set_metaG_include_file() -- good param passed");
	
	#if a bad value is given to the setter it will revert to the old value
	$dummy_param = "bad_value";
	$params_to->set_metaG_include_file($dummy_param);
	is($params_to->get_metaG_include_file(),$old_value, "set_metaG_include_file() -- bad parameter passed");
}

#set_metaG_exclude_file
{
	$params_thref = get_test_href();
	my $old_value = $params_thref->{metaG_exclude_file};
	$params_thref->{metaG_exclude_file} = "filler";
	$params_to = Param_handler->new( {href => $params_thref} );
	$params_to->set_metaG_exclude_file($old_value);
	is($params_to->get_metaG_exclude_file(),$old_value, "set_metaG_exclude_file() -- good param passed");
	
	#if a bad value is given to the setter it will revert to the old value
	$dummy_param = "bad_value";
	$params_to->set_metaG_exclude_file($dummy_param);
	is($params_to->get_metaG_exclude_file(),$old_value, "set_metaG_exclude_file() -- bad parameter passed");
}

#set_count_file_name
{
	$params_thref = get_test_href();
	my $old_value = $params_thref->{count_file_name};
	$params_thref->{count_file_name} = "filler";
	$params_to = Param_handler->new( {href => $params_thref} );
	$params_to->set_count_file_name($old_value);
	is($params_to->get_count_file_name(),$old_value, "set_count_file_name() -- good param passed");
	
	#if a bad value is given to the setter it will revert to the old value
	$dummy_param = "bad_value";
	$params_to->set_count_file_name($dummy_param);
	is($params_to->get_count_file_name(),$old_value, "set_count_file_name() -- bad parameter passed");
}

#set_min_sample_count
{
	$params_thref = get_test_href();
	my $old_value = $params_thref->{min_sample_count};
	$params_thref->{min_sample_count} = "filler";
	$params_to = Param_handler->new( {href => $params_thref} );
	$params_to->set_min_sample_count($old_value);
	is($params_to->get_min_sample_count(),$old_value, "set_min_sample_count() -- good param passed");
	
	#if a bad value is given to the setter it will revert to the old value
	$dummy_param = "-2";
	$params_to->set_min_sample_count($dummy_param);
	is($params_to->get_min_sample_count(),$old_value, "set_min_sample_count() -- bad parameter passed");
}

#set_min_sample_cpm
{
	$params_thref = get_test_href();
	my $old_value = $params_thref->{min_sample_cpm};
	$params_thref->{min_sample_cpm} = "filler";
	$params_to = Param_handler->new( {href => $params_thref} );
	$params_to->set_min_sample_cpm($old_value);
	is($params_to->get_min_sample_cpm(),$old_value, "set_min_sample_cpm() -- good param passed");
	
	#if a bad value is given to the setter it will revert to the old value
	$dummy_param = "-2";
	$params_to->set_min_sample_cpm($dummy_param);
	is($params_to->get_min_sample_cpm(),$old_value, "set_min_sample_cpm() -- bad parameter passed");
}

#set_test_names
{
	$params_thref = get_test_href();
	my $old_value = $params_thref->{test_names};
	$params_thref->{test_names} = "filler";
	$params_to = Param_handler->new( {href => $params_thref} );
	$params_to->set_test_names($old_value);
	is($params_to->get_test_names(),$old_value, "set_test_names() -- good param passed");
	
	#if a bad value is given to the setter it will revert to the old value
	$dummy_param = "bad_value";
	$params_to->set_test_names($dummy_param);
	is($params_to->get_test_names(),$old_value, "set_test_names() -- bad parameter passed");
}

#set_test_col_name
{
	$params_thref = get_test_href();
	my $old_value = $params_thref->{test_col_name};
	$params_thref->{test_col_name} = "filler";
	$params_to = Param_handler->new( {href => $params_thref} );
	$params_to->set_test_col_name($old_value);
	is($params_to->get_test_col_name(),$old_value, "set_test_col_name() -- good param passed");
	
	#if a bad value is given to the setter it will revert to the old value
	$dummy_param = "bad_value";
	$params_to->set_test_col_name($dummy_param);
	is($params_to->get_test_col_name(),$old_value, "set_test_col_name() -- bad parameter passed");
}

#set_grp_meta_file
{
	$params_thref = get_test_href();
	my $old_value = $params_thref->{grp_meta_file};
	$params_thref->{grp_meta_file} = "filler";
	$params_to = Param_handler->new( {href => $params_thref} );
	$params_to->set_grp_meta_file($old_value);
	is($params_to->get_grp_meta_file(),$old_value, "set_grp_meta_file() -- good param passed");
	
	#if a bad value is given to the setter it will revert to the old value
	$dummy_param = "bad_value";
	$params_to->set_grp_meta_file($dummy_param);
	is($params_to->get_grp_meta_file(),$old_value, "set_grp_meta_file() -- bad parameter passed");
}

#set_p3_height
{
	$params_thref = get_test_href();
	my $old_value = $params_thref->{p3_height};
	$params_thref->{p3_height} = "filler";
	$params_to = Param_handler->new( {href => $params_thref} );
	$params_to->set_p3_height($old_value);
	is($params_to->get_p3_height(),$old_value, "set_p3_height() -- good param passed");
}

#set_ref_meta_cols
{
	$params_thref = get_test_href();
	my $old_value = $params_thref->{ref_meta_cols};
	$params_thref->{ref_meta_cols} = "filler";
	$params_to = Param_handler->new( {href => $params_thref} );
	$params_to->set_ref_meta_cols($old_value);
	is($params_to->get_ref_meta_cols(),$old_value, "set_ref_meta_cols() -- good param passed");
	
	#if a bad value is given to the setter it will revert to the old value
	$dummy_param = "bad_value";
	$params_to->set_ref_meta_cols($dummy_param);
	is($params_to->get_ref_meta_cols(),$old_value, "set_ref_meta_cols() -- bad parameter passed");
}


#set_heat_filter
{
	$params_thref = get_test_href();
	my $old_value = $params_thref->{heat_filter};
	$params_thref->{heat_filter} = "filler";
	$params_to = Param_handler->new( {href => $params_thref} );
	$params_to->set_heat_filter($old_value);
	is($params_to->get_heat_filter(),$old_value, "set_heat_filter() -- good param passed");
	
	#if a bad value is given to the setter it will revert to the old value
	$dummy_param = "bad_value";
	$params_to->set_heat_filter($dummy_param);
	is($params_to->get_heat_filter(),$old_value, "set_heat_filter() -- bad parameter passed");
}



### SUBROUTINES ###

#Subroutine to get a test hashref that has test variables inside of it
sub get_test_href {
    my $tempdir = tempdir();

    #put all test data into the hash ref
    my $xml_href = {
       annote_file_name => "all_annote.txt", #Each genome should have one
       grp_genes_by => "kog", # Each annote file should have this column
       gene_id_col => "proteinId", #Check for a column that matches this name
       count_dir => "../test_dir/full_test_dir/DAFE_cnt_results",
       dafe_dir => "../test_dir/full_test_dir/JGI_fungal_db",
       # ^The database directory^
       out_dir => "../test_dir/full_test_dir/test_R_stats",
       ref_meta_file => "../test_dir/full_test_dir/fungi_test.txt", #Need to make sure this annotation file exists
       ref_meta_cols => '["Fraction", "Source", "Label"]', #Columns used in heatmap creation and each should be a column in the metadata file
       ref_include_file => "../test_dir/full_test_dir/names.txt", #The ids of genomes that should be included
       ref_exclude_file => "../test_dir/full_test_dir/empty_file.txt",
       # ^The exclude file should be optional ^
       genome_id_col => "ID", 
       metaG_meta_file => "../test_dir/full_test_dir/metagenome_metadata.txt",
       metaG_include_file => "../test_dir/full_test_dir/sample.names.txt",
       metaG_exclude_file => "../test_dir/full_test_dir/metaG_exclude.txt", #optional parameter
       tree => "../test_dir/full_test_dir/Rooted_test_newick.nwk",
       grp_meta_file => "../test_dir/full_test_dir/kog_metadata.txt", #file containing all the metadata for the features being used
       count_file_name => "gene_counts_id60.txt", #this is the default value
       min_sample_count => 3, #Must be a positive integer
       min_sample_cpm => 0.03, #Must be positive and greater than 0
       test_names => '["BK", "RZ"]', #defines the two sample groups to compare
       test_col_name => "fraction", #where to look at the MetaG meta file to identify which group each experiment is from
       heat_filter => "FALSE", #Do the columns of the heatmap need to be filtered
       p3_height => 8, #Must be a number, but determines the plot height
    };

    return $xml_href;
}

#subroutine that returns a temporary xml file containing all the test info.
sub get_test_file {
    my ($fh, $filename) = tempfile();

    my $xml_href = get_test_href();

    #write the href
    my $xml = XMLout($xml_href);
    print $fh $xml;

    close($fh);
    return $filename;
}
