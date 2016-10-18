#! usr/bin/evn perl

package Param_handler;

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Carp;
use Readonly;
use Path::Class;
use Data::Dumper;
use XML::Simple qw(:strict);
use Class::Std::Utils;
use Scalar::Util::Numeric qw(isneg isint isfloat);
use lib "/proj/cdjones_lab/ncolaian/MyX-Generic-v0.0.2/lib/MyX";
use MyX::Generic;
use YAML::XS qw(LoadFile);


{
  #Obj identifiers with hrefs
  my %attribute_hrefs_ident;

  #Hashes with needed parameters
  Readonly::Hash my %REQUIRED_GLOBAL_TAGS => map { $_ => 1 } qw(
ref_meta_file
ref_include_file
ref_exclude_file
tree
out_dir
							      );

  Readonly::Hash my %REQUIRED_EDGER_TAGS => map { $_ => 1 } qw(
annote_file_name
grp_genes_by
gene_id_col
count_dir
dafe_dir
genome_id_col
metaG_meta_file
metaG_include_file
metaG_exclude_file
count_file_name
min_sample_count
min_sample_cpm
test_names
test_col_name
grp_meta_file
							     );
  
  Readonly::Hash my %REQUIRED_PRINT_TAGS => map { $_ => 1 } qw(
p3_height
ref_meta_cols
heat_filter
							     );
  

  #SUBROUTINES#
  #Global
  sub get_params_href;
  sub get_print_params;
  sub get_global_params;
  sub get_edger_params;
  sub set_params;
  sub check_print_params;
  sub check_edger_params;
  sub spew_xml_file;
	sub spew_out;

  #Internal
  sub _check_run_type;
  sub _check_annote_file_name;
  sub _check_grp_genes_by;
  sub _check_gene_id_col;
  sub _check_count_dir;
  sub _check_dafe_dir;
  sub _check_out_dir;
 # sub _check_ref_out_dir; Been removed but will be hardwired to create a subdirectory within the out directory
  sub _check_ref_meta_file;
  sub _check_ref_meta_cols;
  sub _check_ref_include_file;
  sub _check_ref_exclude_file;
  sub _check_genome_id_col;
  sub _check_metaG_meta_file;
  sub _check_metaG_include_file;
  sub _check_metaG_exclude_file;
  sub _check_tree;
  sub _check_grp_meta_file; #need to add getter, setter, and check
  sub _check_count_file_name;
  sub _check_min_sample_count;
  sub _check_min_sample_cpm;
  sub _check_test_info;
  sub _check_test_col_name;
  sub _check_big_mat;
  sub _check_p3_height;
  sub _check_heat_filter;

  sub _check_global_params;
  sub _are_global_tags_present;
  sub _are_print_tags_present;
  sub _are_edger_tags_present;
  sub _find_unknown_tags;
  
  #############
  #CONSTRUCTOR#
  #############
  
  #constructor that takes a hash of attributes
  sub new {
    my ($class, $arg_href) = @_;

    #Bless a scalar to represet the new object
    my $new_obj = bless \do{my $anon_scalar}, $class;

    #If empty create empty hash to represent object
    if ( !$arg_href->{xml_file} && !$arg_href->{href} && !$arg_href->{yml_file} ) {
      $attribute_hrefs_ident{ident $new_obj} = {};
    }
    # Initialize the attribute
    else {
      $new_obj -> set_params($arg_href);
    }

    return $new_obj;
  }
  #############
  #SUBROUTINES#
  #############
  
  sub get_params_href {
    my ($self) = @_;
    return $attribute_hrefs_ident{ident $self};
  }

  sub set_params {
    my ($self, $arg_href) = @_;

    #Check to see if a href has been passed
    if ( defined $arg_href->{href} ) {
      $attribute_hrefs_ident{ident $self} = $arg_href->{href};
    }
    elsif ( defined $arg_href->{xml_file} ) {
      $attribute_hrefs_ident{ident $self} =
			XML::Simple->new()->XMLin($arg_href->{xml_file}, KeyAttr => [], ForceArray => 0);
    }
		elsif ( defined $arg_href->{yml_file} ) {
			$attribute_hrefs_ident{ident $self} = LoadFile( $arg_href->{yml_file} );
		}
		else {
			croak "Did not pass in attributes correctly. Please pass in a XML file, Hash Reference, or Yaml File with attributes\n";
		}
		

    return 1;
  }
  
  sub check_edger_params {
    my ($self) = @_;
    my $params_href = $self->get_params_href();

     #Check global params and tags
    $self->_check_global_params();
    #check edger tags
    _are_edger_tags_present($params_href);
    #Check print params
    $self->_check_annote_file_name();
    $self->_check_grp_genes_by();
    $self->_check_gene_id_col();
    $self->_check_count_dir();
    $self->_check_dafe_dir();
    $self->_check_genome_id_col();
    $self->_check_metaG_meta_file();
    $self->_check_metaG_include_file();
    $self->_check_metaG_exclude_file();
    $self->_check_count_file_name();
    $self->_check_min_sample_count();
    $self->_check_min_sample_cpm();
    $self->_check_test_names();
    $self->_check_test_col_name();
		$self->_check_grp_meta_file();
  }

  sub check_print_params {
    my ($self) = @_;
    my $param_href = $self->get_params_href();
    
    #Check global params and tags
    $self->_check_global_params();
    #check print tags
    _are_print_tags_present($param_href);
    #Check print params
    $self->_check_p3_height();
    $self->_check_ref_meta_cols();
		$self->_check_heat_filter();
  }
  
  sub _check_global_params {
    my ($self) = @_;
    my $param_href = $self->get_params_href();
    #Check tags
    _are_global_tags_present($param_href);
    #perform all _checks
    $self->_check_ref_meta_file();
    $self->_check_ref_include_file();
    $self->_check_ref_exclude_file();
    $self->_check_tree();
    $self->_check_out_dir();
  }

  sub spew_xml_file {
    my ($self, $out_file) = @_;
    
    open(my $OUT, ">", $out_file) ||
      croak("cannot open file: $out_file\n");

    
    my $xml = XMLout($self->get_params_href(), KeyAttr => []);
    print $OUT $xml;
    close($OUT);
  }
	
	sub spew_out {
		my ($self) = @_;
		my $params_href = $self->get_params_href();
		
		foreach my $key (keys %$params_href) {
			print $key, " => ", $params_href->{$key}, "\n";
		}
	}
  
	sub get_print_params_href {
		my ($self) = @_;
		my $print_params_href = {
			'p3_height'		     => $self->get_p3_height(),
			'ref_meta_cols'    => $self->get_ref_meta_cols(),
			'heat_filter'      => $self->get_heat_filter(),
			'ref_meta_file'    => $self->get_ref_meta_file(),
			'ref_include_file' => $self->getref_include_file(),
			'ref_exclude_file' => $self->get_ref_exclude_file(),
			'tree'             => $self->get_tree(),
			'out_dir'          => $self->get_out_dir(),
		};
		return $print_params_href;
	}
	
	sub get_edger_params_href {
		my ($self) = @_;
		my $edger_params_href = {
			'ref_meta_file'    => $self->get_ref_meta_file(),
			'ref_include_file' => $self->getref_include_file(),
			'ref_exclude_file' => $self->get_ref_exclude_file(),
			'tree'             => $self->get_tree(),
			'out_dir'          => $self->get_out_dir(),
			'annote_file_name' => $self->get_annote_file_name(),
			'grp_genes_by'     => $self->get_grp_genes_by(),
			'gene_id_col'			 => $self->get_gene_id_col(),
			'count_dir'        => $self->get_count_dir(),
			'dafe_dir'         => $self->get_dafe_dir(),
			'genome_id_col'    => $self->get_genome_id_col(),
			'metaG_meta_file'  => $self->get_metaG_meta_file(),
			'metaG_include_file' => $self->get_metaG_include_file(),
			'metaG_exclude_file' => $self->get_metaG_exclude_file(),
			'count_file_name'  => $self->get_count_file_name(),
			'min_sample_count' => $self->get_min_sample_count(),
			'min_sample_cpm'   => $self->get_min_sample_cpm(),
			'test_names'       => $self->get_test_names(),
			'test_col_name'    => $self->get_test_col_name(),
			'grp_meta_file'    => $self->get_grp_meta_file(),
			};
		return $edger_params_href;
	}
	
  sub _are_global_tags_present {
    my ($param_href) = @_;

    foreach my $param ( keys %REQUIRED_GLOBAL_TAGS ) {
      if ( ! defined $param_href->{$param} ) {
	croak "Global parameter $param is not defined\n"
      }
    }
    _find_unknown_tags($param_href);
  }

  sub _are_edger_tags_present {
    my ($param_href) = @_;

    foreach my $param ( keys %REQUIRED_EDGER_TAGS ) {
      if ( ! defined $param_href->{$param} ) {
	croak "EdgeR parameter $param is not defined\n"
      }
    }
  }

  sub _are_print_tags_present {
     my ($param_href) = @_;

    foreach my $param ( keys %REQUIRED_PRINT_TAGS ) {
      if ( ! defined $param_href->{$param} ) {
	croak "Print parameter $param is not defined\n"
      }
    }
   }

  sub _find_unknown_tags {
    my ($param_href) = @_;

    foreach my $tag ( keys %{$param_href} ) {
      if ( ! defined $REQUIRED_GLOBAL_TAGS{$tag} &&
	   ! defined $REQUIRED_EDGER_TAGS{$tag}  &&
	   ! defined $REQUIRED_PRINT_TAGS{$tag}    ) {
	MyX::Generic::Undef::Attribute->throw(
					      error => "Unknown attribute in the parameter hash",
					      att_name => $tag,
					      );
      }#if
    }#foreach
  }#sub

  #Check-Subroutines
	#LOOK AT LOG4PERL
  #Global
  sub _check_ref_meta_file {
    my ($self) = @_;
    my $param = $self->get_ref_meta_file();
    my $href  = $self->get_params_href();
    #Checks if param was passed
    if ( !$href->{ref_meta_file} ) {
      MyX::Generic::Undef::Param->throw(
					error => "Param not set",
					usage => "ref_meta_file",
					);
    }
    #Checks if the file passed exists and isn't empty
    if ( !-e $param || -z $param ) {
      MyX::Generic::DoesNotExist::File->throw(
					      error => "Reference Metadata file does not exist at the current path, or it is empty",
					      file_name => $param,
					      );
    }
    return 1;
  }
  sub _check_ref_include_file {
    my ($self) = @_;
    my $param = $self->get_ref_include_file();
    my $href  = $self->get_params_href();
    #Checks Param
    if ( !$href->{ref_include_file} ) {
      MyX::Generic::Undef::Param->throw(
					error => "Param not set",
					usage => "ref_include_file",
					);
    }
    #Check if file passed exists
    if ( !-e $param || -z $param ) {
      MyX::Generic::DoesNotExist::File->throw(
					      error => "The names to include in the analysis of the statistics wasn't provided",
					      file_name => $param,
					     );
    }

    #Check if references exist in ref_meta_file
    my $meta_file_obj = file( $self->get_ref_meta_file );
    my $param_file_obj = file( $param );
    my $slurped_file = $meta_file_obj->slurp();
    my @include = $param_file_obj->slurp( chomp=>1 );
    foreach my $ref (@include) {
      if ( $slurped_file !~ qr/$ref/i ) {
	MyX::Generic::BadValue->throw(
				      error => "$ref isn't found in metadata file",
				      value => $ref,
				     );
	last;
      }
    }
    return 1;
  }
  sub _check_ref_exclude_file {
    my ($self) = @_;
    my $param = $self->get_ref_exclude_file();
    my $href = $self->get_params_href();
    if ( $href->{ref_exclude_file} ) {
      if ( !-e $param ) {
	MyX::Generic::DoesNotExist::File->throw(
						error => "The optional file containing the references to exclude does not exist",
						file_name => $param,
					       );
      }
      elsif ( !-z $param ) {
	my $include_file_obj = file($self->get_ref_include_file());
	my $exclude_file_obj = file($param);
	my $slurped_include_file = $include_file_obj->slurp();
	my @exclude = $exclude_file_obj->slurp( chomp=>1 );
	foreach my $xcld (@exclude) {
	  if ( $slurped_include_file !~ qr/$xcld/i ) {
	    MyX::Generic::BadValue->throw(
					  error => "$xcld (exclude id) isn't found in include file",
					  value => $xcld,
					 );
	    last;
	  }
	}
      }
    }
    return 1;
  }
  sub _check_tree {
    my ($self) = @_;
    my $param = $self->get_tree();
    my $href  = $self->get_params_href();
    if ( !$href->{tree} ) {
      MyX::Generic::Undef::Param->throw(
					error => "Param not set",
					usage => "tree",
					);
    }
    if ( !-e $param || -z $param ) {
      MyX::Generic::DoesNotExist::File->throw(
					      error => "The tree file does not exist",
					      file_name => $param,
					     );
    }
    return 1;
  }
  sub _check_out_dir {
    my ($self) = @_;
    my $param = $self->get_out_dir();
    my $href  = $self->get_params_href();
    if ( !$href->{out_dir} ) {
      MyX::Generic::Undef::Param->throw(
					error => "Param not set",
					usage => "out_dir",
					);
    }
    elsif (-d $param) {
      rmdir($param);
      mkdir($param);
      #carp("Directory re-created");
    }
    else {
      mkdir($param);
      carp("New directory created");
    }
    return 1;
  }
  #EdgeR
  sub _check_annote_file_name {
    my ($self) = @_;
    my $param = $self->get_annote_file_name();
    my $href  = $self->get_params_href();
    if ( !$href->{annote_file_name} ) {
      MyX::Generic::Undef::Param->throw(
					error => "Param not set",
					usage => "annote_file",
				       );
    }
    else {
      #Check to make sure each genome has an annotation file
      my $dir_path = $self->get_dafe_dir();
      my $include_file_obj = file( $self->get_ref_include_file() );
      my $exclude_file_obj = file( $self->get_ref_exclude_file() );
      my $slurped_excluded = $exclude_file_obj->slurp();
      my @included = $include_file_obj->slurp( chomp=>1 );
      foreach my $genome (@included) {
	if ( $slurped_excluded !~ qr/$genome/ ) {
	  my $full_path = "$dir_path/$genome/$param";
	  if ( !-e $full_path || -z $full_path ) {
	    MyX::Generic::BadValue->throw(
					  error => "Annotation file name is not found in $genome",
					  value => $param,
					 );
	    last;
	  }
	}
      }
    }
    return 1;
  }
  sub _check_grp_genes_by {
    my ($self) = @_;
    my $param = $self->get_grp_genes_by();
    my $href  = $self->get_params_href();
    if ( !$href->{grp_genes_by} ) {
      MyX::Generic::Undef::Param->throw(
					error => "Param not set",
					usage => "grp_genes_by",
					);
    }
    else {
      my $name    = $self->get_annote_file_name();
      my $db_path = $self->get_dafe_dir();
      my $include_file_obj = file( $self->get_ref_include_file() );
      my @all_genome = $include_file_obj->slurp( chomp=>1 );
      my $annote_file_obj = file("$db_path/$all_genome[0]/$name");
      my @lines = $annote_file_obj->slurp( chomp=>1 );
      if ( $lines[0] !~ qr/$param/ ) {
	#need to check an annotation file for the column name specified
	MyX::Generic::BadValue->throw(
				      error => "Column used to group the data does not exist in annotation file",
				      value => $param,
				     );
      }
    }
    return 1;
  }
  sub _check_gene_id_col {
    my ($self) = @_;
    my $param = $self->get_gene_id_col();
    my $href  = $self->get_params_href();
    if ( !$href->{gene_id_col} ) {
      MyX::Generic::Undef::Param->throw(
					error => "Param not set",
					usage => "gene_id_col",
					);
    }
    else {
      my $name    = $self->get_annote_file_name();
      my $db_path = $self->get_dafe_dir();
      my $include_file_obj = file( $self->get_ref_include_file() );
      my @all_genome = $include_file_obj->slurp( chomp=>1 );
      my $annote_file_obj = file("$db_path/$all_genome[0]/$name");
      my @lines = $annote_file_obj->slurp( chomp=>1 );
      if ( $lines[0] !~ qr/$param/ ) {
	#need to check an annotation file for the column name specified
	MyX::Generic::BadValue->throw(
				      error => "$param column used to group the data does not exist in annotation file",
				      value => $param,
				     );
      }
    }
    return 1;
  }
  sub _check_count_dir {
    my ($self) = @_;
    my $param = $self->get_count_dir();
    my $href  = $self->get_params_href();
    if ( !$href->{count_dir} ) {
      MyX::Generic::Undef::Param->throw(
					error => "Param not set",
					usage => "count_dir",
					);
    }
    if (!-d $param) {
      MyX::Generic::DoesNotExist::Dir->throw(
					     error => "Count Data directory does not exist",
					     dir_name => $param,
					     );
    }
    return 1;
  }
  sub _check_dafe_dir {
    my ($self) = @_;
    my $param = $self->get_dafe_dir();
    my $href  = $self->get_params_href();
    if ( !$href->{dafe_dir} ) {
      MyX::Generic::Undef::Param->throw(
					error => "Param not set",
					usage => "dafe_dir",
					);
    }
    if (!-d $param) {
      MyX::Generic::DoesNotExist::Dir->throw(
					     error => "DAFE directory containing scripts does not exist",
					     dir_name => $param,
					     );
    }
    return 1;
  }
  sub _check_genome_id_col {
    my ($self) = @_;
    my $param = $self->get_genome_id_col();
    my $href  = $self->get_params_href();
    if ( !$href->{genome_id_col} ) {
      MyX::Generic::Undef::Param->throw(
					error => "Param not set",
					usage => "genome_id_col",
					);
    }
    else {
      my $metafile_obj = file( $self->get_ref_meta_file() );
      my @slurped_metafile = $metafile_obj->slurp( chomp => 1 );
      if ( $slurped_metafile[0] !~ qr/$param/i ) {
	#need to check the ref_meta_file for the column name specified
	MyX::Generic::BadValue->throw(
				      error => "No column named $param in the ref_meta_file",
				      value => $param,
				     );
      }
    }
    return 1;
  }
  sub _check_metaG_meta_file {
    my ($self) = @_;
    my $param = $self->get_metaG_meta_file();
    my $href  = $self->get_params_href();
    if ( !$href->{metaG_meta_file} ) {
      MyX::Generic::Undef::Param->throw(
					error => "Param not set",
					usage => "metaG_meta_file",
					);
    }
    if ( !-e $param || -z $param ) {
      MyX::Generic::DoesNotExist::File->throw(
					      error => "MetaG Metadata file does not exist",
					      file_name => $param,
					      );
    }
    return 1;
  }
  sub _check_metaG_include_file {
    my ($self) = @_;
    my $param = $self->get_metaG_include_file();
    my $href  = $self->get_params_href();
    if ( !$href->{metaG_include_file} ) {
      MyX::Generic::Undef::Param->throw(
					error => "Param not set",
					usage => "metaG_include_file",
					);
    }
    elsif ( !-e $param || -z $param ) {
      MyX::Generic::DoesNotExist::File->throw(
					      error => "The names to include from the metaG metavfile wasn't provided",
					      file_name => $param,
					     );
    }
    else {
      #Check if all the included id's are in the metaG_meta_file
      my $metaG_meta_file_obj = file( $self->get_metaG_meta_file() );
      my $include_file_obj = file( $param );
      my $slurped_metaG_file = $metaG_meta_file_obj->slurp();
      my @include = $include_file_obj->slurp( chomp=> 1 );
      foreach my $id (@include) {
	if ( $slurped_metaG_file !~ qr/$id/i ) {
	  MyX::Generic::BadValue->throw(
					error => "MetaG inclusion file has value(s) that do not exist in the meta file",
					value => $id,
				       );
	  last;
	}
      }
    }
    return 1;
  }
  sub _check_metaG_exclude_file {
    my ($self) = @_;
    my $param = $self->get_metaG_exclude_file();
    my $href  = $self->get_params_href();
    if ( $href->{metaG_exclude_file} ) {
      if ( !-e $param ) {
	MyX::Generic::DoesNotExist::File->throw(
						error => "The optional metaG exclusion file was passed, but it does not exist",
						file_name => $param,
					       );
      }
      elsif ( !-z $param ) {
	#Check to make sure the id's to exclude are in the inclusion file
	my $include_file_obj = file( $self->get_metaG_include_file() );
	my $exclude_file_obj = file( $param );
	my $slurped_include_file = $include_file_obj->slurp();
	my @exclude = $exclude_file_obj->slurp( chomp=>1 );
	foreach my $id (@exclude) {
	  if ($slurped_include_file !~ qr/$id/i ) {
	    MyX::Generic::BadValue->throw(
					  error => "Exclude file id's do not match inclusion file id's",
					  value => $id,
					 );
	    last;
	  }
	}
      }
    }
    return 1;
  }
  sub _check_count_file_name {
    my ($self) = @_;
    my $param = $self->get_count_file_name();
    my $href  = $self->get_params_href();
    if ( !$href->{count_file_name} ) {
      MyX::Generic::Undef::Param->throw(
					error => "Param not set",
					usage => "count_file",
					);
    }
    else {
      my $include_file_obj = file( $self->get_ref_include_file() );
      my $count_dir = $self->get_count_dir();
      my @include = $include_file_obj->slurp( chomp=>1 );
      if ( !-e "$count_dir/$include[0]/$param" ) {
	#Check if the gene counts file exists in each needed genome
	MyX::Generic::BadValue->throw(
				      error => "Genome doesn't have count file",
				      value => "",
				     );
      }
    }
    return 1;
  }
  sub _check_min_sample_count { #Look at Scalar::Util::Numeric
    my ($self) = @_;
    my $count = $self->get_min_sample_count();
    my $href  = $self->get_params_href();
    if ( !$href->{min_sample_count} ) {
      MyX::Generic::Undef::Param->throw(
					error => "Param not set",
					usage => "min_sample_count",
					);
    }
    if ( isneg $count || !isint $count) {
      MyX::Generic::Digit->throw(
				 error => "Value can not be zero or negative and must be an integer",
					  );
    }
    return 1;
  }
  sub _check_min_sample_cpm {
    my ($self) = @_;
    my $count = $self->get_min_sample_cpm();
    my $href  = $self->get_params_href();
    if ( !$href->{min_sample_cpm} ) {
      MyX::Generic::Undef::Param->throw(
					error => "Param not set",
					usage => "min_sample_cpm",
					);
    }
    if ( isneg $count ) {
      MyX::Generic::Digit::TooSmall->throw(
					   error => "Value cannot be negative",
					   value => $count,
					   MIN => 0,
					  );
    }
    return 1;
  }
  
  sub _check_test_names {
    my ($self) = @_;
    my $test = $self->get_test_names();
    my $href  = $self->get_params_href();
    if ( !$href->{test_names} ) {
      MyX::Generic::Undef::Param->throw(
					error => "Param not set",
					usage => "test",
				       );
    }
    else {
      #Check the test_column in the meta file for test_names
      my $test_col = $self->get_test_col_name();
      my $metaG_meta_file_obj = file($self->get_metaG_meta_file());
      my @slurped_meta_file = $metaG_meta_file_obj->slurp( chomp => 1, split => qr/\t/);
      my @split_test_names = split /\"/, $test;
      my $col_num = 0;
      my $first_line_aref = $slurped_meta_file[0];
      foreach my $col (@$first_line_aref) {
	if ($col eq $test_col) {
	  last;
	}
	$col_num++;
      }
      my $total_names = 0;
      my $found_names = 0;
      foreach my $test_names (@split_test_names) {
	if ( $test_names =~ /[a-z]/i ) {
	  $total_names++;
	  foreach my $line_aref (@slurped_meta_file) {
	    if ( @$line_aref[$col_num] eq $test_names ) {
	      $found_names++;
	      last;
	    }
	  }
	}
      }
      if ($total_names != $found_names) {
	MyX::Generic::BadValue->throw(
				      error => "test names do not match values under the test_col_name in the metaG_meta_file",
				      value => $test,
				     );
      }
    }#elsif
    return 1;
  }
  sub _check_test_col_name {
    my ($self) = @_;
    my $test_col_name = $self->get_test_col_name();
    my $href  = $self->get_params_href();
    if ( !$href->{test_col_name} ) {
      MyX::Generic::Undef::Param->throw(
					error => "Param not set",
					usage => "test_col_name",
					);
    }
    else {
      my $metaG_meta_file_obj = file( $self->get_metaG_meta_file() );
      my @slurped_metaG_meta_file = $metaG_meta_file_obj->slurp( chomp =>1 );
      if ( $slurped_metaG_meta_file[0] !~ qr/$test_col_name/i ) {
	MyX::Generic::BadValue->throw(
				      error => "test_col_name not found in metaG_meta_file",
				      value => $test_col_name,
				     );
      }
    }
    return 1;
  }
  sub _check_grp_meta_file {
    my ($self) = @_;
    my $meta_file = $self->get_grp_meta_file();
    my $href  = $self->get_params_href();
    if ( !$href->{grp_meta_file} ) {
      MyX::Generic::Undef::Param->throw(
					error => "Param not set",
					usage => "grp_meta_file",
					);
    }
    if ( !-e $meta_file || -z $meta_file ) {
      MyX::Generic::DoesNotExist::File->throw(
					      error => "Metadata file containing grouping information does not exist",
					      file_name => $meta_file,
					     );
    }
    return 1;
  }
  
  #Print
  sub _check_p3_height {
    my ($self) = @_;
    my $href  = $self->get_params_href();
    if ( !$href->{p3_height} ) {
      MyX::Generic::Undef::Param->throw(
					error => "Param not set",
					usage => "p3_height",
					);
    }
    return 1;
  }
  sub _check_ref_meta_cols {
    my ($self) = @_;
    my $param = $self->get_ref_meta_cols();
    my $href  = $self->get_params_href();
    if ( !$href->{ref_meta_cols} ) {
      MyX::Generic::Undef::Param->throw(
					error => "Param not set",
					usage => "ref_meta_cols",
					);
    }
    else {
      #Check to make sure the columns exist in the meta data file
      my $meta_file_obj = file( $self->get_ref_meta_file() );
      my @split_param = split /\"/, $param;
      my @slurped_meta = $meta_file_obj->slurp( chomp => 1);
      my $expected_cols = 0;
      my $actual_cols = 0;
      foreach my $col (@split_param) {
	if ( $col =~ /[a-z]/i ) {
	  $expected_cols++;
	  if ( $slurped_meta[0] =~ qr/$col/i ) {
	    $actual_cols++;
	  }
	}
      }
      if ( $expected_cols != $actual_cols ) {
	MyX::Generic::BadValue->throw(
				      error => "Not all the ref_meta_cols exist in the ref_meta_file",
				      value => $param,
				     );
      }
    }
    return 1;
  }

	sub _check_heat_filter {
		my ($self) = @_;
		my $param = $self->get_heat_filter();
		my $href  = $self->get_params_href();
		if ( !$href->{heat_filter} ) {
			MyX::Generic::Undef::Param->throw(
					error => "Param not set",
					usage => "heat_filter",
					);
		}
		else {
			if ($param =~ /[a-z]/i) {
				if ($param =~ qr/t/i || $param =~ qr/true/i) {
					$param = "TRUE";
					$href->{heat_filter} = $param;
				}
				elsif ($param =~ qr/f/i || $param =~ qr/false/i) {
					$param = "FALSE";
					$href->{heat_filter} = $param;
				}
			}
		}
    if ($param ne "TRUE" && $param ne "FALSE") {
			MyX::Generic::BadValue->throw(
				      error => "Need to pass a boolean value for heat_filter",
				      value => $param,
				     );
    }
		return 1;
  }

  #SETTERS#
  #Global
  sub set_ref_meta_file {
    my ($self,$param) = @_;
    my $params_href = $self->get_params_href();
    my $old = $params_href->{ref_meta_file};
    $params_href->{ref_meta_file} = $param;
    eval {$self->_check_ref_meta_file()};
    if ( my $err = Exception::Class->caught() ) {
      $params_href->{ref_meta_file} = $old;
    }
    return 1;
  }
  sub set_ref_include_file {
    my ($self,$param) = @_;
    my $params_href = $self->get_params_href();
    my $old = $params_href->{ref_include_file};
    $params_href->{ref_include_file} = $param;
    eval {$self->_check_ref_include_file()};
    if ( my $err = Exception::Class->caught() ) {
      $params_href->{ref_include_file} = $old;
    }
    return 1;
  }
  sub set_ref_exclude_file {
    my ($self,$param) = @_;
    my $params_href = $self->get_params_href();
    my $old = $params_href->{ref_exclude_file};
    $params_href->{ref_exclude_file} = $param;
    eval {$self->_check_ref_exclude_file()};
    if ( my $err = Exception::Class->caught() ) {
      $params_href->{ref_exclude_file} = $old;
    }
    return 1;
  }
  sub set_tree {
    my ($self,$param) = @_;
    my $params_href = $self->get_params_href();
    my $old = $params_href->{tree};
    $params_href->{tree} = $param;
    eval {$self->_check_tree()};
    if ( my $err = Exception::Class->caught() ) {
      $params_href->{tree} = $old;
    }
    return 1;
  }
  sub set_out_dir {
    my ($self,$param) = @_;
    my $params_href = $self->get_params_href();
    my $old = $params_href->{out_dir};
    $params_href->{out_dir} = $param;
    eval {$self->_check_out_dir()};
    if ( my $err = Exception::Class->caught() ) {
      $params_href->{out_dir} = $old;
    }
    return 1;
  }
  #EdgeR
  sub set_annote_file_name {
    my ($self,$param) = @_;
    my $params_href = $self->get_params_href();
    my $old = $params_href->{annote_file_name};
    $params_href->{annote_file_name} = $param;
    eval {$self->_check_annote_file_name()};
    if ( my $err = Exception::Class->caught() ) {
      $params_href->{annote_file_name} = $old;
    }
    return 1;
  }
  sub set_grp_genes_by {
    my ($self,$param) = @_;
    my $params_href = $self->get_params_href();
    my $old = $params_href->{grp_genes_by};
    $params_href->{grp_genes_by} = $param;
    eval {$self->_check_grp_genes_by()};
    if ( my $err = Exception::Class->caught() ) {
      $params_href->{grp_genes_by} = $old;
    }
    return 1;
    }
  }
  sub set_gene_id_col {
    my ($self,$param) = @_;
    my $params_href = $self->get_params_href();
    my $old = $params_href->{gene_id_col};
    $params_href->{gene_id_col} = $param;
    eval {$self->_check_gene_id_col()};
    if ( my $err = Exception::Class->caught() ) {
      $params_href->{gene_id_col} = $old;
    }
    return 1;
  }
  sub set_count_dir {
    my ($self,$param) = @_;
    my $params_href = $self->get_params_href();
    my $old = $params_href->{count_dir};
    $params_href->{count_dir} = $param;
    eval {$self->_check_count_dir()};
    if ( my $err = Exception::Class->caught() ) {
      $params_href->{count_dir} = $old;
    }
    return 1;
  }
  sub set_dafe_dir {
    my ($self,$param) = @_;
    my $params_href = $self->get_params_href();
    my $old = $params_href->{dafe_dir};
    $params_href->{dafe_dir} = $param;
    eval {$self->_check_dafe_dir()};
    if ( my $err = Exception::Class->caught() ) {
      $params_href->{dafe_dir} = $old;
    }
    return 1;
  }
  sub set_genome_id_col {
    my ($self,$param) = @_;
    my $params_href = $self->get_params_href();
    my $old = $params_href->{genome_id_col};
    $params_href->{genome_id_col} = $param;
    eval {$self->_check_genome_id_col()};
    if ( my $err = Exception::Class->caught() ) {
      $params_href->{genome_id_col} = $old;
    }
    return 1;
  }
  sub set_metaG_meta_file {
   my ($self,$param) = @_;
    my $params_href = $self->get_params_href();
    my $old = $params_href->{metaG_meta_file};
    $params_href->{metaG_meta_file} = $param;
    eval {$self->_check_metaG_meta_file()};
    if ( my $err = Exception::Class->caught() ) {
      $params_href->{metaG_meta_file} = $old;
    }
    return 1;
 }
  sub set_metaG_include_file {
    my ($self,$param) = @_;
    my $params_href = $self->get_params_href();
    my $old = $params_href->{metaG_include_file};
    $params_href->{metaG_include_file} = $param;
    eval {$self->_check_metaG_include_file()};
    if ( my $err = Exception::Class->caught() ) {
      $params_href->{metaG_include_file} = $old;
    }
    return 1;
  }
  sub set_metaG_exclude_file {
    my ($self,$param) = @_;
    my $params_href = $self->get_params_href();
    my $old = $params_href->{metaG_exclude_file};
    $params_href->{metaG_exclude_file} = $param;
    eval {$self->_check_metaG_exclude_file()};
    if ( my $err = Exception::Class->caught() ) {
      $params_href->{metaG_exclude_file} = $old;
    }
    return 1;
  }
  sub set_count_file_name {
    my ($self,$param) = @_;
    my $params_href = $self->get_params_href();
    my $old = $params_href->{count_file_name};
    $params_href->{count_file_name} = $param;
    eval {$self->_check_count_file_name()};
    if ( my $err = Exception::Class->caught() ) {
      $params_href->{count_file_name} = $old;
    }
    return 1;
  }
  sub set_min_sample_count {
    my ($self,$param) = @_;
    my $params_href = $self->get_params_href();
    my $old = $params_href->{min_sample_count};
    $params_href->{min_sample_count} = $param;
    eval {$self->_check_min_sample_count()};
    if ( my $err = Exception::Class->caught() ) {
      $params_href->{min_sample_count} = $old;
    }
    return 1;
  }
  sub set_min_sample_cpm {
    my ($self,$param) = @_;
    my $params_href = $self->get_params_href();
    my $old = $params_href->{min_sample_cpm};
    $params_href->{min_sample_cpm} = $param;
    eval {$self->_check_min_sample_cpm()};
    if ( my $err = Exception::Class->caught() ) {
      $params_href->{min_sample_cpm} = $old;
    }
    return 1;
  }
  sub set_test_names {
    my ($self,$param) = @_;
    my $params_href = $self->get_params_href();
    my $old = $params_href->{test_names};
    $params_href->{test_names} = $param;
    eval {$self->_check_test_names()};
    if ( my $err = Exception::Class->caught() ) {
      $params_href->{test_names} = $old;
    }
    return 1;
  }
  sub set_test_col_name {
    my ($self,$param) = @_;
    my $params_href = $self->get_params_href();
    my $old = $params_href->{test_col_name};
    $params_href->{test_col_name} = $param;
    eval {$self->_check_test_col_name()};
    if ( my $err = Exception::Class->caught() ) {
      $params_href->{test_col_name} = $old;
    }
    return 1;
  }
  sub set_grp_meta_file {
    my ($self,$param) = @_;
    my $params_href = $self->get_params_href();
    my $old = $params_href->{grp_meta_file};
    $params_href->{grp_meta_file} = $param;
    eval {$self->_check_grp_meta_file()};
    if ( my $err = Exception::Class->caught() ) {
      $params_href->{grp_meta_file} = $old;
    }
    return 1;
  }
  
  #Print
  sub set_p3_height {
    my ($self,$param) = @_;
    my $params_href = $self->get_params_href();
    my $old = $params_href->{p3_height};
    $params_href->{p3_height} = $param;
    eval {$self->_check_p3_height()};
    if ( my $err = Exception::Class->caught() ) {
       $params_href->{p3_height} = $old;
     }
    return 1;
  }
  sub set_ref_meta_cols {
    my ($self,$param) = @_;
    my $params_href = $self->get_params_href();
    my $old = $params_href->{ref_meta_cols};
    $params_href->{ref_meta_cols} = $param;
    eval {$self->_check_ref_meta_cols()};
    if ( my $err = Exception::Class->caught() ) {
      $params_href->{ref_meta_cols} = $old;
    }
    return 1;
  }
  sub set_heat_filter {
    my ($self,$param) = @_;
    my $params_href = $self->get_params_href();
    my $old = $params_href->{heat_filter};
    $params_href->{heat_filter} = $param;
    eval {$self->_check_heat_filter()};
    if ( my $err = Exception::Class->caught() ) {
      $params_href->{heat_filter} = $old;
    }
    return 1;
  }

  #GETTERS#
  #Global
  sub get_ref_meta_file {
    my ($self) = @_;
    my $params_href = $self->get_params_href;
    return $params_href->{ref_meta_file};
  }
  sub get_ref_include_file {
    my ($self) = @_;
    my $params_href = $self->get_params_href;
    return $params_href->{ref_include_file};
  }
  sub get_ref_exclude_file {
    my ($self) = @_;
    my $params_href = $self->get_params_href();
    return $params_href->{ref_exclude_file};
  }
  sub get_tree {
    my ($self) = @_;
    my $params_href = $self->get_params_href();
    return $params_href->{tree};
  }
  sub get_out_dir {
    my ($self) = @_;
    my $params_href = $self->get_params_href;
    return $params_href->{out_dir};
  }

  #EdgeR
  sub get_annote_file_name {
    my ($self) = @_;
    my $params_href = $self->get_params_href;
    return $params_href->{annote_file_name};
  }
  sub get_grp_genes_by {
    my ($self) = @_;
    my $params_href = $self->get_params_href;
    return $params_href->{grp_genes_by};
  }
  sub get_gene_id_col {
    my ($self) = @_;
    my $params_href = $self->get_params_href;
    return $params_href->{gene_id_col};
  }
  sub get_count_dir {
    my ($self) = @_;
    my $params_href = $self->get_params_href;
    return $params_href->{count_dir};
  }
  sub get_dafe_dir {
    my ($self) = @_;
    my $params_href = $self->get_params_href;
    return $params_href->{dafe_dir};
  }
  sub get_genome_id_col {
    my ($self) = @_;
    my $params_href = $self->get_params_href;
    return $params_href->{genome_id_col};
  }
  sub get_metaG_meta_file {
    my ($self) = @_;
    my $params_href = $self->get_params_href;
    return $params_href->{metaG_meta_file};
  }
  sub get_metaG_include_file {
    my ($self) = @_;
    my $params_href = $self->get_params_href;
    return $params_href->{metaG_include_file};
  }
  sub get_metaG_exclude_file {
    my ($self) = @_;
    my $params_href = $self->get_params_href;
    return $params_href->{metaG_exclude_file};
  }
  sub get_count_file_name {
    my ($self) = @_;
    my $params_href = $self->get_params_href;
    return $params_href->{count_file_name};
  }
  sub get_min_sample_count {
    my ($self) = @_;
    my $params_href = $self->get_params_href;
    return $params_href->{min_sample_count};
  }
  sub get_min_sample_cpm {
    my ($self) = @_;
    my $params_href = $self->get_params_href;
    return $params_href->{min_sample_cpm};
  }
  sub get_test_names {
    my ($self) = @_;
    my $params_href = $self->get_params_href;
    return $params_href->{test_names};
  }
  sub get_test_col_name {
    my ($self) = @_;
    my $params_href = $self->get_params_href();
    return $params_href->{test_col_name};
  }
  sub get_grp_meta_file {
    my ($self) = @_;
    my $params_href = $self->get_params_href();
    return $params_href->{grp_meta_file};
  }
  
  #Print
  sub get_p3_height {
    my ($self) = @_;
    my $params_href = $self->get_params_href();
    return $params_href->{p3_height};
  }
  sub get_ref_meta_cols {
    my ($self) = @_;
    my $params_href = $self->get_params_href();
    return $params_href->{ref_meta_cols};
  }
  sub get_heat_filter {
    my ($self) = @_;
    my $params_href = $self->get_params_href();
    return $params_href->{heat_filter};
  }

1;
  
__END__

=head1 Paramater Handler Object

This is a paramater object that handles parameters submitted into the DAFE_count
program. It combines the EdgeR and test params.

head1 VERSION

This documentation refers to Justify 0.0.1

=head1 Included Modules

use Getopt::Long;
use Pod::Usage;
use Carp;
use Readonly;
use Path::Class;
use Data::Dumper;
use XML::Simple qw(:strict);
use Class::Std::Utils;
use Scalar::Util::Numeric qw(isneg isint isfloat);
use MyX::Generic;
use YAML::XS qw(LoadFile);

=head1 Inherit
    
    NA
		
=head1 SYNOPSIS

	=head2 Object Creation
	
		use Param_handler;
		my $param_obj = Param_handler->new( {href => $attributes} );
										OR
										Param_handler->new( {xml_file => $attribute_xml} );
										OR
										Param_handler->new();
		$param_obj->set_params({href => $attributes});
		OR
		$param_obj->set_params({xml_file => $attribute_xml});
		$param_obj->check_edger_params();
		$param_obj->check_print_params();
		
		$param_obj->spew_xml_file($out_file_path);
		$param_obj->spew_out();
		
	=head2 Getters
	
		$param_href         = $param_obj->get_params_href();
		$ref_meta_file      = $param_obj->get_ref_meta_file();
		$ref_include_file   = $param_obj->get_ref_include_file();
		$ref_exclude_file   = $param_obj->get_ref_exclude_file();
		$tree               = $param_obj->get_tree();
		$out_dir            = $param_obj->get_out_dir();
		$annote_file_name   = $param_obj->get_annote_file_name();
		$grp_genes_by       = $param_obj->get_grp_genes_by();
		$gene_id_col        = $param_obj->get_gene_id_col();
		$count_dir          = $param_obj->get_count_dir();
		$dafe_dir           = $param_obj->get_dafe_dir();
		$genome_id_col      = $param_obj->get_genome_id_col();
		$metaG_meta_file    = $param_obj->get_metaG_meta_file();
		$metaG_include_file = $param_obj->get_metaG_include_file();
		$metaG_exclude_file = $param_obj->get_metaG_exclude_file();
		$count_file_name    = $param_obj->get_count_file_name();
		$min_sample_count   = $param_obj->get_sample_count();
		$min_sample_cpm     = $param_obj->get_min_sample_cpm();
		$test_names         = $param_obj->get_test_names();
		$test_col_name      = $param_obj->get_test_col_name();
		$grp_meta_file      = $param_obj->get_grp_meta_file();
		$p3_height          = $param_obj->get_p3_height();
		$ref_meta_cols      = $param_obj->get_ref_meta_cols();
		$heat_filter        = $param_obj->get_heat_filter();
		
	=head2 Setters
	
		$param_obj->set_params_href();
		$param_obj->set_ref_meta_file();
		$param_obj->set_ref_include_file();
		$param_obj->set_ref_exclude_file();
		$param_obj->set_tree();
		$param_obj->set_out_dir();
		$param_obj->set_annote_file_name();
		$param_obj->set_grp_genes_by();
		$param_obj->set_gene_id_col();
		$param_obj->set_count_dir();
		$param_obj->set_dafe_dir();
		$param_obj->set_genome_id_col();
		$param_obj->set_metaG_meta_file();
		$param_obj->set_metaG_include_file();
		$param_obj->set_metaG_exclude_file();
		$param_obj->set_count_file_name();
		$param_obj->set_sample_count();
		$param_obj->set_min_sample_cpm();
		$param_obj->set_test_names();
		$param_obj->set_test_col_name();
		$param_obj->set_grp_meta_file();
		$param_obj->set_p3_height();
		$param_obj->set_ref_meta_cols();
		$param_obj->set_heat_filter();
		
=head1  DESCRIPTION

Param_handler is used for keeping track of the parameters needed for full
analysis of the metagenome data with data from genome directories. The parameters
that need to be passed into the object include: A reference meta file, reference
include and exclude file, a tree file, out directory, the annotation file names,
what to group the genes by, the gene id column, the count directory, minimal
sample count, the minimal sample counts per million, the names of the things to
test, the test column names, the group metafile, the p3 height, refrence meta
files, and heat filter.

=head1 METHODS

	new()
	get_params_href()
	set_params()
	check_edger_params()
	check_print_params()
	spew_xml_file()
	spew_out()
	get_print_params_href()
	get_edger_params_href()
	get_params_href();
	get_ref_meta_file();
	get_ref_include_file();
	get_ref_exclude_file();
	get_tree();
	get_out_dir();
	get_annote_file_name();
	get_grp_genes_by();
	get_gene_id_col();
	get_count_dir();
	get_dafe_dir();
	get_genome_id_col();
	get_metaG_meta_file();
	get_metaG_include_file();
	get_metaG_exclude_file();
	get_count_file_name();
	get_sample_count();
	get_min_sample_cpm();
	get_test_names();
	get_test_col_name();
	get_grp_meta_file();
	get_p3_height();
	get_ref_meta_cols();
	get_heat_filter();
	Setters ar the same but with set_ not get_

=head1 METHODS DESCRIPTION

=head2 new()

	Title:		new
	Usage:		Param_handler->new( $arg_href );
	Function:	This creates a new instance of the Param_handler object
	Returns:	Param_handler
	Args:			$arg_href => A hash ref that either has a defined href or xml file.
	Throws:		croak
	Comments:	Make sure you have either a xml file or hash ref in the arg_ref
	See Also: NA
	
=head2 get_params_href()

	Title:		get_params_href
	Usage:		$param_href = $param_obj->get_params_href();
	Function:	This gets the stored hashref containing all the parameters that are
						needed for the metagenome analysis.
	Returns:	Hash Reference
	Args:			- $param_obj => a bless param object
	Throws:		NA
	Comments:	NA
	See Also:	NA
	
=head2 set_params()

	Title:		set_params	
	Usage:		$param_obj->set_params( $arg_href );
	Function:	Will set any parameters that you pass into this method. Can pass in
						a xml file or hash reference
	Returns:	NA
	Args:			- $param_obj => a blessed Param_handler object
						- $arg_href => a hash ref that either has a defined href or xml_file
							key
	Throws:		croak
	Comments:	Make sure to run the checks on the newly set parameters so that you
						are sure everything is ok with the passed parameters
	See Also:	NA
	
=head2  check_edger_params()

	Title:		check_edger_params
	Usage:		$param_obj->check_edger_params();
	Function:	Will check all the edger parameters
	Returns:	NA
	Args:			- $param_obj => a blessed Param_handler object
	Throws:		MyX::Generic::Undef::Param
						MyX::Generic::DoesNotExist::File
						MyX::Generic::BadValue
						MyX::Generic::Digit::TooSmall
						MyX::Generic::Digit
						MyX::Generic::DoesNotExist::Dir
						Croak
	Comments:	Need to call this or check_print_params to make sure the paremeters
						are all good.
	See Also:	NA
	
=head2		check_print_params()

	Title:		check_print_params
	Usage:		$param_obj->check_print_params();
	Function:	This will check all the parameters assocociated with the print
						object
	Returns:	NA
	Args:			-	$param_obj => a blessed Param_handler object
	Throws:		MyX::Generic::Undef::Param
						MyX::Generic::Undef::Attribute
						MyX::Generic::BadValue
						MyX::Generic::DoesNotExist::File
	Comments:	Need to call this or check_edger_params to make sure the parameters
						are all correct
	See Also:	NA
	
=head2		spew_xml_file()

	Title:		spew_xml_file
	Usage:		$param_obj->spew_xml_file( $outfile_path );
	Function:	Creates an xml file of all the parameters stored in the object
	Returns:	NA
	Args:			-	$outfile_path => the path you want the file to be created at
	Throws:		croak
	Comments:	NA
	See Also:	NA

=head2		spew_out()

	Title:		spew_out
	Usage:		$param_obj->spew_out()
	Function: prints the parameters to the terminal window
	Returns:	NA
	Args:			NA
	Throws:		NA
	Comments:	NA
	See Also:	NA
	
=head2		get_print_params_href()

	Title:		get_print_params_href
	Usage:		$print_href = $param_obj->get_print_params_href();
	Function:	This will return a hash reference containing only the parameters
						needed for the print part of analysis
	Returns:	Hash Reference
	Args:			- $param_obj => a blessed Param_handler object
	Throws:		NA
	Comments:	NA
	See Also:	NA
	
=head2		get_edger_params_href()

	Title:		get_edger_params_href
	Usage:		$edger_href = $param_obj->get_edger_param_href()
	Function:	Returns the parameters needed for the edger analysis
	Returns:	Hash Reference
	Args:			-	$param_obj = A blessed Param_handler object
	Throws:		NA	
	Comments:	NA
	See Also:	NA
	
=head2		set_{	ref_meta_file, ref_include_file, ref_exclude_file, tree
								out_dir, annote_file_name, grp_genes_by, gene_id_col,
								count_dir, dafe_dir, genome_id_col, metaG_meta_file,
								metaG_include_file, metaG_exclude_file, count_file_name,
								min_sample_count, min_sample_cpm, test_names, test_col_name,
								grp_meta_file, p3_height, ref_meta_cols, heat_filter }()
								
	Title:		Individual Setters
	Usage:		$param_obj->set_param_name( $param );
	Function:	Will set an individual parameter only if it passes the check for
						that parameter.
	Returns:	NA
	Args:			- $param => the parameter value you want to set the param to
	Throws:		Whatever the parameters check throws
	Comments:	NA
	See Also:	NA
	
=head2		get_{ ref_meta_file, ref_include_file, ref_exclude_file, tree
								out_dir, annote_file_name, grp_genes_by, gene_id_col,
								count_dir, dafe_dir, genome_id_col, metaG_meta_file,
								metaG_include_file, metaG_exclude_file, count_file_name,
								min_sample_count, min_sample_cpm, test_names, test_col_name,
								grp_meta_file, p3_height, ref_meta_cols, heat_filter }()
								
	Title:		Individual Getters
	Usage:		$param = $param_object->set_param_name();
	Function:	Will return the stored value for that param in the Param_handler
						object
	Returns:	String
	Args:			NA
	Throws:		NA
	Comments:	NA
	See Also:	NA
	
=head1 CONFIGURATION AND ENVIRONMENT

    Param_handler requires no configuration files or environment variables.

=head1 DEPENDENCIES

    NA
    
=head1 INCOMPATIBILITIES

    None reported.

=head1 BUGS AND LIMITATIONS

No bugs have been reported.

Please report any bugs or feature requests	
	
=head1 Author

Nicholas Colaianni
contact via C<< <ncolaian@live.unc.edu> >>

=head1 LICENCE AND COPYRIGHT

Copyright (c) 2016, Nicholas Colaianni
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
