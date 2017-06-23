#!/usr/bin/env Rscript

# This function saves a params_obj

# NOTES:
# 1.The test parameter should follow the edgeR specs for the groups parameter
#   in the dge function.  eg ("BK", "RZ").  The base group should go first. In
#   example a positive logFC means that observation was higher in RZ than BK.

create_edgeR_params_obj = function(yaml_file, set_env=T) {
  params = list(working_root_dir = NULL,
                ref_include_file = NULL,
                ref_include_ids = NULL,
                ref_exclude_file = NULL,
                ref_exclude_ids = NULL,
                metaG_meta_file = NULL,
                metaG_include_file = NULL,
                metaG_exclude_file = NULL,
                count_file_name = NULL,
                agg_count_file_prefix = NULL,
                grp_genes_by = NULL,
                min_sample_count = NULL,
                min_sample_cpm = NULL,
                test = NULL,
                test_col_name = NULL,
                test_groups = NULL)
  
  require(yaml)

  # load the params from the yaml input file
  params_yaml = yaml.load_file(yaml_file)
  
  # save the params from the yaml input file
  params$working_root_dir = params_yaml$count_dir
  params$ref_include_file = params_yaml$ref_include_file
  params$ref_exclude_file = params_yaml$ref_exclude_file
  params$metaG_meta_file = params_yaml$metaG_meta_file
  params$metaG_include_file = params_yaml$metaG_include_file
  params$metaG_exclude_file = params_yaml$metaG_exclude_file
  params$count_file_name = params_yaml$count_file_name
  params$grp_genes_by = params_yaml$grp_genes_by
  params$min_sample_cpm = params_yaml$min_sample_cpm
  params$min_sample_count = params_yaml$min_sample_count
  params$test = params_yaml$test
  params$test_col_name = params_yaml$test_col_name
  
   # set the reference ids
  params$ref_include_ids = read_names(params$ref_include_file)
  params$ref_exclude_ids = read_names(params$ref_exclude_file)
  params$ref_include_ids = params$ref_include_ids[
    !params$ref_include_ids %in% params$ref_exclude_ids]

  # set the metaG metadata obj ids
  params$metaG_metadata_tbl = create_metaG_metadata_tbl(
    metadata_file = params$metaG_meta_file, 
    include_file = params$metaG_include_file, 
    exclude_file = params$metaG_exclude_file
    )
  
  # set the test_groups
  params$test_groups = get_test_groups(params)

  # set the count_file_prefix
  params$agg_count_file_prefix = get_agg_count_file_prefix(params)
  
  # make and return the params_obj as a class object
  class(params) = "edgeR_params_obj"
  return(params)
}

get_agg_count_file_prefix = function(params_obj) {
  prefix = gsub(".txt", "", params_obj$count_file_name)
  prefix = paste(prefix, "_", params_obj$grp_genes_by, "_agg", sep="")
  
  return(prefix)
}

get_test_groups = function(params) {
  # check if the test_col_name in the params file is one of the headers
  # of the metaG_metadata_tbl (which is a data frame)
  if ( params$test_col_name %in% colnames(params$metaG_metadata_tbl) ) {
    index = which(colnames(params$metaG_metadata_tbl) == params$test_col_name)
  } else {
    # throw a warning about how the test_col_name does not match a column
    # in the metagenome_metadata file
    message = paste("ERRORa: test_col_name (", 
                    params$test_col_name, 
                    ") is not found in the metaG_meta_file" )
    stop(message)
  }
  
  # get the groups based on the found index from above
  test_groups = params$metaG_metadata_tbl[,index]
  
  # check that each value in the test_col_name column is included in the 
  # test vector
  if ( ! (setequal(unique(test_groups), params$test)) ) {
    message = paste("ERROR: test value params don't match values in test_col_name",
                    paste("test_col_name values:", 
                    paste(unique(test_groups), collapse=",")),
                    paste("test parameter values:",
                    paste(params$test, collapse=",")),
                    sep="\n",
                    collapse="")
    
    stop(message)
  }
  
  return(test_groups)
}


params_obj.write = function(params_obj, file) {
  write(as.yaml(params_obj), file = file)
}

read_names = function(file) {
  # this is used to read the names in the exclude and include file
  
  out = tryCatch({
    names = vector()
    if ( file.info(file)$size > 0 ) {
      names = read.table(file, header=F)$V1
    } else {
      names = NA
    }
  }, error = function(err){
    print(paste("ERROR: ", err$message))
    return(NA)
  })
  
  return(names)
}

is_set = function(val) {
  if ( is.null(val) ) {
    return(FALSE)
  }
  if ( is.na(val) ) {
    return(FALSE)
  }
  if ( val == "" ) {
    return(FALSE)
  }
  if ( val == "NA" ) {
    return(FALSE)
  }
  if ( val == "NULL" ) {
    return(FALSE)
  }
  
  return(TRUE)
}
