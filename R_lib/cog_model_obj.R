

# NOTE: this object also contains information about how to build the
# edgeR tests

create_cog_model_obj = function(cog_annot_tbl, test, params_obj) {
  cog_model_obj = list(mat = NULL,
                   gpc_mat = NULL,
                   cog_annot_tbl = NULL,
                   samples = NULL,
                   groups = NULL,
                   test = NULL)
  
  cog_model_obj$cog_annot_tbl = cog_annot_tbl
  cog_model_obj$samples = row.names(params_obj$metaG_metadata_obj)
  cog_model_obj$groups = get_groups(params_obj)
  cog_model_obj$test = params_obj$test
  
  # create the mat and gpc_mat objects (they are currently empty)
  r = nrow(cog_annot_tbl)
  c = length(params_obj$ref_include_ids)
  n = list(c(), c(params_obj$ref_include_ids))
  cog_model_obj$mat = matrix(, nrow=r, ncol=c, dimnames=n)
  cog_model_obj$gpc_mat = matrix(, nrow=r, ncol=c, dimnames=n)
  
  class(cog_model_obj) = "cog_model_obj"
  return(cog_model_obj)
}

### gets the group vector that is used by edgeR to assign each
# sample to a group for testing.  eg ("BK", "BK", "RZ", "RZ")
# The values in this vector must match one of the values in the
# params$test vector.
get_groups = function(params_obj) {
  # check if the test_col_name in the params file is one of the headers
  # of the metaG_metadata_obj (which is a data frame)
  if ( params_obj$test_col_name %in% colnames(params_obj$metaG_metadata_obj) ) {
    index = which(colnames(params_obj$metaG_metadata_obj) == params_obj$test_col_name)
  } else {
    # throw a warning about how the test_col_name does not match a column
    # in the metagenome_metadata file
    message = paste("ERROR: test_col_name (", 
                    params_obj$test_col_name, 
                    ") is not found in the metaG_meta_file" )
    stop(message)
  }
  
  # get the groups based on the found index from above
  groups = params_obj$metaG_metadata_obj[,index]
  
  # check that each value in the test_col_name column is included in the 
  # test vector
  if ( ! (setequal(unique(groups), params_obj$test)) ) {
    message = paste("ERROR: test value params don't match values in test_col_name",
                    paste("test_col_name values:", 
                          paste(groups, collapse=",")),
                    paste("test parameter values:",
                          paste(params_obj$test, collapse=",")),
                    sep="\n",
                    collapse="")
    
    stop(message)
  }
  
  return(groups)
}

### removes COGs that are non-detectable or not present across all genomes
clean_cogs = function(cog_model_obj) {
  mat = cog_model_obj$mat
  
  # remove cogs that are never present in any genomes
  mat = as.matrix(mat, rownames.force = row.names(mat))
  to_keep = apply(mat, 1, sum) > (-2*ncol(mat))
  
  cog_model_obj$mat = cog_model_obj$mat[to_keep,]
  cog_model_obj$gpc_mat = cog_model_obj$gpc_mat[to_keep,]
  
  return(cog_model_obj)
}

### distinguish between absent and non-detectable COGs
classify_non_detect_cogs = function(cog_model_obj) {
  # There is an important distinction between a COG that is absent and
  # a COG that is not detectable with our data.  This function helps
  # distinguish between these two factors by specifying the absent 
  # COGs with a -3 value in tbl.  So -3 == absent while -2 == not detectable
  
  tbl = cog_model_obj$mat
  gpc = cog_model_obj$gpc_mat
  
  abs_or_non_detec = tbl == -2
  absent = abs_or_non_detec & gpc == 0  # gets bool matrix where absent cogs == TRUE
  tbl[absent] = -3
  
  cog_model_obj$mat = tbl
  
  return(cog_model_obj)
}

# removes genomes (ie columns from a )
remove_genomes = function(cog_model_obj, genomes_to_remove) {
  cog_model_obj$mat = cog_model_obj$mat[, !colnames(cog_model_obj$mat) %in% genomes_to_remove]
  cog_model_obj$gpc_mat = cog_model_obj$gpc_mat[, !colnames(cog_model_obj$gpc_mat) %in% genomes_to_remove]
  return(cog_model_obj)
}

# keep genomes (ie columns from a)
keep_genomes = function(cog_model_obj, genomes_to_keep) {
  cog_model$obj$mat = cog_model_obj$mat[, colnames(cog_model_obj$mat) %in% genomes_to_keep]
  cog_model_obj$gpc_mat = cog_model_obj$gpc_mat[, colnames(cog_model_obj$gpc_mat) %in% genomes_to_keep]
  return(cog_model_obj)
}

