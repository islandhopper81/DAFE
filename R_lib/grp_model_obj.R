

# NOTE: this object also contains information about how to build the
# edgeR tests

###
# grp_metadata_tbl--a table with a row for each group(eg cluster, cog, etc).  
###
create_grp_model_obj = function(grp_metadata_tbl, params_obj) {
  grp_model_obj = list(mat = NULL,
                   gpg_mat = NULL,
                   grp_metadata_tbl = NULL,
                   samples = NULL,
                   test_groups = NULL,
                   test = NULL)
  
  grp_model_obj$grp_metadata_tbl = grp_metadata_tbl
  grp_model_obj$samples = row.names(params_obj$metaG_metadata_obj)
  grp_model_obj$test_groups = get_test_groups(params_obj)
  grp_model_obj$test = params_obj$test
  
  # create the mat and gpg_mat objects (they are currently empty)
  r = nrow(grp_metadata_tbl)
  c = length(params_obj$ref_include_ids)
  n = list(c(), c(params_obj$ref_include_ids))
  grp_model_obj$mat = matrix(, nrow=r, ncol=c, dimnames=n)
  grp_model_obj$gpg_mat = matrix(, nrow=r, ncol=c, dimnames=n)
  
  class(grp_model_obj) = "grp_model_obj"
  return(grp_model_obj)
}

###
# load_grp_model_obj - an alterntative constructor that makes a grp_model_obj
# from a previously generated count table.  The count table should have a row
# for each group and a column for each genome. Each cell should be a DA value
# (ie -3,-2,-1,0,1)
###
load_grp_model_obj = function(count_tbl_file, params_obj) {
  grp_model_obj = list(mat = NULL,
                       gpg_mat = NULL,
                       grp_metadata_tbl = NULL,
                       samples = NULL,
                       test_groups = NULL,
                       test = NULL)
  
  # load the matrix
  count_tbl = read.table(count_tbl_file, header=T, sep="\t", check.names=F)
  
  grp_model_obj$grp_metadata_tbl = as.data.frame(colnames(count_tbl))
  grp_model_obj$samples = row.names(params_obj$metaG_metadata_obj)
  grp_model_obj$test_groups = get_test_groups(params_obj)
  grp_model_obj$test = params_obj$test
  
  # create the mat and object
  r = nrow(count_tbl)
  c = ncol(count_tbl)
  n = list(row.names(count_tbl), colnames(count_tbl))
  grp_model_obj$mat = as.matrix(count_tbl)
  grp_model_obj$gpg_mat = matrix(, nrow=r, ncol=c, dimnames=n) # this one is empty
  
  class(grp_model_obj) = "grp_model_obj"
  return(grp_model_obj)
}

### gets the edgeR group vector that is used by edgeR to assign each
# sample to an edgeR group for statistical testing.  eg ("BK", "BK", "RZ", "RZ")
# The values in this vector must match one of the values in the
# params$test vector.
get_test_groups = function(params_obj) {
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
  test_groups = params_obj$metaG_metadata_obj[,index]
  
  # check that each value in the test_col_name column is included in the 
  # test vector
  if ( ! (setequal(unique(test_groups), params_obj$test)) ) {
    message = paste("ERROR: test value params don't match values in test_col_name",
                    paste("test_col_name values:", 
                          paste(test_groups, collapse=",")),
                    paste("test parameter values:",
                          paste(params_obj$test, collapse=",")),
                    sep="\n",
                    collapse="")
    
    stop(message)
  }
  
  return(test_groups)
}

### removes annotation groups that are non-detectable or not present across all genomes
clean_grps = function(grp_model_obj) {
  mat = grp_model_obj$mat
  
  # remove annotaiton groups that are never present in any genomes
  mat = as.matrix(mat, rownames.force = row.names(mat))
  to_keep = apply(mat, 1, sum) > (-2*ncol(mat))
  
  grp_model_obj$mat = grp_model_obj$mat[to_keep,]
  grp_model_obj$gpg_mat = grp_model_obj$gpg_mat[to_keep,]
  
  return(grp_model_obj)
}

### distinguish between absent and non-detectable annotation groups
classify_non_detect_grps = function(grp_model_obj) {
  # There is an important distinction between an annotation group that is 
  # absent and an annotation group that is not detectable with our data.  
  # This function helps distinguish between these two factors by specifying 
  # the absent groups with a -3 value in tbl.  So -3 == absent while 
  # -2 == not detectable
  
  print("inside function")
  
  tbl = grp_model_obj$mat
  gpg = grp_model_obj$gpg_mat
  
  print("tbl summary")
  print(summary(tbl))
  print("gpg summary")
  print(summary(gpg))
  
  abs_or_non_detec = tbl == -2
  tmp = gpg == 0
  print(paste("abs_or_non_detec dim: ", dim(abs_or_non_detec)))
  print(paste("tmp dim: ", dim(tmp)))
  absent = abs_or_non_detec & gpg == 0  # gets bool matrix where absent grps == TRUE
  tbl[absent] = -3
  
  grp_model_obj$mat = tbl
  
  return(grp_model_obj)
}

# removes genomes (ie columns from a )
remove_genomes = function(grp_model_obj, genomes_to_remove) {
  grp_model_obj$mat = grp_model_obj$mat[, !colnames(grp_model_obj$mat) %in% genomes_to_remove]
  grp_model_obj$gpg_mat = grp_model_obj$gpg_mat[, !colnames(grp_model_obj$gpg_mat) %in% genomes_to_remove]
  return(grp_model_obj)
}

# keep genomes (ie columns from a)
keep_genomes = function(grp_model_obj, genomes_to_keep) {
  grp_model$obj$mat = grp_model_obj$mat[, colnames(grp_model_obj$mat) %in% genomes_to_keep]
  grp_model_obj$gpg_mat = grp_model_obj$gpg_mat[, colnames(grp_model_obj$gpg_mat) %in% genomes_to_keep]
  return(grp_model_obj)
}

