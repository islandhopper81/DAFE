# This file was originally copied from cog_model_obj.  The point is to be able
# to have these functions extend to any type of clustering schema like CD-Hit
# or KO, etc.

create_clstr_model_obj = function(clstr_annot_tbl, test, params_obj) {
  clstr_model_obj = list(mat = NULL,
                       gpc_mat = NULL,
                       clstr_annot_tbl = NULL,
                       samples = NULL,
                       groups = NULL,
                       test = NULL)
  
  clstr_model_obj$clstr_annot_tbl = clstr_annot_tbl
  clstr_model_obj$samples = row.names(params_obj$metaG_metadata_obj)
  clstr_model_obj$groups = params_obj$metaG_metadata_obj$fraction
  clstr_model_obj$test = test
  
  # create the mat and gpc_mat objects (they are currently empty)
  r = nrow(clstr_annot_tbl)
  c = length(params_obj$ref_include_ids)
  n = list(c(), c(params_obj$ref_include_ids))
  clstr_model_obj$mat = matrix(, nrow=r, ncol=c, dimnames=n)
  clstr_model_obj$gpc_mat = matrix(, nrow=r, ncol=c, dimnames=n)
  
  class(clstr_model_obj) = "clstr_model_obj"
  return(clstr_model_obj)
}

### removes clstrs that are non-detectable or not present across all genomes
clean_clstrs = function(clstr_model_obj) {
  mat = clstr_model_obj$mat
  
  # remove clstrs that are never present in any genomes
  mat = as.matrix(mat, rownames.force = row.names(mat))
  to_keep = apply(mat, 1, sum) > (-2*ncol(mat))
  
  clstr_model_obj$mat = clstr_model_obj$mat[to_keep,]
  clstr_model_obj$gpc_mat = clstr_model_obj$gpc_mat[to_keep,]
  
  return(clstr_model_obj)
}

### distinguish between absent and non-detectable clstrs
classify_non_detect_clstrs = function(clstr_model_obj) {
  # There is an important distinction between a clstr that is absent and
  # a clstr that is not detectable with our data.  This function helps
  # distinguish between these two factors by specifying the absent 
  # clstrs with a -3 value in tbl.  So -3 == absent while -2 == not detectable
  
  tbl = clstr_model_obj$mat
  gpc = clstr_model_obj$gpc_mat
  
  abs_or_non_detec = tbl == -2
  absent = abs_or_non_detec & gpc == 0  # gets bool matrix where absent clstrs == TRUE
  tbl[absent] = -3
  
  clstr_model_obj$mat = tbl
  
  return(clstr_model_obj)
}

# removes genomes (ie columns from a )
remove_genomes = function(clstr_model_obj, genomes_to_remove) {
  clstr_model_obj$mat = clstr_model_obj$mat[, !colnames(clstr_model_obj$mat) %in% genomes_to_remove]
  clstr_model_obj$gpc_mat = clstr_model_obj$gpc_mat[, !colnames(clstr_model_obj$gpc_mat) %in% genomes_to_remove]
  return(clstr_model_obj)
}

# keep genomes (ie columns from a)
keep_genomes = function(clstr_model_obj, genomes_to_keep) {
  clstr_model$obj$mat = clstr_model_obj$mat[, colnames(clstr_model_obj$mat) %in% genomes_to_keep]
  clstr_model_obj$gpc_mat = clstr_model_obj$gpc_mat[, colnames(clstr_model_obj$gpc_mat) %in% genomes_to_keep]
  return(clstr_model_obj)
}


