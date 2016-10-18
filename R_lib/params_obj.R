# This function saves a params_obj

# NOTES:
# 1.The test parameter should follow the edgeR specs for the groups parameter
#   in the dge function.  eg ("BK", "RZ").  The base group should go first. In
#   example a positive logFC means that observation was higher in RZ than BK.

# annote_file_name - name of the annotation file in the dafe database (all_annote.txt)
# group_genes_by - name of column in the annotation file for which to group genes
# grp_metadata_file - file with 3 columns -- 1) grpID, 2) higher group, 3) desc
#                   - the most important is column 1 because that is the full list of IDs 
# feature - if this is set to "GENE" it will not group the genes into higher groups
# skip_big_mat -- when this parameter is set to TRUE the DA_test.R script will skip adding
#                 the vector of DA calls to a large full matrix.  This large matrix is
#                 extremely memory and time intensive when there are a large number of
#                 clusters.

create_params_obj = function(yaml_file, set_env=T) {
  params = list(run_name = NULL,
                run_id = NULL,
                feature = NULL,
                annote_file_name = "all_annote.txt",
                group_genes_by = NULL,
                gene_id_col = NULL,
                data_dir = NULL,
                dafe_dir = NULL,
                out_root_dir = NULL,
                ref_out_dir = NULL,
                ref_meta_file = NULL,
                ref_meta_cols = NULL,
                ref_include_file = NULL,
                ref_exclude_file = NULL,
                ref_include_ids = NULL,
                ref_exclude_ids = NULL,
                justify_ref_meta_by = NULL,
                metaG_meta_file = NULL,
                metaG_include_file = NULL,
                metaG_exclude_file = NULL,
                tree_file = NULL,
                grp_metadata_file = NULL,
                grp_model_obj = NULL,
                count_file_name = NULL,
                use_clstr_tbl = FALSE,
                min_sample_count = NULL,
                min_sample_cpm = NULL,
                test = NULL,
                test_col_name = NULL,
                skip_big_mat = FALSE)
  
  # load the params from the yaml input file
  params_yaml = yaml.load_file(yaml_file)
  
  # save the params from the yaml input file
  params$run_name = params_yaml$run_name
  params$run_id = params_yaml$run_id
  params$feature = params_yaml$feature
  params$annote_file_name = params_yaml$annote_file_name
  params$group_genes_by = params_yaml$group_genes_by
  params$gene_id_col = params_yaml$gene_id_col
  params$data_dir = params_yaml$data_dir
  params$dafe_dir = params_yaml$dafe_dir
  params$out_root_dir = params_yaml$out_root_dir
  params$ref_out_dir = params_yaml$ref_out_dir
  params$ref_meta_file = params_yaml$ref_meta_file
  params$ref_meta_cols = params_yaml$ref_meta_cols
  params$ref_include_file = params_yaml$ref_include_file
  params$ref_exclude_file = params_yaml$ref_exclude_file
  params$justify_ref_meta_by = params_yaml$justify_ref_meta_by
  params$metaG_meta_file = params_yaml$metaG_meta_file
  params$metaG_include_file = params_yaml$metaG_include_file
  params$metaG_exclude_file = params_yaml$metaG_exclude_file
  params$tree_file = params_yaml$tree_file
  params$grp_metadata_file = params_yaml$grp_metadata_file
  params$grp_model_obj = params_yaml$grp_model_obj
  params$clstr_annot_tbl_file = params_yaml$clstr_annot_tbl_file
  params$count_file_name = params_yaml$count_file_name
  params$use_clstr_tbl = params_yaml$use_clstr_tbl
  params$min_sample_cpm = params_yaml$min_sample_cpm
  params$min_sample_count = params_yaml$min_sample_count
  params$test = params_yaml$test
  params$test_col_name = params_yaml$test_col_name
  params$skip_big_mat = params_yaml$skip_big_mat
  
  
  # set the output dir
  params$output_dir = paste(params$out_root_dir, "/", params$run_name, "/", sep="")
  
  # set the reference ids
  params$ref_include_ids = read_names(params$ref_include_file)
  params$ref_exclude_ids = read_names(params$ref_exclude_file)
  params$ref_include_ids = params$ref_include_ids[!params$ref_include_ids %in% params$ref_exclude_ids]
  
  # set the metaG metadata obj ids
  params$metaG_metadata_obj = create_metaG_metadata_obj(
    metadata_file = params$metaG_meta_file, 
    include_file = params$metaG_include_file, 
    exclude_file = params$metaG_exclude_file
    )
  
  # set the env
  if (set_env == T) {
    set_env(params)
  }
  
  # do some checks on the parameters
  params$feature = check_type(params$feature)
  
  # make and return the params_obj as a class object
  class(params) = "params_obj"
  return(params)
}

params_obj.write = function(params_obj, file) {
  write(as.yaml(params_obj), file = file)
}

set_env = function(params_obj) {
  # create the output directory if needed
  dir.create(file.path(params_obj$output_dir), showWarnings=F)
}

read_names = function(file) {
  # this is used to read the names in the exclude and include file
  names = vector()
  if ( file.info(file)$size > 0 ) {
    names = read.table(file, header=F)$V1
  } else {
    names = NA
  }
  
  return(names)
}

check_type = function(feature) {
  if ( ! is_set(feature) ) {
    return("COG")
  }
  
  # this parameter must be in upper case
  feature = toupper(feature)
  
  if ( grepl("COG", feature) ) { return(feature) }
  if ( grepl("GENE", feature) ) { return(feature) }
  if ( grepl("CLSTR", feature) ) { return(feature) }
  
  warning("Unknown feature parameter.  Options are COG, GENE, and CLSTR")
  
  return(FALSE)
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
