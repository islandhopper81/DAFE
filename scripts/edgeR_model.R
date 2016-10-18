
# file for running edgeR on the gene count tables of each genome

print("Loading libraries")
library(edgeR)
library(gplots)
library(reshape2)
library(ggplot2)
library(ape) # to read in omri's tree
library(dendextend)  # for manipulating the tree
library(ggdendro) # to create the cog dendrogram
library(grid)  # to create the grid for the dendrogram and DA figure
library(yaml)  # to load the input parameters
library("phytools") # for read.newick function
#library(CompMetaG)
library(R.utils)
library(pryr)

### ARGS
# params_file=
# source_dir=

print("Getting args")
args = commandArgs(trailingOnly=T)

for(i in 1:length(args)){
  eval(parse(text=args[[i]]))
}

### check params
if ( ! exists("params_file") ) {
  stop("Must provide a param.yaml file")
} else {
  print(paste("params_file: ", params_file))
}

if ( ! exists("source_dir") ) {
  source_dir = "/netscr/yourston/compMetaG_R/R"
  print(paste("setting source_dir=", source_dir))
} else {
  print(paste("source_dir: ", source_dir))
}
###


# source necessary files
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    if(trace) cat(nm,":")           
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}
sourceDir(source_dir)

# load parameters
# params_file = "/Users/Scott/Dangl/comparative_metaG/analyses/at_isolates_id60/params.yaml"
# params_file = "/Users/Scott/Dangl/comparative_metaG/analyses/at_pop_id60_20160113/params.yaml"
# params_file = "/Users/Scott/Dangl/comparative_metaG/analyses/caulobacteraceae/params.yaml"
# params_file = "/Users/Scott/Dangl/comparative_metaG/analyses/at_pop_id95_20160205/params.yaml"
# params_file = "/Users/Scott/Dangl/comparative_metaG/analyses/batch1_id60_20160314/params.yaml"
params_obj = create_params_obj(yaml_file = params_file)

# write all the params to the log file
log_file = paste(params_obj$output_dir, "log.yaml", sep="")
params_obj.write(params_obj, file = log_file)

if ( params_obj$feature != "GENE" ) {
  ### read in the COG annotation table
  grp_metadata_tbl = read.table(params_obj$grp_metadata_file, 
                             header=T,row.names=1, sep="\t", quote="", 
                             comment.char="")
  
  
  ### COG ANALYSIS ###
  print("Building the grp_model_obj")
  
  # create the cog model object for storing info about the model inputs and outputs
  grp_model_obj = create_grp_model_obj(grp_metadata_tbl = grp_metadata_tbl,
                                       params_obj = params_obj)
  
  # build the model using edgeR
  # I removed the EC tests and focus only on RZ vs BK
  # note that when params_obj$skip_big_mat is set to TRUE nothing useful will be in
  # the grp_model_obj variable.  The skip_big_mat variable is my way of controling 
  # when I want to make this matrix.  Sometimes it can be so large the program quickly
  # runs out of memeory.
  grp_model_obj = run_edgeR_grps(grp_model_obj = grp_model_obj,
                                 params_obj = params_obj)
  
  # clean -- groups that are always undetectable
  # NOTE this step must come before calling classify_non_detect_grps
  # MOVED TO THE run_edgeR_grps FUNCTION
  #grp_model_obj = clean_grps(grp_model_obj)
  
  # specify distinction between absent and non-detectable groups
  # MOVED TO THE run_edgeR_grps FUNCTION
  #grp_model_obj = classify_non_detect_grps(grp_model_obj)
  
  # reset working dir to root dir
  setwd(params_obj$out_root_dir)
  
  # save the COG matricies if params_obj$skip_big_mat is set to FALSE
  if ( params_obj$skip_big_mat == FALSE ) {
    save(grp_model_obj, file=paste(params_obj$output_dir, "grp_model_obj.RData", sep=""))
  }
} else if ( params_obj$feature == "GENE" ) {
  # run the edgeR tests on the genes
  run_edgeR_genes(params_obj)
}  
