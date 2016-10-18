# this script make the tree/matrix plot(s)

library(gplots)
library(reshape2)
library(ggplot2)
library(ape) # to read in omri's tree
library(dendextend)  # for manipulating the tree
library(ggdendro) # to create the grp dendrogram
library(grid)  # to create the grid for the dendrogram and DA figure
library(yaml)  # to load the input parameters
library("phytools") # for read.newick function
#library(CompMetaG)
library(R.utils)

# the user must provide a saved_grp_model_obj or a count_tbl
# if a count_tbl is provided it will be used to make a 
# grp_model_obj

# parameters
# params_file="..."
# saved_grp_model_obj="..."
# count_tbl_file="..."
# source_dir="..."
# tree_file="..."  (must be newick)
# tree_obj="..."  (must be a saved R obj named tree)
# out_file_name="..."
# legend_name="..."
# legend_values="..."


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

if ( ! exists("saved_grp_model_obj") ) {
  if ( ! exists("count_tbl_file") ) {
    stop("Must provide a saved_grp_model_obj file or count_tbl_file")
  }
}
if ( exists("saved_grp_model_obj") ) {
  print(paste("saved_grp_model_obj: ", saved_grp_model_obj))
}
if ( exists("count_tbl_file") ) {
  print(paste("count_tbl_file: ", count_tbl_file))
}

if ( ! exists("tree_file") ) {
  print("Using tree specified in the params file")
} else {
  print(paste("tree_file: ", tree_file))
}

if ( ! exists("tree_obj") ) {
  tree_obj = NULL
} else {
  print(paste("tree_obj: ", tree_obj))
}

if ( ! exists("out_file_name") ) {
  out_file_name = "RZvsBK_cog_da_figure_with_tree.png"
  print(paste("out_file_name set to ", out_file_name))
} else {
  print(paste("out_file_name: ", out_file_name))
}

if ( exists("legend_name") ) {
  print(paste("legend_name: ", legend_name))
} else {
  legend_name = NULL
}
if ( exists("legend_values") ) {
  print(paste("legend_values:", legend_values, collapse = " "))
} else {
  legend_values = NULL
}

if ( exists("grp_names") ) {
  print(paste("grp_names: ", grp_names, collapse = ", "))
} else {
  grp_names = NULL
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

# load the params object from the given file
params_obj = create_params_obj(yaml_file = params_file)

# load the grp_model_obj
if ( exists("count_tbl_file") ) {
  grp_model_obj = load_grp_model_obj(count_tbl_file, params_obj)
} else {
  load(saved_grp_model_obj)
}

# set a new tree file as necessary
if ( exists("tree_file") ) {
  print("Setting new tree_file in params_obj")
  params_obj$tree_file = tree_file
}

# Print and save the figures
tree_matrix_obj = make_tree_matrix_obj(tree_file = params_obj$tree_file,
                                       meta_file = params_obj$ref_meta_file,
                                       meta.by = params_obj$justify_ref_meta_by,
                                       meta_cols = params_obj$ref_meta_cols,
                                       matrix = grp_model_obj$mat,
                                       tree_obj = tree_obj)
save(tree_matrix_obj, file=paste(params_obj$output_dir, "tree_matrix_obj.RData", sep=""))

# NOTE: this step has parameters that are manually adjusted to make the plot
# look good.  You just have to play with them.  To change the height of the
# plots use the heights paramateers.  To change the positions of the plot 
# use the *x_vp and *y_vp parameters.  P1 is the plot farthest to the left
# and p2 is next and so on.
get_tree_matrix_plot(tree_matrix_obj = tree_matrix_obj, 
                     output_dir = params_obj$output_dir, 
                     output_file=out_file_name, 
                     p3_height=8.87, p3_y_vp=.517, 
                     order_by_cog_groups=FALSE, 
                     grp_metadata_tbl = grp_model_obj$grp_metadata_tbl, 
                     make_legends = TRUE,
                     meta.by = params_obj$justify_ref_meta_by,
                     res = 100,
                     legend.name = legend_name,
                     legend.values = legend_values,
                     grp_names = grp_names
)

# p3_height = 8.87
# p3_y_vp = .517

# params for single grp graphing
# p3_heigth = 10
# p3_y_vp = .45

