# This function saves a tree_matrix_plot_params_obj holding 
# the tree_matrix_plot params

# REQUIRED: grp_names, output_dir, saved_tree_matrix_obj, save_grp_model_obj, source_dir
# OPTIONAL: file_name, save, 
#           gtitle="", legend.limits = c(-3,-2,-1,0,1),
#           legend.values = c("lightgrey", "ivory", "blue", "yellow", "red"),
#           legend.name = paste("COG", "Abundance", sep="\n"),
#           legend.labels = c("Absent", "Undetectable", "RZ < BK", "RZ = BK", "RZ > BK"),
#           p1_x_vp=0.085, p2_x_vp=0.2, p3_x_vp=0.62,
#           p1_y_vp=0.51, p2_y_vp=0.5, p3_y_vp=0.515,
#           p1_mar=c(0,0,0,.1), p2_mar=c(1,0,0,-.3), 
#           p3_mar=c(0,0,0,-.1),
#           unit="inches",
#           plot_height=10, p1_height=9.75, p2_height=10, p3_height=8.9,
#           make_legends=FALSE, res=100,
#           grp_annot_tbl, meta.by="ID",
#           color_family=FALSE, color_phylum=FALSE,
#           add_leaf_node_annot=FALSE,
#           order_by_grp_groups=FALSE)

create_tree_matrix_plot_params_obj = function(yaml_file, set_env=T) {
  # this sets the defaults
  # any set to NULL are required or set in within the function
  obj = list(grp_names = NULL,
              output_dir = NULL,
              output_file = NULL,
              saved_tree_matrix_obj = NULL,
              tree_matrix_obj = NULL,
              saved_grp_model_obj = NULL,
              grp_model_obj = NULL,
              grp_annot_tbl = NULL,
              gtitle = "", 
              legend.limits = c(-4,-3,-2,-1,0,1),
              legend.values = c("darkgrey", "lightgrey", "ivory", "blue", "yellow", "red"),
              legend.name = paste("grp", "Abundance", sep="\n"),
              legend.labels = c("Absent", "Undetectable", "Low Abundance", "RZ < BK", "RZ = BK", "RZ > BK"),
              p1_x_vp=0.085, 
              p2_x_vp=0.2, 
              p3_x_vp=0.62,
              p1_y_vp=0.51, 
              p2_y_vp=0.5, 
              p3_y_vp=0.515,
              p1_mar=c(0,0,0,.1), 
              p2_mar=c(1,0,0,-.3), 
              p3_mar=c(0,0,0,-.1),
              plot_height=10, 
              p1_height=9.75, 
              p2_height=10, 
              p3_height=8.9,
              unit="inches",
              make_legends=FALSE, 
              res=100, 
              meta.by="ID",
              color_family=FALSE, 
              color_phylum=FALSE,
              add_leaf_node_annot=FALSE,
              order_by_grp_groups=FALSE
  )
  
  # load the params from the yaml input file
  params_yaml = yaml.load_file(yaml_file)
  
  # save the params from the yaml input file
  for (p in names(params_yaml)) {
    # skips things that are not supposed to be in this object
    if ( p %in% names(obj) ) {
      #print(paste("p:", p))
      obj[[p]] = params_yaml[[p]]
    }
  }
  
  # set the tree_matrix_obj based on the saved_tree_matrix_obj value
  load(obj$saved_tree_matrix_obj)
  if ( is.na(tree_matrix_obj) ) {
    stop("Variable in saved_tree_matrix_obj must be named tree_matrix_obj")
  }
  obj$tree_matrix_obj = tree_matrix_obj
  
  # set the grp_model_obj based on the saved_grp_model_obj value
  obj$grp_model_obj = set_grp_model_obj(obj)
  
  # set the grp_annot_tbl based on what is in the grp_model_obj
  obj$grp_annot_tbl = obj$grp_model_obj$grp_annot_tbl
  
  # check to make sure the output_dir is given.
  # it is required
  if ( is.null(obj$output_dir) ) {
    stop("Must provide an output_dir")
  }
  
  # set the output file 
  obj$output_file = set_output_file_name(obj)
  
  # make and return the tree_matrix_plot_params_obj as a class object
  class(obj) = "tree_matrix_plot_params_obj"
  return(obj)
}

set_grp_model_obj = function(obj) {
  # try setting the saved_grp_model_obj if needed
  if ( ! is_set(obj$saved_grp_model_obj) ) {
    # this means this value was not passed in with the params file
    # we might still be able to find it based on other values in the params file
    file = find_grp_model_obj(obj)
    if ( file != FALSE ) {
      print(paste("Setting grp_model_obj to:", file))
      obj$saved_grp_model_obj = file
    }
  }
  
  # try loading the saved_grp_model_obj
  load(obj$saved_grp_model_obj)
  if ( is.na(grp_model_obj) ) {
    stop("Variable in saved_grp_model_obj must be named grp_model_obj")
  }
  
  return(grp_model_obj)
}

find_grp_model_obj = function(obj) {
  print("Finding grp_model_obj")
  found_file = FALSE
  if ( is_set(obj$run_name) 
       & is_set(obj$out_root_dir) ) {
    dir = paste(obj$out_root_dir, obj$run_name, sep="/")
    file = paste(dir, "/grp_model_obj.RData", sep="")
    print(paste("Looking for grp_model_obj:", file))
    if ( file.exists(file) ) {
      print(paste("Found:", file))
      found_file = file
    } else {
      print(Paste("Not Found:", file))
    }
  } else {
    # The info needed to look for the grp_model_obj includes
    # run_name and out_root_dir
    print("Not enough info to find grp_model_obj")
  }
  
  return(found_file)
}

set_output_file_name = function(obj) {
  if ( is.null(obj[["output_file"]]) ) {
    if ( length(obj$grp_names) > 1 ) {
      output_file = paste(obj$grp_names[1], "_etc_grp_da_figure_with_tree.png", sep="")
    } else {
      print("here")
      output_file = paste(obj$grp_names, "_grp_da_figure_with_tree.png", sep="")
    }
    print(paste("setting output_file=", output_file))
  } else {
    output_file = obj$output_file
    print(paste("output_file: ", output_file))
  }
  
  return(output_file)
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
