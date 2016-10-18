# this file create an object for counting the number of reads that map
# to each of the genomes


build_abund_matrix = function(params_obj, tree_matrix_obj, 
                              output_dir, save=T,
                              output_file="genome_abund_figure_with_tree.png", 
                              gtitle="",
                              p1_x_vp=0.085, p2_x_vp=0.2, p3_x_vp=0.62,
                              p1_y_vp=0.51, p2_y_vp=0.5, p3_y_vp=0.515,
                              p1_mar=c(0,0,0,.1), p2_mar=c(1,0,0,-.3), 
                              p3_mar=c(0,0,0,-.1),
                              unit="inches",
                              plot_height=10, p1_height=9.75, p2_height=10, p3_height=8.9, 
                              color_phylum=FALSE, color_family=FALSE,
                              add_leaf_node_annot=FALSE,
                              res=100) {
  
  if (tree_matrix_obj$is_justified != TRUE ) {
    stop("ERROR: tree_matrix_obj must have been justified")
  }
  
  #col_names = c("soil_yng", "soil_old", "col_yng_rz", "col_old_rz")
  col_names = row.names(params_obj$metaG_metadata_obj)
  
  seq_totals = matrix(NA, ncol = 4, nrow = length(params_obj$ref_include_ids),
                      dimnames=list(c(), col_names))
  i = 1
  for (ref in params_obj$ref_include_ids) {
    print(ref)
    
    # set and create the output directory
    ref_out_dir = paste(params_obj$out_root_dir, params_obj$run_name, "ref_out", ref, sep="/")
    dir.create(path = ref_out_dir, recursive = T, showWarnings = F)
    setwd(ref_out_dir)
    
    # read in the file
    # the get_count_tbl function is in DA_test.R
    tbl = get_count_tbl(ref, params_obj)
    
    # remove the features that start with "__"
    to_remove = c("__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned", "__alignment_not_unique")
    tbl = tbl[!(row.names(tbl) %in% to_remove),]
    
    # add to the matrix
    seq_totals[i,] = apply(tbl, 2, sum)
    i = i+1
  }
  
  # set the row names
  row.names(seq_totals) = params_obj$ref_include_ids
  
  # reorder based on the tree
  g_order = tree_matrix_obj$meta$ID
  seq_totals = seq_totals[as.character(g_order),]
  
  # melt
  seq_totals_m = melt(seq_totals)
  
  # make the matrix plot
  p3 = ggplot(seq_totals_m, aes(x=factor(Var2), y=factor(Var1), fill=value)) + 
    geom_tile() + 
    scale_y_discrete(limits=g_order) +
    scale_x_discrete("Sample") +
    scale_fill_gradient("Abundance") +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(angle=90, hjust=1),
          axis.title.x = element_text(),
          plot.background = element_rect(color="white"),
          panel.background = element_rect(color="white"),
          plot.margin = unit(p3_mar, unit))
  
  # get the tree plot
  p1 = get_dendextend_ggtree(tree_matrix_obj$tree, tree_matrix_obj$meta, 
                             mar=p1_mar, color_phylum = color_phylum, 
                             color_family = color_family, 
                             add_leaf_node_annot = add_leaf_node_annot)
  
  # get the p2 plot (meta data plot)
  p2 = get_meta_plot(tree_matrix_obj = tree_matrix_obj,
                     mar = p2_mar)
  
  # set some dimentions for the plot
  p1_width = 2
  p2_width = 1
  if ( length(names) == 1 ) {
    p3_width = 4
  } else if ( length(names) > 1 & length(names) <= 4 ) {
    p3_width = 1 * length(names)
  } else {
    p3_width = 8
  }
  
  # create an output file if specified
  if ( save == T) {
    png(paste(output_dir, output_file, sep=""), 
        units="in", width = (p1_width + p2_width + p3_width), 
        height = plot_height, res=res)
  }
  
  grid.newpage()
  print(p1, vp=viewport(x=p1_x_vp, y=p1_y_vp, h=unit(p1_height, unit), w=unit(p1_width, unit)))
  print(p2, vp=viewport(x=p2_x_vp, y=p2_y_vp, h=unit(p2_height, unit), w=unit(p2_width, unit)))
  print(p3, vp=viewport(x=p3_x_vp, y=p3_y_vp, h=unit(p3_height, unit), w=unit(p3_width, unit)))
  if ( save == T ) { dev.off() }
}