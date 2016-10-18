# If you want to use the standard functionality of this object
# first make a matrix_tree_obj by calling the make_tree_matrix_obj
# function.


####
# set some global variables use here
MAX_KELLY = 20
MAX_PRIMARY = 26
SEED = 10




make_tree_matrix_clstr_obj = function(tree_file, meta_file, meta_cols, matrix, meta.by) {
  # the meta_cols are the columns that I want to include metadata for in the plot
  
  tree_matrix_obj = list(tree = NULL, meta = NULL, meta_cols = NULL, matrix = NULL)
  
  # read in the tree to a dendrogram object
  print("load_tree")
  tree = load_tree(tree_file)
  
  # load the meta data -- see function for important notes about this file format
  print("load_meta")
  meta = load_meta(meta_file, meta_cols)
  
  tree_matrix_obj$tree = tree
  tree_matrix_obj$meta = meta
  tree_matrix_obj$meta_cols = meta_cols
  tree_matrix_obj$matrix = matrix
  tree_matrix_obj$is_justified = FALSE
  
  class(tree_matrix_obj) = "tree_matrix_obj"
  
  print("Justify")
  tryCatch({
    tree_matrix_obj = justify_data(tree_matrix_obj, meta.by)
  }, error = function(err) {
    print(paste("ERROR: ", err))
    traceback(2)
  });
  
  
  return(tree_matrix_obj)
}

justify_data = function(tree_matrix_obj, meta.by="Name") {
  # get the IDs corresponding to the names in the tree
  # so this is the list of IDs that are in the tree (and in the meta table)
  
  meta_i = which(colnames(tree_matrix_obj$meta) == meta.by)
  
  #overlapping_labels = tree_matrix_obj$meta$Name %in% labels(tree_matrix_obj$tree)
  overlapping_labels = tree_matrix_obj$meta[,meta_i] %in% labels(tree_matrix_obj$tree)
  tipIDs = tree_matrix_obj$meta[overlapping_labels,]$ID
  
  # intersect all three lists of IDs
  ids_to_keep = Reduce(intersect, list(tipIDs, 
                                       tree_matrix_obj$meta$ID, 
                                       colnames(tree_matrix_obj$matrix)))
  
  # keep only the IDs I want in the meta table
  tree_matrix_obj$meta = tree_matrix_obj$meta[tree_matrix_obj$meta$ID %in% ids_to_keep,]
  
  # keep only the IDs I want in the tbl
  tree_matrix_obj$matrix = subset(tree_matrix_obj$matrix, subset = T, select = ids_to_keep)
  
  # keep only the IDs I want in the tree
  # those IDs now correspond with the names that are in the meta table
  to_drop = labels(tree_matrix_obj$tree)[!labels(tree_matrix_obj$tree) %in% tree_matrix_obj$meta[,meta_i]]
  
  if ( is.dendrogram(tree_matrix_obj$tree) ) {
    tree_matrix_obj$tree = prune(tree_matrix_obj$tree, leaves = to_drop)
  } else if ( is.phylo(tree_matrix_obj$tree ) ) {
    tree_matrix_obj$tree = prune(tree_matrix_obj$tree, tip=to_drop)
  } else {
    stop("ERROR: unknown tree object type")
  }
  
  # reorder the meta to reflect the tree
  tree_matrix_obj$meta = tree_matrix_obj$meta[match(labels(tree_matrix_obj$tree), tree_matrix_obj$meta[,meta_i]), ]
  
  # reorder the tbl to reflect the tree
  tree_matrix_obj$matrix = subset(tree_matrix_obj$matrix, 
                                  select = match(tree_matrix_obj$meta$ID, 
                                                 colnames(tree_matrix_obj$matrix)))
  
  # set the boolean value to indicate that it has been justified
  tree_matrix_obj$is_justified = TRUE
  
  return(tree_matrix_obj)
}


# Function for creating the clstr tree figure
get_tree_matrix_clstr_plot = function(tree_matrix_obj, output_dir, save=T, names=NULL, 
                                output_file="clstr_da_figure_with_tree.png", 
                                gtitle="", legend.limits = c(-3,-2,-1,0,1),
                                legend.values = c("lightgrey", "ivory", "blue", "yellow", "red"),
                                legend.name = paste("Cluster", "Abundance", sep="\n"),
                                legend.labels = c("Absent", "Undetectable", "RZ < BK", "RZ = BK", "RZ > BK"),
                                p1_x_vp=0.085, p2_x_vp=0.2, p3_x_vp=0.62,
                                p1_y_vp=0.51, p2_y_vp=0.5, p3_y_vp=0.515,
                                p1_mar=c(0,0,0,.1), p2_mar=c(1,0,0,-.3), 
                                p3_mar=c(0,0,0,-.1),
                                unit="inches",
                                plot_height=10, p1_height=9.75, p2_height=10, p3_height=8.9,
                                make_legends=FALSE, res=100,
                                clstr_annot_tbl, meta.by="ID",
                                color_family=FALSE, color_phylum=FALSE,
                                add_leaf_node_annot=FALSE,
                                order_by_clstr_groups=FALSE) {
  
  # handle all input cases (ie full matrix, single subset, multiple subset)
  print("get_clstr_df()")
  clstr_df = get_clstr_df(clstr_mat = tree_matrix_obj$matrix, names = names)
  
  # get the tree plot
  print("get_dendextend_ggtree")
  tryCatch({
    p1 = get_dendextend_ggtree(tree = tree_matrix_obj$tree, 
                               meta = tree_matrix_obj$meta, 
                               mar=p1_mar, unit="cm",
                               color_family=FALSE, color_phylum=FALSE,
                               add_leaf_node_annot=FALSE, meta.by = meta.by)
  }, error = function(err) {
    print(err)
    traceback(2)
  })
  
  ### create the clstr dendrogram -- mostly to get the order of the clstrs
  if ( nrow(clstr_df) > 1 ) {
    #clstr.dendro = as.dendrogram(hclust(dist(clstr_df)))
    #clstr.ord = order.dendrogram(clstr.dendro)
    clstr.ord.names = get_clstr_name_ord(clstr_df, 
                                     clstr_annot_tbl = clstr_annot_tbl,
                                     meta.by = meta.by,
                                     order_by_clstr_groups)
    #clstr.ord.names = rownames(clstr_df[clstr.ord,])
  } else {
    clstr.ord.names = rownames(clstr_df)
  }
  
  # get the genome order
  g_order = tree_matrix_obj$meta$ID
  
  # get the p2 plot (meta data plot)
  print("get_meta_plot()")
  p2 = get_meta_plot(tree_matrix_obj = tree_matrix_obj,
                     mar = p2_mar)
  
  # melt the data frame
  clstr.mat.melt = melt(as.matrix(clstr_df))
  
  # create the DA table
  print("get DA matrix")
  p3 = ggplot(clstr.mat.melt, aes(x=factor(Var1), y=factor(Var2), fill=factor(value))) + 
    geom_tile() + 
    scale_x_discrete(limits=clstr.ord.names) + 
    scale_y_discrete(limits=g_order) +
    scale_fill_manual(limits=legend.limits,
                      values=legend.values,
                      name=legend.name,
                      labels=legend.labels) + 
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(angle=90, hjust=1, size=12),
          axis.title = element_blank(),
          legend.title = element_text(size=22),
          legend.text = element_text(size=20),
          plot.margin = unit(p3_mar, unit),
          panel.background = element_rect(fill="white"))
  
  # remove the x axis information if graphing all clstrs
  if (length(names) == 0) {
    p3 = p3 + theme(axis.text.x = element_blank(),
                    axis.ticks.x = element_blank())
  }
  
  # set some dimentions for the plot
  p1_width = 2
  p2_width = 1
  if ( length(names) == 1 ) {
    p3_width = 4
  } else if ( length(names) > 1 & length(names) <= 4 ) {
    p3_width = 2 * length(names)
  } else {
    p3_width = 8
  }
  
  # set up the output file to save the figure to
  if ( save == T) {
    png(paste(output_dir, output_file, sep=""), 
        units="in", width = (p1_width + p2_width + p3_width), 
        height = plot_height, res=res)
  }
  
  # add the plots to a grid
  grid.newpage()
  if ( exists("p1") ) {
    print(p1, vp=viewport(x=p1_x_vp, y=p1_y_vp, h=unit(p1_height, unit), w=unit(p1_width, unit)))
  }
  print(p2, vp=viewport(x=p2_x_vp, y=p2_y_vp, h=unit(p2_height, unit), w=unit(p2_width, unit)))
  print(p3, vp=viewport(x=p3_x_vp, y=p3_y_vp, h=unit(p3_height, unit), w=unit(p3_width, unit)))
  if ( save == T ) { dev.off() }
  
  # print the individual plots
  if ( save == T ) {
    file_base = sub("^([^.]*).*", "\\1", output_file)
    
    # plot1 -- tree
    png(paste(output_dir, file_base, "_p1.png", sep = ""),
        units="in", width=p1_width, height=plot_height, res=res)
    if ( exists("p1") ) {
      print(p1, vp=viewport())
    }
    
    # plot 2 -- meta
    png(paste(output_dir, file_base, "_p2.png", sep = ""),
        units="in", width=p2_width, height=plot_height, res=res)
    print(p2, vp=viewport())
    
    # plot 3 -- matrix heatmap
    png(paste(output_dir, file_base, "_p3.png", sep = ""),
        units="in", width=p3_width, height=plot_height, res=res)
    print(p3, vp=viewport())
    
    dev.off()
  }
  
  # make the legends and save them in their seperate files
  # they will have to be added to the figure manually using PowerPoint or something similar
  if ( make_legends == TRUE ) { 
    print("Make phyla legends")
    make_phyla_legends(tree_matrix_obj$meta, output_dir)
    
    print("Make meta legends")
    make_meta_legends(tree_matrix_obj$meta, output_dir, tree_matrix_obj$meta_cols)
  }
}




### HELPER FUNCTIONS ###
load_tree = function(tree_file) {
  # NOTE: the tree should be in Newick format
  
  # to make this tree I:
{
  # 1. opened Ormi's tree in ITOL.
  # 2. exported the tree in Newick format
  # 3. opened the tree in geneious
  # 4. exported the tree in Newick format
  # 5. substitute the ' characters out of the file
}

tree = read.newick(tree_file)

if ( ! is.ultrametric(tree) ) {
  tree = tryCatch({
    tree = tree %>% collapse.singles %>% chronos
  }, error = function(err) {
    print(paste("ERROR: ", err, sep=""))
    traceback(2)
    tree = tree %>% collapse.singles %>% compute.brtime(force.positive = TRUE)
  });
}

if ( ! is.binary.tree(tree) ) {
  tree = multi2di(tree)
}

tree = tree %>% as.dendrogram

return(tree)
}

load_meta = function(meta_file, meta_cols) {
  #NOTE: the meta_file should have the following columns
  # ID, Soil, Name, contig_count, group1, color1
  # more groups can possibly be added by adding columns "group2, color2"
  
  meta = read.table(meta_file, header=T, sep="\t", comment.char="", quote="")
  
  # substitute the " " characters in the Name columns with "_"
  meta$Name = gsub(" ", "_", meta$Name)
  
  # assign all the color values
  meta = assign_meta_colors(meta, meta_cols)
  
  return(meta)
}

### sets all the color values
# - phylogenty (phylum, family)
# - source values
assign_meta_colors = function(meta, meta_cols) {
  i = 1
  for( col in meta_cols ) {
    meta = assign_col_colors(meta=meta, col_name=col, seed=i)
    i = i + 1
  }
  meta = assign_phylogeny_colors(meta)
  
  # transform false values to NA
  
  return(meta)
}

### sets a color value for each phylum and family
assign_phylogeny_colors = function(meta) {
  # set class and genus colors
  meta$phylum_color = rep("NA", times=nrow(meta))
  meta$family_color = rep("NA", times=nrow(meta))
  phylum_levels = levels(meta$Phylum)
  for (p_i in 1:length(phylum_levels)) {
    colors = get_kelly_colors()  # reset the colors
    
    #paste("phylum: ", phylum_levels[p_i], sep="") %>% print
    indecies = which(meta$Phylum == phylum_levels[p_i])
    
    # set the phylum color
    phylum_color = colors[p_i]
    meta[indecies, "phylum_color"] = phylum_color
    
    # temporarily remove from colors so a phylum can't be the same color as a family
    # I don't need this step anymore because I'm using the primary
    # colors for the family colors
    #colors = colors[which(colors!=phylum_color)]
    family_colors = get_primary_colors()
    
    # set the family colors
    tmp_meta = meta[indecies,]
    family_levels = unique(tmp_meta$Family)
    for (f_i in 1:length(family_levels)) {
      #paste("family: ", family_levels[f_i], sep="") %>% print
      indecies = which(meta$Family == family_levels[f_i])
      meta[indecies, "family_color"] = family_colors[f_i]
    }
  }
  
  return(meta)
}

### For a given column in the metadata tbl assign unique colors to each value
# A column is added to the meta tbl with the assigned colors for the values
# in the specified column.  The name of the new color column is col_name +
# "_color".  For example, if the col_name == source then the new column name
# is source_color.
assign_col_colors = function(meta, col_name, seed=SEED) {
  new_col_name = paste(col_name, "_color", sep="")
  uniq_vals = unique(meta[,col_name])
  
  if ( length(uniq_vals) < MAX_KELLY ) {
    colors = get_kelly_colors(n=length(uniq_vals), rand=T, seed=seed)
  } else if ( length(uniq_vals) < MAX_PRIMARY ) {
    colors = get_primary_colors(length(uniq_vals))
  } else {
    stop(paste("Not enough colors.  Need at least ", length(uniq_vals)))
  }
  
  meta[,new_col_name] = rep("NA", times=nrow(meta))
  for( s_i in 1:length(uniq_vals) ) {
    indecies = which(meta[,col_name] == uniq_vals[s_i])
    meta[,new_col_name][indecies] = colors[s_i]
  }
  
  return(meta)
}


### converts the clstr_mat into a clstr_df
# subsets the matrix if "names" is specified
# "names" is a vector of columns (ie clstrs) to keep
get_clstr_df = function(clstr_mat=NULL, names=NULL) {
  if ( length(names) == 0 ) {
    clstr_df = as.data.frame(clstr_mat)
  } else if ( length(names) == 1 ) {
    clstr_df = t(as.data.frame(clstr_mat[names,]))
    row.names(clstr_df) = names
  } else if ( length(names) > 1 ) {
    clstr_df = as.data.frame(clstr_mat[names,])
  } else {
    print("ERROR: wrong type for names value")
  }
  
  return(clstr_df)
}



### create a tree using dendextend with ggplot functionality
get_dendextend_ggtree = function(tree, meta, mar=c(0,0,0,0), unit="cm",
                                 color_family=FALSE, color_phylum=FALSE,
                                 add_leaf_node_annot=FALSE,
                                 meta.by="ID") {
  #meta = assign_meta_colors(tree_matrix_obj$meta) # wont work anymore
  #tree = tree_matrix_obj$tree
  
  meta_i = which(colnames(meta) == meta.by)
  
  # convert the tree from phylo to dendrogram if needed
  if ( ! is.dendrogram(tree) ) {
    tree = as.dendrogram(tree)
  }
  
  # colors the upper branches of the tree by either phylum or class
  if ( color_phylum == TRUE ) {
    colors = get_kelly_colors()
    phylum_levels = levels(meta$Phylum)
    for (p_i in 1:length(phylum_levels)) {
      print(phylum_levels[p_i])
      tree = set(tree, "by_labels_branches_col", 
                 value = meta[which(meta$Phylum == phylum_levels[p_i]), meta_i],
                 TF_values = c(colors[p_i], Inf))
    }
  }
  
  #class_levels = levels(meta$Class)
  #for (c_i in 1:length(class_levels)) {
  #  print(class_levels[c_i])
  #  tree = set(tree, "by_labels_branches_col", 
  #             value = meta[which(meta$Class == class_levels[c_i]), "Name"],
  #             TF_values = c(colors[c_i], Inf))
  #}
  
  # colors the leaf branches of the tree by family
  if ( color_family == TRUE ) {
    tree = assign_values_to_leaves_edgePar(tree, value=meta$family_color, edgePar = "col")
  }
  
  # color only the leaf branches of the tree by omri's color scheme (color1)
  #tree = assign_values_to_leaves_edgePar(object=tree, value=meta$color1, edgePar = "col")
  
  # set branch width
  tree = set(tree, "branches_lwd", 1)
  
  # set leaf color from the source column
  if ( add_leaf_node_annot == TRUE) {
    tree = tree %>% 
      set("leaves_pch", 19) %>% 
      set("leaves_cex", 2) %>% 
      set("leaves_col", meta$source_color)
  }
  
  # convert to a ggplot acceptable format
  ggdend = as.ggdend(tree)
  
  # I was getting a weird error caused by these being characters instead
  # of integer values.
  ggdend$segments$lwd = as.integer(ggdend$segments$lwd) 
  
  # build the plot
  p = ggplot(ggdend, labels=F, horiz=T)
  
  p = p + theme(plot.margin = unit(mar, unit))
  
  return(p)
}

get_meta_plot = function(tree_matrix_obj, mar=c(0,0,0,0), unit="cm") {
  # start with the meta object
  # lots of this is hard coded.  Sorry but I kind of have to assume that I 
  # know what is coming in via the metadata file.
  
  meta = tree_matrix_obj$meta
  
  # know what attributes I want to save.  Each of these should have a color column
  #attr_to_save = c("phylum", "family", "source")
  attr_to_save = c("phylum", "family", tree_matrix_obj$meta_cols)
  attr_color_to_save = paste(attr_to_save, "_color", sep="")
  
  # create the melted matrix used to plot in ggplot
  ID = rep(meta$ID, times = length(attr_to_save))  # this will end up being the y-axis
  attr = rep(attr_to_save, each = nrow(meta))  # this will be the y-axis
  color = c()  # this will be the color (set in the next few lines)
  for( attr_color in attr_color_to_save ) {
    color = c(color, tree_matrix_obj$meta[,attr_color])
  }
  
  meta_melt = as.data.frame(cbind(ID, attr, color))
  meta_melt$color = as.factor(meta_melt$color)
  
  # make the plot
  meta_plot = ggplot(meta_melt, aes(x=attr, y=ID)) +
    geom_tile(aes(fill=color)) +
    scale_fill_manual(values=levels(meta_melt$color),
                      guide=FALSE) + 
    scale_x_discrete(limits=attr_to_save,
                     labels=capitalize(attr_to_save)) + 
    scale_y_discrete(limits=meta$ID) + 
    theme(axis.title = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size=16, angle=90, hjust=.5),
          plot.margin = unit(mar, unit),
          panel.background = element_rect(color="white"))
  meta_plot
  
  return(meta_plot)
}

# gets the order to graph the clstrs along the y axis
# the default is to order using a hierarchical clustering
get_clstr_name_ord = function(clstr_df, order_by_clstr_groups=FALSE, clstr_annot_tbl, 
                            meta.by) {
  clstr.ord.names = vector() # the return value
  
  if ( order_by_clstr_groups == TRUE ) {
    print("Ordering by clstr groups")
    # clstr_df is a matrix with rows being clstrs and columns being names
    # I want to order the rows by what function they match with.  Then
    # for each function I want to order (within that function) using
    # hierarchical clustering
    
    # merge the clstr_annot_tbl and clstr_df
    merged = merge(x = clstr_df, y = clstr_annot_tbl,
                   all.x=TRUE, all.y = FALSE,
                   by.x = "row.names", by.y = "row.names")
    
    # for functions that are assigned to more than on clstr group randomly pick one
    single_func = sapply(merged$func, 
                         function(x) sample(strsplit(as.character(x), "")[[1]], 1))
    merged$single_func = single_func
    
    # order by func
    merged_ordered = merged[order(merged$single_func),]
    
    # for each unique function, f
    # subset the table for f
    # hierarchical cluster f if f has more than 1 entry
    # save the order of clstrs using the names
    for ( f in unique(merged_ordered$single_func) ) {
      subset = merged_ordered[merged_ordered$single_func == f,]
      if ( nrow(subset) <= 1 ) {
        clstr.ord.names = c(clstr.ord.names, subset$Row.names)
        next
      }
      clstr.dendro = as.dendrogram(hclust(dist(clstr_df[subset$Row.names,])))
      clstr.ord = order.dendrogram(clstr.dendro)
      clstr.ord.names = c(clstr.ord.names, subset[clstr.ord,]$Row.names)
    }
    
    # return final clstr order
    
  } else {
    # order by hierarchical clustering
    clstr.dendro = as.dendrogram(hclust(dist(clstr_df)))
    clstr.ord = order.dendrogram(clstr.dendro)
    clstr.ord.names = rownames(clstr_df[clstr.ord,])
  }
  
  return(clstr.ord.names)
}


# 20 colors
get_kelly_colors = function(n=20, rand=F, seed=10) {
  MAX = 20  # there are 20 colors here
  
  colors = c(rgb(255,179,0,max=255),
             rgb(128,62,117,max=255),
             rgb(255,104,0,max=255),
             rgb(166,189,215,max=255),
             rgb(193,0,32,max=255),
             rgb(206,162,98,max=255),
             rgb(129,112,102,max=255),
             rgb(0,125,52,max=255),
             rgb(246,118,142,max=255),
             rgb(0,83,138,max=255),
             rgb(225,122,92,max=255),
             rgb(83,55,122,max=255),
             rgb(255,142,0,max=255),
             rgb(179,40,81,max=255),
             rgb(244,200,0,max=255),
             rgb(127,24,13,max=255),
             rgb(147,170,0,max=255),
             rgb(89,51,21,max=255),
             rgb(241,58,19,max=255),
             rgb(35,44,22,max=255)
  )
  
  if ( rand == T) {
    set.seed(seed)
    to_keep = sample(1:MAX, n)
  }
  else {
    to_keep = 1:n
  }
  
  return(colors[to_keep])
}

get_primary_colors = function(n=26) {
  require(colorRamps)
  # there are 26 optional colors here
  
  MAX = 26  # the number of primary.colors
  if ( n > MAX) {
    stop("Trying to get more than 26 colors.  n must be < 26")
  }
  
  return(primary.colors()[1:n])
}
