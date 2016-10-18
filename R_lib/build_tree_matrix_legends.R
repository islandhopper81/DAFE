# note: a legend for each phylum is made an stored in a 
# file names <phyla>.png in the given output directory

# note: these functions are a little poorly named.  When I say
# "make_phyla" I mean make a legend for each file that contains
# the color information for each family

# note: the meta parameter should be the tree_matrix_obj$meta
# object after the justify_data() function has been called.
# This will ensure that only the genome that you want in the
# tree are included in the metadata that is used to generate
# the legends.

make_phyla_legends = function(meta, output_dir, res=50) {
  library(gridExtra)
  
  ### Make the legends for the tree showing the phyla and families
  phyla = unique(meta$Phylum)
  
  for (p in phyla) {
    out_file = paste(output_dir, paste(p, ".png", sep=""), sep="/")
    png(filename = out_file, res=res)
    phyla_meta = meta[meta$Phylum == p, c("Family", "family_color")]
    color = unique(meta[meta$Phylum == p, c("phylum_color")])
    phy_legend = get_phylum_legend(phyla_meta, p, title_col = color)
    grid.draw(phy_legend)
    dev.off()
  }
}

g_legend<-function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)
} 

get_phylum_legend = function(tmp_meta, title, title_col) {
  #tmp_meta = meta[,c("Family", "family_color")]
  # convert each column to a factor and conserve the level order as observed in the column
  tmp_meta = lapply(tmp_meta, function(x) factor(x, levels=unique(x)))
  
  # convert to a data frame
  tmp_meta = as.data.frame(tmp_meta)
  
  print("tmp_meta:")
  print(tmp_meta)
  print(paste("title: ", title))
  print(paste("title_col: ", title_col))
  
  print("breaks:")
  print(levels(tmp_meta[,2]))
  print("labels:")
  print(rev(levels(tmp_meta[,1])))
  
  # build the ggplot
  p = ggplot(tmp_meta, aes_string(x=colnames(tmp_meta[1]), y="..count..", fill=colnames(tmp_meta[2]))) + 
    geom_bar() + 
    scale_fill_manual(breaks=levels(tmp_meta[,2]),
                      values=rev(levels(tmp_meta[,2])), 
                      labels=rev(levels(tmp_meta[,1])),
                      guide=guide_legend(title=title)) + 
    theme(legend.title = element_text(color=title_col, size=22),
          legend.text = element_text(size=20))
  
  # get the legend
  legend = g_legend(p)
  
  return(legend)
}

make_meta_legends = function(meta, output_dir, meta_cols, title_col = "black") {
  require(gridExtra)
  require(R.utils)
  
  # convert to a data frame
  meta_df = data.frame(meta)
  
  # build the ggplot
  for( col in meta_cols ) {
    col_color = paste(col, "_color", sep="")  # name of the color column
    vals = unique(meta_df[,col_color])
    print(vals)
    lmts = as.character(unique(meta_df[,col]))
    print(lmts)
#    brks = 
    p = ggplot(meta_df, aes_string(x=col, y="..count..", fill=col)) +
      geom_bar() +
      scale_fill_manual(values=vals,
                        limits=lmts,
                        guide=guide_legend(title=capitalize(col))) +
      theme(legend.title = element_text(color=title_col, size=22),
            legend.text = element_text(size=20))
    
    # get the legend
    legend = g_legend(p)
    
    # save the legend in a file
    out_file = paste(output_dir, paste(col, ".png", sep=""), sep="/")
    png(filename = out_file)
    grid.draw(legend)
    dev.off()
  }
}





### DEPRECIATED - replaced by make_meta_legends which can make a legend for all the
# columns in the metadata plot (except for the family and phylum ones which are 
# made by make_phyla_legends).
make_source_legend = function(meta, output_dir, title="Source", title_col="black") {
  library(gridExtra)
  
  # convert to a data frame
  meta_df = data.frame(meta)
  
  # build the ggplot
  p = ggplot(meta_df, aes(x=source, y=..count.., fill=source)) + 
    geom_bar() + 
    scale_fill_manual(values=levels(as.factor(meta_df$source_color)),
                      guide=guide_legend(title=title)) +
    theme(legend.title = element_text(color=title_col, size=22),
          legend.text = element_text(size=20))
  
  # get the legend
  legend = g_legend(p)
  
  # save the legend in a file
  out_file = paste(output_dir, "source.png", sep="/")
  png(filename = out_file)
  grid.draw(legend)
  dev.off()
  
}
