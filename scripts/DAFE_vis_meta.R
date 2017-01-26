# this script make two types of plots:
# 1. the metadata plot
# 2. the legends

# NOTE: the out_dir parameter is just the prefix to the file
# names.  So if you want you can add a file prefix.  For
# example, -d /Users/Me/id_60_ would create all the output 
# files with that prefix.  However, if you want to print
# to a dir you MUST ensure that the -d parameter ends in 
# a "/" character.

library(ggplot2)
library(Hmisc) # for capitalize
library(grid) # for making legends
library(gridExtra) # for making legends
library(getopt)

### Parameter variables
# The params matrix
# Each row is a parameter with 4 columns: 
#   1) long value, 
#   2) short value, 
#   3) argument type where 0 = no argument, 1 = required , 2 = optional
#   4) data type (logical, integer, double, complex, character)
# When you call this script the args are named like --verbose or -v
params = matrix(c(
  "metadata_file", "m", 1, "character",
  "order_file", "r", 1, "character",
  "Rutils", "R", 1, "character",
  "out_dir", "d", 1, "character"
), byrow=TRUE, ncol=4)
opt = getopt(params)

source(opt$Rutils)



###########
# GLOBALS #
###########
MAX_KELLY = 20
MAX_PRIMARY = 26
SEED = 10

#########################
# 1. Make Metadata Plot #
#########################
'%!in%' <- function(x,y)!('%in%'(x,y))

load_meta = function(meta_file, meta_cols) {
  #NOTE: the meta_file should have the following columns
  # ID, Soil, Name, contig_count, group1, color1
  # more groups can possibly be added by adding columns "group2, color2"
  
  meta = read.table(meta_file, header=T, sep="\t", comment.char="", quote="")
  
  # substitute the " " characters in the Name columns with "_"
  meta$Name = gsub(" ", "_", meta$Name)
  
  # collapse any families that have too many levels
  # too many levels means there are more levels than colors
  meta = collapse_low_abund_fam(meta)
  
  # assign all the color values
  meta = assign_meta_colors(meta, meta_cols)
  
  return(meta)
}
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
      indecies = which(meta$Family == family_levels[f_i] &
                         meta$Phylum == phylum_levels[p_i])
      meta[indecies, "family_color"] = family_colors[f_i]
    }
  }
  
  return(meta)
}

collapse_low_abund_fam = function(meta, min=20) {
  phylum_levels = levels(meta$Phylum)
  
  for (p_i in 1:length(phylum_levels)) {
    #print(paste("phylum: ", phylum_levels[p_i], sep=""))
    
    indecies = which(meta$Phylum == phylum_levels[p_i])
    tmp_meta = meta[indecies,]
    family_levels = unique(tmp_meta$Family)
    
    # return if there is no need to collapse
    # ie there are fewer than min families
    if ( length(family_levels) < min ) {
      next
    }
    
    # keep the highest min-1 families
    # convert the other (ie low) families to "other"
    order = sort(table(tmp_meta[,"Family"]))
    to_keep = names(tail(order, n=(min-1)))
    to_collapse = names(order[1:(length(order) - length(to_keep))])
    
    indecies = which(meta$Phylum == "Proteobacteria" & meta$Family %in% to_collapse)
    levels(meta$Family) = c(levels(meta$Family), "Other")
    meta[indecies, "Family"] = "Other"
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

get_meta_plot = function(meta, attr_to_save = c("phylum", "family", "source"), 
                         mar=c(0,0,0,0), unit="cm") {
  # start with the meta object
  
  # reverse the meta data frame
  #meta = meta[rev(row.names(meta)),]
  
  # know what attributes I want to save.  Each of these should have a color column
  #attr_to_save = c("phylum", "family", "source")
  attr_color_to_save = paste(attr_to_save, "_color", sep="")
  
  # create the melted matrix used to plot in ggplot
  ID = rep(meta$ID, times = length(attr_to_save))  # this will end up being the y-axis
  attr = rep(attr_to_save, each = nrow(meta))  # this will be the y-axis
  color = c()  # this will be the color (set in the next few lines)
  for( attr_color in attr_color_to_save ) {
    color = c(color, meta[,attr_color])
  }
  
  if ( all(color == "NA") ) {
    print("All meta data values are NA")
    return(FALSE)
  }
  
  
  meta_melt = as.data.frame(cbind(ID, attr, color))
  meta_melt$color = as.factor(meta_melt$color)
  
  # if all the values in the meta_melt matrix are NA then return FALSE
  #print("meta_melt: ")
  #print(meta_melt)
  if ( all(is.na(meta_melt)) ) {
    print("All the meta values are false")
    return(FALSE)
  }
  
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
          axis.text.x = element_text(size=16, angle=90, hjust=1, vjust=0.4),
          plot.margin = unit(mar, unit),
          panel.background = element_rect(color="white"))
  #meta_plot
  
  return(meta_plot)
}


# read in the metadata
meta_cols = c("Source", "Fraction",  "Label")
print("Loading the metadata")
meta = load_meta(meta_file = opt$metadata_file, meta_cols = meta_cols)

# order the metadata
ordered_ids = read.table(opt$order_file)
row.names(meta) = meta$ID
meta = meta[as.character(ordered_ids$V1),]

# make metadata plot
print("Making the metadata plot")
p2 = get_meta_plot(meta, c("phylum", "family", meta_cols))

# save the plot
png(paste(opt$out_dir, "metadata.png", sep=""))
plot(p2)
dev.off()



#######################
# 2. Make the legends #
#######################
# note: meta_cols and meta generated in step 1 is used below.

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
  ### Make the legends for the tree showing the phyla and families
  phyla = unique(meta$Phylum)
  
  for (p in phyla) {
    #print(paste("Phylum: ", p, sep=""))
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
  # convert each column to a factor and conserve the level order as observed in the column
  tmp_meta = lapply(tmp_meta, function(x) factor(x, levels=unique(x)))
  
  # convert to a data frame
  tmp_meta = as.data.frame(tmp_meta)
  
  #print("tmp_meta:")
  #print(tmp_meta)
  #print(paste("title: ", title))
  #print(paste("title_col: ", title_col))
  
  #print("breaks:")
  #print(levels(tmp_meta[,2]))
  #print("labels:")
  #print(rev(levels(tmp_meta[,1])))
  
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
    #print(vals)
    lmts = as.character(unique(meta_df[,col]))
    #print(lmts)
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

print("Making the metadata legends")
make_meta_legends(meta, opt$out_dir, meta_cols)
print("Making the phyla legends")
make_phyla_legends(meta, opt$out_dir)





