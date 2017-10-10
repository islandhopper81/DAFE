# This script creates a DAFE heatmap

# There are two parameters
# 1. The matrix of DA values (rows are functions, cols are genomes)
# 2. The png output file name


library(gplots)
library(reshape2)
library(ggplot2)
library(ape) # to read in omri's tree
library(grid)  # to create the grid for the dendrogram and DA figure
library(yaml)  # to load the input parameters
library("phytools") # for read.newick function
library(R.utils)
library(getopt)
library(cluster) # for daisy funciton for clustering COGs


##############
# Parameters #
##############

# The params matrix
# Each row is a parameter with 4 columns: 
#   1) long value, 
#   2) short value, 
#   3) argument type where 0 = no argument, 1 = required , 2 = optional
#   4) data type (logical, integer, double, complex, character)
# When you call this script the args are named like --verbose or -v
params = matrix(c(
  "matrix_file", "m", 1, "character",
  "out_file_pre", "o", 1, "character",
  "show_xlabs", "x", 0, "logical",
	"png", "n", 0, "logical"
), byrow=TRUE, ncol=4)
opt = getopt(params)
matrix_file = opt$matrix_file
out_file_pre = opt$out_file_pre

print(paste("opt$show_xlabs: ", opt$show_xlabs, sep=""))
if ( ! is.null(opt$show_xlabs) ) {
	show_xlabs = T
	print("setting show_xlabs to T (TRUE)")
} else {
	show_xlabs = F
}
print(paste("show_xlabs: ", show_xlabs, sep=""))

if ( ! is.null(opt$png) ) {
	png = T
} else {
	png = F
}

#############
# Functions #
#############

# convert it back to a dataframe
get_grp_df = function(grp_mat=NULL, names=NULL) {
  if ( length(names) == 0 ) {
    grp_df = as.data.frame(grp_mat)
  } else if ( length(names) == 1 ) {
    grp_df = t(as.data.frame(grp_mat[names,]))
    row.names(grp_df) = names
  } else if ( length(names) > 1 ) {
    common_names = names %in% row.names(grp_mat)
    print(paste("common names length: ", length(common_names)))
    print("missing names:")
    print(names[!common_names])
    grp_df = as.data.frame(grp_mat[names,])
  } else {
    print("ERROR: wrong type for names value")
  }
  
  return(grp_df)
}

# gets the order to graph the groups along the x axis
# the default is to order using a hierarchical clustering
get_grp_name_ord = function(grp_df, order_by_grp_groups=FALSE, grp_metadata_tbl, 
                            meta.by) {
  
  ret = list() # the return value
  #grp.ord.names = vector() # the return value
  
  if ( order_by_grp_groups == TRUE ) {
    print("Ordering by higher groups")
    # grp_df is a matrix with rows being groups and columns being names
    # I want to order the rows by what function they match with.  Then
    # for each function I want to order (within that function) using
    # hierarchical clustering
    
    # merge the grp_metadata_tbl and grp_df
    merged = merge(x = grp_df, y = grp_metadata_tbl,
                   all.x=TRUE, all.y = FALSE,
                   by.x = "row.names", by.y = "row.names")
    
    # for functions that are assigned to more than on higher group randomly pick one
    single_func = sapply(merged$func, 
                         function(x) sample(strsplit(as.character(x), "")[[1]], 1))
    merged$single_func = single_func
    
    # order by func
    merged_ordered = merged[order(merged$single_func),]
    
    # for each unique function, f
    # subset the table for f
    # hierarchical cluster f if f has more than 1 entry
    # save the order of groups using the names
    for ( f in unique(merged_ordered$single_func) ) {
      subset = merged_ordered[merged_ordered$single_func == f,]
      if ( nrow(subset) <= 1 ) {
        grp.ord.names = c(grp.ord.names, subset$Row.names)
        next
      }
      grp.dendro = as.dendrogram(hclust(dist(grp_df[subset$Row.names,])))
      grp.ord = order.dendrogram(grp.dendro)
      grp.ord.names = c(grp.ord.names, subset[grp.ord,]$Row.names)
    }
    
    # return final group order
    
  } else {
    # order by hierarchical clustering
    grp_df = t(grp_df)
    grp.dendro = as.dendrogram(hclust(daisy(as.data.frame(apply(grp_df, 2, as.factor)))))
    grp.ord = order.dendrogram(grp.dendro)
    grp.ord.names = rownames(grp_df[grp.ord,])
    ret = list(ord = grp.ord.names,
               dendro = grp.dendro)
  }
  
  return(ret)
}

make_figure = function(mat, file_pre, show_xlabs = F, png=F) {
  
  df = get_grp_df(mat)
  
  # met the df so that I can use ggplot to print it
  grp.mat.melt = melt(as.matrix(df))
  grp.mat.melt$Var1 = as.character(grp.mat.melt$Var1)
  grp.mat.melt$Var2 = as.character(grp.mat.melt$Var2)
  
  
  # set up some variables for plotting
  legend.values = c("lightgrey", "ivory", "blue", "yellow", "red")
  legend.limits = c(-3,-2,-1,0,1)
  legend.name = paste("KOG", "Abundance", sep="\n")
  legend.labels = c("Absent", "Undetectable", "RZ < BK", "RZ = BK", "RZ > BK")
  genome_order = factor(rownames(mat))
  grp_order = get_grp_name_ord(df)$ord
  
  # plot it!
  p3 = ggplot(grp.mat.melt, aes(x=factor(Var2), y=factor(Var1), fill=factor(value))) + 
    geom_tile() +
    scale_fill_manual(limits=legend.limits,
                      values=legend.values,
                      name=legend.name,
                      labels=legend.labels) + 
    scale_x_discrete(limits=grp_order) + 
    scale_y_discrete(limits=genome_order) +
    xlab("KOGs") + 
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          legend.title = element_text(size=22),
          legend.text = element_text(size=20),
          panel.background = element_rect(fill="white"),
		  panel.grid.major = element_blank(),
		  panel.grid.minor = element_blank())

	if ( show_xlabs == F ) {
		p3 = p3 + theme(axis.text.x = element_blank(),
				   axis.ticks.x = element_blank(),
					panel.background = element_rect(file="white"),
					panel.grid.major = element_blank(),
					panel.grid.minor = element_blank())
	}
  
	if ( png == T ) {
		file = paste(file_pre, ".png", sep="")
#		png(file, units = "in", res = 300, width = 5, height = 5)
		png(file, units="in", width=5, height=5, res=300)
	} else {
		file = paste(file_pre, ".pdf", sep="")
		pdf(file, width = 5, height = 5)
	}
	print(p3, vp=viewport())
	dev.off()
}

########
# MAIN #
########

# read in the file as a data frame
tbl = read.table(matrix_file, sep="\t")

# change it to a matirx (so I can get the names right)
mat = as.matrix(tbl[2:nrow(tbl), 2:ncol(tbl)])
rownames(mat) = tbl[2:nrow(tbl),1]
colnames(mat) = tbl[1,2:ncol(tbl)]

# transpose it (so that rows are genomes, cols are functional groups)
mat = t(mat)
make_figure(mat, out_file_pre, show_xlabs, png=png)


