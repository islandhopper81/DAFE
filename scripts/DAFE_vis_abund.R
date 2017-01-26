# make the abundance heatmap.

# Note: This script uses some external code (heatmap.3) that
# it retrieves from the web.  So this script must have access
# to the internet.


library(reshape)
library(ggplot2)
library(RColorBrewer)
library(getopt)
library(devtools)
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")


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
  "meta_file", "m", 1, "character",
  "count_file", "c", 1, "character",
  "out_file", "o", 1, "character"
), byrow=TRUE, ncol=4)
opt = getopt(params)

#############
# Functions #
#############

get_colors = function(meta, palette) {
  my_colors = data.frame(matrix(nrow = nrow(meta)))
  i = 1
  for (name_i in 1:length(names(meta))) {
    vals = unique(meta[,name_i])
    vals = vals[!is.na(vals)]
    my_n = length(vals)
    my_s = i
    name = names(meta)[name_i]
    
    #print(paste("name: ", name, sep=""))
    #print(paste("my_n: ", my_n, sep=""))
    #print(paste("my_s: ", my_s, sep=""))
    
    col_vec = get_color_vec(palette, n=my_n, s=my_s)
    names(col_vec) = vals
    
    merged = merge(col_vec, meta, by.x="row.names", by.y=names(meta)[name_i], 
                   all.y=T, sort=F)
    my_colors[,names(meta)[name_i]] = as.character(merged$x)
    
    i = my_n + i
  }
  
  my_colors = my_colors[2:(ncol(meta) + 1)]
  my_colors_col_names = names(my_colors)
  my_colors = matrix(unlist(my_colors), ncol=ncol(meta), nrow=nrow(meta))
  colnames(my_colors) = my_colors_col_names
  rownames(my_colors) = row.names(meta)
  
  return(my_colors)
}
get_color_vec = function(palette, n, s) {
  # for a special case when n < 3 to avoid an error/warning
  if ( n+s-1 < 3 ){
    color_vec = brewer.pal(3, palette)
  } else {
    color_vec = brewer.pal(n+s-1, palette)
  }
  return(color_vec[s:(n+s-1)])
}
.simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1, 1)), substring(s, 2),
        sep = "", collapse = " ")
}


########
# MAIN #
########

# get the sample metadata
#meta_file = "/Users/Scott/Projects/CompMetaG/experiments/nature_paper_44_metaG/bacterial/metagenome_metadata.txt"
meta = read.table(opt$meta_file, row.names=1, header=T)
meta = meta[23:66,c("genotype", "fraction", "age")]

# get the count data
#count_file = "/Users/Scott/Projects/CompMetaG/experiments/nature_paper_44_metaG/bacterial/figures/id60_abundance/bacteria_counts_by_genome.txt"
counts = read.table(opt$count_file, sep="\t", header=T, row.names=1, check.names = F)


my_mat = data.matrix(matrix(unlist(counts), ncol=ncol(counts), byrow=F))
my_mat = t(my_mat)
mode(my_mat) <- "numeric" 
colnames(my_mat) = row.names(counts)
my_mat = apply(my_mat, 2, rev)
my_mat = log10(my_mat + 1)

cols = get_colors(meta, "Set1")
colnames(cols) = sapply(colnames(cols), .simpleCap)

legend_names = unique(unlist(lapply(meta, unique)))
legend_names = legend_names[!is.na(legend_names)]
legend_cols = unlist(apply(cols, 2, unique))
legend_cols = legend_cols[!is.na(legend_cols)]

png(opt$out_file)

heatmap.3(my_mat, dendrogram = "column",
         Rowv = NA, Colv = T, scale="column",
         ColSideColors = cols,
         labRow=F, labCol=F, xlab="Samples", ylab="Genomes",
         ColSideColorsSize = 3,
         hclustfun = hclust,
         distfun = dist,
         margins=c(2,2))
legend("left",
       legend=legend_names,
       fill = legend_cols,
       border=FALSE, bty="n", y.intersp = 1.2, cex=1)

dev.off()

