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

get_colors = function(meta, pal = "set1", columns) {
  if (length(columns) == 0 ) {
    columns = names(meta)
  }
  uniq_vals = levels(unique(melt(meta, measure.vars = columns)$value))
  col_count = length(uniq_vals)
  col_count
  
  uniq_cols = brewer.pal(col_count, pal)
  
  df = data.frame(val = uniq_vals, col = uniq_cols)
  
  
  col_df = data.frame(matrix(nrow=nrow(meta), ncol=ncol(meta)))
  names(col_df) = names(meta)
  row.names(col_df) = row.names(meta)
  for ( c in colnames(meta) ) {
    cols = droplevels(df$col[match(meta[,c], df$val)])
    col_df[,c] = cols
  }
  cols = col_df
  colnames(cols) = sapply(colnames(cols), .simpleCap)
  
  return(as.matrix(cols))
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

cols = get_colors(meta, "Set1", columns=c("genotype", "fraction", "age"))
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

