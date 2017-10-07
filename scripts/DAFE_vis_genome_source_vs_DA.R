#!/usr/bin/evn Rscript

### Description
# create a figure showing how the DA values match against the source 
# of the genome.

### R Libraries
require("getopt", quietly=T)
require("ggplot2", quietly=T)
require("reshape2", quietly=T)

### Defualt variables
verbose = FALSE

### Parameter variables
# The params matrix
# Each row is a parameter with 4 columns: 
#   1) long value, 
#   2) short value, 
#   3) argument type where 0 = no argument, 1 = required , 2 = optional
#   4) data type (logical, integer, double, complex, character)
# When you call this script the args are named like --verbose or -v
params = matrix(c(
  "da_tbl", "d", 1, "character",
  "genome_meta", "g", 1, "character",
  "out", "o", 1, "character",
  "verbose", "v", 0, "logical"
  ), byrow=TRUE, ncol=4)
opt = getopt(params)

# define parameter specified varaibles
if (! is.null(opt$verbose)) {
  verbose = opt$verbose
}


### Functions
main = function() {
	# read in the genome metadata
	meta = read.table(opt$genome_meta, header=T, row.names=1, sep="\t", quote="", comment.char="")

	# read in the da table
	da = read.table(opt$da_tbl, header=T, row.names=1, sep="\t")

	# count the number of RZ enriched and BK enriched instances for each genome
	rz_up = apply(da, 2, function(x) sum(x==1))
	bk_up = apply(da, 2, function(x) sum(x==-1))

	names = colnames(da)
	names = sapply(names, function(x) substring(x, 2))
	tmp_df = data.frame("genome" = names, "rz_up" = rz_up, "bk_up" = bk_up)

	print(head(tmp_df))

	tmp_df = melt(tmp_df)

	# merge
	final_df = merge(meta, tmp_df, by.x="row.names", by.y="genome")

	# make the figure
	ggplot(final_df, aes(x=Source, y=log10(value), col=variable)) +
		geom_boxplot(outlier.size=0) + 
		geom_point(position=position_jitterdodge(jitter.width=0.5), size=.5) +
		ylab("Number of Enrichmed COGs (log10)") +
		theme(axis.text = element_text(size=16),
				axis.title = element_text(size=16),
				legend.text = element_text(size=16),
				 axis.text.x = element_text(angle = 60, hjust=1))
	ggsave(opt$out)
}

# run the main function to execute the program
main()
