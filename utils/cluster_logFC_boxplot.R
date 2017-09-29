#!/usr/bin/evn Rscript

### Description
# create a boxplot comparing the logFC between RZ and BK
# for a single cluster

### R Libraries
require("getopt", quietly=T)
require("ggplot2", quietly=T)

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
  "data", "d", 1, "character",
  "out_fig", "o", 1, "character",
  "verbose", "v", 0, "logical"
  ), byrow=TRUE, ncol=4)
opt = getopt(params)

# define parameter specified varaibles
if (! is.null(opt$verbose)) {
  verbose = opt$verbose
}


### Functions
main = function() {
	# read in the data
	data = read.table(opt$data, header=T, row.names=1, sep="\t");

	p = ggplot(data, aes(x=Enrichment, y=logFC)) +
		geom_boxplot(outlier.size=NA) +
		geom_jitter() +
		theme(
			text = element_text(size=20))
		ylab("logFC (RZ / BK)")

	ggsave(p, file = opt$out_fig)
}

# run the main function to execute the program
main()
