#!/usr/bin/evn Rscript

### Description
# creates boxplot of a subset of genomes to show their abundance 
# and how it is distributed across the different fractions

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

	p = ggplot(data, aes(x=Genome, y=abund, color=Fraction)) +
		geom_boxplot(outlier.size=NA) +
		geom_point(position = position_jitterdodge()) + 
		theme(
			text = element_text(size=20),
			axis.text.x = element_text(angle = 90, hjust = 1)) + 
		ylab("Marker Read Counts")

	ggsave(p, file = opt$out_fig);

	# a second option where the RZ and BK genomes are all combined
	p = ggplot(data, aes(x=Fraction, y=abund)) +
		geom_boxplot(outlier.size=NA) +
		geom_jitter(aes(color=Genotype)) +
		theme(
			text = element_text(size=20),
			axis.text.x = element_text(angle = 90, hjust = 1)) +
		ylab("Marker Read Counts")

	ggsave(p, file = opt$out_fig)
}

# run the main function to execute the program
main()
