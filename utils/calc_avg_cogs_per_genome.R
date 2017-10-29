#!/usr/bin/evn Rscript

### Description
# creates a boxplot of points where each point shows the fraction
# of the genome that has COG annotations.

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
	data = read.table(opt$data, header=F, sep="\t")

	p = ggplot(data, aes(y=V2, x=1)) + 
		geom_boxplot(outlier.size=NA) +
		geom_jitter(width=0.3) +
		xlab("") +
		ylab("Fraction") +
		ggtitle("Fraction of Genome with COG Annotations") +
		coord_flip() +
		theme(text = element_text(size=20),
				plot.title = element_text(hjust = 0.5),
				axis.ticks.y = element_blank(),
				axis.text.y = element_blank(),
				panel.grid.major.y = element_blank(),
				panel.grid.minor.y = element_blank())

	ggsave(opt$out_fig, p, width=8, height=2.5)
}

# run the main function to execute the program
main()
