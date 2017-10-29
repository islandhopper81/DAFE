#!/usr/bin/evn Rscript

### Description
# looks for a relationship between logFC and the number of RZ enriched
# COG.  I'm doing this because there wasn't a clear relation ship between
# mean RZ abundance and the number of RZ enriched COGs.  However, this 
# could be because some things that are very abundant in the RZ might 
# also be very abundant in the soils which would lessen the chance of
# that genome having many RZ enriched COGs.

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
  "tags", "t", 1, "character",
  "genome_sum", "g", 1, "character",
  "out_fig_prefix", "o", 1, "character",
  "verbose", "v", 0, "logical"
  ), byrow=TRUE, ncol=4)
opt = getopt(params)

# define parameter specified varaibles
if (! is.null(opt$verbose)) {
  verbose = opt$verbose
}


### Functions
main = function() {
	tags = read.table(opt$tags, header=T, row.names=1, sep="\t")
	rz_count = read.table(opt$genome_sum, header=T, row.names=1, sep="\t")

	# combine the two tables
	tags = tags[row.names(rz_count),]
	all_data = rz_count
	all_data[,"logFC"] = tags$logFC

	# create the RZ enrichment plot
	p = ggplot(all_data, aes(x=RZ.Enriched, y=logFC)) +
		geom_point() +
		geom_smooth(method="lm") +
		xlab("Number of RZ Enriched COGs") +
		ylab("Abundance Log Fold Change") +
		theme(text = element_text(size=20))

	file = paste(opt$out_fig_prefix, "_RZ_enr.pdf", sep="")
	ggsave(file, p)
	
	# create the BK enrichment plot
	p = ggplot(all_data, aes(x=BK.Enriched, y=logFC)) +
		geom_point() +
		geom_smooth(method="lm") +
		xlab("Number of BK Enriched COGs") +
		ylab("Abundance Log Fold Change") +
		theme(text = element_text(size=20))

	file = paste(opt$out_fig_prefix, "_BK_enr.pdf", sep="")
	ggsave(file, p)
}

# run the main function to execute the program
main()
