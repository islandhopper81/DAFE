#!/usr/bin/evn Rscript

### Description
# creates the mapping summary figures from the DAFE count output

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
  "count_file", "c", 1, "character",
  "ref_meta_file", "r", 1, "character",
  "out_prefix", "o", 1, "character",
  "title", "t", 1, "character",
  "verbose", "v", 0, "logical"
  ), byrow=TRUE, ncol=4)
opt = getopt(params)

# define parameter specified varaibles
if (! is.null(opt$verbose)) {
  verbose = opt$verbose
}


### Functions
main = function() {
	# read in the counts	
	counts = read.table(opt$count_file, header=T, row.names=1, sep="\t")
	counts = t(counts)  # transpose so genomes are rows
	row.names(counts) = substring(row.names(counts), 2) # remove that stupid X character from names
	sums = apply(counts, 1, sum)

	# read in the meta file
	meta = read.table(opt$ref_meta_file, header=T, row.names=1, sep="\t", quote="", comment.char="");
	meta = meta[names(sums),]
	
	# add the counts
	meta[,"counts"] = sums

	# make the pa/npa/soil plot
	p1 = ggplot(meta, aes(x=Label, y=log10(counts))) +
		geom_boxplot(outlier.size=NA) + 
		geom_point(position=position_jitter(width=0.3)) +
		theme(text = element_text(size=20),
			plot.title = element_text(hjust = 0.5)) +
		scale_x_discrete(limits = c("PA", "Soil", "NPA")) +
		ggtitle(opt$title) +
		ylab("Reads Mapped (log10)") +
		xlab("Plant Association Group")

	out_file = paste(opt$out_prefix, "_pa_npa_soil.pdf", sep="")
	ggsave(out_file, p1) 

	p2 = ggplot(meta, aes(x=Source, y=log10(counts))) +
		geom_boxplot(outlier.size=NA) +
		geom_point(position=position_jitter(width=.2)) +
		theme(text = element_text(size=20),
			plot.title = element_text(hjust = 0.5),
			axis.text.x = element_text(angle = 90, hjust = 1)) +
		ggtitle(opt$title) +
		ylab("Reads Mapped (log10)") +
		xlab("Isolation Environment")

	out_file = paste(opt$out_prefix, "_iso_env.pdf", sep="")
	ggsave(out_file, p2)

		
}

# run the main function to execute the program
main()
