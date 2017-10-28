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
  "sample_meta_file", "m", 1, "character",
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
		geom_point(position=position_jitter(width=0.25)) +
		theme(text = element_text(size=20),
			plot.title = element_text(hjust = 0.5)) +
		scale_x_discrete(limits = c("PA", "Soil", "NPA")) +
		ylab("Reads Mapped (log10)")
		
	  # I removed these so I could make it skinny
	  		#ggtitle(opt$title) +
			#xlab("Plant Association Group")

	out_file = paste(opt$out_prefix, "_pa_npa_soil.pdf", sep="")
	ggsave(out_file, p1, width=3, height=6) 

	# make the plot graphed by genome environment
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

	# make the plot showing the number of reads that map from each 
	# metagenome environment
	sample_sums = apply(counts, 2, sum)
	
	sample_meta = read.table(opt$sample_meta_file, header=T, row.names=1, sep="\t")
	sample_meta = sample_meta[names(sample_sums),]
	sample_meta[,"counts"] = sample_sums

	p3 = ggplot(sample_meta, aes(x=fraction, y=counts)) +
		geom_boxplot(outlier.size=NA) +
		geom_point(position=position_jitter(width=.2)) +
		theme(text = element_text(size=20),
			plot.title = element_text(hjust=0.5)) +
		ylab("Reads Mapped (log10)") +
		xlab("Metagenome Samples")

	out_file = "reads_mapped_per_meta_env.pdf"
	ggsave(out_file, p3)
	
	
	### make a plot similar to p1 but split by fraction
	# remember the genomes are rows in the counts table
	rz_samples = sample_meta$fraction == "RZ"
	bk_samples = sample_meta$fraction == "BK"
	rz_sum = apply(counts[,rz_samples], 1, sum)
	bk_sum = apply(counts[,bk_samples], 1, sum)
	
	# add the rz and bk sums to the meta
	meta[,"rz_counts"] = rz_sum
	meta[,"bk_counts"] = bk_sum
	
	# melt the table to get it into long format for plotting
	molten = melt(meta, id.vars=c("Label"), measure.vars = c("rz_counts", "bk_counts"))
	
	# make the pa/npa/soil plot
	p4 = ggplot(molten, aes(x=Label, y=log10(value), color=variable)) +
		geom_boxplot(outlier.size=NA) + 
		geom_point(position=position_jitterdodge()) +
		theme(text = element_text(size=16),
			plot.title = element_text(hjust = 0.5)) +
		scale_x_discrete(limits = c("PA", "Soil", "NPA")) +
		scale_color_discrete("Fraction", labels = c("RZ", "BK")) +
		ylab("Reads Mapped (log10)")

	out_file = paste(opt$out_prefix, "_pa_npa_soil_by_frac.pdf", sep="")
	ggsave(out_file, p4)
	
	
	### make a plot similar to p2 but split by fraction
	# melt the table to get it into long format for plotting
	molten = melt(meta, id.vars=c("Source"), measure.vars = c("rz_counts", "bk_counts"))
	
	# make the pa/npa/soil plot
	p5 = ggplot(molten, aes(x=Source, y=log10(value), color=variable)) +
		geom_boxplot(outlier.size=NA) + 
		geom_point(position=position_jitterdodge()) +
		theme(text = element_text(size=16),
			plot.title = element_text(hjust = 0.5),
			axis.text.x = element_text(angle = 90, hjust = 1)) +
		scale_color_discrete("Fraction", labels = c("RZ", "BK")) +
		ylab("Reads Mapped (log10)")

	out_file = paste(opt$out_prefix, "_iso_env_by_frac.pdf", sep="")
	ggsave(out_file, p5)
	
	
	# print out the top 5 rhizosphere colonizers
	sorted = meta[order(-meta$rz_count),]
	print("Top 5 RZ colonizers")
	print(sorted[1:5,])
	
	# print out the top 5 rhizosphere colonizers
	sorted = meta[order(-meta$bk_count),]
	print("Top 5 BK colonizers")
	print(sorted[1:5,])
}

# run the main function to execute the program
main()
