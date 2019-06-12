#!/usr/bin/evn Rscript

### Description
# creates a barchart showing the number of instances of each groups enrichment
# for both fractions

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
  "out_fig", "o", 1, "character",
  "up_grp", "u", 1, "character",
  "dn_grp", "n", 1, "character",
  "flip_cols", "f", 0, "logical",
  "xmin", "i", 0, "integer",
  "xmax", "a", 0, "integer",
  "verbose", "v", 0, "logical"
  ), byrow=TRUE, ncol=4)
opt = getopt(params)

# define parameter specified varaibles
if (! is.null(opt$verbose)) {
  verbose = opt$verbose
}

# set the value for flip_cols
if ( is.null(opt$flip_cols) ) {
	flip_cols = F
} else if ( opt$flip_cols == T ) {
	flip_cols = T
} else {
	flip_cols = F
}


### Functions
main = function() {
	# read in the da table
	tbl = read.table(opt$da_tbl, header=T, row.names=1, sep="\t")

	up_counts = apply(tbl, 1, function(x) sum(x == 1))
	dn_counts = apply(tbl, 1, function(x) sum(x == -1))

	# set the full names
	grp_names = c( "RNA processing and modification",
					"Chromatin structure and dynamics",
					"Energy production and conversion",
					"Cell cycle control and mitosis",
					"Amino Acid metabolism and transport",
					"Nucleotide metabolism and transport",
					"Carbohydrate metabolism and transport",
					"Coenzyme metabolism and transport",
					"Lipid metabolism and transport",
					"Translation",
					"Transcription",
					"Replication and repair",
					"Cell wall/membrane/envelope biogenesis",
					"Cell motility",
					"Post-translational modification, protein turnover, chaperones",
					"inorganic ion transport and metabolism",
					"Secondary metabolites biosynthesis, transport and catabolism",
					"General function prediction only",
					"Function unknown",
					"Signal transduction",
					"Intracellular trafficking and secretion",
					"Defense mechanisms",
					"Extracellular structures",
					"Mobilome: prophages, transposons",
					"Cytoskeleton"
				)

	grp_ids = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", 
				 "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Z")
	grp_names = paste(grp_ids, grp_names, sep=" - ")
	names(grp_names) = grp_ids
	
	data = data.frame(up_counts = up_counts,
						dn_counts = dn_counts,
						grp_names = grp_names[names(up_counts)], 
						id = names(up_counts))

	my_melt = melt(data, varnames=c("rz_count", "dn_counts"))

	# get the order of the groups
	# UPDATE: I decided not to do this.  I want the order to be alphabetic.
	order = c("C", "E", "J", "L", "M", "G", "P", "O", "H", "I", "R", "K", "F",
				 "V", "Q", "S", "T", "x", "N", "D", "U", "B", "Z", "A", "W")

	# set the color and ordering 
	my_cols = c("red", "blue")
	my_labels = c(paste(opt$dn_grp, "Enriched"), paste(opt$up_grp, "Enriched"))
	my_breaks = c("dn_counts", "up_counts")

	if ( flip_cols == T ) {
		my_cols = rev(my_cols)
		my_labels = rev(my_labels)
		my_breaks = rev(my_breaks)
	}

	# create the plot
	p = ggplot(my_melt, aes(x=grp_names, y=value, fill=variable)) +
		geom_bar(stat="identity", position=position_dodge(width=.6), width=0.5) +
		scale_fill_manual(name = "",
							values=my_cols,
							labels = my_labels, 
							breaks = my_breaks) +
		ylab("Number of Enrichments\nAcross All Genomes") +
		xlab("") +
		theme_bw() +
		theme(legend.position="bottom",
				legend.direction="vertical",
				panel.grid.major.y = element_blank(),
				panel.grid.minor.y = element_blank(),
				legend.text = element_text(size=20),
				axis.title.x = element_text(size=18),
				axis.text.x = element_text(size=16)) +
		scale_x_discrete(position = "top", limits = rev(factor(data$grp_names))) +
		scale_y_reverse() +
		coord_flip()

	print(opt$xmin)
	if ( opt$xmax > 0 ) {
		print("here!!!!")
		p = p + scale_y_reverse(limits = c(opt$xmax, opt$xmin))
	}

	ggsave(opt$out_fig, p)

}

# run the main function to execute the program
main()
