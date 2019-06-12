#!/usr/bin/evn Rscript

### Description
# compares the number of reads that map to the database of genomes
# to the number of reads that map to the de novo assembly

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
  "db", "d", 1, "character",
	"assembly", "a", 1, "character",
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

		#Read in the files
		db <- read.delim(opt$db, header=TRUE, sep="\t")
		asm <- read.delim(opt$assembly, header=TRUE, sep="\t")
		colnames(db) <- colnames(asm)

		# run a t test to see if they are significnatly different
		my_test = t.test(db$Perc_map, asm$Perc_map)
		print(my_test)

		#merge the two datasets for plotting
		db$Map_type <- "Database"
		asm$Map_type <- "Assembly"
		combined_stat <- rbind(db,asm)

		#make into data frame
		stat_df <- as.data.frame(combined_stat)

		stat_plot <- ggplot(stat_df, aes( x=Map_type, y=Perc_map) )+
			geom_boxplot(outlier.size = NA)+
			geom_jitter(width=.2) +
			theme(	text = element_text(size = 12),
					axis.title.y = element_text(size=10),
					plot.title = element_text(hjust = .5),
					axis.text.x = element_text(angle = 90, hjust = 1),
					legend.position = "none",
					panel.grid.minor = element_blank()
				)+
			labs(x = "", y="% of Reads Mapped")

		ggsave(opt$out_fig, plot=stat_plot, width = 1.95, height=2.5)
}

# run the main function to execute the program
main()
