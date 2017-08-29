#!/usr/bin/evn Rscript

### Description
# This R script prints the files for the top X cogs that are enriched in 
# either the rhizosphere or bulk soil (or other comparison groups).  It 
# creates the file that is used as an input to create the evolview files
# for printing a evolview tree with the top examples.

### R Libraries
require("getopt", quietly=T)

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
  "da_tbl", "t", 1, "character",
  "count", "c", 1, "integer",
  "direction", "d", 1, "character",
  "out_dir", "o", 1, "character",
  "verbose", "v", 0, "logical"
  ), byrow=TRUE, ncol=4)
opt = getopt(params)

# define parameter specified varaibles
if (! is.null(opt$verbose)) {
  verbose = opt$verbose
}


### Functions
main = function() {
	# read in the da_table
	da_tbl = read.table(opt$da_tbl, header=T, row.names=1, sep="\t")

	# set the direction
	if ( opt$direction == "up" ) {
		direction = 1
	} else if ( opt$direction == "dn" ) {
		direction = -1
	} else {
		stop("--direction must be \"up\" or \"dn\"")
	}

	# find the top COGs
	# first look at in how many genomes a COG is called as DA
	da_counts = apply(da_tbl == direction, 1, sum)

	# get the top X COGs that are called enriched in the most genomes
	top = sort(da_counts, decreasing = T)[0:opt$count]
	print(top[0:opt$count])

	# print the vector of calls for each genome for each of the top COGs
	for ( t in names(top) ) {
		print(t)
		file = paste(opt$out_dir, "/", t, ".txt", sep="")
		vec = t(da_tbl[t,])
		write.table(vec, file = file, sep="\t", row.names=T, col.names=F, quote=F)
	}
	
}

# run the main function to execute the program
main()
