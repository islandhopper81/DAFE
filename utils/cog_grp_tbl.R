#!/usr/bin/evn Rscript

### Description
# creates a table of COG group info that will be submitted in the paper

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
	# read in the data
	data = read.table(opt$da_tbl, header=T, row.names=1, sep="\t")

	my_list = apply(as.matrix(data), 1, function(x) table(x))

#	list$A["-1"] = 0
#	list$B["-1"] = 0
#	list$B["1"] = 1

	df = data.frame("-2" = integer(), "-1" = integer(), "0" = integer(), "1" = integer())
	df = data.frame(matrix(,ncol=5, nrow=length(names(my_list))))
	row.names(df) = names(my_list)
	colnames(df) = c("-3", "-2", "-1", "0", "1")

	for (name in names(my_list)) {
		df[name,"-3"] = my_list[[name]]["-3"]
		df[name,"-2"] = my_list[[name]]["-2"]
		df[name,"-1"] = my_list[[name]]["-1"]
		df[name,"0"] = my_list[[name]]["0"]
		df[name,"1"] = my_list[[name]]["1"]
	}

	# convert NA values to 0
	df[is.na(df)] = 0

	# reorder the columns
	df = df[,c("1", "-1", "0", "-2", "-3")]

	# rename the columns
	colnames(df) = c("RZ Enriched", "BK Enriched", "No Difference", "Undetectable", "Absent")
	

	print(df)
	write.table(df, file=opt$out, sep="\t", quote=F)
}

# run the main function to execute the program
main()
