#!/usr/bin/evn Rscript

### Description
# Prints summary numbers and info about a da table

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
  "da_tbl", "d", 1, "character",
  "verbose", "v", 0, "logical"
  ), byrow=TRUE, ncol=4)
opt = getopt(params)

# define parameter specified varaibles
if (! is.null(opt$verbose)) {
  verbose = opt$verbose
}


### Functions
main = function() {
	# read in the da table
	tbl = read.table(opt$da_tbl, header=T, row.names=1, sep="\t")

	# print the da talbe file
	print(paste("DA table file:", opt$da_tbl))

	# print the current workind dir
	print(paste("Current working directory:", getwd()))

	# print the date
	print(paste("Date:", Sys.Date()))

	# print the number of tests (ie all cells with >= -2
	tests = table(tbl > -2)["TRUE"]
	print(paste("Total number of tests (ie > -2):", tests))

	# print the number of 1 values (enriched in target group)
	val_1_count = table(tbl == 1)["TRUE"]
	print(paste("Number of enriched target group (ie 1 values):", val_1_count))

	# print the number of -1 values (enriched in reference group)
	val_m1_count = table(tbl == -1)["TRUE"]
	print(paste("Number of enriched reference group (ie -1 values):", val_m1_count))

	# print number of genomes that have at least one enriched value
	tmp = tbl
	tmp[tmp == -1] = 1
	genome_with_enr = sum(apply(tmp, 2, function(x) any(x == 1)))
	print(paste("Number of genomes with at least one enriched value:", genome_with_enr))

	# print number of genomes that have at least one "1" value
	genome_with_1 = sum(apply(tbl, 2, function(x) any(x==1)))
	print(paste("Number of genomes with at least one target enrichment (ie 1):", genome_with_1))

	# print number of genomes that have at least one "-1" value
	genome_with_m1 = sum(apply(tbl, 2, function(x) any(x==-1)))
	print(paste("Number of genomes with at least one referene enrichment (ie -1):", genome_with_m1))

	# print number of genoems that have at least one "1" AND at least one "-1" value
	genome_with_1 = apply(tbl, 2, function(x) any(x==1))
	genome_with_m1 = apply(tbl, 2, function(x) any(x ==-1))
	genome_with_both = sum(genome_with_1 & genome_with_m1)
	print(paste("Number of genomes with at least one of both enrichments:", genome_with_both))

	# print number of COGs that had at least one genome where there was a target enrichment
	cog_with_1 = sum(apply(tbl, 1, function(x) any(x==1)))
	print(paste("Number of COGs with at least one genome where the COG was target enriched:", cog_with_1))

	# print number of COGs that had at least one genome where there was a reference enrichment
	cog_with_m1 = sum(apply(tbl, 1, function(x) any(x==-1)))
	print(paste("Number of COGs with at least one genome where the COG was reference enriched:", cog_with_m1))

	# print number of COGs that had at least one genome where there was either a target or reference enrichment
	cog_with_either = sum(apply(tmp, 1, function(x) any(x==1)))
	print(paste("Number of COGs with at least one goenme where the COG was either target or reference enriched:", cog_with_either))

	# print number of COGs that had one of each enrichment
	cog_with_1 = apply(tbl, 1, function(x) any(x==1))
	cog_with_m1 = apply(tbl, 1, function(x) any(x==-1))
	cog_with_both = sum(cog_with_1 & cog_with_m1)
	print(paste("Number of COGs with both a target and reference enrichment:", cog_with_both))

	# print the top 10 most frequently DA functional elements for the 1 group
	top_1 = sort(apply(tbl, 1, function(x) sum(which(x == 1))), decreasing = T)
	print("Top 10 most frequently DA functional elements in the target group (ie 1):")
	print(top_1[1:10])

	# print the top 10 most frequently DA functional elements for the -1 group
	top_m1 = sort(apply(tbl, 1, function(x) sum(which(x == -1))), decreasing = T)
	print("Top 10 most frequently DA functional elements in the reference group (ie -1):")
	print(top_m1[1:10])
	
	# print the top 10 genomes with the most 1 group
	top_1_genomes = sort(apply(tbl, 2, function(x) sum(which(x == 1))), decreasing = T)
	print("Top 10 most genomes with the most DA calls for the target group (ie 1):")
	print(top_1_genomes[1:10])

	# print the top 10 genomes with the most -1 group
	top_m1_genomes = sort(apply(tbl, 2, function(x) sum(which(x == -1))), decreasing = T)
	print("Top 10 most genomes with the most DA calls for the reference group (ie -1):")
	print(top_m1_genomes[1:10])

	# print the top 10 genomes with the most total DA functions
	top_both_genomes = sort(apply(tmp, 2, function(x) sum(which(x == 1))), decreasing = T)
	print("Top 10 most genomes with the most DA calls for both group (ie 1 and -1):")
	print(top_both_genomes[1:10])

	# other ?
}

# run the main function to execute the program
main()
