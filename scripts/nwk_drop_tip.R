#!/usr/bin/evn Rscript

### Description
# script for removing tip nodes from a phylogenetic tree

### R Libraries
require("getopt", quietly=T)
require(ape)

### Parameter variables
# The params matrix
# Each row is a parameter with 4 columns: 
#   1) long value, 
#   2) short value, 
#   3) argument type where 0 = no argument, 1 = required , 2 = optional
#   4) data type (logical, integer, double, complex, character)
# When you call this script the args are named like --verbose or -v
params = matrix(c(
  "tree", "t", 1, "character",
  "remove", "r", 1, "character",
  "out", "o", 1, "character"
), byrow=TRUE, ncol=4)
opt = getopt(params)

# opt$verbose


### Functions
main = function() {
  # read in the tips to remove
  to_remove = read.table(file = opt$remove, header=F)
	to_remove$V1 = as.character(to_remove$V1)
	print(to_remove)
  
  # read in the tree
  tree = read.tree(opt$tree)
  
  # remove the nodes
  new_tree = drop.tip(tree, to_remove$V1)
  
  # output the new tree
  write.tree(new_tree, file = opt$out)
}

# run the main function to execute the program
main()
