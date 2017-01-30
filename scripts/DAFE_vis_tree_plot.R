# this script creates the tree plot

library(ggplot2)
library(ggtree)
library(ape)
library(phytools) # for read.newick
library(getopt)

##############
# Parameters #
##############

# The params matrix
# Each row is a parameter with 4 columns: 
#   1) long value, 
#   2) short value, 
#   3) argument type where 0 = no argument, 1 = required , 2 = optional
#   4) data type (logical, integer, double, complex, character)
# When you call this script the args are named like --verbose or -v
params = matrix(c(
  "nwk_file", "n", 1, "character",
  "out_file", "o", 1, "character"
), byrow=TRUE, ncol=4)
opt = getopt(params)

# read in the tree file
tree = read.newick(opt$nwk_file)
tree = collapse.singles(tree)

# don't use ggtree!!!  I flips the tips around so they do not match
# the correct order!!!
# CORRECTION: you can use the laddreize=F parameter to fix the flipping
#ggtree(tree)

png(opt$out_file)
plot(tree, show.tip.label=F)
dev.off()