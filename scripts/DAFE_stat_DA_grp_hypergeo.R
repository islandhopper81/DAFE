# this script contains the functions for running a hypergeometric 
# test on a give set of genomes

# x = number of times a COG group is DA in a genome for a given direction
# m = number of times a COG group is observed in a genome
# n = number of times all COGs groups besides those counted in "m"
# k = number of DA COGs in a genome

# I need a double for loop
# 1) loop over all the genomes in the matrix
# 2) loop over all the COGs and run the test

library("stringi") # for function stri_detect
require("getopt", quietly=T)


###
# Inputs
# cog grp metadata file (for list of all pissible cog grps)
#cog_grp_metadata_file = "/Users/Scott/temp/temp/id60_hypergeo/cog_grp_metadata.txt"
#cog_grps = read.table(cog_grp_metadata_file, header=T)
#cog_group_names = cog_grps$COG_GRP

# DA matrix
#mat_file = "/Users/Scott/temp/temp/id60_hypergeo/full_da_tbl_18254.txt"
#mat = read.table(mat_file, header=T, row.names=1, sep="\t")

# COG metadata file
#cog_metadata_file = "/Users/Scott/temp/temp/id60_hypergeo/cognames2003-2014.tab"
#cog_meta = read.table(cog_metadata_file, header=T, sep="\t", na.strings = "", 
#                      comment.char = "", quote="")

# output file
#out_file = "/Users/Scott/temp/temp/id60_hypergeo/out.txt"


### Parameter variables
# The params matrix
# Each row is a parameter with 4 columns: 
#   1) long value, 
#   2) short value, 
#   3) argument type where 0 = no argument, 1 = required , 2 = optional
#   4) data type (logical, integer, double, complex, character)
# When you call this script the args are named like --verbose or -v
params = matrix(c(
  "mat_file", "m", 1, "character",
  "cog_grp_metadata_file", "g", 1, "character",
  "cog_metadata_file", "c", 1, "character", 
  "out_file", "o", 1, "character"
), byrow=TRUE, ncol=4)
opt = getopt(params)


# The main function that runs everything
main = function(mat_file, cog_grp_metadata_file, cog_metadata_file, out_file) {
  mat = read.table(mat_file, header=T, row.names=1, sep="\t")
  
  cog_grps = read.table(cog_grp_metadata_file, header=T)
  cog_group_names = cog_grps$COG_GRP
  
  cog_meta = read.table(cog_metadata_file, header=T, sep="\t", na.strings = "", 
                        comment.char = "", quote="")
  
  
  up_mat = main_hyper(cog_group_names, mat, cog_meta, "up")
  dn_mat = main_hyper(cog_group_names, mat, cog_meta, "down")
  
  # convert the 1's to 2's so I can add the two matrices and 
  # still distinguish between the up and down values
  #hypergeo_dn_mat_obj_tmp$binary_mat[hypergeo_dn_mat_obj_tmp$binary == 1] = 2
  dn_mat[dn_mat == 1] = 2
  
  # add the two matrices
  #both_binary_matrix = hypergeo_dn_mat_obj_tmp$binary_mat + hypergeo_up_mat_obj$binary_mat
  combined_mat = dn_mat + up_mat
  
  write.table(combined_mat, out_file, quote=F, row.names=T, col.names=T, sep="\t")
}

# function for running all the preprocessing and hypergeometric test
main_hyper = function(cog_group_names, mat, cog_meta, direction) {
  # create a matrix to store which COG groups are enriched in the DA COGs
  x_mat = matrix(nrow = length(cog_group_names),
                 ncol = ncol(mat), 
                 dimnames = list(cog_group_names, colnames(mat)))
  m_mat = matrix(nrow = length(cog_group_names),
                 ncol = ncol(mat), 
                 dimnames = list(cog_group_names, colnames(mat)))
  
  # set the direction vector
  direction_vec = get_direction_vec(direction)
  
  x_mat = get_x_mat(mat, x_mat, 
                    cog_meta, direction_vec,
                    cog_group_names)
  
  # get m matrix
  m_mat = get_x_mat(mat, m_mat, 
                    cog_meta, c(-1, 0, 1),
                    cog_group_names)
  
  # run hyper genometric test for each COG group
  hyper_mat = run_hypergeo(x_mat, m_mat)
  
  # get the binary values (ie convert into significant and not significant classes)
  bin_hyper_mat = convert_to_binary(hyper_mat, 0.05)
  
  return(bin_hyper_mat)
  
  # print the hypergeo matrix
  #write.table(bin_hyper_mat, out_file, quote=F, row.names=T, col.names=T, sep="\t")
}

convert_to_binary = function(hyper_mat, pval) {
  hyper_mat[hyper_mat > pval] = 2
  hyper_mat[hyper_mat <= pval] = 1
  hyper_mat[hyper_mat == 2] = 0
  
  return(hyper_mat)
}

run_hypergeo = function(x_mat, m_mat) {
  # make sure the dimensions of the matries are the same
  if ( nrow(x_mat) != nrow(m_mat) ) {
    stop("ERROR: nrow are not equal for x_mat and m_mat")
  }
  if ( ncol(x_mat) != ncol(m_mat) ) {
    stop("ERROR: ncol are not equal for x_mat and m_mat")
  }
  
  # build a matrix that will hold 0 if a COG group is not 
  # significant and 1 if it is significant for each COG group
  # for each genome
  out_mat = matrix(nrow = nrow(x_mat), ncol = ncol(x_mat),
                   dimnames = list(row.names(x_mat), colnames(x_mat)))
  
  # x = number of times a COG group is DA in a genome for a given direction
  # m = number of times a COG group is observed in a genome
  # n = number of times all COGs groups besides those counted in "m"
  # k = number of DA COGs in a genome
  
  # for each genome test each COG group to see if it is 
  # significant
  # Remember the columns are genomes and the rows are COG groups
  for (r in 1:nrow(x_mat)) {
    for (c in 1:ncol(x_mat)) {
      x = x_mat[r,c]
      m = m_mat[r,c]
      n = sum(m_mat[,c]) - m
      k = sum(x_mat[,c])
      out_mat[r,c] = phyper(q = x, m = m, n = n, k = k)
    }
  }
  
  return(out_mat)
}

get_x_mat = function(mat, x_mat, cog_annot_tbl, 
                     codes_to_keep, cog_group_names) {
  # for each genome in the cog_model_obj$mat object
  # get a COG group vector for the COGs that are DA (ie >1)
  # and add it to the x_mat
  for (i in 1:ncol(mat)) {
    vec = get_COG_group_vec(mat[,i], cog_annot_tbl, codes_to_keep)
    vec = add_missing_COG_groups(vec, cog_group_names)
    x_mat[,i] = vec$x
  }
  return(x_mat)
}

# this function generates a vector (which is actually just a column
# in the returned table -- sorry the names are bad) that has the number
# of times a COG group is observed.  The "codes_to_keep" parameter is
# a vector of DA codes (ie -2, -1, 0, 1) to count as present towards a 
# COG group
get_COG_group_vec = function(cog_da_vec, cog_annot_tbl, codes_to_keep) {
  cog_da_vec = cog_da_vec[which(cog_da_vec %in% codes_to_keep)]
  
  # check to see if there is anything in cog_da_vec
  # sometimes there are no COGs that are ==1 in which
  # case return the empty table.
  if ( length(cog_da_vec) == 0 ) {
    return(data.frame(func=factor(), x=integer()))
  }
  
  # generate a table that is a combination of the COG, its group,
  # the annotation (ie name), and the count (which will be 1)
  tbl = merge(cog_annot_tbl, as.matrix(cog_da_vec),
              by.x="row.names", by.y="row.names")
  
  new_tbl = split_multi_groups(tbl)
  
  # change all the values in the V1 column to 1 so that the sum function
  # when I aggregate will be appropriate.  
  new_tbl$V1 = rep(1, length(new_tbl$V1))
  
  vec = aggregate(as.integer(new_tbl$V1), by=list(func=new_tbl$func), 
                  FUN=sum, na.rm=T)
  
  return(vec)
}

# This function splits the groups that have two COG groups codes (ie AB)
# into two different groups.  Right now I am assuming that when a COG is in
# two groups it is part of BOTH groups.  So it gets counted as in the "A"
# group and in the "B" group.
split_multi_groups = function(tbl) {
  lengths = sapply(as.character(tbl$func), nchar)
  double_entries = tbl[lengths > 1,]
  
  # if there are no double entries then simply return the original tbl
  if ( length(double_entries$func) == 0) {
    return(tbl)
  }
  
  # set the row names of the double entries as increasing integers
  row.names(double_entries) = 1:nrow(double_entries)
  
  # make a matrix to hold the lines that have double entries when converted to 
  # a single entry
  nrow = sum(lengths[lengths>1])
  split_mat = matrix(nrow=nrow, ncol=ncol(tbl), 
                     dimnames=list(c(), colnames(double_entries)))
  
  # create a table where the double entries have been converted
  # to single entries
  j = 1 # current index of split_mat when adding a new single entry
  # i = current index in the double_entries matrix
  for (i in 1:nrow(double_entries)) {
    values = strsplit(as.character(double_entries$func[i]), split="")
    for (v in values[[1]]) {
      new_single_entry = double_entries[i,]
      new_single_entry$func = v
      split_mat[j,] = as.matrix(new_single_entry)
      j = j + 1
    }
  }
  
  new_tbl = rbind(tbl[lengths<2,], split_mat)
  
  return(new_tbl)
}

# this function adds the COG groups that are missing from a 
# genome.  They are added to the vec object with a count of 0
# where the count column is really named "x".
add_missing_COG_groups = function(vec, cog_group_names) {
  missing_groups_names = cog_group_names[!(cog_group_names %in% vec$func)]
  missing_groups_df = data.frame(func=as.factor(missing_groups_names),
                                 x=rep.int(0, times=length(missing_groups_names)))
  
  new_vec = rbind(vec, missing_groups_df)
  
  # sort the new_vec by the func column
  new_vec = new_vec[with(new_vec, order(func)),]
  
  return(new_vec)
}


# gets the vector of number to look for depending on the direction
# parameter value
get_direction_vec = function(direction="up") {
  direction = tolower(direction)
  if (direction == "up") {
    vec = c(1)
  } else if (direction == "down") {
    vec = c(-1)
  } else {
    stop(paste("Unrecognized direction value: ", direction, sep=""))
  }
  
  return(vec)
}





# run the main function now
main(opt$mat_file, opt$cog_grp_metadata_file, 
     opt$cog_metadata_file, opt$out_file)

