# this object contains the functions for running a hypergeometric 
# test on a give set of genomes

# x = number of times a COG group is DA in a genome for a given direction
# m = number of times a COG group is observed in a genome
# n = number of times all COGs groups besides those counted in "m"
# k = number of DA COGs in a genome

# I need a double for loop
# 1) loop over all the genomes in the matrix
# 2) loop over all the COGs and run the test

library("stringi") # for function stri_detect

# direction -- "up" test for COG groups that are enriched in the up COGs
#              "down" tet for COG groups that are enriched in the down COGs
hypergeo = function(cog_model_obj, direction, output_dir, genomes=NULL) {
  # subset the matrix by the given genomes
  if ( ! is.null(genomes)) {
    cog_model_obj = keep_genomes(cog_model_obj, genomes_to_keep = genomes)
  }
  
  # get a vector of all possible COG groups from the cog annotation file
  group_lengths = sapply(levels(cog_model_obj$cog_annot_tb$func), nchar)
  cog_group_names = names(which(group_lengths == 1))
  
  # create a matrix to store which COG groups are enriched in the DA COGs
  x_mat = matrix(nrow = length(cog_group_names),
                 ncol = ncol(cog_model_obj$mat), 
                 dimnames = list(cog_group_names, colnames(cog_model_obj$mat)))
  m_mat = matrix(nrow = length(cog_group_names),
                 ncol = ncol(cog_model_obj$mat), 
                 dimnames = list(cog_group_names, colnames(cog_model_obj$mat)))
  
  # set the direction vector
  direction_vec = get_direction_vec(direction)
  
  # x = number of times a COG group is DA in a genome for a given direction
  # m = number of times a COG group is observed in a genome
  # n = number of times all COGs groups besides those counted in "m"
  # k = number of DA COGs in a genome
  
  # get x matrix
  x_mat = get_x_mat(cog_model_obj$mat, x_mat, 
                    cog_model_obj$cog_annot_tbl, direction_vec,
                    cog_group_names)
  
  # get m matrix
  m_mat = get_x_mat(cog_model_obj$mat, m_mat, 
                    cog_model_obj$cog_annot_tbl, c(-1, 0, 1),
                    cog_group_names)
  
  # NOTE: when I run the hypergeometric test I can caluculate n and k
  # using the x and m matricies
  
  # run hyper genometric test for each COG group
  hyper_mat = run_hypergeo(x_mat, m_mat)
  
  # get the binary values (ie convert into significant and not significant classes)
  bin_hyper_mat = convert_to_binary(hyper_mat, 0.05)
  
  # now that all the elements of the object are calculated make the actual object
  hypergeo_mat_obj = list(binary_mat = NULL, pval_mat = NULL, 
                          m_mat = NULL, x_mat = NULL)
  
  hypergeo_mat_obj$binary_mat = bin_hyper_mat
  hypergeo_mat_obj$pval_mat = hyper_mat
  hypergeo_mat_obj$x_mat = x_mat
  hypergeo_mat_obj$m_mat = m_mat
  
  class(hypergeo_mat_obj) = "hypergeo_mat_obj"
  
  return(hypergeo_mat_obj)
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



