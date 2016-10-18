# CPM Note:
# edgeR uses the lib.size varaible when calculating cpm.  In these data
# 1 read is ~0.01; 18 reads is ~ 0.02; 100 reads is ~0.4

# DOCUMENTATION #
# These functions are for doing differential abundance testings on the 
# output from DAFE_count.pl.


run_edgeR_grps = function(grp_model_obj, params_obj) {
  # NOTE: there are global variable that are used in this function
  
  # get the empty matricies that will be populated
  mat = grp_model_obj$mat
  gpg_mat = grp_model_obj$gpg_mat
  
  # a vector of any genomes that should be removed from mat and gpg_mat at the end
  genomes_to_remove = vector()
  
  ### the exact test for only the RZ and BK subset
  for (ref in params_obj$ref_include_ids) {
    print(ref)
    
    # set and create the output directory
    ref_out_path = paste(params_obj$out_root_dir, params_obj$run_name,
                         params_obj$ref_out_dir, ref, sep="/")
    dir.create(path = ref_out_path, recursive = T, showWarnings = F)
    setwd(ref_out_path)
    
    # read in the file
    tbl = get_count_tbl(ref, params_obj)
    
    # remove the features that start with "__"
    to_remove = c("__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned", "__alignment_not_unique")
    tbl = tbl[!(row.names(tbl) %in% to_remove),]
    
    # combine the technical replicates 
    # (neccessary for differential expression; we only care about biological variation)
    # I should use this function when I have tech reps
    #tbl = combine_technical_reps(tbl)
    
    # convert genes count table to higher groups (ie COG, KO, custers, etc.)
    skip_ref=F
    out = tryCatch({
      annote_file = paste(params_obj$dafe_dir, ref, params_obj$annote_file_name, sep="/")
      result = genes_to_grp(annote_file, tbl, params_obj$group_genes_by, params_obj$gene_id_col)
      tbl = result$grp_tbl
      gpg_tbl = result$gene_counts
    }, error = function(err) {
      print(paste("ERROR: ", err$message, sep=""))
      traceback(2)
      skip_ref=T
    })
    if ( skip_ref == T ) {
      print("SKIPPING -- Cannot calculate group table")
      genomes_to_remove = c(genomes_to_remove, ref)
      next
    }
    
    # create the DGEList object
    dge = DGEList(counts=tbl, group=grp_model_obj$test_groups)  # exact test with count table
    
    # set the library size for normalization
    dge$samples$lib.size = params_obj$metaG_metadata_obj[row.names(dge$samples), "reads"]
    
    # calculte the normalization factors
    dge = calcNormFactors(dge)
    
    # filter -- must have min_sample_count samples with greater than min_sample_cpm
    keep = rowSums(cpm(dge) > params_obj$min_sample_cpm) > params_obj$min_sample_count
    if ( length(keep) == 0 | all(is.na(keep)) ) {
      # no genes pass the filter
      print("FINISH EARLY -- Didn't pass filter")
      
      # finish early because this genome didn't pass the filter
      if ( params_obj$skip_big_mat == FALSE ) {
        mat[,as.character(ref)] = make_grp_vector(grp_model_obj$grp_metadata_tbl, data.frame(), -2, 
                                                  paste(ref_out_path, "da_grps.txt", sep="/"))
        gpg_mat[,as.character(ref)] = make_grp_vector(grp_model_obj$grp_metadata_tbl, data.frame(), 0, 
                                                      paste(ref_out_path, "genes_per_grp.txt", sep="/"))
      }
      next
    } else {
      dge.filt = dge[keep,]
    }
    
    # estimate dispersion parameters
    out = tryCatch({
      dge.filt = estimateCommonDisp(dge.filt)
      #dge.filt = estimateTrendedDisp(dge.filt)
      dge.filt = estimateTagwiseDisp(dge.filt)
      TRUE  # return TRUE.  This is to avoid an unnecessary warning
    }, error = function(err){
      print(paste(err$message))
      return(NA)
    })
    if (is.na(out)) {
      print("FINISH EARLY -- Cannot estimate dispersion")
      
      # finish early because this genome the dispersion could not be estimated
      if ( params_obj$skip_big_mat == FALSE ) {
        mat[,as.character(ref)] = make_grp_vector(grp_model_obj$grp_metadata_tbl, data.frame(), -2, 
                                                  paste(ref_out_path, "da_grps.txt", sep="/"))
        gpg_mat[,as.character(ref)] = make_grp_vector(grp_model_obj$grp_metadata_tbl, data.frame(), 0, 
                                                      paste(ref_out_path, "genes_per_grp.txt", sep="/"))
      }
      next
    }
    
    # make the BCV plot
    tiff(paste(ref_out_path, "BCV.tiff", sep="/"))
    plotBCV(dge.filt)
    dev.off()
    
    # run the exact binomial test
    exact = exactTest(dge.filt, pair=grp_model_obj$test)
    
    # get table of all tags.
    tags = topTags(exact, n=nrow(tbl))
    
    # save the tables
    write.table(tags, paste(ref_out_path, "tags_data.txt", sep="/"), quote=F)
    
    # summarize differential abundace directionality
    da = summerize_da_direction(exact, ref_out_path)
    
    # Make smear plot
    make_smear_plot(exact, da, dge.filt, ref_out_path)
    
    # some debugging code
#     {
#       print("===== Debug block 1 =====")
#       print(paste("da type: ", typeof(da)))
#       print(paste("da size: ", object_size(da)))
#       print(paste("mat type: ", typeof(mat)))
#       print(paste("mat size: ", object_size(mat)))
#       print(paste("mem used: ", mem_used()))
#       gc()
#       print(paste("mem used after gc: ", mem_used()))
#     }
    
    # everything after this point is related to adding the DA values to the large matrix
    # So I can skip everything if the skip_big_mat value is FALSE
    if ( params_obj$skip_big_mat == TRUE ) {
      next
    }
    
    # each genome should have a vector where rows are groups and 
    # values are a factor (up, down, nd, np)
    # where nd == not different and np == not present
    # add to the mat matrix
    mat[,as.character(ref)] = make_grp_vector(grp_model_obj$grp_metadata_tbl, da, -2, 
                                              paste(ref_out_path, "da_grps.txt", sep="/"))
  
    # debugging code
    {
#       print("===== Debug block 2 =====")
#       print(paste("mat type: ", typeof(mat)))
#       print(paste("mat size: ", object_size(mat)))
#       print(paste("mem used: ", mem_used()))
#       gc()
#       print(paste("mem used after gc: ", mem_used()))
    }
    
    # doing the same opertion for the genes per group matrix
    # add to the gpg matrix
    gpg_mat[,as.character(ref)] = make_grp_vector(grp_model_obj$grp_metadata_tbl, gpg_tbl, 0, 
                                                  paste(ref_out_path, "genes_per_grp.txt", sep="/"))
  }
  
  # when I get to this point I have finished with all the tests
  # for the individual references (ie genomes).  Now I just have
  # to do some clean up stuff
  
  # Note that when params_obj$skip_big_mat is set to TRUE the grp_model_obj
  # values are basically meaningless.  I added the skip_big_mat parameter to
  # allow for an analysis involving many rows in the matrix.  I found that the
  # matrix quickly sucked up more memory than was available.  So the solution
  # is to calculate the final matrix using a perl script later (ie after this one
  # completes) using the tags_data.txt file
  if ( params_obj$skip_big_mat == FALSE ) {
    # set the row names in the matricies
    row.names(mat) = row.names(grp_model_obj$grp_metadata_tbl)
    row.names(gpg_mat) = row.names(grp_model_obj$grp_metadata_tbl)
    
    # reset the grp_model_obj values
    grp_model_obj$mat = mat
    grp_model_obj$gpg_mat = gpg_mat
    
    # remove any genomes that were specified to remove in the loop
    grp_model_obj = remove_genomes(grp_model_obj, genomes_to_remove)
    
    # clean -- groups that are always undetectable
    # NOTE this step must come before calling classify_non_detect_grps
    # the clean_grps function is defined in the grp_model_obj
    grp_model_obj = clean_grps(grp_model_obj)
    
    # specify distinction between absent and non-detectable groups
    # the classify_non_detect_grps function is defined in the grp_model_obj
    grp_model_obj = classify_non_detect_grps(grp_model_obj)
  }
    
  return(grp_model_obj)
}


### HELPER FUNCTIONs ###



### combine_technical_reps
# merges the numbers for the techincal reps
combine_technical_reps = function(tbl) {
  # combine (ie add together) the technical reps for the FACS samples
  tbl$Col_CL_FACS = tbl$Col_CL_FACS + tbl$Col_CL_FACS_r2
  tbl$Col_MF_FACS = tbl$Col_MF_FACS + tbl$Col_MF_FACS_r2
  tbl$Cvi_MF_FACS = tbl$Cvi_MF_FACS + tbl$Cvi_MF_FACS_r2
  tbl$Cvi_CL_FACS = tbl$Cvi_CL_FACS + tbl$Cvi_CL_FACS_r2
  tbl = tbl[,!(names(tbl) %in% c("Col_CL_FACS_r2",
                                 "Col_MF_FACS_r2",
                                 "Cvi_CL_FACS_r2",
                                 "Cvi_MF_FACS_r2"))]
  
  return(tbl)
}

### Convert a gene table to a group table and get corresponding gene counts per group table
genes_to_grp = function(annote_file, gene_tbl, group_genes_by, gene_id_col) {
  # read in the annotaiton table file
  annote_tbl = read.table(annote_file, header=T, sep="\t", quote="", comment.char="")
  
  # make sure group_genes_by is a column in the annote_tbl
  if ( !(group_genes_by %in% colnames(annote_tbl)) ) {
    msg = paste("Group by column (", group_genes_by, ") ",
                "not a header in the annotation file: ", annote_file, sep="")
    stop(msg)
  }
  
  # subset the table to get only the gene_names and grouping column
  annote_tbl = annote_tbl[,c(gene_id_col, group_genes_by)]
  
  # remove any NA values in the group_genes_by column
  na_entries = is.na(annote_tbl[,c(group_genes_by)])
  annote_tbl = annote_tbl[!na_entries,]
  
  # remove duplicates -- see notes above remove_duplicate_genes function
  annote_tbl = remove_duplicate_genes(annote_tbl)
  
  # get the gene counts in each group (a vector)
  gene_counts = apply(table(annote_tbl), 2, sum)
  
  # merge the tables
  tbl_merge = merge(x=gene_tbl, y=annote_tbl, by.x = "row.names", by.y = gene_id_col)
  
  # there are some genes that are removed from the tbl_merge table because they
  # were not assigned a cluster.  These need to be removed from the gene_tbl
  gene_tbl = gene_tbl[row.names(gene_tbl) %in% annote_tbl[,c(gene_id_col)],]
  
  # aggregate
  grp_tbl = aggregate(gene_tbl, by=list(tbl_merge[,c(group_genes_by)]), FUN=sum, na.rm=T)
  
  row.names(grp_tbl) = grp_tbl$Group.1
  grp_tbl = grp_tbl[,!names(grp_tbl) %in% c("Group.1")]
  
  counts = as.data.frame(gene_counts)
  row.names(counts) = names(gene_counts)
  
  # both are returned as data frames
  return(list(grp_tbl = grp_tbl, gene_counts = counts))
}

# remove duplicated genes
# There are some genes that are duplicated in the annotation
# file. This is because they are found multiple times in the
# genbank file.  I only looked at one case where this happened
# and it was because one genes was completely (and exactly, 
# I think) encased in the other other (kind of like having
# multiple start sites.)
# NOTE: This function simply removes the second instance of 
# the duplicates.  In the single case that I looked at the
# second gene was the smaller one.
remove_duplicate_genes = function (annotation_tbl) {
  duplicates = duplicated(annotation_tbl$gene_name)
  
  if ( sum(duplicates) < 1) {
    # there are no duplicates
    return(annotation_tbl)
  }
  
  return(annotation_tbl[(!duplicates),])
}

### get a vector of with values for all groups where the NA 
# ones are replaced by empty_val
make_grp_vector = function(grp_metadata_tbl, da_tbl, empty_val="-2",
                           out_file="da_grps.txt", write_bool=T) {
  # grp_metadata_tbl => the 1st column in this tbl gives me all possible clusters.
  # So if a gene does not have a cluster I can put the empty_val in that cell.
  
  if ( typeof(da_tbl) == "list" ) {
    # convert it to a matrix
    da_tbl = data.matrix(da_tbl, rownames.force=TRUE)
  }
  
  save(da_tbl, file="da_tbl.RData")
  save(grp_metadata_tbl, file="grp_metadata_tbl.RData")

  if ( nrow(da_tbl) == 0 ) {
    # this is the case where there are no DA grps.
    merged_vec = rep(NA, nrow(grp_metadata_tbl))
  }
  else {
    # this is the standard case
    merged_tbl = merge(grp_metadata_tbl, da_tbl, by.x="row.names", by.y="row.names", 
                      all.x = T, sort=T)
    
    # this is not really a vector but it looks like one.
    merged_vec = merged_tbl[,ncol(merged_tbl)]
    #print(paste("merged_vec type: ", typeof(merged_vec)))
    
    merged_vec[is.na(merged_vec)] = empty_val
  }
  
  names(merged_vec) = row.names(grp_metadata_tbl)

  if (write_bool == T) {
    write.table(merged_vec, out_file, quote=F, col.names=F)
  }
  
  return(merged_vec)
}

summerize_da_direction = function(exact_test, output_dir, output_file = "da_direction.txt", 
                                  write_bool = T) {
  da.exact = decideTestsDGE(exact_test, p=0.05, adjust="BH")
  row.names(da.exact) = row.names(exact_test)
  da.direction = summary(da.exact)
  if (write_bool == T) {
    file = paste(output_dir, output_file, sep="/")
    write.table(da.direction, file, quote = F, col.names = F)
  }
  
  return(da.exact)
}

### used in edgeR_model.R
# makes a smear plot from the edgeR output
# edgeR module must be sourced to use the function plotSmear
make_smear_plot = function(exact_test, da_exact, dge, output_dir, output_file = "smear.tiff",
                           write_bool = T) {
  datags_exact = rownames(dge)[as.logical(da_exact)]
  if ( write_bool == T ) { tiff(paste(output_dir, output_file, sep="/")) }
  plotSmear(exact_test, de.tags=datags_exact)
  if ( write_bool == T ) { dev.off() }
}

### Gets the count table from each group 
# only keeps the columns that are associated with the samples in 
# the metaG_metadata_obj
get_count_tbl = function(ref, params_obj) {
  ## OLD
#   tables = list()
  #   for (g in unique(params_obj$metaG_metadata_obj$group)) {
#     dir = paste(params_obj$count_data_dir, g, ref, sep="/")
#     file = paste(dir, "/", params_obj$count_file_name, sep="")
#     
#     # read the data
#     data = read.table(file, header=T, row.names=1, sep="\t")
#     
#     # store the data in a list to be concatenated later
#     tables[[g]] = data
#   }
#   
#   # concatenate all the data frames
#   tbl = data.frame(row.names = row.names(tables[[1]]))
#   for (t in tables) {
#     tbl = cbind(tbl, t)
#   }
  ##
  
  dir = paste(params_obj$data_dir, ref, sep="/")
  file = paste(dir, "/", params_obj$count_file_name, sep="")
  
  # read in the data
  tbl = read.table(file, header=T, row.names=1, sep="\t")
  
  # keep only the columns in that correspond to the metaG_include IDs
  tbl = tbl[,row.names(params_obj$metaG_metadata_obj)]
  
  
  return(tbl)
}


# Function for building the model using genes (not groups)
run_edgeR_genes = function(params_obj) {
  # NOTE: there are global variable that are used in this function
  
  ### the exact test for only the RZ and BK subset
  for (ref in params_obj$ref_include_ids) {
    print(ref)
    
    # set and create the output directory
    ref_out_path = paste(params_obj$out_root_dir, params_obj$run_name, 
                        params_obj$ref_out_dir, ref, sep="/")
    dir.create(path = ref_out_path, recursive = T, showWarnings = F)
    setwd(ref_out_path)
    
    # read in the file
    tbl = get_count_tbl(ref, params_obj)
    
    # remove the features that start with "__"
    to_remove = c("__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned", "__alignment_not_unique")
    tbl = tbl[!(row.names(tbl) %in% to_remove),]
    
    # create the DGEList object
    # exact test with gene count table
    dge = DGEList(counts=tbl, group=params_obj$metaG_metadata_obj$fraction)
    
    # set the library size for normalization
    dge$samples$lib.size = params_obj$metaG_metadata_obj[row.names(dge$samples), "reads"]
    
    # calculte the normalization factors
    dge = calcNormFactors(dge)
    
    # filter -- right now this also only considers BK and EC samples
    keep = rowSums(cpm(dge) > params_obj$min_sample_cpm) > params_obj$min_sample_count
    if ( length(keep) == 0 | all(is.na(keep)) ) {
      # no genes pass the filter
      print("SKIPPING -- Didn't pass filter")
      next
    } else {
      dge.filt = dge[keep,]
    }
    
    # estimate dispersion parameters
    out = tryCatch({
      dge.filt = estimateCommonDisp(dge.filt)
      #dge.filt = estimateTrendedDisp(dge.filt)
      dge.filt = estimateTagwiseDisp(dge.filt)
    }, error = function(err){
      print(paste(err$message))
      dge.filt = NA
    })
    if (is.na(out)) {
      print("SKIPPING -- Cannot estimate dispersion")
      next
    }
    
    # make the BCV plot
    tiff(paste(ref_out_path, "BCV_genes.tiff", sep="/"))
    plotBCV(dge.filt)
    dev.off()
    
    # run the exact binomial test
    # remember the group order is 1=>soil, 2=>EC, 3=>RZ
    exact = exactTest(dge.filt, pair=params_obj$test)
    
    # get table of all tags.
    tags = topTags(exact, n=nrow(tbl))
    
    # save the tables
    write.table(tags, paste(ref_out_path, "tags_genes.RData", sep="/"), quote=F)
    
    # summarize differential abundace directionality
    da = summerize_da_direction(exact, ref_out_path, output_file = "da_direction_genes.txt")
    
    # Make smear plot
    make_smear_plot(exact, da, dge.filt, ref_out_path, output_file = "smear_genes.tiff")
    
    # write the da genes
    # NOTE: This is ONLY mostly complete. The problem is that there are some genes that 
    # are not tested.  Ie they don't show up in the htseq mapped counts.
    # I need to add the -2 value to represent these genes....maybe.
    write.table(da, file = paste(ref_out_path, "da_genes.txt", sep="/"), quote = F, col.names=F)
    write.table(tags, file = paste(ref_out_path, "tags_genes.txt", sep="/"), quote=F, col.names=T)
  }
  # right now there is nothing that is returned from this function
  # there are two output files that are made in the 
}

