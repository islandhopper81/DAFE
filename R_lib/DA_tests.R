# CPM Note:
# edgeR uses the lib.size varaible when calculating cpm.  In these data
# 1 read is ~0.01; 18 reads is ~ 0.02; 100 reads is ~0.4

# DOCUMENTATION #
# These functions are for doing differential abundance testings on the 
# output from DAFE_count.pl.

run_edgeR_grps = function(params_obj) {
  require(edgeR)
  main_dir = getwd()
  #print(paste("main_dir: ", main_dir))
  
  for ( ref in params_obj$ref_include_ids ) {
    print(ref)
    
    # the input and output dirs are the same
    # set them as working_dir
    setwd(main_dir)
    working_dir = paste(params_obj$working_root_dir, ref, sep="/")
    setwd(working_dir)
    #print(paste("getwd(): ", getwd()))
    
    # read in the counts
    tbl = get_count_tbl(ref, params_obj)

    # create the empty output vector
    da_vec = as.vector(rep(-2,nrow(tbl)))
    names(da_vec) = row.names(tbl)
    
    # create the DGEList object (exact test)
    dge = DGEList(counts=tbl, group=params_obj$test_groups)  
    
    # set the library size for normalization
    dge[["samples"]][["lib.size"]] = params_obj$metaG_metadata_tbl[row.names(dge[["samples"]]), "reads"]
    
    # calculte the normalization factors
    dge = calcNormFactors(dge)
    
    # filter -- by min_sample_count and min_sample_cpm
    keep = rowSums(cpm(dge) > params_obj$min_sample_cpm) > params_obj$min_sample_count
    if ( length(keep) == 0 | all(is.na(keep)) ) {
      # no genes in this genome pass the filter
      print("FINISH EARLY -- Didn't pass filter")
      
      # print the da direction file (they will all be -2)
      out_f = paste(params_obj$agg_count_file_prefix, "_da_vec.txt", sep="")
      write.table(da_vec, out_f, quote=F, col.names=F)
      
      # move on to the next genome
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
      # dispersion parameters could not be calculated for this genome
      print("FINISH EARLY -- Cannot estimate dispersion")
      
      # print the da direction file (they will all be -2)
      out_f = paste(params_obj$agg_count_file_prefix, "_da_vec.txt", sep="")
      write.table(da_vec, out_f, quote=F, col.names=F)
      
      # move on to the next genome
      next
    }
    
    # make the BCV plot
    out_f = paste(params_obj$agg_count_file_prefix, "_BCV.tiff", sep="")
    tiff(out_f)
    plotBCV(dge.filt)
    dev.off()
    
    # run the exact binomial test
    exact = exactTest(dge.filt, pair=params_obj$test)
    
    # get table of all tags.
    tags = topTags(exact, n=nrow(tbl))
    
    # save the tables
    out_f = paste(params_obj$agg_count_file_prefix, "_tags_data.txt", sep="")
    write.table(tags, out_f, quote=F)
    
    # summarize differential abundace directionality
    out_f = paste(params_obj$agg_count_file_prefix, "_da_tbl.txt", sep="")
    da = summerize_da_direction(exact, output_file = out_f)
    
    # write the da_vec file
    out_f = paste(params_obj$agg_count_file_prefix, "_da_vec.txt", sep="")
    da_vec = get_da_vec(tags, da_vec)
    write.table(da_vec, out_f, quote=F, col.names=F)
    
    # Make smear plot
    out_f = paste(params_obj$agg_count_file_prefix, "_smear.tiff", sep="")
    make_smear_plot(exact, da, dge.filt, output_file = out_f)
  }
}



### HELPER FUNCTIONs ###

# Makes a da_vec by combining the empty_da_vec
# with the calls from the tags object.  A vector
# of DA codes is returned (ie -2,-1,0,1)
get_da_vec = function(tags, empty_da_vec) {
  # add a da column to the tags table
  da_col = apply(tags$table, 1, function(x) get_da_code(x[1], x[4]))
  
  # merge the da column in tags with empty_da_vec
  da_vec = empty_da_vec
  da_vec[row.names(tags)] = da_col
  
  return(da_vec)
}

get_da_code = function(logFC, fdr) {
  da_code = -2
  
  if ( fdr < 0.05 ) {
    if ( logFC > 0 ) {
      da_code = 1
    } else {
      da_code = -1
    }
  } else {
    da_code = 0
  }
  
  return(da_code)
}

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

summerize_da_direction = function(exact_test, da_vec, 
                                  output_file = "da_tbl.txt", 
                                  write_bool = T) {
  # multiple testing adjustment
  da.exact = decideTestsDGE(exact_test, p=0.05, adjust="BH")
  
  row.names(da.exact) = row.names(exact_test)
  da.direction = summary(da.exact)
  
  if (write_bool == T) {
    write.table(da.direction, output_file, quote = F, col.names = F)
  }
  
  return(da.exact)
}

### Make a smear plot
# makes a smear plot from the edgeR output
# edgeR module must be sourced to use the function plotSmear
make_smear_plot = function(exact_test, da_exact, dge, output_file = "smear.tiff",
                           write_bool = T) {
  datags_exact = rownames(dge)[as.logical(da_exact)]
  if ( write_bool == T ) { tiff(output_file) }
  plotSmear(exact_test, de.tags=datags_exact)
  if ( write_bool == T ) { dev.off() }
}

### Gets the count table from each group 
# only keeps the columns that are associated with the samples in 
# the metaG_metadata_tbl
get_count_tbl = function(count_file, params_obj) {
  # NOTE: the working dir is where the count files are already located.

  file = paste(params_obj$agg_count_file_prefix, ".txt", sep="")
    
  # read in the data
  tbl = read.table(file, header=T, row.names=1, sep="\t", check.names=F)
  
  # keep only the columns in that correspond to the metaG_include IDs
  tbl = tbl[,row.names(params_obj$metaG_metadata_tbl)]
  
  return(tbl)
}



