# CPM Note:
# edgeR uses the lib.size varaible when calculating cpm.  In these data
# 1 read is ~0.01; 18 reads is ~ 0.02; 100 reads is ~0.4

# DOCUMENTATION #
# These functions are for doing differential abundance testings on the 
# output from DAFE_count.pl.

### GLOBALS ###
# these are the coding for DA values
ENR = 1   # Enriched
NS = 0    # Not significantly different
DEP = -1  # Depleted
LOW = -2  # Low abundance (too low to run stats on)
UMS = -3  # Unmeasurable -- there are no reads that map at all
ABS = -4  # Abset from the genome (probably not used in this script)

###

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

	# create the output file prefix
	outfile_prefix = get_out_file_prefix(params_obj)
    
    # read in the counts -- features are row and samples are cols
    tbl = get_count_tbl(ref, params_obj)

    # create the empty output vector -- set all the value to UMS (unmeasurable)
    da_vec = as.vector(rep(UMS,nrow(tbl)))
    names(da_vec) = row.names(tbl)

	# filter the count table
	filt_res = filter(tbl = tbl,
					metadata = params_obj$metaG_metadata_tbl, 
					test_col = params_obj$test_col_name,
					t1 = params_obj$test[1],
					t2 = params_obj$test[2],
					min_sample = params_obj$min_sample_count,
					min_count = params_obj$min_sample_cpm)
	tbl = filt_res$tbl

	# any features (eg COGs) that are not in the table at all
	# designate as UMS.  Use LOW to indicate features that we 
	# can measure but are not statistically testable due to 
	# low reads counts
	da_vec[filt_res$low_names] = LOW
	da_vec[filt_res$ums_or_abs_names] = UMS  # set these as UMS for now. handled later in Decouple.pm
	
	# if the table is empty then everything was removed
	# in this case I want to write out some files and move
	# to the next genome
	if ( nrow(tbl) < 1 ) {
		print("FINISH EARLY -- Didn't pass filter")

		# print the da direction file (they will all be -2)
		out_f = paste(outfile_prefix, "_da_vec.txt", sep="")
		write.table(da_vec, out_f, quote=F, col.names=F)
		
		# move to the next genome
		next
	}
    
    # create the DGEList object (exact test)
    dge = DGEList(counts=tbl, group=params_obj$test_groups)  
    
    # set the library size for normalization
    dge[["samples"]][["lib.size"]] = params_obj$metaG_metadata_tbl[row.names(dge[["samples"]]), "reads"]
    
    # calculte the normalization factors
    dge = calcNormFactors(dge)
    
    # estimate dispersion parameters
    out = tryCatch({
      dge = estimateCommonDisp(dge)
      #dge = estimateTrendedDisp(dge)
      dge = estimateTagwiseDisp(dge)
      TRUE  # return TRUE.  This is to avoid an unnecessary warning
    }, error = function(err){
      print(paste(err$message))
      return(NA)
    })
    if (is.na(out)) {
      # dispersion parameters could not be calculated for this genome
      print("FINISH EARLY -- Cannot estimate dispersion")
      
      # print the da direction file (they will all be UMS or LOW (unmeasurable or low abundance))
      out_f = paste(outfile_prefix, "_da_vec.txt", sep="")
      write.table(da_vec, out_f, quote=F, col.names=F)
      
      # move on to the next genome
      next
    }
    
    # make the BCV plot
    out_f = paste(outfile_prefix, "_BCV.tiff", sep="")
    tiff(out_f)
    plotBCV(dge)
    dev.off()
    
    # run the exact binomial test
    exact = exactTest(dge, pair=params_obj$test)
    
    # get table of all tags.
    tags = topTags(exact, n=nrow(tbl))
    
    # save the tables
    out_f = paste(outfile_prefix, "_tags_data.txt", sep="")
    write.table(tags, out_f, quote=F)
    
    # summarize differential abundace directionality
    out_f = paste(outfile_prefix, "_da_tbl.txt", sep="")
    da = summerize_da_direction(exact, output_file = out_f)
    
    # write the da_vec file
    out_f = paste(outfile_prefix, "_da_vec.txt", sep="")
    da_vec = get_da_vec(tags, da_vec)
    write.table(da_vec, out_f, quote=F, col.names=F)
    
    # Make smear plot
    #out_f = paste(outfile_prefix, "_smear.tiff", sep="")
    #make_smear_plot(exact, da, dge, output_file = out_f)
  }
}



### HELPER FUNCTIONs ###

## FILTER ##
# removes features (e.g. COGs) that have fewer than do not have min_count
# reads in at least min_sample samples.
filter = function(tbl, metadata, test_col, t1, t2, min_sample, min_count) {
	# this function will need to be updated if I ever update the code
	# to use multiple test groups (eg col, cvi, and oy). Rigth now I have 
	# assumes that there will only be two test groups considered at one time.
	# This implies that I can't do things like interaction terms.  If I
	# ever stop using the exact test this function will have to change

	# the tbl contains features that are ENR, NS, DEP, LOW, UMS, or ABS.
	# it is easy to tell of something is ENR, NS, or DEP as those are set
	# when I run the edgeR statistics. Things that are removed from this
	# table in this function are things that are LOW, UMS, or ABS.  It is
	# also easy to tell what is LOW because it will have a few reads. I
	# use Nick's Decouple.pm library to differentiate between UMS and ABS
	# later on in the pipeline.


	# get the sample names for each group
	t1_names = row.names(metadata[metadata[,test_col] == t1,])
	t2_names = row.names(metadata[metadata[,test_col] == t2,])

	# make sure the sample names are in the table for all these test
	t1_names = t1_names[t1_names %in% colnames(tbl)]
	t2_names = t2_names[t2_names %in% colnames(tbl)]

	# check for the feature in each group
	t1_keep = apply(tbl[,t1_names], 1, function(x) sum(x > min_count)) > min_sample
	t2_keep = apply(tbl[,t2_names], 1, function(x) sum(x > min_count)) > min_sample

	# combine the keeps
	keep = t1_keep & t2_keep

	# get a list of names that are LOW
	has_counts = apply(tbl, 1, sum) > 0
	low_names = row.names(tbl[has_counts & !keep,])

	# get a list of names that are UMS or ABS -- things that have no counts
	ums_or_abs_names = row.names(tbl[!has_counts,])

	# filter the table by removing features that do not pass the filter
	# also return 
	res = list(
			tbl = tbl[keep,], 
			low_names = low_names,
			ums_or_abs_names = ums_or_abs_names
	);
	return(res)
}


# Makes a da_vec by combining the empty_da_vec
# with the calls from the tags object.  A vector
# of DA codes is returned (ie -3,-2,-1,0,1)
get_da_vec = function(tags, empty_da_vec) {
  # add a da column to the tags table
  da_col = apply(tags$table, 1, function(x) get_da_code(x[1], x[4]))
  
  # merge the da column in tags with empty_da_vec
  da_vec = empty_da_vec
  da_vec[row.names(tags)] = da_col
  
  return(da_vec)
}

get_da_code = function(logFC, fdr) {
  da_code = UMS  # set as default to unmeasurable
  
  if ( fdr < 0.05 ) {
    if ( logFC > 0 ) {
      da_code = ENR # enriched
    } else {
      da_code = DEP # depleted
    }
  } else {
    da_code = NS # not significant
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
  tryCatch({
  	datags_exact = rownames(dge)[as.logical(da_exact)]
  	if ( write_bool == T ) { tiff(output_file) }
  	plotSmear(exact_test, de.tags=datags_exact)
  	if ( write_bool == T ) { dev.off() }
  }, error = function(err) {
	print("Cannot make smear plot")
	print(paste(err$message))
	return(NA)
  })
}

### Gets the count table from each group 
# only keeps the columns that are associated with the samples in 
# the metaG_metadata_tbl
get_count_tbl = function(count_file, params_obj) {
  # NOTE: the working dir is where the count files are already located.

  file = paste(params_obj$agg_count_file_prefix, ".txt", sep="")
    
  # read in the data
  tbl = read.table(file, header=T, row.names=1, sep="\t", check.names=F, na.strings = "")
  
  # keep only the columns in that correspond to the metaG_include IDs
  tbl = tbl[,row.names(params_obj$metaG_metadata_tbl)]
  
  return(tbl)
}

### Creates the output file prefix
# this prefix is used to create all the edgeR output files
# it is the aggregated count file prefix and the test information
# combined. This allows for multiple tests to be ran on the same
# aggregated count file without overwriting the last test's output
get_out_file_prefix = function(params_obj) {
      out_f = paste(
		params_obj$agg_count_file_prefix, 
		"_",
		params_obj$test[1],
		"_v_",
		params_obj$test[2],
		sep="")
	
	return(out_f)
}

