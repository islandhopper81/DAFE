# this is an alternative algorithm than the one in the run_edgeR_grps
# function above.  It runs the statistical tests on a single large
# table where all the counts have been combined as opposed to doing
# the tests on each genome one at a time.

# NOTE that the params object features are interpreted a little differently
# in this function than in the run_edgeR_grps function.

# count_file_name => the full path to the single count file

run_edgeR_grps_cbn = function(params_obj) {
		
	require(edgeR)
	main_dir = getwd()

	# read in the table
	tbl = read.table(params_obj$count_file_name, header=T, row.names=1,
						sep="\t", check.names=F, na.strings="")

	# keep only the columns that correspond to the metaG_include IDs
	tbl = tbl[,row.names(params_obj$metaG_metadata_tbl)]

	# create the DGEList object (exact test)
	dge = DGEList(counts=tbl, group=params_obj$test_groups)

	# set the library size for normalization
	dge[["samples"]][["lib.size"]] = params_obj$metaG_metadata_tbl[row.names(dge[["samples"]]), "reads"]

	# calculte the normalization factors
	dge = calcNormFactors(dge)

	# filter -- by min_sample_count and min_sample_cpm
	keep = rowSums(cpm(dge) > params_obj$min_sample_cpm) > params_obj$min_sample_count
	if ( length(keep) == 0 | all(is.na(keep)) | sum(keep) == 0) {
	  	# no genes in this genome pass the filter
	  	print("FINISH EARLY -- Didn't pass filter")
		quit()
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
	if (is.na(out) ) {
		# dispersion parameters could not be calculated for this genome
		print("FINISH EARLY -- Cannot estimate dispersion")
		quit()
	}

	# run the exact binomial test
	exact = exactTest(dge.filt, pair=params_obj$test)

	# get the table of all tags
	tags = topTags(exact, n=nrow(tbl))

	# save the table
	out_f = "tags.txt"
	write.table(tags, out_f, quote=F, sep="\t")

	# save only the cases that are significant
	sig_tags = topTags(exact, p.value=0.05, n=nrow(tbl))
	write.table(sig_tags, file="tags_sig.txt", sep="\t", quote=F)
}
