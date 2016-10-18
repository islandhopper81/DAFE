

# CPM Note:
# edgeR uses the lib.size varaible when calculating cpm.  In these data
# 1 read is ~0.01; 18 reads is ~ 0.02; 100 reads is ~0.4

# Runs edgeR using CD-hit tables.  Think of this file as inheritting
# from DA_test.R.  There are several functions from that file that 
# are called here (get_count_tbl, )

run_edgeR_cdhit = function(clstr_model_obj, params_obj) {
  # NOTE: there are global variable that are used in this function
  
  # get the empty matricies that will be populated
  mat = clstr_model_obj$mat
  gpc_mat = clstr_model_obj$gpc_mat
  
  # a vector of any genomes that should be removed from mat and gpc_mat at the end
  genomes_to_remove = vector()
  
  ### the exact test for only the RZ and BK subset
  for (ref in params_obj$ref_include_ids) {
    print(ref)
    
    # set and create the output directory
    ref_out_dir = paste(params_obj$out_root_dir, params_obj$run_name, "ref_out", ref, sep="/")
    dir.create(path = ref_out_dir, recursive = T, showWarnings = F)
    setwd(ref_out_dir)
    
    # read in the file
    tbl = get_count_tbl(ref, params_obj)
    
    # remove the features that start with "__"
    to_remove = c("__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned", "__alignment_not_unique")
    tbl = tbl[!(row.names(tbl) %in% to_remove),]
    
    # combine the technical replicates 
    # (neccessary for differential expression; we only care about biological variation)
    # I should use this function when I have tech reps
    #tbl = combine_technical_reps(tbl)
    
    # convert any genes that have CD-Hit clusters into clusters
    # This give me a genes + CD-Hit clusters table
    # I leave the genes that are not in clusters in the table to help
    # edgeR estimate the dispersion.  I anticipate that there will be
    # very few clusters that have multiple genes (from the same genome).
    # That is important as the dispersion parameters will be estimated
    # from single genes.
    if (params_obj$use_clstr_tbl) { 
      skip_ref=F
      out = tryCatch({
        annote_file_path = paste(params_obj$dafe_dir, ref, "cdhit_annot.txt", sep="/")
        res = genes_to_clstrs(annote_file_path, tbl)
        tbl = res$clstr_tbl
        gpc_tbl = res$gene_counts
        clstr_ids = res$clstr_ids
      }, error = function(err) {
        print(paste(err$message))
        skip_ref=T
      })
      if ( skip_ref == T ) {
        print("SKIPPING -- Cannot calculate COG table")
        genomes_to_remove = c(genomes_to_remove, ref)
        next
      }
    }
    
    # create the DGEList object
    dge = DGEList(counts=tbl, group=clstr_model_obj$groups)  # exact test with gene count table
    
    # set the library size for normalization
    dge$samples$lib.size = params_obj$metaG_metadata_obj[row.names(dge$samples), "reads"]
    
    # calculte the normalization factors
    dge = calcNormFactors(dge)
    
    # filter -- right now this also only considers BK and EC samples
    keep = rowSums(cpm(dge) > params_obj$min_sample_cpm) > params_obj$min_sample_count
    if ( length(keep) == 0 | all(is.na(keep)) ) {
      # no genes pass the filter
      print("SKIPPING -- Didn't pass filter")
      
      # remember which genomes I need to remove for after the loop
      genomes_to_remove = c(genomes_to_remove, ref)
      next
    } else {
      dge.filt = dge[keep,]
    }
    
    # estimate dispersion parameters
    errFlag=F
    out = tryCatch({
      dge.filt = estimateCommonDisp(dge.filt)
      #dge.filt = estimateTrendedDisp(dge.filt)
      dge.filt = estimateTagwiseDisp(dge.filt)
      TRUE  # this is to avoid a warning
    }, error = function(err){
      print(paste(err$message))
      return(NA)
    })
    if (is.na(out)) {
      print("SKIPPING -- Cannot estimate dispersion")
      genomes_to_remove = c(genomes_to_remove, ref)
      next
    }
    
    # make the BCV plot
    tiff(paste(ref_out_dir, "BCV.tiff", sep="/"))
    plotBCV(dge.filt)
    dev.off()
    
    # run the exact binomial test
    # remember the group order is 1=>soil, 2=>EC, 3=>RZ
    exact = exactTest(dge.filt, pair=clstr_model_obj$test)
    
    # get table of all tags.
    tags = topTags(exact, n=nrow(exact$table))
    
    # save the tables
    write.table(tags, paste(ref_out_dir, "tags.RData", sep="/"), quote=F)
    
    # summarize differential abundace directionality
    da = summerize_da_direction(exact, ref_out_dir)
    
    # Make smear plot
    make_smear_plot(exact, da, dge.filt, ref_out_dir)
    
    # get only the resuts of the clusters
    # Note: there are some clusters that are filtered out.  When I get the 
    # clusters in da I have to get only the ones that were not filtered.
    print(da)
    da = da[row.names(as.data.frame(da)) %in% clstr_ids,]
    print(da)
    
    if (params_obj$use_clstr_tbl) {
      # each genome should have a vector where rows are cogs and 
      # values are a factor (up, down, nd, np)
      # where nd == not different and np == not present
      da_vect = make_clstr_vector(clstr_model_obj$clstr_annot_tbl, da, -2, 
                                paste(ref_out_dir, "da_clstr.txt", sep="/"))
      
      # add the da.cogs vector to the cog.matrix
      mat[,as.character(ref)] = da_vect
      
      # doing the same opertion for the genes per cog matrix
      gpc_vect = make_clstr_vector(clstr_model_obj$clstr_annot_tbl, gpc_tbl, 0, 
                                   paste(ref_out_dir, "genes_per_clstr.txt", sep="/"))
      
      # add the gpc_vect to the gpc matrix
      gpc_mat[,as.character(ref)] = gpc_vect
    }
    
  }
  
  # set the row names in the matricies
  row.names(mat) = row.names(clstr_model_obj$cog_annot_tbl)
  row.names(gpc_mat) = row.names(clstr_model_obj$cog_annot_tbl)
  
  # reset the clstr_model_obj values
  clstr_model_obj$mat = mat
  clstr_model_obj$gpc_mat = gpc_mat
  
  # remove any genomes that were specified to remove in the loop
  clstr_model_obj = remove_genomes(clstr_model_obj, genomes_to_remove)
  
  return(clstr_model_obj)
}

genes_to_clstrs = function(annot_file, gene_tbl) {
  # annot_file => file with annotations for each gene.  In this case
  #               it will have a line for any genes that are assigned
  #               a cluster and another column for which cluster they
  #               belong to.
  annote_tbl = read.table(annot_file, header=T, sep="\t", quote="", comment.char="")
  annote_tbl = annote_tbl[,c("gene_name", "clstr_id")]
  clstr_ids = annote_tbl$clstr_id
  
  # remove duplicates -- see notes above function
  annote_tbl = remove_duplicate_genes(annote_tbl)
  
  # get the gene counts in each cluster (a vector)
  gene_counts = apply(table(annote_tbl), 2, sum)
  
  # merge the tables
  tbl.merge = merge(x=gene_tbl, y=annote_tbl, by.x = "row.names", by.y = "gene_name")
  last_col_index = ncol(gene_tbl) + 1  # actually the last count column
  clstr_tbl = aggregate(tbl.merge[2:last_col_index], 
                        by=list(tbl.merge$clstr_id), 
                        FUN=sum, na.rm=T)
  row.names(clstr_tbl) = clstr_tbl$Group.1
  clstr_tbl = clstr_tbl[,!names(clstr_tbl) %in% c("Group.1")]
  
  # now remove the genes that have counts in the clstr_tbl from the gene_tbl
  new_gene_tbl = gene_tbl[!(row.names(gene_tbl) %in% annote_tbl$gene_name),]
  
  # now combine the gene_tbl and clstr_tbl
  clstr_tbl = rbind(new_gene_tbl, clstr_tbl)
  
  return(list(clstr_tbl = clstr_tbl, 
              gene_counts = gene_counts, 
              clstr_ids = clstr_ids)
         )
}

# remove duplicated genes
# There are some genes that are duplicated in the annotation
# file. This is because they are found multiple times in the
# genbank file.  I only looked at one case where this happend
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

### get a vector of with values for all clusters
# where the NA ones are replaced by empty_val
make_clstr_vector = function(clstr_annote_tbl, da_tbl, empty_val="-2", 
                             out_file="da_clstrs.txt", write_bool=T) {
  # annot_tbl => give me all possible clusters.  So if a gene does not
  # have a cluster I can put the empty_val in that cell.
  
  merged_tbl = merge(clstr_annote_tbl, da_tbl, by.x="row.names", 
                     by.y="row.names",all.x = T)
  
  # this is not really a vector but it looks like one.
  merged_vec = merged_tbl[,ncol(merged_tbl)
                          ]
  merged_vec[is.na(merged_vec)] = empty_val
  names(merged_vec) = row.names(clstr_annote_tbl)
  if (write_bool == T) {
    write.table(merged_vec, out_file, quote=F, col.names=F)
  }
  
  return(merged_vec)
}