### DEPRECIATED
# The main functionality of this method was to combine the data from all three
# summary files (EC, RZ, BK) into a single object.  These are now manually combined
# in the metadata file.  The second main purpose of this function was to combined 
# the summary numbers for the technical reps in the FACS data.  I am no longer
# using the FACS data.  If I do end up using it in the future or have other technical
# reps I should update this function and use it to combine the summary numbers for
# the technical reps.  This should be done using the information in the metadata
# file.  -sy 20151121


### get_normalization_numbers
# I am normalizing by the number of non_contaminant reads that are in each sample
# I don't think I can normalize based on the number of reads that map from each sample
get_normalization_numbers = function(metaG_meta_file, metaG_include_ids) {
  
  # read in the table
  tbl = read.table(metaG_meta_file, header=T, row.names="sample_names", sep="\t")
  
  # keep only the necessary samples
  tbl = tbl[metaG_include_ids,]
  
  facs_summ = read.table(facs_summary_file, header=T, sep="\t")
  facs_summ = facs_summ[which(facs_summ$type == "non_contaminants"),]  # don't consider contaminant seqs
  facs_summ = facs_summ[which(facs_summ$sort_method != "No_Sort"),] # don't consider the no sort samples
  facs_summ_tot = facs_summ$reads
  names(facs_summ_tot) = paste(facs_summ$geno, facs_summ$soil, facs_summ$sort_method, facs_summ$rep, sep="_")
  
  # read in the summary table for RZ to get the total reads counts for normalization
  rz_summ = read.table(rz_summary_file, header=T, sep="\t")
  rz_summ_tot = rz_summ$reads
  names(rz_summ_tot) = rz_summ$sample_name
  
  # read in the summary table for soil to get the total reads counts for normalization
  soil_summary_file = "/Users/Scott/Dangl/comparitive_metaG/soil_metadata.txt"
  soil_summ = read.table(soil_summary_file, header=T, sep="\t")
  soil_summ_tot = soil_summ$reads
  names(soil_summ_tot) = soil_summ$sample_name
  
  # combine the summary tables
  full_summ = c(facs_summ_tot, rz_summ_tot, soil_summ_tot)
  full_summ['Col_CL_FACS'] = full_summ['Col_CL_FACS_r1'] + full_summ['Col_CL_FACS_r2']
  full_summ['Col_MF_FACS'] = full_summ['Col_MF_FACS_r1'] + full_summ['Col_MF_FACS_r2']
  full_summ['Cvi_CL_FACS'] = full_summ['Cvi_CL_FACS_r1'] + full_summ['Cvi_CL_FACS_r2']
  full_summ['Cvi_MF_FACS'] = full_summ['Cvi_MF_FACS_r1'] + full_summ['Cvi_MF_FACS_r2']
  full_summ = full_summ[!(names(full_summ) %in% c("Col_CL_FACS_r1",
                                                  "Col_CL_FACS_r2",
                                                  "Col_MF_FACS_r1",
                                                  "Col_MF_FACS_r2",
                                                  "Cvi_CL_FACS_r1",
                                                  "Cvi_CL_FACS_r2",
                                                  "Cvi_MF_FACS_r1",
                                                  "Cvi_MF_FACS_r2"))]
  
  # return all the summ data in a dataframe
  ret.list = list("full_summ" = full_summ, 
                  "facs_summ" = facs_summ,
                  "rz_summ" = rz_summ,
                  "soil_summ" = soil_summ)
  
  return(ret.list)
}

