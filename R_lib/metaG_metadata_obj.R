

create_metaG_metadata_obj = function(metadata_file, include_file, exclude_file) {
  # NOTE: this is actually returned as a dataframe.  Most R object are returned
  # as a list.
  
  tbl = read.table(metadata_file, header=T, sep="\t")
  
  include = read_names(include_file) # read_names() in is params_obj.R
  tbl = tbl[tbl$sample_name %in% include,]
  
  exclude = read_names(exclude_file)
  tbl = tbl[!(tbl$sample_name %in% exclude),]
  
  tbl = set_names(tbl)
  
  return(tbl)
}

set_names = function(tbl) {
  # combine tech reps -- I need to finish!
  #tbl = combine_tech_reps(tbl)
  
  row.names(tbl) = tbl$sample_name
  
  return(tbl)
}

# I still need to finish this
# but I'm not actually using it right now so I can do it later.
combined_tech_reps = function(tbl) {
  tech_rep_names = tbl[grep("_r\\d+", tbl$sample_name, perl = T), "sample_name"]
  
  for (t in tech_rep_names) {
    
  }
}