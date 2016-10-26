
# file for running edgeR on the gene count tables of each genome

print("Loading libraries")
library(edgeR)
library(yaml)  # to load the input parameters

### ARGS
# params_file=
# source_dir=

print("Getting args")
args = commandArgs(trailingOnly=T)

for(i in 1:length(args)){
  eval(parse(text=args[[i]]))
}

### check params
if ( ! exists("params_file") ) {
  stop("Must provide a param.yaml file")
} else {
  print(paste("params_file: ", params_file))
}

if ( ! exists("source_dir") ) {
  source_dir = "R_lib"
  print(paste("setting source_dir=", source_dir))
} else {
  print(paste("source_dir: ", source_dir))
}
###


# source necessary files
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    if(trace) cat(nm,":")           
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}
sourceDir(source_dir)

# create the parameters object
params_obj = create_edgeR_params_obj(yaml_file = params_file)


### GRP ANALYSIS ###
# a grp is a feature e.g. KOG, COG, cluster, etc.
print("Run edgeR")
run_edgeR_grps(params_obj)

