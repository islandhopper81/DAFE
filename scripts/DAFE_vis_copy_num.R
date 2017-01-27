library(ggplot2)
library(reshape2)

# an analysis on KOG0001 to see if 
# more instanance with having KOG0001 correlate
# with being DA up.

# uncomment these for testing when I don't run via Rscript
opt = list()
opt$grp_name = "KOG0001"
opt$grp_file = "/Users/Scott/temp/KOG0001/KOG0001_da_calls_final.txt"
opt$copy_num_file = "/Users/Scott/temp/KOG0001/KOG0001_in_fungal_genomes_sorted.txt"

### Parameter variables
# The params matrix
# Each row is a parameter with 4 columns: 
#   1) long value, 
#   2) short value, 
#   3) argument type where 0 = no argument, 1 = required , 2 = optional
#   4) data type (logical, integer, double, complex, character)
# When you call this script the args are named like --verbose or -v
params = matrix(c(
  "da_file", "d", 1, "character",
  "copy_num_file", "c", 1, "character",
  "grp_name", "n", 1, "character"
), byrow=TRUE, ncol=4)
opt = getopt(params)

# read in the data
da_data = read.table(opt$grp_file, 
                     header=F, sep="\t", row.names=1)
copy_num_data = read.table(opt$copy_num_file,
                           header=F, sep="\t", row.names=1)

# for the rows that are in the da_data but not in the 
# copy_num_data that means they are absent from those
# genomes.  The copy number should be 0
copy_num_zero = row.names(da_data)[!row.names(da_data) %in% row.names(copy_num_data)]
to_add = data.frame("V2" = rep(0, length(copy_num_zero)), row.names = copy_num_zero)
copy_num_data = rbind(copy_num_data, to_add)



# merge the two tables.
# there are rows in the copy_num_data that are not in the 
# da data because those genomes were not able to be incorperated
# into the final tree used to run DAFE_count.  These extra rows
# can be ignored
merged = merge(da_data, copy_num_data, by.x = "row.names",
               by.y = "row.names")

# prepare the x-axis labels and ordering
my_labels = c("-3" = "Undetectable", "-2" = "Absent", "-1" = "DA Down", 
           "0" = "DA Equal", "1" = "DA up")
my_limits = c("-3", "-2", "-1", "0", "1")


# make the plot
ggplot(merged, aes(x=as.character(V2.x), y=V2.y, group=V2.x)) + 
  geom_boxplot(outlier.size = NA) +
  geom_point(position = position_jitter(width = 0.2)) + 
  scale_x_discrete("DA Category", labels = my_labels,
                   limits = my_limits) +
  ylab("Copy Number") + 
  ggtitle(paste(opt$grp_name, " Genomic Copy Number vs DA Category", sep="")) +
  theme(plot.title = element_text(hjust = .5, size=20),
        axis.text = element_text(size=16),
        axis.title = element_text(size=18))


# run a t-test to see if there is a difference between DA up (ie 1) and DA equal (ie 0)
t_test = t.test(merged[merged$V2.x == 0, "V2.y"], merged[merged$V2.x == 1, "V2.y"])



### make the plot for all the data (ie full DA table and copy num table)

# read in the copy num tbl
cpy_num = read.table("/Users/Scott/Projects/CompMetaG/experiments/nature_paper_44_metaG/fungal/copy_number/kog_count_tbl_rheat.txt",
                     header=T, sep="\t")

# read in the full DA tbl
da_tbl = read.table("/Users/Scott/Projects/CompMetaG/experiments/nature_paper_44_metaG/fungal/copy_number/full_da_tbl_29763.txt",
                    header=T, row.names=1, sep="\t", check.names=F)
da_tbl_melt = melt(as.matrix(da_tbl))
colnames(da_tbl_melt) = c("kogID", "genomeID", "DA Value")

# merge and order the two tables
# there were some genomes that were left out of this analysis
# when the DA table was created.  I need to remove those genomes
# from the copy number table and I need everything in the same order
cpy_num$tmp_name = paste(cpy_num$kogID, cpy_num$genomeID, sep="_")
da_tbl_melt$tmp_name = paste(da_tbl_melt$kogID, da_tbl_melt$genomeID, sep="_")

write.table(cpy_num, file="/Users/Scott/Projects/CompMetaG/experiments/nature_paper_44_metaG/fungal/copy_number/cpy_num_tmp.txt",
            quote=F, sep="\t", row.names=F, col.names=F)

write.table(da_tbl_melt, file="/Users/Scott/Projects/CompMetaG/experiments/nature_paper_44_metaG/fungal/copy_number/da_tbl_tmp.txt",
            quote=F, sep="\t", row.names=F, col.names=F)

# i had to merge and order the tables using perl
new_da_tbl = read.table("/Users/Scott/Projects/CompMetaG/experiments/nature_paper_44_metaG/fungal/copy_number/new_da_tbl.txt",
                        header=F, sep="\t")
colnames(new_da_tbl) = c("kogID", "genomeID", "DA_value", "tmp")

new_cpy_num_tbl = read.table("/Users/Scott/Projects/CompMetaG/experiments/nature_paper_44_metaG/fungal/copy_number/new_cpy_num.txt",
                        header=F, sep="\t")
colnames(new_cpy_num_tbl) = c("kogID", "genomeID", "Copy_Num", "tmp")

tbl = new_cpy_num_tbl[,1:3]
tbl$DA_value = new_da_tbl$DA_value


# make a box plot figure for each KOG of interest
kogs_of_interest = c("KOG0519", "KOG2161", "KOG3076", "KOG0001", "KOG0052",
                     "KOG1376", "KOG2571", "KOG1745", "KOG0399", "KOG2004",
                     "KOG0524", "KOG1317")
figs = list()
for ( k in kogs_of_interest ) {
  tmp = tbl[tbl$kogID == k,]
  p = ggplot(tmp, aes(x=as.character(DA_value), y=Copy_Num, group=DA_value)) + 
    geom_boxplot(outlier.size = NA) +
    geom_point(position = position_jitter(width = 0.2)) + 
    scale_x_discrete("DA Category", labels = my_labels,
                     limits = my_limits) +
    ylab("Copy Number") + 
    ggtitle(paste(k, " Genomic Copy Number vs DA Category", sep="")) +
    theme(plot.title = element_text(hjust = .5, size=20),
          axis.text = element_text(size=16),
          axis.title = element_text(size=18))
    
  figs[[k]] = p
}

library(grid)
library(gridExtra)

multiplot(figs)

ggplot(tbl, aes(x=as.character(DA_value), y=Copy_Num, group=DA_value)) + 
  geom_boxplot(outlier.size = NA) +
  geom_point(position = position_jitter(width = 0.2)) + 
  scale_x_discrete("DA Category", labels = my_labels,
                   limits = my_limits) +
  ylab("Copy Number") + 
  ggtitle(paste(opt$grp_name, " Genomic Copy Number vs DA Category", sep="")) +
  theme(plot.title = element_text(hjust = .5, size=20),
        axis.text = element_text(size=16),
        axis.title = element_text(size=18))
