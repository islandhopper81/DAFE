#!/usr/bin/evn Rscript

library(ggplot2)
library(reshape2)
library("getopt")
library("plotly")

# Each row is a parameter with 4 columns: 
#   1) long value, 
#   2) short value, 
#   3) argument type where 0 = no argument, 1 = required , 2 = optional
#   4) data type (logical, integer, double, complex, character)
# When you call this script the args are named like --verbose or -v
params = matrix(c(
  "con_scores", "c", 1, "character",
  "da_tbl", "d", 1, "character",
  "out", "o", 1, "character",
  "verbose", "v", 0, "logical"
  ), byrow=TRUE, ncol=4)
opt = getopt(params)

# set the ggplot theme
theme_update(plot.title = element_text(hjust=.5, size=22),
axis.text = element_text(size=16),
axis.title = element_text(size=16),
legend.title = element_text(size=18),
legend.text = element_text(size=16)
)

# read in the con scores
print("Reading in con scores")
con_scores = read.table(opt$con_scores, header=T, sep="\t")

# read in the DA table
print("Reading in DA table")
da_tbl = read.table(opt$da_tbl, header=T, sep="\t")

# melt the DA table
print("Melting the DA table")
da_tbl_melt = melt(as.matrix(da_tbl))
print(head(da_tbl_melt))

# filtering
print("filtering")
da_tbl_melt2 = da_tbl_melt[da_tbl_melt$value > -2,]
write.table(da_tbl_melt2, file = "melted.txt")

# merge the melted da table and con scores table
print("Merging")
print(summary(con_scores))
print(summary(da_tbl_melt2))
merged = merge(con_scores, da_tbl_melt2, by.x = "trait", by.y = "Var1")
print(summary(merged))
save(merged, file=paste(opt$out, "_tbl.RData", sep=""))

# make the FN mean plot
# this is the plot that I care most about because it is the only one
# that will go in my paper at this time
print("Making FN Mean plot")
bk_label = paste("Bulk Soil", "Enriched", sep="\n")
eq_label = "No Difference"
rz_label = paste("Rhizosphere", "Enriched", sep="\n")
ggplot(merged, aes(x=factor(value), y=FN_score, color=factor(value))) + 
  geom_boxplot(outlier.size=NA) +
  geom_jitter(width = 0.2) +
  theme(axis.text = element_text(size=10),
		axis.title = element_text(size=10)) +
  xlab("Enrichment Call") +
  ylab("Conservation Score") +
  scale_x_discrete(labels = c("-2" = bk_label, "-1" = bk_label, "0" = eq_label, "1" = rz_label)) + 
  scale_color_manual(values = c("-1" = "blue", "0" = "yellow", "1" = "red"),
	guide=F)

ggsave(paste(opt$out, "_FN_mean.png", sep=""), width=3.5, height=3)

exit()


# make the FN median plot
print("Making the FN Median Plot")
ggplot(merged, aes(x=factor(value), y=FN_median)) +
	geom_boxplot(outlier.size=NA) +
	xlab("Enrichment Call") +
	ylab("Conservation Score (Median Furthest Neighbor)") +
	scale_x_discrete(labels = c("-2" = "UK", "-1" = "BK Up", "0" = "Eq", "1" = "RZ Up")) + 
 	geom_jitter(width = 0.2)

ggsave(paste(opt$out, "_FN_median.png", sep=""))

# make NN mean plot
print("Making NN mean plot")
ggplot(merged, aes(x=factor(value), y=NN_score)) +
  geom_boxplot(outlier.size=NA) +
  xlab("Enrichment Call") +
  ylab("Conservation Score (Mean Nearest Neighbor)") +
  scale_x_discrete(labels = c("-2" = "UK", "-1" = "BK Up", "0" = "Eq", "1" = "RZ up")) +
  geom_jitter(width = 0.2)

ggsave(paste(opt$out, "_NN_mean.png", sep=""))

# make NN medain plot
print("Making NN medain plot")
ggplot(merged, aes(x=factor(value), y=NN_median)) +
  geom_boxplot(outlier.size=NA) +
  xlab("Enrichment Call") +
  ylab("Conservation Score (Median Nearest Neighbor)") +
  scale_x_discrete(labels = c("-2" = "UK", "-1" = "BK Up", "0" = "Eq", "1" = "RZ up")) +
  geom_jitter(width = 0.2)

ggsave(paste(opt$out, "_NN_median.png", sep=""))

# FN mean t.test Eq vs RZ up
print("FN mean t-test: Eq vs RZ up")
v1 = merged[merged$value == 0,"FN_score"]
v2 = merged[merged$value == 1, "FN_score"]
my_t = t.test(v1, v2)
print(my_t)

# FN median t.test Eq vs RZ up
print("FN median t-test: Eq vs RZ up")
v1 = merged[merged$value == 0,"FN_median"]
v2 = merged[merged$value == 1, "FN_median"]
my_t = t.test(v1, v2)
print(my_t)

# FN mean t.test Eq vs BK up
print("FN mean t-test: Eq vs BK up")
v1 = merged[merged$value == 0, "FN_score"]
v2 = merged[merged$value == -1, "FN_score"]
my_t = t.test(v1, v2)
print(my_t)

# FN median t.test Eq vs BK up
print("FN median t-test: Eq vs BK up")
v1 = merged[merged$value == 0, "FN_median"]
v2 = merged[merged$value == -1, "FN_median"]
my_t = t.test(v1, v2)
print(my_t)

# FN mean t.test BK up vs RZ up
print("FN mean t-test: BK up vs RZ up")
v1 = merged[merged$value == -1, "FN_score"]
v2 = merged[merged$value == 1, "FN_score"]
my_t = t.test(v1, v2)
print(my_t)

# FN median t.test BK up vs RZ up
print("FN median t-test: BK up vs RZ up")
v1 = merged[merged$value == -1, "FN_median"]
v2 = merged[merged$value == 1, "FN_median"]
my_t = t.test(v1, v2)
print(my_t)

# NN mean t.test Eq vs RZ up
print("NN mean t-test: Eq vs RZ up")
v1 = merged[merged$value == 0,"NN_score"]
v2 = merged[merged$value == 1, "NN_score"]
my_t = t.test(v1, v2)
print(my_t)

# NN median t.test Eq vs RZ up
print("NN median t-test: Eq vs RZ up")
v1 = merged[merged$value == 0,"NN_median"]
v2 = merged[merged$value == 1, "NN_median"]
my_t = t.test(v1, v2)
print(my_t)

# NN mean t.test Eq vs BK up
print("NN mean t-test: Eq vs BK up")
v1 = merged[merged$value == 0, "NN_score"]
v2 = merged[merged$value == -1, "NN_score"]
my_t = t.test(v1, v2)
print(my_t)

# NN median t.test Eq vs BK up
print("NN median t-test: Eq vs BK up")
v1 = merged[merged$value == 0, "NN_median"]
v2 = merged[merged$value == -1, "NN_median"]
my_t = t.test(v1, v2)
print(my_t)

# NN mean t.test BK up vs RZ up
print("NN mean t-test: BK up vs RZ up")
v1 = merged[merged$value == -1, "NN_score"]
v2 = merged[merged$value == 1, "NN_score"]
my_t = t.test(v1, v2)
print(my_t)

# NN median t.test BK up vs RZ up
print("NN median t-test: BK up vs RZ up")
v1 = merged[merged$value == -1, "NN_median"]
v2 = merged[merged$value == 1, "NN_median"]
my_t = t.test(v1, v2)
print(my_t)


# make figure comparing the NN and FN for the 1 group (ie RZ up)
merged_up = merged[merged$value == 1,]
p = ggplot(merged_up, aes(x=NN_score, y=FN_score)) + 
	geom_point() +
	ggtitle("RZ Up Clusters")
ggsave(paste(opt$out, "_NN_vs_FN.png", sep=""), plot=p)

#gg = ggplotly(p)
#save(gg, file=paste(opt$out, "_NN_vs_FN.RData", sep=""))
