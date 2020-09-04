###
### Down-sampling for perMANOVA tests
###

library("caret"); packageVersion("caret")
library("phyloseq")
library("vegan")
library("cba")
library("usedist")

load(file = "./data/ps_rar8663.RData")
ps

# full data
sampledf <- data.frame(sample_data(ps))
ps_bc <- phyloseq::distance(ps, method = "bray")
adonis(ps_bc ~ bleach*type, data = sampledf)

# Permutational down-sampling
set.seed(300)
nperm = 1000
vals = 5

output <- matrix(ncol = vals, nrow = nperm)

for (perm in 1:nperm) {
  
sub_meta <- downSample(sampledf, sampledf$bleach.type)
ps_bc_sub <- dist_subset(ps_bc, sub_meta$sampleid)
x <- adonis(ps_bc_sub ~ bleach*type, data = sub_meta)
output[perm,] <- x$aov.tab$`Pr(>F)`

}

output_df <- as.data.frame(output)
sum(output_df$V1 > 0.05)
sum(output_df$V2 > 0.05)
sum(output_df$V3 > 0.05)

#PERMDISP
set.seed(301)
nperm = 1000
vals = 2

output <- matrix(ncol = vals, nrow = nperm)

# repeat with sampledf$bleach and sampledf$type
for (perm in 1:nperm) {
  
  sub_meta <- downSample(sampledf, sampledf$bleach.type)
  ps_bc_sub <- dist_subset(ps_bc, sub_meta$sampleid)
  x <- anova(betadisper(ps_bc, sampledf$bleach.type, bias.adjust = TRUE))
  output[perm,] <- x$`Pr(>F)`
  
}

output_df <- as.data.frame(output)
sum(output_df$V1 > 0.05)
sum(output_df$V2 > 0.05)

###########################################################
# Repeat analysis with Aquarickettsia sequences removed
load(file = "./data/ps_rar823.RData")
ps

set.seed(302)
nperm = 1000
vals = 5

output <- matrix(ncol = vals, nrow = nperm)

for (perm in 1:nperm) {
  
  sub_meta <- downSample(sampledf, sampledf$bleach.type)
  ps_bc_sub <- dist_subset(ps_bc, sub_meta$sampleid)
  x <- adonis(ps_bc_sub ~ bleach*type, data = sub_meta)
  output[perm,] <- x$aov.tab$`Pr(>F)`
  
}

output_df <- as.data.frame(output)
sum(output_df$V1 > 0.05)
sum(output_df$V2 > 0.05)
sum(output_df$V3 > 0.05)

#PERMDISP
set.seed(303)
nperm = 1000
vals = 2

output <- matrix(ncol = vals, nrow = nperm)

# repeat with sampledf$bleach and sampledf$type
for (perm in 1:nperm) {
  
  sub_meta <- downSample(sampledf, sampledf$bleach.type)
  ps_bc_sub <- dist_subset(ps_bc, sub_meta$sampleid)
  x <- anova(betadisper(ps_bc, sampledf$type, bias.adjust = TRUE))
  output[perm,] <- x$`Pr(>F)`
  
}

output_df <- as.data.frame(output)
sum(output_df$V1 > 0.05)
sum(output_df$V2 > 0.05)
