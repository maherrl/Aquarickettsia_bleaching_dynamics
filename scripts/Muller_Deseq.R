####################################################
### DESeq2 
####################################################

rm(list=ls())

library("phyloseq"); packageVersion("phyloseq")
library("DESeq2")
packageVersion("DESeq2")
library("ggplot2")

load(file = "./data/ps.RData")
ps
# agglomerate to genus
ps <- subset_taxa(ps, Genus != "NA")
ps = filter_taxa(ps, function(x) sum(x > 10) > (0.2*length(x)), TRUE)
ps <- tax_glom(ps, "Genus")


####################################################
# First model
diagdds <- phyloseq_to_deseq2(ps, ~1)
design(diagdds) <- ~bleach*type

# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

# investigate test results table
resultsNames(diagdds)
res = results(diagdds)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.05
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
head(sigtab)

# specific contrasts
res1 <- results(diagdds, name = "bleach_Sep_vs_Aug", cooksCutoff = FALSE)
res1 = res1[which(res1$padj < alpha),]
res1 = cbind(as(res1, "data.frame"), as(tax_table(ps)[rownames(res1), ], "matrix"))
res1$comparison <- "Bleached vs. Pre-Bleached"
res2 <- results(diagdds, name = "type_susceptible_vs_resistant", cooksCutoff = FALSE)
res2 = res2[which(res2$padj < alpha),]
res2 = cbind(as(res2, "data.frame"), as(tax_table(ps)[rownames(res2), ], "matrix"))
res2$comparison <- "Susceptible vs. Resistant"
res3 <- results(diagdds, name = "bleachSep.typesusceptible", cooksCutoff = FALSE)
res3 = res3[which(res3$padj < alpha),]
res3 = cbind(as(res3, "data.frame"), as(tax_table(ps)[rownames(res3), ], "matrix"))
res3$comparison <- "Bleaching:Resistance"
res4 <- results(diagdds, name = "Intercept", cooksCutoff = FALSE)
res4 = res4[which(res4$padj < alpha),]
res4 = cbind(as(res4, "data.frame"), as(tax_table(ps)[rownames(res4), ], "matrix"))
res4$comparison <- "Intercept"

res_tab <- rbind(res1, res2, res3, res4)


# plot results
theme_set(theme_bw())
sigtabgen = subset(sigtab, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels=names(x))
ggplot(sigtab, aes(y=Genus, x=log2FoldChange, color=Phylum)) +
#  facet_grid(~comparison) +
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

####################################################
# 2nd model
diagdds <- phyloseq_to_deseq2(ps, ~1)
design(diagdds) <- ~bleach.type

# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

# investigate test results table
resultsNames(diagdds)
res = results(diagdds)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.05
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
head(sigtab)

# results by contrasts
res1 <- results(diagdds, contrast = c("bleach.type", "Sepsusceptible","Sepresistant"))
res1 = res1[which(res1$padj < alpha),]
res1 = cbind(as(res1, "data.frame"), as(tax_table(ps)[rownames(res1), ], "matrix"))
res1$comparison <- "Bleached Susceptible vs Resistant"
res2 <- results(diagdds, contrast = c("bleach.type", "Sepsusceptible","Augsusceptible"))
res2 = res2[which(res2$padj < alpha),]
res2 = cbind(as(res2, "data.frame"), as(tax_table(ps)[rownames(res2), ], "matrix"))
res2$comparison <- "Susceptible Bleached vs Pre-Bleach"
res3 <- results(diagdds, contrast = c("bleach.type", "Augsusceptible","Augresistant"))
res3 = res3[which(res3$padj < alpha),]
res3 = cbind(as(res3, "data.frame"), as(tax_table(ps)[rownames(res3), ], "matrix"))
res3$comparison <- "Pre-Bleached Susceptible vs Resistant"
res4 <- results(diagdds, contrast = c("bleach.type", "Sepresistant","Augresistant"))
res4 <- results(diagdds, name = "Intercept", cooksCutoff = FALSE)
res4 = res4[which(res4$padj < alpha),]
res4 = cbind(as(res4, "data.frame"), as(tax_table(ps)[rownames(res4), ], "matrix"))
res4$comparison <- "Resistant Bleached vs Pre-Bleach"

res_tab <- rbind(res1, res2, res3, res4)
# plot results
theme_set(theme_bw())
# Phylum order
x = tapply(res_tab$log2FoldChange, res_tab$Order, function(x) max(x))
x = sort(x, TRUE)
res_tab$Order = factor(as.character(res_tab$Order), levels=names(x))
# Genus order
x = tapply(res_tab$log2FoldChange, res_tab$Genus, function(x) max(x))
x = sort(x, TRUE)
res_tab$Genus = factor(as.character(res_tab$Genus), levels=names(x))
myCol <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4',
           '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', 
           '#008080', '#e6beff', '#9a6324', 
           '#800000', '#aaffc3', '#808000', '#ffd8b1', 
           '#000075', '#808080', '#ffffff', '#000000')
ggplot(res_tab, aes(y=Genus, x=log2FoldChange, color=Order)) +
  facet_grid(~comparison) +
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=4) + 
  scale_color_manual(values = myCol) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

