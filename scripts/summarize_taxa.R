## Simple script for summarizing relative abundance by genus and groups
rm(list=ls())

library(phyloseq)
library(dplyr)
library(microbiome)

# find most abundant genera
load(file = "./data/ps_rar8663.RData")
ps

# functions
sderr <- function(x) {sd(x)/sqrt(length(x))}

ps_rel_genus_melt <- ps %>%
  transform_sample_counts(function(x) x / sum(x)) %>%
  tax_glom(taxrank = "Genus", bad_empty = c(NA, "", " ", "\t")) %>%
  psmelt()

genus_bleach_type <- ps_rel_genus_melt %>% group_by(Genus, bleach, type) %>% summarise(average = mean(Abundance), sem = sderr(Abundance))
genus_bleach <- ps_rel_genus_melt %>% group_by(Genus, bleach) %>% summarise(average = mean(Abundance), sem = sderr(Abundance))
