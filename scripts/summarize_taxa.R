## Simple script for summarizing relative abundance by genus and groups

# find most abundant genera
file = "./data/ps_rar8663.RData"

# functions
sderr <- function(x) {sd(x)/sqrt(length(x))}
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sderr(x[[col]]), na.rm=TRUE)
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

# code
ps <- tax_glom(ps, taxrank = "Genus", bad_empty = c(NA, "", " ", "\t"))
ps

ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
ps_melt <- psmelt(ps_rel)

head(sort(taxa_sums(ps_rel), decreasing = TRUE))

# summarize
genus <- ps_melt %>% group_by(Genus, bleach.type) %>% summarise(average = mean(Abundance), sd = sderr(Abundance))
genus_bleach <- ps_melt %>% group_by(Genus, bleach) %>% summarise(average = mean(Abundance), sd = sderr(Abundance))


genus_sd <- ps_melt %>% group_by(Genus, bleach.type) %>% summarise(sd = sderr(Abundance))

