#################################################
## Ricketsiales
library(dplyr)
library(emmeans)

asinTransform <- function(p) { asin(sqrt(p)) }

normality.plots <- function(x) {
  par(mfrow=c(2,2))
  hist(residuals(x), main = "Histogram", xlab = "Values")
  boxplot(residuals(x), main = "Boxplot", ylab = "Values")
  qqnorm(residuals(x))
  qqline(residuals(x))
  plot(density(residuals(x)), main = "Kernel Density Estimate")
}

load(file = "./data/ps_rar8692.RData")
ps

ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
sample_sums(ps_rel)
sample_sums(ps)

# only keep Rickettsiales
rick_rel <- subset_taxa(ps_rel, Genus=="MD3-55")
rick_rel
rick <- subset_taxa(ps, Genus=="MD3-55")
rick
# calculate frequency of each Rickettsiales ASV per sample
rick_freq <- psmelt(rick_rel)
rick_freq$bleach.type <- paste(rick_freq$bleach, rick_freq$type, sep = "_")
head(rick_freq)

# Take melted data and sum up the relative abundance of taxa from the order Rickettsiales by sample
rick_freq_summary<-rick_freq %>%  dplyr::group_by(Sample) %>% summarise(SUM=sum(Abundance))
head(rick_freq_summary)

rick_meta<-rick_freq %>% select(Sample, bleach, geno.num, type, bleach.type, diseaseAug, diseaseSep) %>% dplyr::distinct()
head(rick_meta)
rick_freq_summary<- dplyr::full_join(rick_freq_summary, rick_meta, by="Sample")
head(rick_freq_summary)

## PLOTTING
#levels(rick_freq_summary$geno) <- c("G1", "G3", "G4","G5","G7","G9","G10","G13","G20","G41","G44","G46","G47","G50","G57","G58")
#levels(rick_freq_summary$bleach) <- c("PreBleach", "Bleached")
myCol <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', 
           '#f032e6', '#bcf60c', '#008080',
           '#9a6324', '#800000', '#808000', 
           '#000075', '#808080', '#000000')

bleach.labs <- c("Pre Bleach", "Bleached")
bleach_labeller <- function(variable,value){
  return(bleach.labs[value])
}

ggplot(rick_freq_summary, aes(x=geno.num, y=SUM, color = geno.num)) +
  geom_point() +
  scale_colour_manual(values = myCol) +
  facet_grid(. ~ rick_freq_summary$bleach, labeller = bleach_labeller) +
  theme_bw() + 
  ylab("Relative abundance Rickettsiales Genus MD3-55") +
  xlab("Genotype") +
  labs(color = "Genotype") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))

#############
# relative abundance of Rickettsiales is not normally distributed.
hist(rick_freq_summary$SUM, n=50)
hist(asinTransform(rick_freq_summary$SUM))
shapiro.test(asinTransform(rick_freq_summary$SUM))     
qqnorm(rick_freq_summary$SUM)

kruskal.test(SUM ~ bleach, rick_freq_summary)
kruskal.test(SUM ~ geno.num, rick_freq_summary)
kruskal.test(SUM ~ type, rick_freq_summary)

lm1 <- lm(asinTransform(SUM) ~ bleach * geno.num, data = rick_freq_summary)
lm3 <- lm(asinTransform(SUM) ~ bleach, data = rick_freq_summary)
lm2 <- lm(asinTransform(SUM) ~ bleach * type, data = rick_freq_summary)
qqnorm(resid(lm1))
normality.plots(lm2)
anova(lm3)
plot(SUM ~ time, data = rick_freq_summary)
lm1e <- emmeans(lm1, list(pairwise ~ bleach), adjust = "tukey")
