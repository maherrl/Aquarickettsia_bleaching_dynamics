## Alpha diversity plots for Muller corals

rm(list=ls())
library("ggplot2")
library("ggthemes")

alphadiv <- read.csv("./data/alphadiv8692.csv")
head(alphadiv)

levels(alphadiv$geno) <- c("G1", "G3", "G4","G5","G7","G9","G10","G13","G20","G41","G44","G46","G47","G50","G57","G58")
names(alphadiv$bleach) <- gsub(x = names(alphadiv$bleach), pattern = "Aug", replacement = "Pre Bleach")
names(alphadiv$bleach) <- gsub(x = names(alphadiv$bleach), pattern = "Sep", replacement = "Bleaching")
levels(alphadiv$bleach)

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
  data_sum <- plyr::rename(data_sum, c("mean" = varname))
  return(data_sum)
}

myCol <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', 
           '#f032e6', '#bcf60c', '#008080',
           '#9a6324', '#800000', '#808000', 
          '#000075', '#808080', '#000000')

# alpha diversity boxplot
D <- ggplot(alphadiv, aes(x=bleach, y=faithPD)) +
  geom_boxplot(outlier.shape = NA, color = "gray35") +
  geom_point(aes(color = geno), position = position_jitter(width = .25, height = 0)) +
  scale_x_discrete(breaks=breakss, labels=labelss) +
  ylab("Faith's Phylogenetic Diversity") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10)) +
  scale_colour_manual(values = myCol)
D
# alpha diversity scatter plot
bleach.labs <- c("Pre Bleach", "Bleached")
bleach_labeller <- function(variable,value){
  return(bleach.labs[value])
}

ggplot(alphadiv, aes(x=geno, y=shannon, color = geno)) +
  geom_point() +
  scale_colour_manual(values = myCol) +
  facet_grid(. ~ alphadiv$bleach, labeller = bleach_labeller) +
  theme_bw() + 
  ylab("Shannon diversity index") +
  xlab("Genotype") +
  theme(axis.text.x = element_text(angle = 90))
  
