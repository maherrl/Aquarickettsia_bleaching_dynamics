## Alpha diversity plots for Muller corals

rm(list=ls())
library("ggplot2")
library("ggthemes")


alphadiv <- read.csv("./data/alphadiv8692.csv")
head(alphadiv)
alphadiv$geno.num <- as.factor(alphadiv$geno.num)


levels(alphadiv$geno) <- c("G1", "G3", "G4","G5","G7","G9","G10","G13","G20","G41","G44","G46","G47","G50","G57","G58")
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

# Distinct colors from https://sashat.me/2017/01/11/list-of-20-simple-distinct-colors/
myCol <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', 
           '#f032e6', '#bcf60c', '#008080',
           '#9a6324', '#800000', '#808000', 
          '#000075', '#808080', '#000000')
breakss <- c("Aug","Sep")
labelss <- c("Pre Bleach", "Bleached")

# alpha diversity boxplot
A <- ggplot(alphadiv, aes(x=bleach, y=faithPD)) +
  geom_boxplot(outlier.shape = NA, color = "gray35") +
  geom_point(aes(color = geno.num), 
             position = position_jitter(width = .25, height = 0)) +
  scale_x_discrete(breaks=breakss, labels=labelss) +
  ylab("Faith's phylogenetic diversity") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10)) +
  scale_colour_manual(values = myCol)
A
ggsave(A, "./plots/observed_.pdf")

D <- ggplot(alphadiv, aes(x=bleach, y=faithPD)) +
  geom_boxplot(outlier.shape = NA, color = "gray35") +
  geom_point(aes(color = geno.num), position = position_jitter(width = .25, height = 0)) +
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

A <- ggplot(alphadiv, aes(x=geno.num, y=faithPD, color = geno.num)) +
  geom_point() +
  scale_colour_manual(values = myCol) +
  facet_grid(. ~ alphadiv$bleach, labeller = bleach_labeller) +
#  facet_grid(. ~ alphadiv$bleach) +
  theme_bw() + 
  ylab("Faith's phylogenetic diversity") +
  xlab("Genotype") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))
ggsave(A, filename = "./plots/faithPD.pdf")  
