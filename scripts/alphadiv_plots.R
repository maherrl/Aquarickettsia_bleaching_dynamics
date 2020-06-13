## Alpha diversity plots for Muller corals

rm(list=ls())
library("ggplot2")
library("ggthemes")
library("cowplot")

# load alpha diversity datasets for with and without aquarickettsia
alphadiv <- read.csv("./data/alphadiv_8663.csv")
alphadiv2 <- read.csv("./data/alphadiv_823.csv")
alphadiv2 <- alphadiv2[,1:5]
colnames(alphadiv2) <- c("X","observeda", "chao1a","simpsona","shannona")
alphadiv <- cbind(alphadiv, alphadiv2)
alphadiv <- alphadiv[,-17]
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
levels(alphadiv$bleach) <- c("Pre-Bleach", "Bleached")
levels(alphadiv$type) <- c("Resistant","Susceptible")

myCol <- c("#E69F00","#56B4E9")

# alpha diversity boxplot
A <- ggplot(alphadiv, aes(x=type, y=shannon)) +
  facet_wrap(~bleach)+
  geom_boxplot(outlier.shape = NA, color = "gray35") +
  geom_point(aes(color = type), 
             position = position_jitter(width = .25, height = 0)) +
  scale_x_discrete(breaks=breakss, labels=labelss) +
  ylab("Shannon's diversity index") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10),
        legend.position = "none") +
  scale_colour_manual(values = myCol) +
  ylim(2,5.5) +
  ggtitle("With Aquarickettsia")
A
B <- ggplot(alphadiv, aes(x=type, y=shannona)) +
  facet_wrap(~bleach)+
  geom_boxplot(outlier.shape = NA, color = "gray35") +
  geom_point(aes(color = type), 
             position = position_jitter(width = .25, height = 0)) +
  scale_x_discrete(breaks=breakss, labels=labelss) +
  ylab("Shannon's diversity index") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10),
        legend.position = "none") +
  scale_colour_manual(values = myCol) +
  ylim(2,5.5) +
  ggtitle("Without Aquarickettsia")
B

plot_grid(A,B, labels = c("A","B"))

C <- ggplot(alphadiv, aes(x=type, y=chao1)) +
  facet_wrap(~bleach)+
  geom_boxplot(outlier.shape = NA, color = "gray35") +
  geom_point(aes(color = type), 
             position = position_jitter(width = .25, height = 0)) +
  scale_x_discrete(breaks=breakss, labels=labelss) +
  ylab("Chao1 index") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10),
        legend.position = "none") +
  scale_colour_manual(values = myCol) +
  ylim(25,200) +
  ggtitle("With Aquarickettsia")
C
D <- ggplot(alphadiv, aes(x=type, y=chao1a)) +
  facet_wrap(~bleach)+
  geom_boxplot(outlier.shape = NA, color = "gray35") +
  geom_point(aes(color = type), 
             position = position_jitter(width = .25, height = 0)) +
  scale_x_discrete(breaks=breakss, labels=labelss) +
  ylab("Shannon's diversity index") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10),
        legend.position = "none") +
  scale_colour_manual(values = myCol) +
  ylim(25,200) +
  ggtitle("Without Aquarickettsia")
D

plot_grid(A,B,C,D, nrow = 2, labels = c("A","B","C","D"))

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
