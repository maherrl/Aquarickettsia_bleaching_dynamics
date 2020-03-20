#####################################################################
## Figures for the Acropora bleaching Rickettsiales manuscript
#####################################################################


rm(list=ls())

# libraries
library("phyloseq")
library("data.table")
library("plyr")
library("dplyr")
library("tidyr")
library("ggplot2")
library("reshape2")
library("indicspecies")
library("ggnetwork")
library("ape")
library("microbiome")
library("ggthemes")
library(cowplot)
library("ggsignif")

#functions
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

## Figure 1
#####################################################################
load(file = "./data/ps_rar8663.RData")
# subset to pre-bleached corals
ps <- subset_samples(ps, bleach == "Aug")
# summaryize to genus
ps <- tax_glom(ps, taxrank = "Genus", bad_empty = c(NA, "", " ", "\t"))
# relative abundance
ps = transform_sample_counts(ps, function(x) x/sum(x)) # turn into relative abundances
# melt the data
ps_melt <- psmelt(ps)

# summarize by Genus, sort, get top taxa
genus <- ps_melt %>% group_by(Genus) %>% summarise(average = mean(Abundance))
genus_sort <- arrange(genus, desc(average))
top <- genus_sort[which(genus_sort$average > 0.01),]
top <- top$Genus
ps_melt_top <- ps_melt[ps_melt$Genus %in% top,]
top[1] <- "Aquarickettsia"

# plot
ps_melt_top$Genus <- factor(ps_melt_top$Genus, levels = c("MD3-55", "Corynebacterium", "Exiguobacterium",
                                                          "Staphylococcus", "Cloacibacterium", "Acinetobacter"))
colorblind_pallette = c("#CC3300", "#E69F00", "#CC79A7", "#56B4E9", "#ffffff", "#009E73", "#0072B2", "#660066", )
colorblind_pallette = c("#009E73","#E69F00","#ffffff","#CC3300","#660066","#FFFF00")

supp.labs <- c("Resistant", "Susceptible")
names(supp.labs) <- c("resistant","susceptible")

p <- ggplot(ps_melt_top, aes(x = Sample, y = reorder(OTU, Abundance))) +
  geom_point(aes(size=Abundance, fill = Genus), shape = 21) +
  facet_grid(. ~ type, scales = "free", space = "free", switch = "y", 
             labeller = labeller(type = supp.labs)) +
  theme_facet() +
  scale_fill_manual(values = colorblind_pallette, breaks = c("MD3-55", "Corynebacterium", "Exiguobacterium",
                                                             "Staphylococcus", "Cloacibacterium", "Acinetobacter"),
                    labels = c("Aquarickettsia", "Corynebacterium", "Exiguobacterium",
                               "Staphylococcus", "Cloacibacterium", "Acinetobacter")) +
  guides(fill = guide_legend(override.aes = list(size=4))) +
  labs(fill = "Genus", size = "Relative \nAbundance")
p



# Figure 2 - Taxa plot
#######################################################
# load the rarefied data table with all species
load(file = "./data/ps_rar8663.RData")
ps
# transform to relative abundance
ps_rel <- transform(ps, "compositional")
# melt the data at the Genus level
ps_rel_genus_melt <- ps_rel %>%
  tax_glom(taxrank = "Genus") %>%
  psmelt()
head(ps_rel_genus_melt)

# Get genera with mean realtive abundance >0.01 across all samples 
genus_sum <- ps_rel_genus_melt %>% group_by(Genus) %>% dplyr::summarise(Aver = mean(Abundance))
genus_sub <- genus_sum[which(genus_sum$Aver > 0.01),]
names <- genus_sub$Genus
# Replace genera with <0.01 abundance with "NA"
ps_rel_genus_melt$genus <- ps_rel_genus_melt$Genus

ps_rel_genus_melt$genus[ps_rel_genus_melt$genus != "MD3-55" & 
                          ps_rel_genus_melt$genus != "Vibrio" &
                          ps_rel_genus_melt$genus != "Acinetobacter" &
                          ps_rel_genus_melt$genus != "Corynebacterium" &
                          ps_rel_genus_melt$genus != "Pseudoalteromonas" &
                          ps_rel_genus_melt$genus != "Cloacibacterium" &
                          ps_rel_genus_melt$genus != "Staphylococcus" &
                          ps_rel_genus_melt$genus != "Alteromonas" &
                          ps_rel_genus_melt$genus != "Aestuariibacter"] <- NA

levels(ps_rel_genus_melt$bleach) <- c("Pre-Bleach", "Bleached")
levels(ps_rel_genus_melt$type) <- c("Resistant","Susceptible")
levels(ps_rel_genus_melt$genus) <- c("Pseudoalteromonas","Acinetobacter","Staphylococcus","Aestuariibacter",
                                     "Cloacibacterium","Corynebacterium","Vibrio","Alteromonas","MD3-55")
# plot
bar_species = ggplot(ps_rel_genus_melt, aes(x = reorder(geno, geno.num), y=Abundance)) + 
  geom_bar(stat="identity", position="fill", aes(fill = reorder(genus, Abundance))) +
  scale_fill_manual(values=c("#999999", "#FFFF00",  "#CC3300","#0072B2","#660066", "#E69F00", "#56B4E9","#CC79A7",  "#009E73"), 
                    breaks = c("Pseudoalteromonas","Acinetobacter","Staphylococcus","Aestuariibacter",
                               "Cloacibacterium","Corynebacterium","Vibrio","Alteromonas","MD3-55"),
                    labels = c("Pseudoalteromonas","Acinetobacter","Staphylococcus","Aestuariibacter",
                               "Cloacibacterium","Corynebacterium","Vibrio","Alteromonas","Aquarickettsia")) +
  facet_grid(bleach~type, scales = "free_x", space = "free_x") +
  theme_bw() +
#  guides(fill = guide_legend(keywidth = 0.3, , keyheight =.40, ncol=1)) +
  ylab("Relative Abundance") +
  xlab("Genotype") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(fill = "Genus", size = "Genus")
bar_species


# Figure 3 - Alpha div
#######################################################
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


levels(alphadiv$bleach) <- c("Pre-Bleach", "Bleached")
levels(alphadiv$type) <- c("Resistant","Susceptible")

myCol <- c("#E69F00","#56B4E9")

# alpha diversity boxplot
A <- ggplot(alphadiv, aes(x=type, y=simpson)) +
  facet_wrap(~bleach)+
  geom_boxplot(outlier.shape = NA, color = "gray35") +
  geom_point(aes(color = type), 
             position = position_jitter(width = .25, height = 0)) +
  scale_x_discrete(breaks=breakss, labels=labelss) +
  ylab("Simpson's diversity index") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10),
        legend.position = c(0.8,0.2),
        legend.title = element_blank()) +
  scale_colour_manual(values = myCol) +
  ylim(0.85,1) +
  ggtitle("With Aquarickettsia")
A
B <- ggplot(alphadiv, aes(x=type, y=simpsona)) +
  facet_wrap(~bleach)+
  geom_boxplot(outlier.shape = NA, color = "gray35") +
  geom_point(aes(color = type), 
             position = position_jitter(width = .25, height = 0)) +
  scale_x_discrete(breaks=breakss, labels=labelss) +
  ylab("Simpson's diversity index") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10),
        legend.position = "none") +
  scale_colour_manual(values = myCol) +
  ylim(0.85,1) +
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
  ylab("Chao1 index") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10),
        legend.position = "none") +
  scale_colour_manual(values = myCol) +
  ylim(25,200) +
  ggtitle("Without Aquarickettsia")
D

plot_grid(A,B,C,D, nrow = 2, labels = c("A","B","C","D"))
