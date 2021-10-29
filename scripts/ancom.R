########################################################
# This script is for performing an ANCOM (Analysis of
# Composition of Microbiomes (ANCOM) to detect 
# differentially abundant taxa in microbial surveys.
# Created by Rebecca Maher
# Using Muller-Rickettsiales data
########################################################

# clear workspace-----------------------------
rm(list=ls())

# load libraries
library(exactRankTests)
library(nlme)
library(stats)
library(ggplot2)
library(dplyr)
library(phyloseq)


#########################################################################
## Trying new code from https://github.com/FrederickHuangLin/ANCOM

rm(list=ls())

library(exactRankTests)
library(nlme)
library(dplyr)
library(ggplot2)
library(compositions)
library(phyloseq)
library(robustbase)

source("scripts/ancom_v2.1.R")

# load data from phyloseq for differential abundance analysis
load(file = "./data/ps.RData")
ps = subset_samples(ps, geno.num != 20)
# agglomerate to genus
ps <- subset_taxa(ps, Genus != "NA")
ps = filter_taxa(ps, function(x) sum(x > 10) > (0.2*length(x)), TRUE)
ps <- tax_glom(ps, "Genus")
                 
# Here I am subsetting the data to run several ANCOM tests with different contrasts
# ancom will be run 6 independent times for each of these contrasts.
# For bleach == "Aug", ancom will contrast August resistant versus August susceptible samples,
                 # for type == "resistant", ancom will contrast resistant August versus resistant September samples
                 # for bleach.type == ..., ancom will comapre august susceptible versus september reseistant and vice versa
ps <- subset_samples(ps, bleach == "Aug")
ps <- subset_samples(ps, bleach == "Sep")
ps <- subset_samples(ps, type == "resistant")
ps <- subset_samples(ps, type == "susceptible")
ps <- subset_samples(ps, bleach.type == "Augresistant" | bleach.type == "Sepsusceptible")
ps <- subset_samples(ps, bleach.type == "Augsusceptible" | bleach.type == "Sepresistant")

# OTU data or taxa data: This should be a data frame with each
# sample in rows and OTUs (or taxa) in columns. The first 
# column should be the sample identifier with column name
# "Sample.ID"
# the rest of the code should be repeated for each ps variable (each contrast)
OTUdf <- as.data.frame(t(otu_table(ps)))
metadf <- read.csv(file = "./data/map.csv")
colnames(metadf)[1] <- "Sample.ID"
#levels(metadf$bleach.type) <- c("Sepsusceptible", "Sepresistant", "Augsusceptible", "Augresistant")

# Step 1: Data preprocessing

feature_table = OTUdf; meta_data = metadf; sample_var = "Sample.ID"; group_var = NULL
out_cut = 0.05; zero_cut = 0.90; lib_cut = 1000; neg_lb = FALSE
prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var, 
                                   out_cut, zero_cut, lib_cut, neg_lb)
feature_table = prepro$feature_table # Preprocessed feature table
meta_data = prepro$meta_data # Preprocessed metadata
struc_zero = prepro$structure_zeros # Structural zero info


# Step 2: ANCOM

main_var = "bleach.type"; p_adj_method = "fdr"; alpha = 0.05
adj_formula = NULL; rand_formula = NULL

res = ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, 
            alpha, adj_formula, rand_formula)

resdf <- as.data.frame(res$out)

# compiling results from each contrast for the figure
figdf_sus <- as.data.frame(res$fig$data)
figdf_res <- as.data.frame(res$fig$data)
figdf_Sep <- as.data.frame(res$fig$data)
figdf_Aug <- as.data.frame(res$fig$data)
figdf_b.t <- as.data.frame(res$fig$data)
figdf_t.b <- as.data.frame(res$fig$data)

figdf_sus$contrast <- "sus"
figdf_res$contrast <- "res"
figdf_Sep$contrast <- "Sep"
figdf_Aug$contrast <- "Aug"
figdf_b.t$contrast <- "ArSs"
figdf_t.b$contrast <- "AsSr"
                 
# add taxonomy
tax<-as(tax_table(ps),"matrix")
tax_cols <- c("Kingdom", "Phylum", "Class", "Order","Family","Genus", "Species")
tax<-as.data.frame(tax)
colnames(tax)
tax$taxa_id <- rownames(tax)
rownames(tax) <- NULL
dim(tax)
tax <- tax[,c(8,1:7)]
                 
                 
figdf <- rbind(figdf_b.t, figdf_t.b)
figdf <- rbind(figdf_Aug, figdf_Sep, figdf_res, figdf_sus)
figdf <- merge(figdf,tax, by = "taxa_id")
write_csv(figdf, "data/ancomv2_all_fig.csv")

resdf_sus <- resdf
resdf_res <- resdf
resdf_Sep <- resdf
resdf_Aug <- resdf
resdf_b.t <- resdf
resdf_t.b <- resdf

resdf_sus$contrast <- "sus"
resdf_res$contrast <- "res"
resdf_Sep$contrast <- "Sep"
resdf_Aug$contrast <- "Aug"
resdf_b.t$contrast <- "ArSs"
resdf_t.b$contrast <- "AsSr"

resdf <- rbind(resdf_b.t, resdf_t.b)
resdf <- rbind(resdf_Aug, resdf_Sep, resdf_res, resdf_sus)
resdf <- merge(resdf,tax, by = "taxa_id")
write_csv(resdf, "data/ancomv2_all_res.csv")



# Step 3: Volcano Plot

# Number of taxa except structural zeros
n_taxa = ifelse(is.null(struc_zero), nrow(feature_table), sum(apply(struc_zero, 1, sum) == 0))
# Cutoff values for declaring differentially abundant taxa
cut_off = c(0.9 * (n_taxa -1), 0.8 * (n_taxa -1), 0.7 * (n_taxa -1), 0.6 * (n_taxa -1))
names(cut_off) = c("detected_0.9", "detected_0.8", "detected_0.7", "detected_0.6")

# Annotation data
dat_ann = data.frame(x = min(figdf$x), y = cut_off["detected_0.6"], label = "W[0.6]")

# fig = res$fig +  
#   geom_hline(yintercept = cut_off["detected_0.6"], linetype = "dashed") + 
#   geom_text(data = dat_ann, aes(x = x, y = y, label = label), 
#             size = 4, vjust = -0.5, hjust = 0, color = "orange", parse = TRUE)
# fig  
# 
# # specialized plot
# figdf <- as.data.frame(fig$data)
# figdf <- cbind(figdf,tax, by = "taxa.id")
# figdf <- figdf[,-6]
# figdf$group <- as.factor(figdf$group)
# levels(figdf$group) <- c("Pre-Bleach Susceptible", "Bleached Resistant", "Bleached Susceptible")
# levels(figdf$group) <- c("Pre-Bleach Resistant", "Pre-Bleached Susceptible", "Bleached Resistant")

figdf <- read.csv(file = "./data/ancomv2_all_fig.csv")
# Replace name
figdf <- figdf %>% 
  mutate(Genus = as.character(Genus)) %>% 
  mutate(Genus = replace(Genus, Genus == 'MD3-55', 'Aquarickettsia'))
# order genus
x = tapply(figdf$y, figdf$Genus, function(x) max(x))
x = sort(x, TRUE)
figdf$Genus = factor(as.character(figdf$Genus), levels=names(x))
figdf$col_genus <- figdf$Genus

figdf$col_genus[figdf$col_genus != "Aestuariibacter" & 
                  figdf$col_genus != "Aquarickettsia" &
                  figdf$col_genus != "Exiguobacterium" &
                  figdf$col_genus != "Marivita" &
                  figdf$col_genus != "HIMB11" &
                  figdf$col_genus != "Pseudoalteromonas" &
                  figdf$col_genus != "Staphylococcus" &
                  figdf$col_genus != "Alteromonas"] <- NA
levels(figdf$col_genus)
# add new factor
figdf$col_genus <- factor(figdf$col_genus, levels = c(levels(figdf$col_genus), "Other"))
# convert NAs to other
figdf$col_genus[is.na(figdf$col_genus)] = "Other"

# change facet names
figdf$contrast <- as.factor(figdf$contrast)
levels(figdf$contrast) <- c("Aug","Sep","res","sus")
figdf$contrast <- factor(figdf$contrast, c("Aug","Sep","sus","res"))

levels(figdf$contrast) <- c("Apparently Healthy Resistant vs\nBleached Susceptible", 
                            "Apparently Healthy Susceptible vs\nBleached Resistant")

levels(figdf$contrast) <- c("Apparently Healthy\nSusceptible vs Resistant","Bleached\nSusceptible vs Resistant", 
                            "Susceptible Bleached\nvs Apparently Healthy","Resistant Bleached\nvs Apparently Healthy")

kelly_colors = c('#F3C300',  '#008856','#875692', '#F38400', '#A1CAF1', '#BE0032', 
                '#C2B280',  '#222222','#848482',  '#E68FAC', '#0067A5', 
                '#F99379', '#604E97', '#F6A600', '#B3446C', '#DCD300', 
                '#882D17', '#8DB600', '#654522', '#E25822', '#2B3D26', 
                '#222222', '#F3C300', '#875692', '#F38400', '#A1CAF1', '#BE0032', 
                '#C2B280', '#848482', '#008856', '#E68FAC', '#0067A5')

ggfig <- ggplot(figdf, aes(x = x, y = y, color = col_genus)) +
  geom_vline(xintercept = 0, color = "grey") +
  geom_point() +
  facet_grid(~contrast) +
  ylab("W statistic") +
  xlab("CLR mean difference") +
  scale_color_manual(name = "Genus", values = kelly_colors) +
  geom_hline(yintercept = 18, linetype = "dashed") + 
  geom_text(data = dat_ann, aes(x = x, y = y, label = label), 
            size = 4, vjust = -0.5, hjust = 0, color = "orange", parse = TRUE)


ggfig
