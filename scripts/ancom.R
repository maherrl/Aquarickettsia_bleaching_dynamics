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

# function
ANCOM.main = function(OTUdat,Vardat,
                      adjusted,repeated,
                      main.var,adj.formula,
                      repeat.var,longitudinal,
                      random.formula,
                      multcorr,sig,
                      prev.cut){
  
  p.zeroes=apply(OTUdat[,-1],2,function(x){
    s=length(which(x==0))/length(x)
  })
  
  zeroes.dist=data.frame(colnames(OTUdat)[-1],p.zeroes,row.names=NULL)
  colnames(zeroes.dist)=c("Taxon","Proportion_zero")
  
  zero.plot = ggplot(zeroes.dist, aes(x=Proportion_zero)) + 
    geom_histogram(binwidth=0.1,colour="black",fill="white") + 
    xlab("Proportion of zeroes") + ylab("Number of taxa") +
    theme_bw()
  
  #print(zero.plot)
  
  OTUdat.thinned=OTUdat
  OTUdat.thinned=OTUdat.thinned[,c(1,1+which(p.zeroes<prev.cut))]
  
  otu.names=colnames(OTUdat.thinned)[-1]
  
  W.detected   <- ancom.W(OTUdat.thinned,Vardat,
                          adjusted,repeated,
                          main.var,adj.formula,
                          repeat.var,longitudinal,random.formula,
                          multcorr,sig)
  
  W_stat       <- W.detected
  
  
  ### Bubble plot
  
  W_frame = data.frame(otu.names,W_stat,row.names=NULL)
  W_frame = W_frame[order(-W_frame$W_stat),]
  
  W_frame$detected_0.9=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.8=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.7=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.6=rep(FALSE,dim(W_frame)[1])
  
  W_frame$detected_0.9[which(W_frame$W_stat>0.9*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.8[which(W_frame$W_stat>0.8*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.7[which(W_frame$W_stat>0.7*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.6[which(W_frame$W_stat>0.6*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  
  final_results=list(W_frame,zero.plot)
  names(final_results)=c("W.taxa","PLot.zeroes")
  return(final_results)
}

# Ancom function
ancom.W = function(otu_data,var_data,
                   adjusted,repeated,
                   main.var,adj.formula,
                   repeat.var,long,rand.formula,
                   multcorr,sig){
  
  n_otu=dim(otu_data)[2]-1
  
  otu_ids=colnames(otu_data)[-1]
  
  if(repeated==F){
    data_comp=data.frame(merge(otu_data,var_data,by="Sample.ID",all.y=T),row.names=NULL)
    #data_comp=data.frame(merge(otu_data,var_data[,c("Sample.ID",main.var)],by="Sample.ID",all.y=T),row.names=NULL)
  }else if(repeated==T){
    data_comp=data.frame(merge(otu_data,var_data,by="Sample.ID"),row.names=NULL)
    # data_comp=data.frame(merge(otu_data,var_data[,c("Sample.ID",main.var,repeat.var)],by="Sample.ID"),row.names=NULL)
  }
  
  base.formula = paste0("lr ~ ",main.var)
  if(repeated==T){
    repeat.formula = paste0(base.formula," | ", repeat.var)
  }
  if(adjusted==T){
    adjusted.formula = paste0(base.formula," + ", adj.formula)
  }
  
  if( adjusted == F & repeated == F ){
    fformula  <- formula(base.formula)
  } else if( adjusted == F & repeated == T & long == T ){
    fformula  <- formula(base.formula)   
  }else if( adjusted == F & repeated == T & long == F ){
    fformula  <- formula(repeat.formula)   
  }else if( adjusted == T & repeated == F  ){
    fformula  <- formula(adjusted.formula)   
  }else if( adjusted == T & repeated == T  ){
    fformula  <- formula(adjusted.formula)   
  }else{
    stop("Problem with data. Dataset should contain OTU abundances, groups, 
         and optionally an ID for repeated measures.")
  }
  
  
  
  if( repeated==FALSE & adjusted == FALSE){
    if( length(unique(data_comp[,which(colnames(data_comp)==main.var)]))==2 ){
      tfun <- exactRankTests::wilcox.exact
    } else{
      tfun <- stats::kruskal.test
    }
  }else if( repeated==FALSE & adjusted == TRUE){
    tfun <- stats::aov
  }else if( repeated== TRUE & adjusted == FALSE & long == FALSE){
    tfun <- stats::friedman.test
  }else if( repeated== TRUE & adjusted == FALSE & long == TRUE){
    tfun <- nlme::lme
  }else if( repeated== TRUE & adjusted == TRUE){
    tfun <- nlme::lme
  }
  
  logratio.mat <- matrix(NA, nrow=n_otu, ncol=n_otu)
  for(ii in 1:(n_otu-1)){
    for(jj in (ii+1):n_otu){
      data.pair <- data_comp[,which(colnames(data_comp)%in%otu_ids[c(ii,jj)])]
      lr <- log((1+as.numeric(data.pair[,1]))/(1+as.numeric(data.pair[,2])))
      
      lr_dat <- data.frame( lr=lr, data_comp,row.names=NULL )
      
      if(adjusted==FALSE&repeated==FALSE){  ## Wilcox, Kruskal Wallis
        logratio.mat[ii,jj] <- tfun( formula=fformula, data = lr_dat)$p.value
      }else if(adjusted==FALSE&repeated==TRUE&long==FALSE){ ## Friedman's 
        logratio.mat[ii,jj] <- tfun( formula=fformula, data = lr_dat)$p.value
      }else if(adjusted==TRUE&repeated==FALSE){ ## ANOVA
        model=tfun(formula=fformula, data = lr_dat,na.action=na.omit)   
        picker=which(gsub(" ","",row.names(summary(model)[[1]]))==main.var)  
        logratio.mat[ii,jj] <- summary(model)[[1]][["Pr(>F)"]][picker]
      }else if(repeated==TRUE&long==TRUE){ ## GEE
        model=tfun(fixed=fformula,data = lr_dat,
                   random = formula(rand.formula),
                   correlation=corAR1(),
                   na.action=na.omit)   
        picker=which(gsub(" ","",row.names(anova(model)))==main.var)
        logratio.mat[ii,jj] <- anova(model)[["p-value"]][picker]
      }
      
    }
  } 
  
  ind <- lower.tri(logratio.mat)
  logratio.mat[ind] <- t(logratio.mat)[ind]
  
  
  logratio.mat[which(is.finite(logratio.mat)==FALSE)] <- 1
  
  mc.pval <- t(apply(logratio.mat,1,function(x){
    s <- p.adjust(x, method = "BH")
    return(s)
  }))
  
  a <- logratio.mat[upper.tri(logratio.mat,diag=FALSE)==TRUE]
  
  b <- matrix(0,ncol=n_otu,nrow=n_otu)
  b[upper.tri(b)==T] <- p.adjust(a, method = "BH")
  diag(b)  <- NA
  ind.1    <- lower.tri(b)
  b[ind.1] <- t(b)[ind.1]
  
  #########################################
  ### Code to extract surrogate p-value
  surr.pval <- apply(mc.pval,1,function(x){
    s0=quantile(x[which(as.numeric(as.character(x))<sig)],0.95)
    # s0=max(x[which(as.numeric(as.character(x))<alpha)])
    return(s0)
  })
  #########################################
  ### Conservative
  if(multcorr==1){
    W <- apply(b,1,function(x){
      subp <- length(which(x<sig))
    })
    ### Moderate
  } else if(multcorr==2){
    W <- apply(mc.pval,1,function(x){
      subp <- length(which(x<sig))
    })
    ### No correction
  } else if(multcorr==3){
    W <- apply(logratio.mat,1,function(x){
      subp <- length(which(x<sig))
    })
  }
  
  return(W)
  }


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

# load data for differential abundance analysis
load(file = "./data/ps.RData")
ps = subset_samples(ps, geno.num != 20)
# agglomerate to genus
ps <- subset_taxa(ps, Genus != "NA")
ps = filter_taxa(ps, function(x) sum(x > 10) > (0.2*length(x)), TRUE)
ps <- tax_glom(ps, "Genus")
# subset contrasts
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

figdf <- rbind(figdf_b.t, figdf_t.b)
figdf <- rbind(figdf_Aug, figdf_Sep, figdf_res, figdf_sus)
figdf <- merge(figdf,tax, by = "taxa_id")
write_csv(figdf, "data/ancomv2_all_fig.csv")

# add taxonomy
tax<-as(tax_table(ps),"matrix")
tax_cols <- c("Kingdom", "Phylum", "Class", "Order","Family","Genus", "Species")
tax<-as.data.frame(tax)
colnames(tax)
tax$taxa_id <- rownames(tax)
rownames(tax) <- NULL
dim(tax)
tax <- tax[,c(8,1:7)]

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
