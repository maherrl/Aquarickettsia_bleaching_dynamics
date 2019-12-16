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


# load data for differential abundance analysis
load("./data/ps.RData")
ps

# Agglomerate samples to family ignoring, but still including taxa not identified to family
ps <- tax_glom(ps, taxrank = "Family", bad_empty = c(NA, "", " ", "\t"))
#ps <- tax_glom(ps, taxrank = "Genus", bad_empty = c(NA, "", " ", "\t"))

# OTU data or taxa data: This should be a data frame with each
# sample in rows and OTUs (or taxa) in columns. The first 
# column should be the sample identifier with column name
# "Sample.ID"
OTUdf1 <- as.data.frame(otu_table(ps))
OTUdf <- cbind(Sample.ID = rownames(OTUdf1), OTUdf1)
rownames(OTUdf) <- NULL


# Metadata: Dataframe with the first columns being the sample 
# identifier with column name "Sample.ID"

Vardat2 <- read.csv(file = "./data/map.csv")
colnames(Vardat2)[1] <- "Sample.ID"


# ANCOM test
# First must run the function ANCOM.main from the ANCOM_updated_code.R files

comp_test = ANCOM.main(OTUdat = OTUdf, 
                                    Vardat = Vardat2, 
                                    adjusted = T,
                                    repeated = F,
                                    main.var = "bleach:geno",
                                    adj.formula = "bleach + geno",
                                    repeat.var = NULL,
                                    longitudinal = F,
                                    random.formula = NULL,
                                    multcorr = 2,
                                    sig=0.05,
                                    prev.cut = 0.90)

comp_test$W.taxa # no significant results

comp_test2 = ANCOM.main(OTUdat = OTUdf, 
                       Vardat = Vardat2, 
                       adjusted = T,
                       repeated = F,
                       main.var = "bleach:type",
                       adj.formula = "bleach + type",
                       repeat.var = NULL,
                       longitudinal = F,
                       random.formula = NULL,
                       multcorr = 2,
                       sig=0.05,
                       prev.cut = 0.90)

comp_test2$W.taxa

comp_test3 = ANCOM.main(OTUdat = OTUdf, 
                        Vardat = Vardat2, 
                        adjusted = F,
                        repeated = F,
                        main.var = "bleach",
                        adj.formula = NULL,
                        repeat.var = NULL,
                        longitudinal = F,
                        random.formula = NULL,
                        multcorr = 2,
                        sig=0.05,
                        prev.cut = 0.90)

comp_test3$W.taxa


# Output results and join with taxonomy
res <- as.data.frame(comp_test3$W.taxa)
head(res)
colnames(res)
dim(res) # 59 6
tax<-as(tax_table(ps),"matrix")
tax_cols <- c("Kingdom", "Phylum", "Class", "Order","Family","Genus", "Species")
tax<-as.data.frame(tax)
colnames(tax)
tax$otu.names <- rownames(tax)
rownames(tax) <- NULL
dim(tax) #181 8
tax <- tax[,c(8,1:7)]

joined <- left_join(res,tax, by = "otu.names")
colnames(joined)
write.csv(res, file = "./data/ancom_bleach.csv")

# prepare for heat map plotting
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
sample_sums(ps_rel)
ps_rel
res_T <- res[which(res$detected_0.6 =="TRUE"),]
ps_rel_T <- prune_taxa(as.vector(res_T$otu.names), ps_rel)
ps_rel_T
plot_heatmap(ps_rel_T, sample.label = "bleach", taxa.order = "Family", 
             taxa.label = "Family", sample.order = "bleach",
             low = "grey", high = "darkblue", na.value = "white")

