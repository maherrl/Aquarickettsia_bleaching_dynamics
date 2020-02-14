#################################################################
## Combine output from Dada2 with a phylogenetic tree from qiime
## and add metadata. Also includes exploratory data.
#################################################################

rm(list=ls())

library("dada2")
library("seqinr")
library("biomformat")
library("data.table")
library("ggplot2")
library("ggthemes")
library("ampvis2")
library("cowplot")
library("phyloseq")
library("data.table")

setwd("/Users/Becca/Box Sync/Muller Bleaching and Rickettsiales/Muller-Acropora/")

# functions
fast_melt = function(physeq){
  # supports "naked" otu_table as `physeq` input.
  otutab = as(otu_table(physeq), "matrix")
  if(!taxa_are_rows(physeq)){otutab <- t(otutab)}
  otudt = data.table(otutab, keep.rownames = TRUE)
  setnames(otudt, "rn", "taxaID")
  # Enforce character taxaID key
  otudt[, taxaIDchar := as.character(taxaID)]
  otudt[, taxaID := NULL]
  setnames(otudt, "taxaIDchar", "taxaID")
  # Melt count table
  mdt = melt.data.table(otudt, 
                        id.vars = "taxaID",
                        variable.name = "SampleID",
                        value.name = "count")
  # Remove zeroes, NAs
  mdt <- mdt[count > 0][!is.na(count)]
  # Calculate relative abundance
  mdt[, RelativeAbundance := count / sum(count), by = SampleID]
  if(!is.null(tax_table(physeq, errorIfNULL = FALSE))){
    # If there is a tax_table, join with it. Otherwise, skip this join.
    taxdt = data.table(as(tax_table(physeq, errorIfNULL = TRUE), "matrix"), keep.rownames = TRUE)
    setnames(taxdt, "rn", "taxaID")
    # Enforce character taxaID key
    taxdt[, taxaIDchar := as.character(taxaID)]
    taxdt[, taxaID := NULL]
    setnames(taxdt, "taxaIDchar", "taxaID")
    # Join with tax table
    setkey(taxdt, "taxaID")
    setkey(mdt, "taxaID")
    mdt <- taxdt[mdt]
  }
  return(mdt)
}

# phyloseq object output from Dada2
ps <- readRDS("./data/ps_object_sp_new.rds")
# sequences object made from Dada2
# seqtab <- readRDS("./data/seqtab_strict.rds") #old strict ps object

# exporting sequences to a fasta file for import into qiime
uniqueSeqs <- as.list(colnames(seqtab))
#write.fasta(uniqueSeqs, uniqueSeqs, "./data/uniqueSeqs.fasta")

# phylogenetic tree made from qiime phylogeny align-to-tree-mafft-fasttree
tree = read_tree("./data/tree.nwk")
# import metadata and merge into phyloseq object
mapfile = "./data/map.txt"
map = import_qiime_sample_data(mapfile)
sample_data(ps) <- map
ps <- merge_phyloseq(ps, tree)

# export .tsv asv_table to import into qiime2
otu<-t(as(otu_table(ps),"matrix"))
otu_biom<-make_biom(data=otu)
write_biom(otu_biom,"./qiime/ps.biom")

# export taxonomy to import into qiime2
tax<-as(tax_table(ps),"matrix")
tax_cols <- c("Kingdom", "Phylum", "Class", "Order","Family","Genus", "Species")
tax<-as.data.frame(tax)
tax$taxonomy<-do.call(paste, c(tax[tax_cols], sep=";"))
for(co in tax_cols) tax[co]<-NULL
write.table(tax, "./qiime/taxonomy.txt", quote=FALSE, col.names=FALSE, sep="\t")

# summary of data
ps
summary(sample_data(ps))

ntaxa(ps)
nsamples(ps)
rank_names(ps)
sample_names(ps)[1:5]
sample_variables(ps)
phy_tree(ps)

# remove mitochondria and chloroplasts, is.na important becuase if not included
# this command will also remove all Family = NA or Order = NA
ps = subset_taxa(ps, (Family!="Mitochondria") | is.na(Family))
ps = subset_taxa(ps, (Order!="Chloroplast") | is.na(Order))

# Filter on prevalence or total counts
pst = fast_melt(ps)
prevdt = pst[, list(Prevalence = sum(count > 0), 
                    TotalCounts = sum(count)),
             by = taxaID]
keepTaxa = prevdt[(Prevalence >=0 & TotalCounts >29), taxaID]
ps = prune_taxa(keepTaxa,ps)
ps
sample_sums(ps)
min(sample_sums(ps))
#save(ps, file = "./data/ps.RData")

sample_data(ps)$geno.num=factor(get_variable(ps, "geno.num"))
# 28 is the bottom quartile frequency per feature (asv)

# exploratory code
observed <- estimate_richness(ps, measures = c('Observed'))
explore.df <- cbind(observed, sample_sums(ps), sample_data(ps)$bleach)
colnames(explore.df) <- c('Observed', 'Sample_Sums', 'Month')
observed_mean <- mean(explore.df$Observed)
sample_sum_mean <- mean(explore.df$Sample_Sums)
ggplot(data = explore.df, aes(x = Sample_Sums, y = Observed, color = Month)) + 
  geom_point() +
  geom_smooth(method="auto", se=TRUE, fullrange=FALSE, level=0.95, 
              inherit.aes = F, mapping = aes(Sample_Sums, Observed),
              data = explore.df) +
  ylab("Observed OTUs") +
  scale_colour_colorblind()

# Let's use ampvis2 again so we can easily make a rarefaction curve

# Need to convert from phyloseq to ampvis
av2_otutable <- data.frame(OTU = rownames(t(phyloseq::otu_table(ps)@.Data)),
                           t(phyloseq::otu_table(ps)@.Data),
                           phyloseq::tax_table(ps)@.Data,
                           check.names = F
)

#Extract metadata from the phyloseq object:
av2_metadata <- data.frame(phyloseq::sample_data(ps), 
                           check.names = F
)

av2_metadata <- cbind(rownames(av2_metadata), av2_metadata)

#Load the data with amp_load:
av2_obj <- amp_load(av2_otutable, av2_metadata)

# RARE CURVE
rare_plot_amp <- amp_rarecurve(data = av2_obj, color_by = "bleach")
rare_curve_plot <- rare_plot_amp + ylab('Observed ASVs (count)') + 
  geom_vline(xintercept=min(sample_sums(ps)), linetype='dashed') +
  scale_colour_colorblind() +
  xlim(c(0, 35000))
rare_curve_plot

summary(explore.df$Sample_Sums)
length(explore.df$Sample_Sums)

# plot alpha diversity
plot_richness(ps, measures = c('Shannon', 'Simpson'))

#################################################################################
# rarefy
num <- min(sample_sums(ps))
ps <- rarefy_even_depth(ps, sample.size = 8663, rngseed = 999) 
ps
sum(sample_sums(ps))
sample_sums(ps)
# plot alpha diversity
richness.rare <- cbind(estimate_richness(ps_rar, 
                                         measures = c('Shannon', 'Simpson')),
                       sample_data(ps_rar)$bleach)
colnames(richness.rare) <- c('Shannon', 'Simpson', 'Bleach')
richness.rare$Labels <- rownames(richness.rare)

ggplot(data = richness.rare, aes(x = Shannon, y = Simpson)) + 
  geom_point()

ggplot(data = richness.rare, aes(x = Shannon, y = (Simpson-1)*-1, color = Bleach)) + 
  geom_point() +
  #  geom_smooth(method="auto", se=TRUE, fullrange=FALSE, level=0.95, 
  #              inherit.aes = F, mapping = aes(Shannon, Simpson),
  #              data = richness.rare) +
  scale_color_colorblind() +
  geom_text(aes(label=ifelse(Shannon<1, Labels, ""), hjust=-0.1),
            show.legend = F) +
  xlim(c(0,max(richness.rare$Shannon))) +
  theme_cowplot()

# save ps object
save(ps, file = "./data/ps_rar8663.RData")
