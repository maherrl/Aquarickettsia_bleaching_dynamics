#################################################################
## Combine output from Dada2 with a phylogenetic tree from qiime
## and add metadata
#################################################################

library("dada2")
library("seqinr")
library("biomformat")

# phyloseq object output from Dada2
ps <- readRDS("./data/ps_object_sp_strict.rds")
# sequences object made from Dada2
seqtab <- readRDS("./data/seqtab_strict.rds")

# exporting sequences to a fasta file for import into qiime
uniqueSeqs <- as.list(colnames(seqtab))
write.fasta(uniqueSeqs, uniqueSeqs, "./data/uniqueSeqs.fasta")

# phylogenetic tree made from qiime phylogeny align-to-tree-mafft-fasttree
tree = read_tree("./qiime/tree.nwk")
# import metadata and merge into phyloseq object
mapfile = "./data/map.txt"
map = import_qiime_sample_data(mapfile)
sample_data(ps) <- map
ps <- merge_phyloseq(ps, tree)

# export .tsv asv_table to import into qiime2
otu<-t(as(otu_table(ps),"matrix"))
otu_biom<-make_biom(data=otu)
write_biom(otu_biom,"./qiime/ps_table.biom")

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
