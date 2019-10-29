#######################################################
## Dada2 script for paired-end sequence processing,
## ASV assignment, and taxonomic assignment
## From Grace Klinges
## Added October 29, 2018
######################################################

# Libraries
alibrary(dada2); packageVersion("dada2")

path <- "~/Bleaching_Rickettsiales/"
# Forward and reverse fastq filenames have format: SAMPLENAME_R1.fq and SAMPLENAME_R2.fq
fnFs <- sort(list.files(path, pattern="_R1.fq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(260,210),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

# File parsing
filtpathF <- "~/Bleaching_Rickettsiales/fastq/filteredF" # CHANGE ME to the directory containing your filtered forward fastqs
filtpathR <- "~/Bleaching_Rickettsiales/fastq/filteredR" # CHANGE ME ...
filtFs <- list.files(filtpathF, pattern="fastq.gz", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="fastq.gz", full.names = TRUE)
sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names
set.seed(100)
# Learn forward error rates
errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE)
# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE)

# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)

seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, "path/to/run1/output/seqtab.rds")
st.all <- readRDS("~/Bleaching_Rickettsiales/seqtab.rds")
seqtab <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE)

tax_silva <- assignTaxonomy(seqtab, "~/Bleaching_Rickettsiales/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
silva_sp <- addSpecies(tax_silva, "~/Bleaching_Rickettsiales/silva_species_assignment_v132.fa.gz")
ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), 
               tax_table(tax_silva))

