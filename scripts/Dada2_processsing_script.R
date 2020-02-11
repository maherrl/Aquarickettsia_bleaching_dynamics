library(dada2); packageVersion("dada2")

pathF <- "~/Bleaching_Rickettsiales/fastq_filter/raw/F"
pathR <- "~/Bleaching_Rickettsiales/fastq_filter/raw/R"

# Forward and reverse fastq filenames have format: SAMPLENAME_R1.fq and SAMPLENAME_R2.fq
fnFs <- sort(list.files(pathF, pattern="_R1.fq", full.names = TRUE))
fnRs <- sort(list.files(pathR, pattern="_R2.fq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(pathF, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(pathR, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(260,210),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)

# File parsing
filtpathR <- "~/Bleaching_Rickettsiales/fastq_filter/raw/R/filtered"
filtpathF <- "~/Bleaching_Rickettsiales/fastq_filter/raw/F/filtered" 
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

#note to self: if this loop doesn't work, you've done something wrong. Stop trying to run it without the for loop. You've probably made an error earlier in the code, probably with the location of your forward and reverse reads and it's generated duplicate sample names
#normal to get warning message that there are duplicate sequences
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
}

st.all <- makeSequenceTable(mergers) #normal to get warning message saying the sequences being tabled vary in length
seqtab <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE)
saveRDS(seqtab, "~/Bleaching_Rickettsiales/seqtab.rds")

tax_silva <- assignTaxonomy(seqtab, "~/Bleaching_Rickettsiales/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
silva_sp <- addSpecies(tax_silva, "~/Bleaching_Rickettsiales/silva_species_assignment_v132.fa.gz")
ps_object <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), 
               tax_table(tax_silva))
ps_object_sp <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), 
                             tax_table(silva_sp))
