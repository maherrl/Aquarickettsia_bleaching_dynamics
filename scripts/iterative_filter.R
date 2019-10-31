# Assessing read loss with iterative filtration steps. How much are we losing with each filtration/trimming step
filtFs_2 <- file.path(path, "filtered2", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs_2 <- file.path(path, "filtered2", paste0(sample.names, "_R_filt.fastq.gz"))
out <- filterAndTrim(fnFs, filtFs_2, fnRs, filtRs_2, truncLen=c(260,210),
                     compress=TRUE, multithread=TRUE)
write.table(out, file = "truncLen_out.txt", sep="\t")

filtFs_3 <- file.path(path, "filtered3", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs_3 <- file.path(path, "filtered3", paste0(sample.names, "_R_filt.fastq.gz"))
out <- filterAndTrim(filtFs_2, filtFs_3, filtRs_2, filtRs_3,
                     maxN=0, compress=TRUE, multithread=TRUE)
write.table(out, file = "maxN.txt", sep="\t")
#NO READS LOST, i guess there were no Ns

filtFs_4 <- file.path(path, "filtered4", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs_4 <- file.path(path, "filtered4", paste0(sample.names, "_R_filt.fastq.gz"))
out <- filterAndTrim(filtFs_3, filtFs_4, filtRs_3, filtRs_4, maxEE=c(2,2), compress=TRUE, multithread=TRUE)
write.table(out, file = "maxEE_out.txt", sep="\t")

filtFs_5 <- file.path(path, "filtered5", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs_5 <- file.path(path, "filtered5", paste0(sample.names, "_R_filt.fastq.gz"))
out <- filterAndTrim(filtFs_4, filtFs_5, filtRs_4, filtRs_5, truncQ=2, compress=TRUE, multithread=TRUE)
write.table(out, file = "truncQ_out.txt", sep="\t")

filtFs_6 <- file.path(path, "filtered6", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs_6 <- file.path(path, "filtered6", paste0(sample.names, "_R_filt.fastq.gz"))
out <- filterAndTrim(filtFs_5, filtFs_6, filtRs_5, filtRs_6, rm.phix=TRUE, compress=TRUE, multithread=TRUE)
write.table(out, file = "phiX_out.txt", sep="\t")
