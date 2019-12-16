myCol <- c('#e6194b', '#3cb44b', '#ffe119', 
           '#4363d8', '#f58231', '#911eb4', '#46f0f0', 
           '#f032e6', '#bcf60c', '#008080',
           '#9a6324', '#800000', '#808000', 
           '#000075', '#808080', '#000000')


nmds_bc <- metaMDSiter(ps_bc, k=2, trymax = 1000, maxit = 1000, autotransform=FALSE)
nmds_wu <- metaMDSiter(ps_wu, k=2, trymax = 1000, maxit = 1000, autotransform=FALSE)

meta <- as.data.frame(sample_data(ps))
head(meta)

G1 <- rownames(meta[which(meta[,3] == "G1"),])
G3 <- rownames(meta[which(meta[,3] == "G3"),])
G4 <- rownames(meta[which(meta[,3] == "G4"),])
G5 <- rownames(meta[which(meta[,3] == "G5"),])
G7 <- rownames(meta[which(meta[,3] == "G7"),])
G9 <- rownames(meta[which(meta[,3] == "G9"),])
G10 <- rownames(meta[which(meta[,3] == "G10"),])
G13 <- rownames(meta[which(meta[,3] == "G13"),])
G20 <- rownames(meta[which(meta[,3] == "G20"),])
G41 <- rownames(meta[which(meta[,3] == "G41"),])
G44 <- rownames(meta[which(meta[,3] == "G44"),])
G46 <- rownames(meta[which(meta[,3] == "G46"),])
G47 <- rownames(meta[which(meta[,3] == "G47"),])
G50 <- rownames(meta[which(meta[,3] == "G50"),])
G57 <- rownames(meta[which(meta[,3] == "G57"),])
G58 <- rownames(meta[which(meta[,3] == "G58"),])

mds.fig <- ordiplot(nmds_bc, display = "sites", type = "none")
points(mds.fig, "sites", pch = 19, col = "#e6194b", select = G1)



a <- plot_ordination(ps, nmds_bc, type = "samples", color = "geno.num", shape = "type") + 
  geom_point(aes(size = "bleach")) + 
  theme_classic() +
  theme(legend.position = "left") +
  labs(color = "Genotype") +
  scale_colour_manual(values = myCol)
a 

b <- plot_ordination(ps, nmds_bc, color = "bleach") + 
  geom_point(size = 2) + 
  theme_classic() +
  theme(legend.position = "left")
b + stat_ellipse(type = "norm", mapping = aes(fill = bleach))
