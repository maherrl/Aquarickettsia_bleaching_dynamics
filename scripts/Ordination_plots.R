## Plotting ordination
library("vegan")

# functions
ordicenter <- function (ord, groups, display = "sites", w = weights(ord, display), 
                        show.groups, ...) 
{
  weights.default <- function(object, ...) NULL
  pts <- scores(ord, display = display, ...)
  w <- eval(w)
  if (length(w) == 1) 
    w <- rep(1, nrow(pts))
  if (is.null(w)) 
    w <- rep(1, nrow(pts))
  if (!missing(show.groups)) 
  {
    take <- groups %in% show.groups
    pts <- pts[take, , drop = FALSE]
    groups <- groups[take]
    w <- w[take]
  }
  out <- seq(along = groups)
  inds <- names(table(groups))
  for (is in inds) 
  {
    gr <- out[groups == is]
    if (length(gr) > 1)
    {
      X <- pts[gr, ]
      W <- w[gr]
      ave <- apply(X, 2, weighted.mean, w = W)
      vegan:::ordiArgAbsorber(ave[1], ave[2], labels = is, FUN = text, ...)
    }
    if (length(gr) == 1)
    {
      X <- pts[gr, ]
      W <- w[gr]
      vegan:::ordiArgAbsorber(X[1], X[2], labels = is, FUN = text, ...)
    }
  }
  invisible()
}

load(file = "./data/ps_rar8663.RData")
ps
sample_sums(ps)
ps_bc <- phyloseq::distance(ps, method = "bray")
nmds_bc <- metaMDSiter(ps_bc, k=2, trymax = 1000, maxit = 1000, autotransform=FALSE)
pcoa_bc <- ordinate(ps, dist = ps_bc, "PCoA")

plot_ordination(ps, object, color = "type", shape = "bleach", axes = 2:3) + 
  geom_point(size = 2) + 
  theme_classic() +
  theme(legend.position = "left") +
  labs(color = "Resistance") +
  scale_colour_manual(values = myCol) 

# Prepare data for plotting
meta <- as.data.frame(sample_data(ps))
table <- as.data.frame(otu_table(ps))
head(meta)

Aug <- rownames(meta[which(meta[,2] == "Aug"),])
Sep <- rownames(meta[which(meta[,2] == "Sep"),])
Aug.res <- rownames(meta[which(meta[,2] == "Aug" & meta[,7] =="resistant"),])
Aug.sus <- rownames(meta[which(meta[,2] == "Aug" & meta[,7] =="susceptible"),])
Sep.res <- rownames(meta[which(meta[,2] == "Sep" & meta[,7] =="resistant"),])
Sep.sus <- rownames(meta[which(meta[,2] == "Sep" & meta[,7] =="susceptible"),])


dims <- c(1,2)
ellp.kind <- "ehull"

# For Weighted Unifra
object <- metaMDSiter(ps_bc, k=2, trymax = 1000, maxit = 1000, autotransform=FALSE)
meta$bleach.type <- paste(meta$bleach, meta$type, sep = "_")
bleach.type <- meta$bleach.type

mds.fig <- ordiplot(object, xlim = c(-1, 1), display = "sites", type = "none", choices = dims)
#ordispider(object, groups, col = "gray")
points(mds.fig, "sites", pch = 19, col = "#E69F00", select = Aug.res)
points(mds.fig, "sites", pch = 19, col = "#56B4E9", select = Aug.sus)
points(mds.fig, "sites", pch = 17, col = "#E69F00", select = Sep.res)
points(mds.fig, "sites", pch = 17, col = "#56B4E9", select = Sep.sus)
ordiellipse(object, bleach.type, conf = 0.95, label = FALSE, choices = dims, kind = ellp.kind, col = "#E69F00", lwd = 2, show.groups = "Sep_resistant")
ordiellipse(object, bleach.type, conf = 0.95, label = FALSE, choices = dims, kind = ellp.kind, col = "#E69F00", lwd = 2, show.groups = "Aug_resistant")
ordiellipse(object, bleach.type, conf = 0.95, label = FALSE, choices = dims, kind = ellp.kind, col = "#56B4E9", lwd = 2, show.groups = "Sep_susceptible")
ordiellipse(object, bleach.type, conf = 0.95, label = FALSE, choices = dims, kind = ellp.kind, col = "#56B4E9", lwd = 2, show.groups = "Aug_susceptible")

ordiellipse(object, groups, conf = 0.95, label = FALSE, choices = dims, kind = ellp.kind, col = "gray", lwd = 2, show.groups = "Aug")
ordiellipse(nmds_bc, groups, conf = 0.95, label = FALSE, choices = dims, kind = ellp.kind, col = "gray", lwd = 2, show.groups = "Sep")
ordicenter(nmds_bc, groups, pch = 4, col = "red", label = FALSE)


load(file = "./data/ps_rar8692.RData")

# extract distance to centroid
sampledf <- data.frame(sample_data(ps))
disp <- betadisper(ps_bc, sampledf$bleach, bias.adjust = TRUE)
dispd <- as.data.frame(disp$distances)
dispd <- cbind(dispd, sample_data(ps))
colnames(dispd)[1] <- "distance"