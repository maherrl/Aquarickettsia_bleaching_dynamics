## Plotting ordination
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

load(file = "../data/ps_rar8692.RData")
ps
sample_sums(ps)
ps_bc <- phyloseq::distance(ps, method = "bray")
nmds_bc <- metaMDSiter(ps_bc, k=2, trymax = 1000, maxit = 1000, autotransform=FALSE)

plot_ordination(ps, nmds_bc, color = "type", shape = "bleach") + 
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



Aug.df <- table[,Aug]
nitrate.df <- table[,nitrate]
ammon.df <- table[,ammon]

test.df <- cbind(control.df, nitrate.df, ammon.df)
classes <- c(rep("Aug", length(Aug)), rep("Sep", length(Sep)))
groups <- classes
classes <- c(rep("Aug.res", length(Aug.res)), rep("Aug.sus", length(Aug.sus)),
             rep("Sep.res", length(Sep.res)), rep("Sep.sus", length(Sep.sus)))
df <- test.df
dims <- c(1,2)
ellp.kind <- "ehull"

# For Weighted Unifra
ps_bc <- phyloseq::distance(ps, method = "wunifrac")
object <- metaMDS(ps_bc, k=2, binary = FALSE)
groups_bc <- meta$nutrient

mds.fig <- ordiplot(nmds_bc, display = "sites", type = "none", choices = dims, xlim = c(-1,1))
ordispider(nmds_bc, groups, col = "gray")
points(mds.fig, "sites", pch = 19, col = "#E69F00", select = Aug.res)
points(mds.fig, "sites", pch = 19, col = "#56B4E9", select = Aug.sus)
points(mds.fig, "sites", pch = 17, col = "#E69F00", select = Sep.res)
points(mds.fig, "sites", pch = 17, col = "#56B4E9", select = Sep.sus)
ordiellipse(nmds_bc, groups, conf = 0.95, label = FALSE, choices = dims, kind = ellp.kind, col = "gray", lwd = 2, show.groups = "Aug")
ordiellipse(nmds_bc, groups, conf = 0.95, label = FALSE, choices = dims, kind = ellp.kind, col = "gray", lwd = 2, show.groups = "Sep")
ordicenter(nmds_bc, groups, pch = 4, col = "red", label = FALSE)
