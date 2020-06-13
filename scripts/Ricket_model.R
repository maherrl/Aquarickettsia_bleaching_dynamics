#################################################
## Rickettsiales
library(dplyr)
library(emmeans)
library(mvabund)
library(cowplot)

asinTransform <- function(p) { asin(sqrt(p)) }

normality.plots <- function(x) {
  par(mfrow=c(2,2))
  hist(residuals(x), main = "Histogram", xlab = "Values")
  boxplot(residuals(x), main = "Boxplot", ylab = "Values")
  qqnorm(residuals(x))
  qqline(residuals(x))
  plot(density(residuals(x)), main = "Kernel Density Estimate")
}

load(file = "./data/ps_rar8663.RData")
map = import_qiime_sample_data(mapfile)
sample_data(ps) <- map
ps

ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
ps_rel
sample_sums(ps_rel)
sample_sums(ps)

# only keep Rickettsiales
rick_rel <- subset_taxa(ps_rel, Genus=="MD3-55")
rick_rel
rick <- subset_taxa(ps, Genus=="MD3-55")
rick
# calculate frequency of each Rickettsiales ASV per sample
rick_freq <- psmelt(rick_rel)
rick_freq$bleach.type <- paste(rick_freq$bleach, rick_freq$type, sep = "_")
head(rick_freq)

rick_freq$geno.num <- as.factor(rick_freq$geno.num)

# Take melted data and sum up the relative abundance of taxa from the order Rickettsiales by sample
rick_freq_summary<-rick_freq %>%  
  dplyr::group_by(Sample, geno.num, bleach, type) %>% 
  summarise(SUM=sum(Abundance))
head(rick_freq_summary)

## PLOTTING
#levels(rick_freq_summary$geno) <- c("G1", "G3", "G4","G5","G7","G9","G10","G13","G20","G41","G44","G46","G47","G50","G57","G58")
#levels(rick_freq_summary$bleach) <- c("PreBleach", "Bleached")
myCol <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', 
           '#f032e6', '#bcf60c', '#008080',
           '#9a6324', '#800000', '#808000', 
           '#000075', '#808080', '#000000')

levels(rick_freq_summary$bleach) <- c("Apparently Healthy", "Bleached")
levels(rick_freq_summary$type) <- c("Resistant","Susceptible")
rick_freq_summary$geno.num <- factor(rick_freq_summary$geno.num, c("3","7","1","4", "5", "9", "10", "13", "20", "41", "44", 
                                                                   "46", "47", "50", "57", "58"))
myCol <- c("#0072B2", "#CC3300")

A <- ggplot(rick_freq_summary, aes(x=geno.num, y=SUM, color = type)) +
  facet_grid(~bleach) +
  geom_point(size = 3) +
  scale_colour_manual(values = myCol) +
  theme_bw() + 
  theme(legend.position = c(0.9,0.8)) +
  xlab("Genotype") +
  ylab("Relative abundance Aquarickettsia") +
  scale_shape_manual(values = c(1,2))

B <- ggplot(rick_freq_summary, aes(x=type, y=SUM, col = type)) +
  facet_grid(~bleach) +
  geom_boxplot() +
  scale_colour_manual(values = myCol) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") +
  xlab("Genotype") +
  ylab("Relative abundance of Aquarickettsia") +
  scale_shape_manual(values = c(1,2))

plot_grid(A, B, rel_widths = c(1, 0.5), labels = c("A","B"))

#############
# relative abundance of Rickettsiales is not normally distributed.
hist(rick_freq_summary$SUM, n=50)
hist(asinTransform(rick_freq_summary$SUM))
shapiro.test(asinTransform(rick_freq_summary$SUM))     
qqnorm(rick_freq_summary$SUM)

kruskal.test(SUM ~ bleach, rick_freq_summary)
kruskal.test(SUM ~ geno.num, rick_freq_summary)
kruskal.test(SUM ~ type, rick_freq_summary)

lm1 <- lm(asinTransform(SUM) ~ bleach * geno.num, data = rick_freq_summary)
lm3 <- lm(asinTransform(SUM) ~ bleach, data = rick_freq_summary)
lm2 <- lm(asinTransform(SUM) ~ bleach * type, data = rick_freq_summary)
qqnorm(resid(lm1))
normality.plots(lm2)
anova(lm3)
plot(SUM ~ time, data = rick_freq_summary)
lm1e <- emmeans(lm1, list(pairwise ~ bleach), adjust = "tukey")

head(rick_freq_summary)

model <- glm(disease ~ SUM*bleach, family = binomial, data = rick_freq_summary)
summary(model)
plot(model)
plot(rick_freq_summary$SUM, 
     rick_freq_summary$disease, 
     ylab="Probability of disease", 
     xlab= "Relative abundance of Aquarickettsia", pch =16, 
     col = rick_freq_summary$bleach)

legend(0.1, 0.6, legend = c("Pre-Bleach", "Bleached"), lty = rep(1,2), col = c("black","red"))

xv <- seq(0,1,0.01)
vn <- rep("Pre-Bleach", length(xv))
yv <- predict(model, list(bleach=factor(vn), SUM = xv), type = "response")
lines(xv,yv,col = "black")

xv <- seq(0,1,0.01)
vn <- rep("Bleached", length(xv))
yv <- predict(model, list(bleach=factor(vn), SUM = xv), type = "response")
lines(xv,yv,col = "red")

legend(0.2, 0., legend = c("Pre-Bleach", "Bleached"), lty = rep(1,2), col = c("black","red"))

disease_Aug <- rick_freq_summary[which(rick_freq_summary$bleach == "Pre-Bleach"),]
model <- glm(disease ~ SUM, family = binomial, data = disease_Aug)
summary(model)
xv <- seq(0,1,0.01)
yv <- predict(model, list(SUM=exp(xv)), type = "response")
plot(disease_Aug$SUM, disease_Aug$disease, ylab="Probability of disease", xlab= "Relative abundance of Aquarickettsia", pch =16, 
     col = "blue")
lines(xv,yv,col = "red")

disease_Sep <- rick_freq_summary[which(rick_freq_summary$bleach == "Bleached"),]
model <- glm(disease ~ SUM, family = binomial, data = disease_Sep)
summary(model)
xv <- seq(0,1,0.01)
yv <- predict(model, list(SUM=exp(xv)), type = "response")
plot(disease_Sep$SUM, disease_Sep$disease, ylab="Probability of disease", xlab= "Relative abundance of Aquarickettsia", pch =16, 
     col = "blue")
lines(xv,yv,col = "red")
plot(model)

