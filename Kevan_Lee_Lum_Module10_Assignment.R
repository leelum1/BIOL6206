# BIOL 6206 Module 10 BIOENV              15/15
# Kevan Lee Lum

library(vegan)

# 1. Import the data
env <- read.table("env.txt", row.names = 1, header=TRUE)
spe <- read.table("spe.txt", row.names = 1, header=TRUE)
View(spe)
View(env)


#2. Double square root transform
trans <- sqrt(sqrt(spe))
View(trans)


#3. NMDS
spe.nmds <- metaMDS(trans, distance="bray")
plot(spe.nmds, type="t", display='sites', main=paste("NMDS/Percentage difference - Stress =", round(spe.nmds$stress,3)))
# The points appear to be well spread.


#4. Ward's cluster analysis
spe.bray <- vegdist(trans, method="bray") 
spe.ward <- hclust(spe.bray, method="ward.D")
dev.new(title="Ward's cluster analysis")
plot(spe.ward)
# At a height of 0.6, there appear to be 3 distinct groups
# The S6 branch seems to be just below the cutoff

# Combination with NMDS result
spe.groups <- cutree(spe.ward, k=3)
grp.lev <- levels(factor(spe.groups))
sit.sc <- scores(spe.nmds)
dev.new(title="NMDS plot with cluster colors")
p <- ordiplot(sit.sc, type="n", main="NMDS/Percentage difference + clusters Ward/Percentage difference")
for (i in 1:length(grp.lev))
{
  points(sit.sc[spe.groups==i,], pch=(14+i), cex=2, col=i+1)
}
text(sit.sc, row.names(spe), pos=4, cex=0.7)
# Add the dendrogram
ordicluster(p, spe.ward, col="dark grey")
# Add a legend interactively
legend(locator(1), paste("Group",c(1:length(grp.lev))), 
       pch=14+c(1:length(grp.lev)), col=1+c(1:length(grp.lev)), pt.cex=2)


#5. Create new env variables
attach(env)
env_trans <- env[c("DEPTH.M.")]
env_trans$log_copper <- log10(COPPER)
env_trans$log_manganese <- log10(MANGANESE)
env_trans$log_cobalt <- log10(COBALT)
env_trans$log_nickel <- log10(NICKEL)
env_trans$log_zinc <- log10(ZINC)
env_trans$log_cadmium <- log10(CADMIUM + 0.1) #To account for 0 value
env_trans$log_lead <- log10(LEAD)
env_trans$log_chromium <- log10(CHROMIUM)
env_trans$sq_carbon <- sqrt(X.CARBON)
env_trans$sq_nitrogen <- sqrt(X.NITROGEN)
detach(env)
attach(env_trans)
View(env_trans)

#Multicollinearity
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
pairs(~+log_copper+log_cadmium+log_chromium+log_cobalt+log_lead+log_manganese+log_nickel+log_zinc+DEPTH.M.+sq_carbon+sq_nitrogen, data=env_trans,
      lower.panel=panel.smooth, upper.panel=panel.cor,
      pch=20, main="Species Ind. Variable Matrix")
# Based on the matrix, we should be concerned about copper, zinc, and lead due to the high correlation values > 0.7


#6. Create new dataset from new variables
newenv <- env_trans[,c(1, 3, 4, 5, 7, 9, 10, 11)]
View(newenv)


#7. PCA
pca <- rda(newenv, scale=TRUE) #redundancy analysis, scale species to unit variance
summary(pca)
biplot(pca, main="PCA - scaling 2") 
# Approximately 84.8% is explained by the first two axes
# The sites are also well distributed


#8. Ward's cluster analysis
env.ward <- hclust(dist(scale(newenv)), "ward.D")
dev.new(title="Ward's cluster analysis")
plot(env.ward)
# At a height of 4, there appears to be 3 distinct groups
# This is very different from the previous cluster analysis.
# 4, 5, and 6 are in a similar group
# 1 and 2 are together in the pca whereas 1 and 8 are together in the nmds

# Tie to pca analysis
# Cut the dendrogram to yield 3 groups
gr <- cutree(env.ward, k=3)
grl <- levels(factor(gr))

# Extract the site scores, scaling 1
sit.sc1 <- scores(pca, display="wa", scaling=1)

# Plot the sites with cluster symbols and colours (scaling 1)
dev.new(title="Ordination and clustering")
p <- plot(pca, display="wa", scaling=1, type="n", 
          main="PCA correlation + clusters")
abline(v=0, lty="dotted")
abline(h=0, lty="dotted")
for (i in 1:length(grl))
{
  points(sit.sc1[gr==i,], pch=(14+i), cex=2, col=i+1)
}
text(sit.sc1, row.names(env), cex=0.7, pos=3)
# Add the dendrogram
ordicluster(p, env.ward, col="dark grey")
# Add legend interactively
legend(locator(1), paste("Cluster", c(1:length(grl))), pch=14+c(1:length(grl)), 
       col=1+c(1:length(grl)), pt.cex=2)


#9. BIOENV
??bioenv
??wisconsin #standardization
sol <- bioenv(trans, newenv)
sol
summary(sol)
# It appears that cadmium best describes the distribution in space of the species. 


#10. ANOSIM
??anosim
attach(env)
newenv$group <- GROUP
detach(env)
View(newenv)
ano <- anosim(spe.bray, newenv$group)
summary(ano)
plot(ano)
# The R statistic is 0.7823 and the significance is 0.006
# This indicates that there are significant differences between groups since the significance is less than 0.05


#11. SIMPER
??simper
(sim <- with(newenv, simper(trans, newenv$group)))
summary(sim)
# The species that are most responsible are:
# A & B: X18
# B & C: X50
# A & C: X50

