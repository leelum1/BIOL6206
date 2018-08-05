# Kevan Lee Lum                                              36/40 = 8/10
# BIOL 6206 Module 8 and 9 Assignment

#####################################
#Cluster Analysis
#####################################

#1. Read the data
env2 <- read.table("env2.txt", row.names = 1, header=TRUE)
View(env2)
str(env2)
# The dataset seems to be a set of abiotic data taken from 8 sample sites
# The row.names parameters tells R that the row header is the site column the txt file, and is not a variable


#2. Euclidean distance, 3.3.3
??dist
env2.de <- dist(scale(env2)) # scale function to standardize
dev.new(width=10, height=5) #opens new window
library(gclus)
source("coldiss.R") #from Chapter 3 Numerical Ecology with R
coldiss(env2.de, nc=16, diag=TRUE)
# "The use of the Euclidean distance on raw data is restricted to datasets that are dimensionally homogeneous"
# The data must be standardized since the Euclidean distance value is strongly influenced by the scale of each descriptor
# Since the variables have different scales, the data must be standardized to a uniform scale
# The numbers represent the sample sites.
                                           #what do the numbers represent? -1

#3. Single linkage agglomerative clustering analysis, 4.3.1
env2.ch.single <- hclust(env2.de, method="single")


#4. Cluster dendrogram
dev.new(width=12, height=8)
plot(env2.ch.single)


#5. Ward's minimum variance cluster analysis, 4.5
env2.ch.ward <- hclust(env2.de, method="ward")
dev.new(width=12, height=8)
plot(env2.ch.ward)
# The two dendrograms are similar in terms of the site position
# The single linkage exaggerates distances between clusters more
                                       #Not really. Site 3 is associated with sites 5,4, and 6, in the single inkage agglomerative case, while it is associated with sites 7 and 8 in the Ward's -1

#6. Bray Curtis matrix, 3.3.1
library(vegan)
??vegdist
env2.db <- vegdist(env2, method="bray") #no standardization needed
env2.db.single <- hclust(env2.db, method="single")
dev.new(width=12, height=8)
plot(env2.db.single)
# The relative positions of the sites have changed significantly from the euclidean distance
# Sites 1 and 2 are more separated, and site 3 is in a different position

                                        #The relative positons haven't changed -1


##############################
#Principle Components Analysis
##############################
#7.Variable transformation
attach(env2)
env2_trans <- env2[c("DEPTH.M.")]
env2_trans$log_copper <- log10(COPPER)
env2_trans$log_manganese <- log10(MANGANESE)
env2_trans$log_cobalt <- log10(COBALT)
env2_trans$log_nickel <- log10(NICKEL)
env2_trans$log_zinc <- log10(ZINC)
env2_trans$log_cadmium <- log10(CADMIUM + 0.1) #To account for 0 value
env2_trans$log_lead <- log10(LEAD)
env2_trans$log_chromium <- log10(CHROMIUM)
env2_trans$sq_carbon <- sqrt(X.CARBON)
env2_trans$sq_nitrogen <- sqrt(X.NITROGEN)
detach(env2)
View(env2_trans)


#8. PCA on transformed data, 5.3.2.2
env2.pca <- rda(env2_trans, scale=TRUE) #redundancy analysis, scale species to unit variance
env2.pca
# The data needs to be scaled since they have different units of measurement
# A uniform scale is needed for an accurate evaluation


#9. PCA summary
summary(env2.pca)
# Inertia is the sum of all correlations of the variables with themselves, which is the number of variables
# Inertia is the maximum correlation possible
# The inertia is 11, meaning that 11 is the maximum correlation possible
# The results are unconstrained since PCA is unconstrained
# 6.1
# Unconstrained analysis involves one data matrix, as opposed to multiple, 
# and reveals its major structure in a graph constructed from a reduced set of orthogonal axes, which is the case here

                  #unconstrained means it isn't constrained by any abiotic parameter -1


#10. Kaiser-Guttman and Broken Stick
source("evplot.R")
ev <- env2.pca$CA$eig # Get eigenvalues
dev.new(width=12, height=8)
evplot(ev)
# The Kaiser-Guttman model suggests 2 axes (0.88 variance)
# The broken stick suggests 1 axis (0.69 variance, which is below the standard 0.75)


#11. Draw PCA plot
dev.new(width=12, height=8)
biplot(env2.pca, main="PCA - scaling 2") 
#correlation biplot: each eigenvector is scaled to the square root of its eigenvalue


#12. PCA biplot
# Percentage nitrogen does not appear to be correlated with depth since at almost orthogonal angle
# The most different site appears to be site 7
# The most different site in the single linkage cluster analysis is site 3
# This result is expected to be different since the transformed data is being compared in PCA.
# PCA is an eigenvalue-based approach whereas Cluster analysis is distance-based.
# Therefore different values are being compared


#13. Ward's minimum variance cluster with PCA biplot
# Cut the previous dendrogram to yield 3 groups
gr <- cutree(env2.ch.ward, k=3)
grl <- levels(factor(gr))

# Extract the site scores, scaling 1
sit.sc1 <- scores(env2.pca, display="wa", scaling=1)

# Plot the sites with cluster symbols and colours (scaling 1)
dev.new(title="Ordination and clustering")
p <- plot(env2.pca, display="wa", scaling=1, type="n", 
          main="PCA correlation + clusters")
abline(v=0, lty="dotted")
abline(h=0, lty="dotted")
for (i in 1:length(grl))
{
  points(sit.sc1[gr==i,], pch=(14+i), cex=2, col=i+1)
}
text(sit.sc1, row.names(env2), cex=0.7, pos=3)
# Add the dendrogram
ordicluster(p, env2.ch.ward, col="dark grey")
# Add legend interactively
legend(locator(1), paste("Cluster", c(1:length(grl))), pch=14+c(1:length(grl)), 
       col=1+c(1:length(grl)), pt.cex=2)

