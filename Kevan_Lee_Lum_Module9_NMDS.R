# BIOL 6206 Module 9: NMDS Assignment                                          9.5/10
# Kevan Lee Lum

#1. Bring in data
moll <- read.table("moll.txt", row.names = 1, header=TRUE)
View(moll)


#2. Import libraries
library(ade4)
library(vegan)
library(gclus)
library(ape)
library(FactoMineR)
library(vegan3d)


#3. NMDS with 2 dimensions using the Bray-Curtis index
moll.nmds <- metaMDS(moll, distance="bray")
moll.nmds
moll.nmds$stress
# The stress value is 0.163
# This means that the rank order of the dissimilarities is being preserved by the inter point distances
# We have a useful representation 


#4. NMDS biplot
dev.new(title="NMDS on moll")
plot(moll.nmds, type="t", display='sites', main=paste("NMDS/Percentage difference - Stress =", round(moll.nmds$stress,3)))
# There is no discernable pattern in the plot
# It is not easy to see any groups 


#5. Shepard plot and goodness of fit
dev.new(title="NMDS - Shepard plot", width=12, height=6)
par(mfrow=c(1,2))
stressplot(moll.nmds, main="Shepard plot")
gof <- goodness(moll.nmds)
plot(moll.nmds, type="t", display='sites', main="Goodness of fit")
points(moll.nmds, display="sites", cex=gof*300)
# The Shepard's plot shows most points are around the monotonically increasing line, which indicates a good fit
# High R squared values indicate a good fit.

                        #But what exactly is the Shepard's plot plotting? -1
                        
# The goodness of fit plot shows that all points are influencing the stress
# Based on the goodness of fit plot, the Trubiska site has the most influence since it has the largest circle
# Cudrak and Klubina are also similar in size


#6. Tie the NMDS biplot to a cluster analysis
moll.bray <- vegdist(moll, method="bray") 
moll.ward <- hclust(moll.bray, method="ward.D")
plot(moll.ward)
# There are 4 groups at a dissimilarity of 2


#7. Combine cluster and NMDS biplot 
# Cut the previous dendrogram to yield 4 groups
moll.groups <- cutree(moll.ward, k=4)
grp.lev <- levels(factor(moll.groups))

# Combination with NMDS result
sit.sc <- scores(moll.nmds)
p <- ordiplot(sit.sc, type="n", main="NMDS/Percentage difference + clusters Ward/Percentage difference")
for (i in 1:length(grp.lev))
{
  points(sit.sc[moll.groups==i,], pch=(14+i), cex=2, col=i+1)
}
text(sit.sc, row.names(moll), pos=4, cex=0.7)
# Add the dendrogram
ordicluster(p, moll.ward, col="dark grey")
# Add a legend interactively
legend(locator(1), paste("Group",c(1:length(grp.lev))), 
       pch=14+c(1:length(grp.lev)), col=1+c(1:length(grp.lev)), pt.cex=2)
# The NMDS now shows better groupings by assigning different colours


#8. 3D NMDS
moll.nmds3d <- metaMDS(moll, k = 3, distance = 'bray')
moll.nmds3d$stress 
ordiplot3d(moll.nmds3d)
# 3D plot created with the vegan3d package
# A dynamic 3D graph could also have been plotted but looked like a lot of struggle
# The stress value is 0.121, which is lower than previous
# This happens because the third dimension adds a new dimension to better represent the rank dissimlarties to the inter point distances
# More space to spread out