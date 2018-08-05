# BIOL6206 Management and Analysis of Environmental Data          40/40 =5/5
# Module 5 - Species Richness and Diversity - DUE October 23rd by 6pm
# Kevan Lee Lum - 816003573

#import necessary packages
library('vegan')
library('BiodiversityR')

#1. How many species were recorded in the sample?
mammals = read.csv("Mammals.csv")
View(mammals)
str(mammals)
ncol(mammals)
#The number of species is 139


#2. Which is the most species rich habitat of all sites pooled? 
environ = read.csv("MammEnv.csv")
View(environ)
str(environ)
??diversitycomp
diversitycomp(x = mammals, y = environ,
              factor1 = "Habitat", factor2 = NULL,
              index = c("richness"),
              method = c("pooled"),
              sortit = TRUE, digits = 8)
#The most species rich habitat is the Forest (101), which is expected, 
#followed by Secondary_forest (95) and Scrub (30)


#3. Create a box and whisker plot of the number of species by ???Use??? type
boxplot(Ind_Abundance ~ Use, data=environ, 
        main="Number of Species by Use", xlab="Use", ylab="Number of Species")
# Is there are significant differences in the number of species by Use?
# The protected area and wildlife reserve seem to have significantly more species than the others


#4. In which habitat is the plot with the largest sample size?
# Habitat is column 3 in the environ dataset
environ[which.max(environ$Ind_Abundance), 3]
# Forest has the largest sample size


#5. Which species is the most abundant across all habitat types?
??which
which.max(colSums(mammals))
#The EASO species is the most abundant


#6. How many individuals of the ???NOCA??? species are there in the total sample? 
sum(mammals$NOCA)
#There are 12 individuals of the NOCA species


#7. How many singletons are in the total sample?
length(which(colSums(mammals) == 1))
#There are 7 singletons in the total sample


#8. Using a sample-based method, plot the species accumulation curve with 95% confidence intervals 
??specaccum
sample_spa <- specaccum(mammals)
plot(sample_spa, xlab = "Plot", ylab = "Number of species", 
     ci = 0.95, ci.type = c("line"), ci.col = "red", ci.lty = 1)
# Add an individual-based method (???rarefaction???) curve
individual_spa <- specaccum(mammals, method = "rarefaction")
plot(individual_spa, add=TRUE)
# From comparing these curves, do individuals from the same species appear to be clustered
# The rarefaction curve rises faster than the sample curve, which suggests that individuals
# from the same species appear to be clustered.


#9. Plot an accumulation curve for all ???Use??? types
??accumcomp
accumcomp(x = mammals, y = environ, factor = "Use", scale = "Ind_Abundance",
          plotit=T, labelit=F, legend=T, rainbow=T, type="l", 
          xlab="Abundance", ylab="Plot")
# Which ???Use??? type requires more samples?
# The Forest Reserve has the lowest abundance and so requires more samples
# In which ???Use??? type did species accumulate fastest?
# Species accumulated fastest in the Protected Area


#10.Plot the overall rank abundance curve, scaled by proportion.
??rankabundance
rankabunplot(rankabundance(mammals), scale="proportion")
# Are most of the species common or rare? 
# Most of the species are rare as the curve is initially very steep


#11. Plot the rank abundance curve of ???Habitat??? types, scaled by log abundance
??rankabundcomp
rankabuncomp(x = mammals, y = environ, factor = "Habitat", scale="logabun",
             scaledx=F, type="l", rainbow=T, legend=T)
# Which Habitat has the lowest evenness? 
# The Scrub has the lowest eveness as it is the steepest

#12. Construct a rank-abundance model of species abundance for the Private_reserve.
??as.matrix
??radfit
pr_rows <- which(environ$Use == "Private_reserve")
pr_data <- t(as.matrix(colSums(mammals[pr_rows])))
plot(radfit(pr_data))
# Which are the two best-fitting models? 
# The two best fitting models are the Lognormal and the Mandelbrot


#13. Which ???Habitat??? type has lower Shannon-weaver diversity than the mean Shannon-weaver
# diversity for all the data?
??diversitycomp
diversitycomp(x = mammals, y = environ,
              factor1 = "Habitat", factor2 = NULL,
              index = "Shannon",
              method = "mean",
              sortit = TRUE, digits = 3)
# The Scrub has the lower Shannon-weaver diversity

# What is the standard deviation?
diversityresult(x = mammals, y = environ, factor = "Habitat", level = "Scrub",
                index = c("Shannon"), method = c("sd"), digits = 3)
# The standard deviation is 0.20


#14.Calculate the Simpson index for ???Habitat??? type. 
diversitycomp(x = mammals, y = environ,
              factor1 = "Habitat", factor2 = NULL,
              index = c("Simpson"),
              method = c("mean"),
              sortit = TRUE, digits = 3)
# Is the habitat with the lowest Simpson index the same as your answer for Q13? 
# No, the Forest has the lowest Simpson index
# What is the index value of the lowest Habitat type (to 2 d.p)?
# The index value of the lowest Habitat type is 0.821


#15.Calculate the Chao index for the ???Management??? factor and compare this to the order for raw
# species richness for Management
diversitycomp(x = mammals, y = environ,
              factor1 = "Management", factor2 = NULL,
              index = c("chao"),
              method = c("mean"),
              sortit = TRUE, digits = 8)
diversitycomp(x = mammals, y = environ,
              factor1 = "Management", factor2 = NULL,
              index = c("richness"),
              method = c("mean"),
              sortit = TRUE, digits = 8)
# Is the order of expected species richness the same as the order for actual species richness for Management
# No, the order is not the same