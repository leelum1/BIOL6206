#BIOL 6206 Management and Analysis of Environmental Data
#Kevan Lee Lum Project
#Use of Botryllus schlosseri as an indicator species

#######################################################
                    #Read in Data
#######################################################

spe <- read.csv("spe.csv", row.names=1, header=TRUE) 
env <- read.csv("env.csv", row.names=1, header=TRUE) 
View(spe)
View(env)

#######################################################
                #Load necessary packages
#######################################################

library(car) # to calculate vifs
library(vegan) #for PCA
source("evplot.R") #to plot PCA


#######################################################
              #Exploratory Data Analysis
#######################################################
attach(env)

#qqplots
qqnorm(salinity)
qqnorm(nitrate)
qqnorm(oxygen)
qqnorm(phosphate)
#All variables seem to be normally distributed

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
pairs(~latitude+salinity+nitrate++phosphate+oxygen, data=env,
      lower.panel=panel.smooth, upper.panel=panel.cor,
      pch=20, main="Species Ind. Variable Matrix")
#All variables have low multicollinearity, so good to go


#######################################################
                  #Logistic Regression
#######################################################
# Th objective is to determine what variables affect the presence of the species and derive a model to predict their presence
# The pres column in the env dataset is species presence/absence of Botryllus Schlosseri
# Logistic regression is preformed with latitude, salinity, nitrate, phosphate, and oxygen
# variables do not have to be normalized for logistic regression

mylogit <- glm(pres ~ latitude + longitude + salinity + nitrate + phosphate + oxygen, family = "binomial", data = env)
summary(mylogit)
vif(mylogit)
#all variance inflation factors are acceptable
#phosphate seems to be the only significant variable where p = 0.0232
#therefore rerun model using on phosphate

#But first lets try these three variables
mylogit2 <- glm(pres ~ nitrate + phosphate + oxygen, family = "binomial", data = env)
summary(mylogit2)
#nope not any better


mylogit3 <- glm(pres ~ phosphate, family = "binomial", data = env)
summary(mylogit3)
#Can't do vifs with just one variable
#For every one unit change in phosphate, the log odds of Pres increases by 6.613

exp(coef(mylogit3))
# For every one unit change in phosphate, the odds of Pres increases by a factor of 744.5
# Seems unreasonable. Needs more data and further research?


#Model Significance
with(mylogit3, pchisq(null.deviance - deviance, df.null - df.residual, lower.tail = FALSE))
# Model is significant with p-value is 0.017


#Influence measures
influence.measures(mylogit3)
# Samoa Pacific is highlighted


#Residuals
residuals(mylogit3, "deviance")
# Point 19 has the largest residual, as expected from the influence measures
hist(residuals(mylogit3, "deviance"))
# Most residuals appear to be positive
residuals(mylogit3, "pearson")
# Point 19 again has the highest residual
hist(residuals(mylogit3, "pearson"))
#Similar to the deviance residuals


#######################################################
            #Principle Components Analysis
#######################################################

#Data transformation to comparable units
#Create abiotic dataset with transformed variables
newenv <- env[,c("latitude", "longitude", "salinity")]
newenv$nitrate <- log10(nitrate)
newenv$phosphate <- log10(phosphate)
newenv$oxygen <- log10(oxygen)
View(newenv)

#Redundancy analysis, scale species to unit variance
env.pca <- rda(newenv, scale=TRUE) 
env.pca
summary(env.pca)
dev.new(width=10, height=5)
biplot(env.pca, display='species', main="PCA - scaling 2") 
#5 axes represent greater than 75% of the variance. 


#Kaiser-Guttman and Broken Stick models
ev <- env.pca$CA$eig # Get eigenvalues
dev.new(width=10, height=5)
evplot(ev)
#The Kaiser-Guttman models suggests 2 axes
#Broken-stick suggests 6 axes.


#Try again without the positional arguments
#Create abiotic dataset with transformed variables
newenv2 <- env[c("salinity")]
newenv2$nitrate <- log10(nitrate)
newenv2$phosphate <- log10(phosphate)
newenv2$oxygen <- log10(oxygen)
View(newenv)

#Redundancy analysis, scale species to unit variance
env.pca2 <- rda(newenv2, scale=TRUE) 
env.pca2
summary(env.pca2)
dev.new(width=10, height=5)
biplot(env.pca2, display='species', main="PCA - scaling 2") 
#This looks better and probably more useful
#Three axes are required 


#Kaiser-Guttman and Broken Stick models
ev2 <- env.pca2$CA$eig # Get eigenvalues
dev.new(width=10, height=5)
evplot(ev2)
#The Kaiser-Guttman models suggests 2 axes
#Broken-stick suggests 4 axes.


#######################################################
                        #BIOENV
#######################################################
bio <- bioenv(spe, newenv)
bio
summary(bio)
# The bioenv shows that latitude best describes the distribution in space of the species. 
# A strange result so let's try again without positioning
bio <- bioenv(spe, newenv2)
bio
summary(bio)
# Phosphate and oxygen best describe the distribution, which is more expected.


#######################################################
                        #ANOSIM
#######################################################
spe.bray = vegdist(spe)
spe.bray
ano <- anosim(spe.bray, env$group)
summary(ano)
plot(ano)
# The R statistic is 0.1902 and the significance is 0.001
# This indicates that there are significant differences between groups since the significance is less than 0.05


#######################################################
                        #SIMPER
#######################################################
sim <- with(newenv, simper(spe, env$group))
summary(sim)
# The species that are most responsible for the differences between A and B is:
# Botrylloides violaceus
# This is species is the most abundant in the dataset.
