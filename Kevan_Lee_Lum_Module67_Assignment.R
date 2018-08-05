#Kevan Lee Lum                                88/90 = 14.7/15
#Module 6 and 7 R Assignment
#Saved excel as csv to make things easier

plants <- read.table('plantabund2017.txt', header=TRUE, sep=",")
View(plants)
str(plants)

#1. Produce scatter plots of the dependent variable (C3) against the two independent variables (LAT and LONG)
plot(x=plants$C3, y=plants$LAT, 
     main='Latitude vs C3', xlab = 'C3', ylab = 'Latitude')
#The abundance of C3 appears to increase with latitude
plot(x=plants$C3, y=plants$LONG, 
     main='Longitude vs C3', xlab = 'C3', ylab = 'Longitude') 
#There doesn't seem to be any shape as the points are well spread.
#Since Latitude seems to have a better shape than longitude, 
#Latitude might have more of an effect on the relative abundance of C3 plants

                     #C3 is the dependent variable and as such ought to be on the Y axis -1

#2. Produce histograms of each of the variables, both dependent and independent 
hist(plants$C3, main="Histogram for C3", xlab="C3")
hist(plants$LONG, main="Histogram for Longitude", xlab="Longitude")
hist(plants$LAT, main="Histogram for Latitude", xlab="Latitude")
# The latitude variable illustrates a normal distribution
library(moments)
kurtosis(plants$C3)
skewness(plants$C3)
# For C3, the skewness is 0.9308(moderately skewed to the right) 
# the kurtosis is 2.9413 (leptokurtic)
kurtosis(plants$LONG)
skewness(plants$LONG)
# For LONG, the skewness is -0.1935(apprximately symmetric)
# the kurtosis is 1.9213 (no conclusion)
kurtosis(plants$LAT)
skewness(plants$LAT)
# For LAT, the skewness is 0.0488(apprximately symmetric)
# The kurtosis is 2.9116 (leptokurtic)

                         #Values close to 3 indicate a mesokurtic variable -1


#3.Transforming the data
# The C3 data needs to be tranformed because it is has a positive skew.
# This can be easily seen on the histrogram where most of the value are to the left of the plot
plants$LC3 <- log10(plants$C3)
View(plants)
# There is a problem with the zero C3 values since log10(0) is undefined


#4.Create LC3
plants$LC3 <- log10(plants$C3 + 0.1)
View(plants)
# No more errors


#5.Create LatLong
plants$LatLong <- (plants$LAT * plants$LONG)
View(plants)


#6. Correlation Analysis
attach(plants)
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
pairs(~+LAT+LONG+LatLong, data=plants,
      lower.panel=panel.smooth, upper.panel=panel.cor,
      pch=20, main="Plants Ind. Variable Matrix")
# Lat and Long correlation = 0.05528
# Long and LatLong correlation = 0.47601
# Lat and LatLong correlation = 0.90326
# We should be concerned about Lat and LatLong since the value is > 0.7 


#7. Regression Analysis
A = lm(LC3~LAT+LONG+LatLong)
summary(A)
# y = 6.6654 - 0.1753*Lat - 0.0864*Long + 0.002106*LatLong
# At a level of 0.05, none of the variables are significant, although LONG and LatLong are flagged as possibly significant
# The adjusted R-squared value is 0.4424, which suggests that less than half is explained by the independant variables
# Therefore the model is poor at explaining the variance in the relative abundance of C3 plants
# The model has a p-value of 4.263e-07, which is less than 0.05, therefore model is significant
library(car)
vif(A)
# The vifs are very large for each variable, which indicates high multicollinearity
# Much larger than the cutoff value of 10. Therefore we have a problem.
# High vifs suggest the the model is significantly impacted by multicollinearity
par(mfrow=c(2,2))
plot(A)
# The Residuals vs Fitted: No underlying trend in the data as shown by horizontal line
# The Normal Q-Q: Some large residuals on both ends of plot. Should be reduced
# The Scale-Location: No pattern. Three points flagged as possibly influential
# The Residuals vs Leverage: No points highly influential (past Cook's distance)
AIC(A) # 4.0513


#8.Anova comparison
anova(A)
# The anova table gives an idea of how the model performs without a particular variable
# The anova table shows that the LAT variable is the most significant
# This is different from the coefficients table since each variable is examined separately
# High multicollinearity would affect the actual significance of each independent vairables
# LatLong is somewhat significant with p-value of 0.0918
# LONG is not significant and can be dropped


#9. Regression revisited
C = lm(LC3~LAT+LatLong)
anova(C)
# Anova table suggests that LatLong is not significant, so try again with only LAT
D = lm(LC3~LAT)
summary(D)
# LC3 = 0.04149*LAT - 2.2342
# LC3 increase as LAT increases, which was expected from the scatterplot
# The model is significant, p-value = 6.733e-07
# The adjusted R-squared is 0.421, which is not great, lower then the previous model, but we know that the previous model is inaccurate
par(mfrow=c(2,2))
plot(D)
# The Residuals vs Fitted: No trends in the data still
# The Normal Q-Q: Does not look any better than befoe
# The Scale-Location: May be influenced by a few points. Point 43 in particular is flagged in all four plots
# The Residuals vs Leverage: Nobody is beyond Cook's distance
anova(D)
AIC(D) #4.203, similar to previous model


#10. Model Prediction
y = 10 ^ (-2.2342 + 0.041491 * 100) - 0.1
print(y)  # y = 82.105
# Regression equation should not be used to extrapolate date, only interpolate
# Any relationships beyond the limits of the original regression are unknown
# Therefore, prone to errors
detach(plants)


#11. Logistic Regression
species <- read.table('species.csv', header=TRUE, sep=",")
View(species)
str(species)
par(mfrow=c(2,2))
qqnorm(species$Pcov)
qqnorm(species$El)
qqnorm(species$dist)
qqnorm(species$Pres)
# The Pres data does not give a useful plot as the data is 1's and 0's


#12. Multicollinearity
attach(species)
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
pairs(~+Pcov+El+dist, data=species,
      lower.panel=panel.smooth, upper.panel=panel.cor,
      pch=20, main="Species Ind. Variable Matrix")
# There seems to be a problem with the Pcov and dist variables. Correlation = 0.8


#13. Variance inflation factors
mylogit <- glm(Pres ~ Pcov + El + dist, data = species, family = "binomial")
vif(mylogit)
# The vifs are low, Pcov = 2.049, El - 1.118, dist = 1.971
# All values are under 10, therefore acceptable. No large changes in confidence intervals


#14.Log Odds
summary(mylogit)
# the Pcov variable is significant, while the p-values of El and dist are > 0.05.
# For every one unit change in Pcov, the log odds of Pres increases by 0.0959
# For every one unit change in El, the log odds of Pres increases by 0.000309
# For every one unit change in dist, the log odds of Pres increases by 0.0250
# The AIC value is 27.358


#15.Odds Ratio
exp(coef(mylogit))
# For every one unit change in Pcov, the odds of Pres increases by a factor of 1.1006
# For every one unit change in El, the odds of Pres increases by a factor of 1.0003
# For every one unit change in dist, the odds of Pres increases by a factor of 1.0253


#16.Model Significance
with(mylogit, null.deviance - deviance)
with(mylogit, df.null - df.residual)
with(mylogit, pchisq(null.deviance - deviance, df.null - df.residual, lower.tail = FALSE))
# Model is significant since 0.00161 is less than 0.05


#17.Influence measures
influence.measures(mylogit)
# We should be concerned about point 9 and 19


#18.Residuals
residuals(mylogit, "deviance")
# Point 19 has the largest residual, as expected from the influence measures
hist(residuals(mylogit, "deviance"))
# There appears to be skew in the histogram, with the outlier clearly visible to the right
residuals(mylogit, "pearson")
# Point 19 again has the highest residual
hist(residuals(mylogit, "pearson"))
# The residuals have more of a normal distribution, with the exception of the maximum residual, which is point 19


#19. Remove unnecessary variables 
F <- glm(Pres ~ Pcov, data = species, family = "binomial")
summary(F)
# For every one unit change in Pcov, the log odds of Pres increases by 0.0766
# The AIC value is 24.049, which is lower than the previous model, and thus suggests a better fit model

exp(coef(F))
# For every one unit change in Pcov, the odds of Pres increases by a factor of 1.0796

with(F, pchisq(null.deviance - deviance, df.null - df.residual, lower.tail = FALSE))
# Model is significant since 0.000135 is much less than 0.05

influence.measures(F)
# Again, point 19 is flagged, while point 9 is not.

residuals(F, "deviance")
# Point 19 has the largest residual, as expected.
hist(residuals(F, "deviance"))
# The residuals have close to a normal distribution, with a larger number of negative residuals than expected, and point 19 as an outlier
residuals(F, "pearson")
# Point 19 has a very large residual compared to the rest of the points
hist(residuals(F, "pearson"))
# The residuals have a clear normal distribution if point 19 is not considered
detach(species)
