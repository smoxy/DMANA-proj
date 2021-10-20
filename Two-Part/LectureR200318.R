library(lme4) # for the analysis
library(haven) # to load the SPSS .sav file
library(tidyverse) # needed for data manipulation.
library(RColorBrewer) # needed for some extra colours in one of the graphs
library(lmerTest)# to get p-value estimations that are not part of the standard lme4 packages
#
popular2data<- read_sav(file ="https://github.com/MultiLevelAnalysis/Datasets-third-edition-Multilevel-book/blob/master/chapter%202/popularity/SPSS/popular2.sav?raw=true")
popular2data<-subset(popular2data, select=c(pupil, class, extrav, sex, texp, popular)) # we select just the variables we will use
head(popular2data) # we have a look at the first 6 observations
library(ggplot2)
ggplot(data=popular2data, aes(x=extrav, y=popular, group=as.factor(class)))+
  geom_point(size=.7, col="grey", alpha=.8, position = "jitter")+
  geom_smooth(method=lm,se=FALSE, col="black", size=.5, alpha=.8)+theme_minimal()+labs(title="Linear Relationship Between Popularity and Extraversion for 100 Classes", subtitle="The linear relationship between the two is not the same for all classes")
# We can also plot the difference in effect for the two genders. 
# We see there is likely an average difference between genders, 
# but no difference in slopes (regression coefficients).
ggplot(data=popular2data, aes(x=extrav, y=popular, col=as.factor(sex)))+
  geom_point(size=1, alpha=.7, position = "jitter")+
  geom_smooth(method=lm,se=T, size=1.5, linetype=1, alpha=.7)+theme_minimal()+labs(title="Linear Relationship Between Popularity and Extraversion for the 2 Genders", subtitle="The linear relationship between the two is similar for both genders, with a clear intercept difference")+
  scale_color_manual(name ="Gender", labels=c("Boys", "Girls"), values=c("lightblue", "pink"))
#
# Intercept only model The first model that we replicate is the intercept only model.
interceptonlymodel<-lmer(popular~1 + (1|class), data=popular2data) #to run the model
#If we look at the different inputs for the LMER function we:
  
# have "popular", which indicates the dependent variable we want to predict.
# a "~", that we use to indicate that we now give the other variables of interest.
# a "1" in the formula the function indicates the intercept.
# since this is an intercept only model, we do not have any other independent variables here.
# between brackets we have the random effects/slopes. 
# Again the value 1 is to indicate the intercept and 
# the variables right of the vertical "|" bar are use to indicate grouping variables. 
# In this case the class ID. So the dependent variable 'popular' is predicted by a intercept 
# and a random error term for the intercept.
# Finally we specify which dataset we want to use after the "data=" command
summary(interceptonlymodel) #to get paramater estimates.

# If we look at the summary output we see under the Random Effects that 0.7021 
# is the residual variance on the class level 
# and 1.2218 is the residual variance on the first level (pupil level). 
# This means that the intraclass correlation (ICC) is 0.7021/(1.2218+0.7021)=.36. 
# Under Fixed Effects the estimate of the intercept is stated, which is 5.078.
#
#Now we can now first add first (student) level predictors. The first level predictors are sex and extraversion. For now we just add them as fixed effects and not yet as random slopes.
model1<-lmer(popular~1 + sex + extrav + (1|class), data=popular2data)
summary(model1)
# As default the lmer function does only give test statistics and estimates, but no p-values. 
# However, because we use the lmerTest package we do get P-values. 
# The intercept is now 2.14, the regression coefficient for sex is 1.25, 
# and the regression coefficient for extraversion 0.44. 
# In the last column of the Fixed effects table of the output we see the P-values, 
# which indicate all regression coefficients are significantly different from 0.
#
# We now also (in addition to the level 1 variables that were both significant) add a predictor variable 
# on the second level (teacher experience).
model2<-lmer(popular~1 + sex + extrav + texp+(1|class), data=popular2data)
summary(model2)
# Now we also want to include random slopes.
model3<-lmer(popular~1 + sex + extrav + texp+(1+sex+extrav |class), data=popular2data)
summary(model3)
#We can see that all the fixed regression slopes are still significant. 
#However, no significance test for the Random effects are given, 
#but we do see that the error term (Variance) for the slope 
# of the variable sex is estimated to be very small (0.0024). #
#This likely means that there is no slope variation of the SEX variable between classes
# and therefore the random slope estimation can be dropped from the next analyses. 
#Since there is no direct significance test for this Variance we can use the ranova #
#function of the lmerTest package, which will give us an ANOVA-like table for random effects. #
#It checks whether the model becomes significantly worse if a certain random effect is dropped 
#(formally known as likelihood ratio tests), if this is not the case, the random effect is not significant.

ranova(model3)
# We see that the random effect of sex is indeed not significant (P=0.6792) 
# and the random effect of extraversion is significant (P<.0001).
#
# We continue after omitting the random slope of sex.
model4<-lmer(popular~1 + sex + extrav + texp+(1+extrav |class), data=popular2data)
summary(model4)
# As a final step, we can add a cross-level interaction between teacher experience and extraversion 
# (since this had a significant random effect, that we might be able to explain).
model5<-lmer(popular~1 + sex + extrav + texp+ extrav:texp+(1+extrav|class), data=popular2data)
summary(model5)
#In a plot we can also clearly see that years of teacher experience has influence 
# on both the intercept and the regression coefficient of extraversion on popularity.

ggplot(data=popular2data, aes(x=extrav, y=popular, col=as.factor(texp)))+
  viridis::scale_color_viridis(discrete = TRUE)+
  geom_point(size=.7, alpha=.8, position = "jitter")+
  geom_smooth(method=lm,se=FALSE, size=2, alpha=.8)+theme_minimal()+
  labs(title="Linear Relationship for Different Years of Teacher Experience as Observed", subtitle="The linear relationship between the two is not the same for all classes", col= "Years of\nTeacher\nExperience")
