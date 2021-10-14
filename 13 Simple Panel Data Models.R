# Simple Panel Data Models

# Outline:
#   Difference-in-differences model
#   Panel data model with first differences

# Data files: 
#   KIELMC.csv
#   wagepan.csv

# setup
rm(list = ls()) 
directory <- "C:/Econometrics/DataR/"

# Install packages
PackageNames <- c("tidyverse", "stargazer", "magrittr", "haven")
for(i in PackageNames){
  if(!require(i, character.only = T)){
    install.packages(i, dependencies = T)
    require(i, character.only = T)
  }
}

# Keep only two years of data
wagepan %<>% filter(year == 1980 | year == 1981)

# Generate dummy for year and interaction term
wagepan %<>% mutate(d1981 = (year == 1981),
                    d1981hours = d1981*hours)

# Get wage from log(wage)
wagepan %<>% mutate(wage = 10^lwage)

# List, describe, and summarize data
select(wagepan, nr, year, wage, hours, educ, exper) %>% head(10)
select(wagepan, nr, year, wage, hours, educ, exper) %>% str
select(wagepan, nr, year, wage, hours, educ, exper) %>% 
  as.data.frame %>%
  stargazer(type = "text")

# Panel data where nr is the cross sectional dimension and year is the time dimension
# Describe and summarize as panel data

# Number of observations
wagepan %>% group_by(nr) %>% summarize(n_nr = n()) # number of obs by 'nr'
wagepan %>% group_by(year) %>% summarize(n_year = n()) # number of obs by 'year'

# Overall variations
wagepan %>% 
  select(wage, hours, educ, exper) %>% 
  mutate_all(function(x) {x - mean(x)}) %>% # variable - overall mean
  as.data.frame %>% 
  stargazer(type = "text", omit.summary.stat = "mean")

# Between variations
wagepan %>% group_by(nr) %>%
  select(wage, hours, educ, exper) %>% 
  summarize_all(mean) %>% 
  as.data.frame %>% 
  select(-nr) %>%
  stargazer(type = "text")

# Within variations
wagepan %>% group_by(nr) %>% 
  select(wage, hours, educ, exper) %>% 
  mutate_all(function(x) {x - mean(x)}) %>% # demean
  as.data.frame %>% 
  select(-nr) %>%
  stargazer(type = "text", omit.summary.stat = "mean")


# Regression model with both years (ignoring that it is panel data set)
model6 <- lm(wage ~ hours + educ + exper, wagepan)
summary(model6)
model7 <- lm(wage ~ hours, wagepan)
summary(model7)
model8 <- lm(wage ~ hours + educ + exper + d1981, wagepan)
summary(model8)

# Regression models for each year
model9 <- update(model7, subset = year == 1980)
summary(model9)
model10 <- update(model7, subset = year == 1981)
summary(model10)

# Regression model with different intercept and slope for both years
model11 <- lm(wage ~ d1981 + hours + d1981hours, wagepan)
summary(model11)

# Generate first differences
diff <- function(x) {x - dplyr::lag(x)}
wagepan %<>% group_by(nr) %>%
  mutate(dwage = diff(wage), 
         deduc = diff(educ), 
         dhours = diff(hours),
         dexper = diff(exper))

wagepan %>% 
  select(nr, year, wage, hours, educ, exper, 
         dwage, dhours, deduc, dexper) %>% 
  head(10)

# Panel data model with first differences
# Cannot be estimated due to perfect collinearity of deduc and dexper
lm(dwage ~ deduc + dexper + dhours, wagepan)

# Panel data model with first differences
lm(dwage ~ dhours, wagepan)
