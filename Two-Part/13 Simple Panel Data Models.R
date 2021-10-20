# Simple Panel Data Models

# Outline:
#   Difference-in-differences model
#   Panel data model with first differences

# Data files: 
#   KIELMC.csv
#   wagepan.csv

# setup


# Install packages
for(i in c("tidyverse", "stargazer", "magrittr", "haven")){
  if(!require(i, character.only = T)){
    install.packages(i, dependencies = T)
    require(i, character.only = T)
  }
}

wagepan <- read_stata("Z://DesktopC//LUMSA//2//Data Mining//wagepan.dta")

wagepan$nr <- as.factor(wagepan$nr)
wagepan$year <- as.factor(wagepan$year)

# Keep only two years of data
#wagepan %<>% filter(year == 1980 | year == 1981)

# Generate dummy for year and interaction term
wagepan %<>% mutate(d81hours = d81*hours,
                    d82hours = d82*hours,
                    d83hours = d83*hours,
                    d84hours = d84*hours,
                    d85hours = d85*hours,
                    d86hours = d86*hours,
                    d87hours = d87*hours)

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
model8 <- lm(wage ~ hours + educ + exper + d81, wagepan)
summary(model8)

# Regression models for each year
model9 <- update(model7, subset = year == 1980)
summary(model9)
model10 <- update(model7, subset = year == 1981)
summary(model10)

# Regression model with different intercept and slope for both years
model11 <- lm(wage ~ hours + d81+d82+d83+d84+d85+d86+d87+d81hours+d82hours+d83hours+d84hours+d85hours+d86hours+d87hours, wagepan)
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
  head(16)

# Panel data model with first differences
# Cannot be estimated due to perfect collinearity of deduc and dexper
lm(dwage ~ deduc + dexper + dhours, wagepan)

# Panel data model with first differences
lm(dwage ~ dhours, wagepan)
