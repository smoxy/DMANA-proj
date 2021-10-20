pks <- c("haven","ggplot2","scales","dplyr","tidyverse","magrittr","gridExtra","pscl","countreg")
for(i in pks){
  if(!require(i, character.only = T)){
    install.packages(i, dependencies = T)
    require(i, character.only = T)
  }
}

custom <- read_stata("Z://DesktopC//LUMSA//2//Data Mining//customers.dta")

custom %<>% filter(monthnumb < 37)

custom$gender <- as.factor(custom$gender)
custom$married <- as.factor(custom$married)
custom$catalogpromo <- round(custom$catalogpromo)
custom$stdMonth <- round(scale(custom$monthnumb),2)
custom$winter <- ifelse(custom$stdMonth==-0.72,1,
                        ifelse(custom$stdMonth==-0.63,1,
                               ifelse(custom$stdMonth==0.43,1,
                                      ifelse(custom$stdMonth==0.53,1,
                                             ifelse(custom$stdMonth==1.59,1,
                                                    ifelse(custom$stdMonth==1.68,1,0))))))
custom$winter <- as.factor(custom$winter)

custom %>% group_by(hh_key) %>% summarize(n_ppl = n())

keyGroup <-custom %>% group_by(hh_key) %>% summarize(items = sum(item),
                                                catalogs = sum(catalogpromo),
                                                retails = sum(retailpromo),
                                                prices = sum(pricepromo),
                                                itemPromo = sum(item[which(pricepromo>=1&item>0)]))

monthGroup <-custom %>% group_by(monthnumb) %>% summarize(avgItem = mean(item),
                                avgCatal = mean(catalogpromo),
                                avgRetail = mean(retailpromo)) %>% 
  mutate(year = ifelse(monthnumb < 13, 1, ifelse(monthnumb < 25, 2, ifelse(monthnumb < 37, 3, "   "))),
  fixedMon = c(seq(1,12,by=1),seq(1,12,by=1),seq(1,12,by=1)),
  winter = ifelse(fixedMon%%12==0, 1, ifelse(fixedMon%%11==0, 1, 0)))

p1 <- ggplot(data = monthGroup, aes(x = monthnumb, y = avgItem, group=1)) + 
      geom_line() +
      scale_x_continuous(name = "time [years]",breaks = c(0,12,24,36), labels = c(0,1,2,3)) +
      theme(plot.title = element_text(hjust = 0.5)) +
      labs(title = "Purchases",
       y = "av. number")

p2 <- ggplot(data = monthGroup, aes(x = monthnumb, y = avgCatal, group=1)) + 
      geom_line() +
      scale_x_continuous(name = "time [years]",breaks = c(0,12,24,36), labels = c(0,1,2,3)) +
      theme(plot.title = element_text(hjust = 0.5)) +
      labs(title = "Catalog promotions",
       y = "av. number")

p3 <- ggplot(data = monthGroup, aes(x = monthnumb, y = avgRetail, group=1)) + 
      geom_line() +
      scale_x_continuous(name = "time [years]",breaks = c(0,12,24,36), labels = c(0,1,2,3)) +
      theme(plot.title = element_text(hjust = 0.5)) +
      labs(title = "Retail promotions",
           y = "av. number")

grid.arrange(p1, p2, p3, ncol=1)


ggplot(data = monthGroup, aes(x = monthnumb, y = avgItem, group=as.factor(winter), color=as.factor(winter))) + 
  geom_boxplot() +
  scale_x_continuous(name = "Month",breaks = c(13,26), labels = c("Non Winter","Winter")) +
  theme() +
  labs(color = "Winter")



ggplot(data = custom, aes(x = monthnumb, y = item, group=as.factor(pricepromo), color=as.factor(pricepromo))) + 
  geom_boxplot()+
  scale_x_continuous(name = "time [years]",breaks = c(0,12,24,36), labels = c(0,1,2,3))+
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title="Purchased items per price promotion over time",y="n. items",color = "N° of price promotions")

as.matrix(tapply(custom$item, list(gender,married), mean))
################################################################
formula=item~pricepromo+retailpromo+catalogpromo+gender+married+winter
fm_pois <- glm(formula, data = custom, family = poisson)
fm_qpois <- glm(formula, data = custom, family = quasipoisson)
fm_nbin <- glm.nb(formula, data = custom)
fm_hurdle <- hurdle(formula, data = custom, dist = "negbin")
fm_zinb <- zeroinfl(formula, data = custom, dist = "negbin")
################################################################
summaryModels <- function(numb){
  fm <- list("ML-Pois" = fm_pois, "Quasi-Pois" = fm_qpois, "NB" = fm_nbin,
             "Hurdle-NB" = fm_hurdle, "ZINB" = fm_zinb)
  invisible(readline(prompt="Press [enter] for summary of models"))
  print(sapply(fm, function(x) coef(x)[1:12]))
  invisible(readline(prompt="Press [enter] for estimated standard errors of models"))
  print(cbind("ML-Pois" = sqrt(diag(vcov(fm_pois))),"Adj-Pois" = sqrt(diag(sandwich(fm_pois))),sapply(fm[-1], function(x) sqrt(diag(vcov(x)))[1:12])))
  invisible(readline(prompt="Press [enter] for logLik and Df analysis"))
  print(rbind(logLik = sapply(fm, function(x) round(logLik(x), digits = 0)),Df = sapply(fm, function(x) attr(logLik(x), "df"))))
  invisible(readline(prompt="Press [enter] for real vs fitted of numb"))
  print(round(c("Obs" = sum(item < numb+1),"ML-Pois" = sum(dpois(numb, fitted(fm_pois))),"NB" = sum(dnbinom(numb, mu = fitted(fm_nbin), size = fm_nbin$theta)),"NB-Hurdle" = sum(predict(fm_hurdle, type = "prob")[,numb+1]),"ZINB" = sum(predict(fm_zinb, type = "prob")[,numb+1]))))
  invisible(readline(prompt="Press [enter] for zero-augmented models (zero-part model)"))
  print(t(sapply(fm[4:5], function(x) round(x$coefficients$zero, digits = 3))))
}


rootogram(fm_hurdle, max = 30)
plot(factor(item==0)~retailpromo+catalogpromo+gender+married+winter,data=custom,main="Zero")
plot(factor(item>0)~retailpromo+catalogpromo+gender+married+winter,data=custom,main="Zero")

# Overall variations
custom %>% 
  select(item, retailpromo, catalogpromo) %>% 
  mutate_all(function(x) {x - mean(x)}) %>% # variable - overall mean
  as.data.frame %>% 
  stargazer(type = "text", omit.summary.stat = "mean")

# Between variations
custom %>% group_by(hh_key) %>%
  select(item, retailpromo, catalogpromo) %>% 
  summarize_all(mean) %>% 
  as.data.frame %>% 
  select(-hh_key) %>%
  stargazer(type = "text")

# Within variations
custom %>% group_by(hh_key) %>% 
  select(item, retailpromo, catalogpromo) %>% 
  mutate_all(function(x) {x - mean(x)}) %>% # demean
  as.data.frame %>% 
  select(-hh_key) %>%
  stargazer(type = "text", omit.summary.stat = "mean")

# Generate first differences
diff <- function(x) {x - dplyr::lag(x)}
fd <- custom %>% group_by(hh_key) %>%
  mutate(ditem = diff(item), 
         dretailpromo = diff(retailpromo), 
         dcatalogpromo = diff(catalogpromo))%>%
  select(hh_key, monthnumb, item, retailpromo, catalogpromo, 
         ditem, dretailpromo, dcatalogpromo)%>% 
  .[order(.$hh_key, .$monthnumb),]

summary(lm(item ~ ., data = fd))
  

fm_hurdle1 <- mixed_model(item ~ retailpromo+catalogpromo+gender+married+winter, 
                          random = ~ 1 | hh_key, 
                          data = custom, 
                          family = hurdle.negative.binomial(), 
                          zi_fixed = ~retailpromo+catalogpromo+gender+married+winter,
                          zi_random = ~ 1 | hh_key)
