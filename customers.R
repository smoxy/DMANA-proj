pks <- c("haven","ggplot2","scales","dplyr","tidyverse","magrittr","gridExtra","pscl","countreg")
for(i in pks){
  if(!require(i, character.only = T)){
    install.packages(i, dependencies = T)
    require(i, character.only = T)
  }
}

df <- read_stata("Z://DesktopC//LUMSA//2//Data Mining//customers.dta")

df$gender <- as.factor(df$gender)
df$married <- as.factor(df$married)
df$catalogpromo <- round(df$catalogpromo)

df %<>% filter(monthnumb < 25)

df %>% group_by(hh_key) %>% summarize(n_ppl = n())

keyGroup<-df %>% group_by(hh_key) %>% summarize(items = sum(item),
                                                catalogs = sum(catalogpromo),
                                                retails = sum(retailpromo),
                                                prices = sum(pricepromo))

df %<>% group_by(monthnumb) %>% 
  summarize(avgItem = mean(item),
            avgCatal = mean(catalogpromo),
            avgRetail = mean(retailpromo)) %>% 
  mutate(year = ifelse(monthnumb < 13, 1, ifelse(monthnumb < 25, 2, ifelse(monthnumb < 37, 3, "   "))),
         fixedMon = c(seq(1,12,by=1),seq(1,12,by=1),seq(1,12,by=1)),
         winter = ifelse(fixedMon%%12==0, 1, ifelse(fixedMon%%11==0, 1, 0)))

p1 <- ggplot(data = df, aes(x = monthnumb, y = avgItem, group=1)) + 
      geom_line() +
      scale_x_continuous(name = "time [years]",breaks = c(0,12,24,36), labels = c(0,1,2,3)) +
      theme(plot.title = element_text(hjust = 0.5)) +
      labs(title = "Purchases",
       y = "av. number")

p2 <- ggplot(data = df, aes(x = monthnumb, y = avgCatal, group=1)) + 
      geom_line() +
      scale_x_continuous(name = "time [years]",breaks = c(0,12,24,36), labels = c(0,1,2,3)) +
      theme(plot.title = element_text(hjust = 0.5)) +
      labs(title = "Catalog promotions",
       y = "av. number")

p3 <- ggplot(data = df, aes(x = monthnumb, y = avgRetail, group=1)) + 
      geom_line() +
      scale_x_continuous(name = "time [years]",breaks = c(0,12,24,36), labels = c(0,1,2,3)) +
      theme(plot.title = element_text(hjust = 0.5)) +
      labs(title = "Retail promotions",
           y = "av. number")

grid.arrange(p1, p2, p3, ncol=1)


ggplot(data = df, aes(x = monthnumb, y = avgItem, group=as.factor(winter), color=as.factor(winter))) + 
  geom_boxplot() +
  scale_x_continuous(name = "Month",breaks = c(13,26), labels = c("Non Winter","Winter")) +
  theme() +
  labs(color = "Winter")

model5 <- hurdle(item ~ as.factor(monthnumb)+as.factor(retailpromo)+catalogpromo+as.factor(gender)+as.factor(married)+as.factor(income), data = customers, dist = "negbin")
AIC(model,model2,model3,model4,model5)
BIC(model,model2,model3,model4,model5)

evaluateModel <- function(model1, model2, data){
  TZero <- sum(data$item == 0)
  Pzero1 <- sum(predict(model1, type = "prob")[,1])
  Pzero2 <- sum(predict(model2, type = "prob")[,1])
  Pones1 <- sum(predict(model1, type = "prob")[,2])
  Pones2 <- sum(predict(model2, type = "prob")[,2])
  Pduos1 <- sum(predict(model1, type = "prob")[,3])
  Pduos2 <- sum(predict(model2, type = "prob")[,3])
  res = paste("True zeros:",as.numeric(TZero),
        "Predicted zeros from model-1:",as.numeric(Pzero1),
          "Predicted zeros from model-2 ",as.numeric(Pzero2),
           "Predicted ones from model-1:",as.numeric(Pones1),
           "Predicted ones from model-2:",as.numeric(Pones2),
           "Predicted duos from model-1:",as.numeric(Pduos1),
           "Predicted duos from model-2:",as.numeric(Pduos2),sep = "\n")
  return(cat(res))
}

expCoef <- exp(coef((model6)))
expCoef <- matrix(expCoef, ncol = 2)
rownames(expCoef) <- names(coef(model7))
colnames(expCoef) <- c("Count_model","Zero_hurdle_model")
expCoef


rootogram(model8, max = 30)
