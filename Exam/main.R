load(url("https://github.com/smoxy/DMANA-proj/blob/main/Exam/DatiPresentazione2021.RData?raw=true"))
source("https://raw.githubusercontent.com/MatteoFasulo/Rare-Earth/main/R/util/coreFunctions.R")
loadPackages(c('ggplot2','survival','mice','magrittr','tidyverse','gamlss'))
rm(BancaItalia)
rm(data_criminal_sim)
rm(data_SRHS_long)
rm(film90)
rm(ibmspko)
rm(PSIDlong)
rm(qgdp)
rm(tenstocks)
setwd('Z:/DesktopC/LUMSA/2/Data Mining/Exam/')
################################################################################
source("https://raw.githubusercontent.com/Silvia-Pand/HMContMiss/main/lmbasic.cont.MISS.R")
source("https://raw.githubusercontent.com/Silvia-Pand/HMContMiss/main/lmcovlatent.cont.MISS.R")
source("https://raw.githubusercontent.com/Silvia-Pand/HMContMiss/main/bootstrap.MISS.R")
source("https://raw.githubusercontent.com/Silvia-Pand/HMContMiss/main/lk_comp_cont_MISS.R")
source("https://raw.githubusercontent.com/Silvia-Pand/HMContMiss/main/lk_comp_latent_cont_MISS.R")
source("https://raw.githubusercontent.com/Silvia-Pand/HMContMiss/main/lk_obs_latent_cont_MISS.R")
source("https://raw.githubusercontent.com/Silvia-Pand/HMContMiss/main/prob_multilogit.R")
source("https://raw.githubusercontent.com/Silvia-Pand/HMContMiss/main/est_multilogit.R")
source("https://raw.githubusercontent.com/Silvia-Pand/HMContMiss/main/prob_post_cov_cont.R")
############################## CLEANING ########################################
pbcseq %<>%
            mutate(id = as.integer(id),
                   age = as.integer(round(age)),
                   trt = ifelse(trt==1,'D-penicillmain','placebo'),
                   ascites = as.factor(ascites),
                   hepato = as.factor(hepato),
                   spiders = as.factor(spiders),
                   edema = ifelse(edema==1,'serious',ifelse(edema==0.5,'evidence','no')))



apply(is.na(pbcseq), 2, which)
################################### MICE #######################################
md.pattern(pbcseq)
numVar <- c('chol','alk.phos','platelet')
binVar <-  c('ascites','hepato','spiders')
df1 <- mice(pbcseq[,-c(8,9,10)], m=10, maxit = 50, method = 'pmm', seed = 1234)
df2 <- mice(pbcseq[,c(8,9,10)], m=10, maxit = 50, method = 'logreg', seed = 1234)

densityplot(df1)

pbcseq2 <- complete(df1,1)
pbcseq2 <- cbind(pbcseq2,complete(df2,1))

pbcseq2 %<>%
  mutate(lbili = log(bili),
         lalbumin = log(albumin),
         lalk.phos = log(alk.phos),
         lchol = log(chol),
         lsgot = log(ast),
         lplatelet = log(platelet),
         lprotime = log(protime))

################################### CLEANED DATA ###############################

pbcseq3 <- pbcseq2 %>% group_by(id) %>% filter(n()>= 5) %>% do(head(., n = 5)) %>% ungroup()
apply(is.na(pbcseq3), 2, which)

################################################################################

#pbcseq %>% group_by(id) %>% mutate('AXA' = log(albumin)) %>% slice(1) %>% summary(.)

#pbcseq6 <- pbcseq5
#pbcseq6['id'] <- rep(c(1:length(unique(pbcseq5$id))), each = 5)

################################################################################
# Create the start-stop-event triplet needed for coxph
first <- with(pbcseq5, c(TRUE, diff(id) !=0)) #first id for each subject
last  <- c(first[-1], TRUE)  #last id

time1 <- with(pbcseq5, ifelse(first, 0, day))
time2 <- with(pbcseq5, ifelse(last,  futime, c(day[-1], 0)))
event <- with(pbcseq5, ifelse(last,  status, 0))

survival <- Surv(time1, time2, event)
################################################################################
#pbcseq5 <- pdata.frame(pbcseq4, index = c("id","day"))
#pbcseq6 <- make.pbalanced(pbcseq5, balance.type = 'shared.times')
################################################################################
pbcseq4['transplant'] <- rep(0, nrow(pbcseq4))
pbcseq4[which(pbcseq4$status==1),'transplant'] <- 1
pbcseq5['dead'] <- rep(0, nrow(pbcseq5))
pbcseq5[which(pbcseq5$status==2),'dead'] <- 1

pbcseq4[which(pbcseq4$status==1),'status'] <- 0
pbcseq4[which(pbcseq4$status==2),'status'] <- 1
################################################################################
histDist(status, family = ZIP(), data=pbcseq5)
histDist(pbcseq5$status, family = ZIP())
model <- gamlss(status ~ ., family=PO(), data = pbcseq5)
summary(model)
plot(model)
model <- gamlss(dead ~ trt + sex + age + edema + bili + alk.phos+ albumin + platelet, family=BI(), data = pbcseq5)
plot(model)
summary(model)

fit <- coxph(Surv(futime,dead) ~ log(age) + hepato + strata(sex) + strata(chol>300), data=pbcseq5)
fitnew <- survfit(fit, conf.type = "plain")
ggsurvplot(fitnew, data = newdata,
           conf.int = TRUE, censor= FALSE,
           ggtheme = theme_minimal())
