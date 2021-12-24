library(LMest)
pbcseq4 <- pbcseq3 %>%
  group_by(id) %>%
  mutate(t = rep(c(1:5), times = length(unique(id)))) %>%
  ungroup()
as.data.frame(pbcseq4)

test<-lmest(responsesFormula = status ~ age + sex,
      latentFormula = NULL,
      data=as.data.frame(pbcseq8),
      index=c('id','t'),
      k=1:5,
      start=1,
      modSel="BIC",
      seed=1234)


mod1 <- lmestData(data = pbcseq8, id = "id", time="t",
          responsesFormula = lchol + lalbumin + last ~ age + sex)

modc <- lmestCont(responsesFormula = lbili+lalbumin+lalk.phos+lchol+lsgot+lplatelet+lprotime ~ NULL,
                  latentFormula = ~ age + sex | age + sex,
                  index = c("id", "t"),
                  data = as.data.frame(pbcseq4),
                  k = 1:5,
                  start=0,
                  modSel = "BIC",
                  output = TRUE,
                  out_se = TRUE,
                  seed=1234)

modc1 <- lmestCont(responsesFormula = lbili+lalbumin+lalk.phos+lchol+lsgot+lplatelet+lprotime ~ NULL,
                  latentFormula = ~ age + sex | age + sex,
                  index = c("id", "t"),
                  data = as.data.frame(pbcseq4),
                  k = 1:5,
                  start=1,
                  modSel = "BIC",
                  output = TRUE,
                  out_se = TRUE,
                  seed=1234)



plot(modc,what="modSel")
plot(modc, what="density")
plot(modc,what="density",components=c(1,2))
plot(modc,what="transitions")

plot(modc1,what="modSel")
plot(modc1,what="transitions")

semodc <- se(modc1)
TabBe <-cbind(modc1$Be, semodc$seBe, modc1$Be/semodc$seBe)
colnames(TabBe) <- c("estBe.1","estBe.2","estBe.3", "s.e.Be.1","s.e.Be.2","s.e.Be.3","t-test.1","t-test.2","t-test.3")
round(TabBe,3)

out <- lmestSearch(responsesFormula =  lbili+lalbumin+lalk.phos+lchol+lsgot+lplatelet+lprotime ~ NULL,
                   latentFormula = ~ age + sex | age + sex,
                   index = c("id","t"),
                   data = as.data.frame(pbcseq4),
                   version ="continuous",
                   k = 1:7,
                   nrep = 10,
                   out_se = TRUE,
                   seed = 1234)




mod6 <- out$out.single[[6]]
plot(mod6, what="marginal")
plot(mod6,what="transitions")

dec <- lmestDecoding(mod6)

pbcseq %<>%
  mutate(id = as.integer(id),
         age = as.double(age+(day/365.25)),
         ascites = as.factor(ascites),
         hepato = as.factor(hepato),
         spiders = as.factor(spiders),
         edema = ifelse(edema==1,'serious',ifelse(edema==0.5,'evidence','no')))
df1 <- mice(pbcseq[,-c(8,9,10)], m=10, maxit = 50, method = 'pmm', seed = 1234)
df2 <- mice(pbcseq[,c(8,9,10)], m=10, maxit = 50, method = 'logreg', seed = 1234)
pbcseq2 <- complete(df1,1)
pbcseq2 <- cbind(pbcseq2,complete(df2,1))
pbcseq2 %<>%
  group_by(id) %>%
  mutate(sex = ifelse(sex=='f',1,0)) %>%
  filter(n()>= 5) %>%
  do(head(., n = 5)) %>%
  ungroup() %>%
  mutate(t = rep(c(1:5), times = length(unique(id))),
         lbili = log(bili),
         lalbumin = log(albumin),
         lalk.phos = log(alk.phos),
         lchol = log(chol),
         lsgot = log(ast),
         lplatelet = log(platelet),
         lprotime = log(protime)) %>%
  as.data.frame(.)


longMatrix <- with(pbcseq4, long2matrices(id = id, X = cbind(sex,
                                                             age),
                                                             Y = status))
S <- longMatrix$YY
X <- longMatrix$XX
X1 <- X[, 1, ]
TT <- 5
X2 <- X[, 2:TT, ]
colnames(X1) <- c("sex", "age")
dimnames(X2)[[3]] <- c("sex", "age")
modCustom <- est_lm_cov_latent(S = S, X1 = X1, X2 = X2, k = 5, start = 0,
                          param = "multilogit", fort = TRUE, output = TRUE, out_se=TRUE)

modCustom <- lmestData(data = pbcseq4, id = "id", time="t",
                  responsesFormula = status ~ lbili+lalbumin+lalk.phos+lchol+lsgot+lplatelet+lprotime,
                  latentFormula = ~ age + sex | age + sex,
                  check.names=TRUE)


modCustomFit <- lmestCont(responsesFormula = lbili+lalbumin+lalk.phos+lchol+lsgot+lplatelet+lprotime ~ NULL,
                  latentFormula = ~ age + sex | age + sex,
                  index = c("id", "t"),
                  data = pbcseq2,
                  k = 1:8,
                  start=1,
                  modSel = "BIC",
                  modBasic = 0,
                  paramLatent = "multilogit",
                  output = TRUE,
                  out_se = TRUE,
                  fort = TRUE,
                  seed=1234)
plot(modCustomFit)


withNa <- pbcseq %>%
  mutate(id = as.integer(id),
         age = as.double(age+(day/365.25)),
         ascites = as.factor(ascites),
         hepato = as.factor(hepato),
         spiders = as.factor(spiders),
         edema = ifelse(edema==1,'serious',ifelse(edema==0.5,'evidence','no')))
withNa %<>%
  group_by(id) %>%
  mutate(sex = ifelse(sex=='f',1,0)) %>%
  filter(n()>= 5) %>%
  do(head(., n = 5)) %>%
  ungroup() %>%
  mutate(t = rep(c(1:5), times = length(unique(id))),
         lbili = log(bili),
         lalbumin = log(albumin),
         lalk.phos = log(alk.phos),
         lchol = log(chol),
         lsgot = log(ast),
         lplatelet = log(platelet),
         lprotime = log(protime)) %>%
  as.data.frame(.)


################################################
out <- lmestSearch(responsesFormula = lbili+lalbumin+lalk.phos+lchol+lsgot+lplatelet+lprotime ~ NULL,
                   index = c("id","times"),
                   data = pbcseq10,
                   version ="continuous",
                   k = 2:8,
                   nrep = 50,
                   out_se = TRUE,
                   seed = 1234) #6

summary(out) # bic e aic di tutti i modelli con la ricerca mista

pbcseq10 <- pbcseq %>% group_by(id) %>% mutate(times = row_number(),
                                              id = as.integer(id),
                                              age = as.double(age+(day/365.25)),
                                              sex = ifelse(sex=='f',1,0),
                                              lbili = log(bili),
                                              lalbumin = log(albumin),
                                              lalk.phos = log(alk.phos),
                                              lchol = log(chol),
                                              lsgot = log(ast),
                                              lplatelet = log(platelet),
                                              lprotime = log(protime),
                                              trt = as.factor(trt))  %>%
                                              filter(n()>= 5) %>%
                                              do(head(., n = 5)) %>%
                                              ungroup() %>%
                                              as.data.frame(.)

# modello migliore deterministico
mod.with.na.det <- lmestCont(responsesFormula = lbili+lalbumin+lalk.phos+lchol+lsgot+lplatelet+lprotime ~ NULL,
                          latentFormula = ~ sex + age + trt | sex + age + trt,
                          index = c("id", "times"),
                          data = pbcseq10,
                          k = 5,
                          start=0,
                          modSel = "BIC",
                          modBasic = 0,
                          paramLatent = "multilogit",
                          output = TRUE,
                          out_se = TRUE,
                          fort = TRUE,
                          tol = 10^-10,
                          seed=1234)

# modello migliore stocastico
mod.with.na.stoc <- lmestCont(responsesFormula = lbili+lalbumin+lalk.phos+lchol+lsgot+lplatelet+lprotime ~ NULL,
                              latentFormula = ~ age + sex + trt | age + sex + trt,
                              index = c("id", "times"),
                              data = pbcseq10,
                              k = 6,
                              start=1,
                              modSel = "BIC",
                              modBasic = 0,
                              paramLatent = "multilogit",
                              output = TRUE,
                              out_se = TRUE,
                              fort = TRUE,
                              tol = 10^-10,
                              seed=1234)

mod.with.na.stoc[["Mu"]] #medie condizionate
mod.with.na.stoc[["Si"]] #var e cov
colMeans(mod.with.na.stoc[["Piv"]],dims=1) #prob iniziali

out <- with(withNa,long2matrices(id = id,X = cbind(sex,age,trt),Y = status))
X <- out$XX
TT <- mod.with.na.stoc$k
X1 <- X[, 1, ]
colnames(X1) <- c('sex','age','trt')
ind1 <- (X1[, 'sex'] == 1 | X1[, 'sex'] == 0)
PI1 <- round(apply(mod.with.na.stoc$PI[ , , ind1, 2:TT], c(1, 2), mean), 4) #transition probs
PI1

boot.det<-bootstrap(mod.with.na.stoc,B=100,seed = 1234)

mod.with.na.stoc[["Be"]]
mod.with.na.stoc[["Ga"]]

ind1 <- (X1[, 'sex'] == 1)
piv1 <- round(colMeans(mod.with.na.stoc$Piv[ind1, ]), 4) #init probs
piv1
PI1 <- round(apply(mod.with.na.stoc$PI[ , , ind1, 2:TT], c(1, 2), mean), 4) #transition probs
PI1

ind1 <- (X1[, 'sex'] == 0)
piv1 <- round(colMeans(mod.with.na.stoc$Piv[ind1, ]), 4) #init probs
piv1
PI1 <- round(apply(mod.with.na.stoc$PI[ , , ind1, 2:TT], c(1, 2), mean), 4) #transition probs
PI1

ind1 <- (X1[, 'age'] > 52)
piv1 <- round(colMeans(mod.with.na.stoc$Piv[ind1, ]), 4) #init probs
piv1
PI1 <- round(apply(mod.with.na.stoc$PI[ , , ind1, 2:TT], c(1, 2), mean), 4) #transition probs
PI1

ind1 <- (X1[, 'sex'] <= 52)
piv1 <- round(colMeans(mod.with.na.stoc$Piv[ind1, ]), 4) #init probs
piv1
PI1 <- round(apply(mod.with.na.stoc$PI[ , , ind1, 2:TT], c(1, 2), mean), 4) #transition probs
PI1

S <- out$YY
X <- out$XX
X1 <- X[, 1, ]
TT <- 5
X2 <- X[, 2:TT, ]
out <- lmestDecoding(test.mod)
head(out$Ug)



