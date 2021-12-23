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

out <- lmestSearch(responsesFormula =  lbili+lalbumin+lalk.phos+lchol+lsgot+lplatelet+lprotime ~ NULL,
                   latentFormula = ~ age + sex | age + sex,
                   index = c("id","t"),
                   data = as.data.frame(pbcseq4),
                   version ="continuous",
                   k = 1:5,
                   nrep = 50,
                   out_se = TRUE,
                   seed = 1234)
plot(out$out.single[[5]],what="transitions")
