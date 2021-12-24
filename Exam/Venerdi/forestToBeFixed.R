splitrule <- c("logrank", "bs.gradient",
               "logrankscore")

nrep <- 100
ntree <- 500
err.rate <- matrix(0, 3, nrep)
names(err.rate) <- splitrule
for (j in 1:3) {
  for (k in 1:nrep) {
    err.rate[j,k] <- rfsrc(Surv(day,status) ~ age+sex,
                           binaryStatus,
                           ntree=ntree,
                           splitrule=splitrule[j])$err.rate[ntree]
  }
}
err.rate <- rbind(
  mean=apply(err.rate, 1, mean),
  std=apply(err.rate, 1, sd))

colnames(err.rate) <- splitrule
print(round(err.rate,3))


binaryStatus %<>%
  mutate(id = as.integer(id),
         age = as.integer(round(age)),
         trt = as.factor(trt),
         ascites = as.factor(ascites),
         hepato = as.factor(hepato),
         spiders = as.factor(spiders),
         edema = as.factor(edema),
         sex = as.factor(sex),
         status = ifelse(status==1,0,ifelse(status==2,1,0)))

#pbc.obj2 <- rfsrc(Surv(day, status) ~ ., pbcseq,
#                  nsplit = 10, na.action = "na.impute", nimpute = 10)

if (library("survival", logical.return = TRUE)
    & library("pec", logical.return = TRUE)
    & library("prodlim", logical.return = TRUE))

{
  ##prediction function required for pec
  predictSurvProb.rfsrc <- function(object, newdata, times, ...){
    ptemp <- predict(object,newdata=newdata,...)$survival
    pos <- sindex(jump.times = object$time.interest, eval.times = times)
    p <- cbind(1,ptemp)[, pos + 1]
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
      stop("Prediction failed")
    p
  }

  require(randomForestSRC)

  ## data, formula specifications
  surv.f <- as.formula(Surv(day, status) ~ age+sex)
  pec.f <- as.formula(Hist(day,status) ~ age+sex)

  ## run cox/rfsrc models
  ## for illustration we use a small number of trees
  cox.obj <- coxph(surv.f, data = binaryStatus, x = TRUE)
  rfsrc.obj <- rfsrc(surv.f, binaryStatus, ntree = 150)

  ## compute bootstrap cross-validation estimate of expected Brier score
  ## see Mogensen, Ishwaran and Gerds (2012) Journal of Statistical Software
  set.seed(17743)
  prederror.pbc <- pec(list(cox.obj,rfsrc.obj), data = binaryStatus, formula = pec.f,
                       splitMethod = "bootcv", B = 50)
  print(prederror.pbc)
  plot(prederror.pbc)

  ## compute out-of-bag C-index for cox regression and compare to rfsrc
  rfsrc.obj <- rfsrc(surv.f, binaryStatus)
  cat("out-of-bag Cox Analysis ...", "\n")
  cox.err <- sapply(1:100, function(b) {
    if (b%%10 == 0) cat("cox bootstrap:", b, "\n")
    train <- sample(1:nrow(binaryStatus), nrow(binaryStatus), replace = TRUE)
    cox.obj <- tryCatch({coxph(surv.f, binaryStatus[train, ])}, error=function(ex){NULL})
    if (!is.null(cox.obj)) {
      get.cindex(binaryStatus$day[-train], binaryStatus$status[-train], predict(cox.obj, binaryStatus[-train, ]))
    } else NA
  })
  cat("\n\tOOB error rates\n\n")
  cat("\tRSF            : ", rfsrc.obj$err.rate[rfsrc.obj$ntree], "\n")
  cat("\tCox regression : ", mean(cox.err, na.rm = TRUE), "\n")
}

################################################################################
gg_dta <- gg_survival(interval = "years",
                      censor = "status",
                      by = "trt",
                      data = as.data.frame(pbcseq.forest.train),
                      conf.int = .95)
plot(gg_dta) +
  labs(y = "Survival Probability",
       x = "Observation Time (years)",
       color = "trt", fill = "trt")+
  theme(legend.position = c(.2,.2))+
  coord_cartesian(y = c(0,1.01))



# C'è un outlier, proviamo a rimuoverlo
which.max(pbcseq.forest.train$day)

toBeRemoved = pbcseq.forest.train[which.max(pbcseq.forest.train$day),'id']
pbcseq.forest.train <- pbcseq.forest.train[-c(which(pbcseq.forest.train$id==toBeRemoved)),]

gg_dta <- gg_survival(interval = "years",
                      censor = "status",
                      by = "trt",
                      data = as.data.frame(pbcseq.forest.train),
                      conf.int = .95)
plot(gg_dta) +
  labs(y = "Survival Probability",
       x = "Observation Time (years)",
       color = "trt", fill = "trt")+
  theme(legend.position = c(.2,.2))+
  coord_cartesian(y = c(0,1.01))

plot(gg_dta, type="cum_haz") +
  labs(y = "Cumulative Hazard",
       x = "Observation Time (years)",
       color = "trt", fill = "trt")+
  theme(legend.position = c(.2,.8))

pbcseq.forest.train$bili_grp <- cut(pbcseq.forest.train$bili,breaks = c(0, .8, 1.3, 3.4,max(pbcseq.forest.train$bili)))
# plot the gg_survival object directly
plot(gg_survival(interval = "years",censor = "status",
                 by = "bili_grp", data = pbcseq.forest.train),
     error = "none") +
  labs(y = "Survival Probability",
       x = "Observation Time (years)",
       color = "Bilirubin")

cls <- sapply(pbcseq.forest.train, class)
labels <- c("ID",
            "AllTime",
            "Status (0 = censor, 1 = death)",
            "Treament (1 = DPCA, 0 = Placebo)",
            "Age (years)",
            "Sex (f = Female, m = Male)",
            "Edema (0, 0.5, 1)",
            "Serum Bilirubin (mg/dl)",
            "Serum Cholesterol (mg/dl)",
            "Albumin (gm/dl)",
            "Alkaline Phosphatase (U/liter)",
            "SGOT (U/ml)",
            "Platelets per cubic ml/1000",
            "Prothrombin time (sec)",
            "Histologic Stage",
            "Presence of Ascites",
            "Presence of Hepatomegaly",
            "Presence of Spiders",
            "lbili",
            "lalbumin",
            "lalk.phos",
            "lchol",
            "lsgot",
            "lplatelet",
            "lprotime",
            "Time (years)",
            "BiliGRP",
            "NA",
            "NA")

dta.labs <- data.frame(cbind(names = colnames(pbcseq.forest.train), label = labels, type = cls))
# Put the "years" variable on top.
dta.labs <- rbind(dta.labs[nrow(dta.labs),], dta.labs[-nrow(dta.labs),])
st.labs <- as.character(dta.labs$label)

pbcseq.forest.train$years50 <- cut(pbcseq.forest.train$age,breaks = c(min(pbcseq.forest.train$age),50,max(pbcseq.forest.train$age)))
# Grow and store the random survival forest
# Use random splitting (nsplit = 10) and impute
rfsrc_pbc <- rfsrc(Surv(years, status) ~ age+sex+edema+spiders+hepato+lbili+lalbumin+lalk.phos+lchol+lsgot+lplatelet+lprotime+trt,
                   data = pbcseq.forest.train,
                   ntree = 600,
                   importance=TRUE,
                   seed=1234)

# Print the forest summary
rfsrc_pbc
plot(gg_rfsrc(rfsrc_pbc,oob=TRUE,by='years50'))
plot(gg_rfsrc(rfsrc_pbc,oob=TRUE,by='trt'))


rfsrc_pbc_test <- predict(rfsrc_pbc,
                          newdata = pbcseq.forest.test)
plot(gg_rfsrc(rfsrc_pbc_test,oob=TRUE,by='sex'))


gg_dta <- vimp(rfsrc_pbc)
plot(gg_dta, lbls = st.labs)

varsel_pbc <- var.select(rfsrc_pbc)
plot(gg_minimal_depth(varsel_pbc))
plot(gg_minimal_vimp(varsel_pbc))

xvar <- varsel_pbc$topvars
ggrf <- gg_variable(rfsrc_pbc, time = c(1, 3),
                    time.labels = c("1 Year", "3 Years"))
plot(ggrf, xvar = "lbili", alpha = .3) +
  labs(y = "Survival", x = st.labs["lbili"]) +
  theme(legend.position = "right") +
  coord_cartesian(y = c(-.01,1.01))

xvar.cat <- c("edema", "stage")
xvar <- xvar[-which(xvar %in% xvar.cat)]

# plot the next 5 continuous variable dependence plots.
plot(ggrf, xvar = xvar[2:6], panel = TRUE,
     se = FALSE, alpha = .3,
     method = "glm", formula = y~poly(x,2)) +
  labs(y = "Survival") +
  theme(legend.position = "none") +
  coord_cartesian(y = c(-.01,1.01))

interaction_pbc <- find.interaction(rfsrc_pbc)
ggint <- gg_interaction(interaction_pbc)
plot(ggint, xvar = xvar) +
  labs(y = "Interactive Minimal Depth") +
  theme(legend.position = "none")






