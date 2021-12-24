pbcseq.forest <- pbcseq  %>%
                  mutate(id = as.integer(id),
                         age = as.double(day+(age*365.25)),
                         status = ifelse(status==0,0,status-1),
                         ascites = as.factor(ascites),
                         hepato = as.factor(hepato),
                         spiders = as.factor(spiders),
                         edema = as.factor(ifelse(edema==1,'serious',ifelse(edema==0.5,'evidence','no')))) %>%
                  group_by(id) %>%
                  do(tail(., n = 1)) %>%
                  ungroup() %>%
                  as.data.frame(.)

pbc.obj <- rfsrc(Surv(day, status) ~ ., pbcseq.forest, na.action = "na.impute", nimpute = 10)
print(pbc.obj)

surv.f <- as.formula(Surv(day, status) ~ sex + age + trt + bili)
pec.f <- as.formula(Hist(day,status) ~ sex + age + trt + bili)
pbcseq.forest.na <- na.omit(pbcseq.forest)

cox.obj <- coxph(surv.f, data = pbcseq.forest.na, x = TRUE)
rfsrc.obj <- rfsrc(surv.f,
                   pbcseq.forest.na,
                   ntree = 500,
                   splitrule = 'logrank',
                   nsplit = 14,
                   nodesize = 5)

prederror.pbc <- pec(list(cox.obj,rfsrc.obj), data = pbcseq.forest.na, formula = pec.f,
                     splitMethod = "bootcv", B = 50)
print(prederror.pbc)
plot(prederror.pbc)

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
  surv.f <- as.formula(Surv(day, status) ~ bili+ascites+protime+edema+albumin+ast+alk.phos+age+chol+platelet+spiders+sex+hepato+trt)
  pec.f <- as.formula(Hist(day,status) ~ bili+ascites+protime+edema+albumin+ast+alk.phos+age+chol+platelet+spiders+sex+hepato+trt)

  ## run cox/rfsrc models
  ## for illustration we use a small number of trees
  cox.obj <- coxph(surv.f, data = pbcseq.forest.na, x = TRUE)
  rfsrc.obj <- rfsrc(surv.f,
                     pbcseq.forest.na,
                     ntree = 500,
                     nsplit = 14,
                     nodesize = 5,
                     splitrule = 'logrank')

  ## compute bootstrap cross-validation estimate of expected Brier score
  ## see Mogensen, Ishwaran and Gerds (2012) Journal of Statistical Software
  set.seed(1234)
  prederror.pbc <- pec(list(cox.obj,rfsrc.obj), data = pbcseq.forest.na, formula = pec.f,
                       splitMethod = "bootcv", B = 50)
  print(prederror.pbc)
  plot(prederror.pbc)

  ## compute out-of-bag C-index for cox regression and compare to rfsrc
  rfsrc.obj <- rfsrc(surv.f,pbcseq.forest.na)

  cat("out-of-bag Cox Analysis ...", "\n")
  cox.err <- sapply(1:100, function(b) {
    if (b%%10 == 0) cat("cox bootstrap:", b, "\n")
    train <- sample(1:nrow(pbcseq.forest.na), nrow(pbcseq.forest.na), replace = TRUE)
    cox.obj <- tryCatch({coxph(surv.f, pbcseq.forest.na[train, ])}, error=function(ex){NULL})
    if (!is.null(cox.obj)) {
      get.cindex(pbcseq.forest.na$day[-train], pbcseq.forest.na$status[-train], predict(cox.obj, pbcseq.forest.na[-train, ]))
    } else NA
  })
  cat("\n\tOOB error rates\n\n")
  cat("\tRSF            : ", rfsrc.obj$err.rate[rfsrc.obj$ntree], "\n")
  cat("\tCox regression : ", mean(cox.err, na.rm = TRUE), "\n")
}




splitrule <- c("logrank", "bs.gradient","logrankscore")
nrep <- 100
ntree <- 500
err.rate <- matrix(0, 3, nrep)
names(err.rate) <- splitrule
for (j in 1:3) {
  for (k in 1:nrep) {
    err.rate[j,k] <- rfsrc(Surv(day,status) ~
                             bili+ascites+protime+edema+albumin+ast+alk.phos+age+chol+platelet+spiders+sex+hepato+trt,
                             pbcseq.forest,
                             ntree=ntree,
                             splitrule=splitrule[j])$err.rate[ntree]
  }
}
err.rate <- rbind(
  mean=apply(err.rate, 1, mean),
  std=apply(err.rate, 1, sd))

colnames(err.rate) <- splitrule
print(round(err.rate,3))



if (library("survival", logical.return = TRUE)) {
  ## use the pbc data from the survival package
  ## events are transplant (1) and death (2)
  pbcseq.forest.na$id <- NULL

  ## modified Gray's weighted log-rank splitting
  ## (equivalent to cause=c(1,1) and splitrule="logrankCR")
  pbc.cr <- rfsrc(Surv(day, status) ~
                    bili+ascites+protime+edema+albumin+ast+alk.phos+age+chol+platelet+spiders+sex+hepato+trt,
                  pbcseq.forest.na,
                  nsplit = 14,
                  nodesize = 5)

  ## log-rank cause-1 specific splitting and targeted VIMP for cause 1
  pbc.log1 <- rfsrc(Surv(day, status) ~
                      bili+ascites+protime+edema+albumin+ast+alk.phos+age+chol+platelet+spiders+sex+hepato+trt,
                    pbcseq.forest.na,
                    nsplit = 14,
                    nodesize = 5,
                    splitrule = "logrank", cause = c(1,0), importance = TRUE)

  ## log-rank cause-2 specific splitting and targeted VIMP for cause 2
  pbc.log2 <- rfsrc(Surv(day, status) ~
                      bili+ascites+protime+edema+albumin+ast+alk.phos+age+chol+platelet+spiders+sex+hepato+trt,
                    pbcseq.forest.na,
                    nsplit = 14,
                    nodesize = 5,
                    splitrule = "logrank", cause = c(0,1), importance = TRUE)

  ## extract VIMP from the log-rank forests: event-specific
  ## extract minimal depth from the Gray log-rank forest: non-event specific
  var.perf <- data.frame(md = max.subtree(pbc.cr)$order[, 1],
                         vimp1 = 100 * pbc.log1$importance,
                         vimp2 = 100 * pbc.log2$importance)
  print(var.perf[order(var.perf$md), ], digits = 3)
}
plot.survival(pbc.cr, cens.model = 'km', collapse = F)
plot(pbc.log1, verbose = TRUE, plots.one.page = T)




