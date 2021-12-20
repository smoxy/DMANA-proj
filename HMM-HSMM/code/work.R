###################################LOAD#########################################
load(url("https://github.com/smoxy/DMANA-proj/blob/main/HMM-HSMM/code/Examples_L31.RData?raw=true"))
source("https://raw.githubusercontent.com/MatteoFasulo/Rare-Earth/main/R/util/coreFunctions.R")
loadPackages(c('ggplot2','gamlss','tidyverse','tidyquant','magrittr',
               'tseries','MVN','mhsmm','doParallel','pacman','gridExtra'))
rm(returns)
rm(pollution)
rm(stock.names)
################################################################################

histDist(test, family = "BCT", nbins=50) #?
################################################################################
StockReturns %<>%
  mutate(Date=as.Date(StockReturns$Date,format="%d/%m/%Y"))

StockReturns2 <- StockReturns %>%
                  select(SP500,NASDAQ,ESTX50,FTSE,Date) %>%
                  filter(Date > "2003/01/01" & Date < "2016/06/23")

dfStruct <- function(dataframe){
  dates <- dataframe$Date
  SP500 <- dataframe$SP500
  NASDAQ <- dataframe$NASDAQ
  FTSE <- dataframe$FTSE
  ESTX50 <- dataframe$ESTX50

  name <- rep(c("SP500","NASDAQ","FTSE","ESTX50"),times=c(nrow(dataframe),nrow(dataframe),nrow(dataframe),nrow(dataframe)))
  dates <- rep(dates, times=4)
  value <- rep(c(SP500,NASDAQ,FTSE,ESTX50))

  df <- data.frame("name"=name,
                       "date"=as.Date(dates,format="%d/%m/%Y"),
                       "value"=value)
  return(df)
}

Stocks <- dfStruct(StockReturns2)

Stocks %>%
  ggplot(., aes(x=date, y=value, group=name))+
    geom_line(size=.7, col="gray13") +
    scale_x_date(date_labels = "%m-%y", date_breaks="4 months") +
    facet_wrap(~ name, ncol = 1, scale = "free_y") +
    labs(title = "ESTX50, FTSE, NASDAQ & SP500 Chart",
         subtitle = "Multiple Stocks",
         y = "Return",
         x = "") +
    theme_tq()

Stocks %>%
  group_by(name) %>%
  summarise("Mean"=round(mean(value),3),
            "Std. dev"=sd(value),
            "Skewness"=skewness(value),
            "Kurtosis"=PerformanceAnalytics::kurtosis(value, method = 'excess'),
            "Jarque–Bera test (p-value)"=paste(round(as.double(jarque.bera.test(value)$statistic),1)," (",as.double(jarque.bera.test(value)$p.value),")",sep=""))

qqnorm(StockReturns2$NASDAQ)
qqline(StockReturns2$NASDAQ)
qqnorm(StockReturns2$ESTX50)
qqline(StockReturns2$ESTX50)
mvn(StockReturns2[,-c(1,4,5)],mvnTest="mardia",multivariatePlot="qq")
###################################HMM##########################################
newDF <- StockReturns2[,-c(1,4,5)]
# initialization <- kmeans(newDF,3)
# kmeans.means <- initialization$centers
# sigma1 <- cov(newDF[initialization$cluster==1,])
# sigma2 <- cov(newDF[initialization$cluster==2,])
# sigma3 <- cov(newDF[initialization$cluster==3,])
#
# K <- 3
# start.val <- hmmspec(init = rep(1/K, K),
#                      trans = matrix(1/K, nrow = K, ncol = K),
#                      parms.emis = list(mu=list(kmeans.means[1,],kmeans.means[2,],kmeans.means[3,]),
#                                        sigma=list(sigma1,sigma2,sigma3)),
#                      dens.emis = dmvnorm.hsmm)
# mod.hmm.k2.v1 <- hmmfit(matrix(unlist(newDF),ncol=2),
#                         start.val, mstep = mstep.mvnorm)
# plot(newDF,col=mod.hmm.k2.v1$yhat)
#
# B <- 1000 # replicates
# mu.boot <- matrix(NA,6,B)
# sigma.boot <- matrix(NA,12,B)
# for (b in 1:B)
# {
#   true.par <- hmmspec(init = mod.hmm.k2.v1$model$init,
#                       trans = mod.hmm.k2.v1$model$transition,
#                       parms.emis = list(mu = mod.hmm.k2.v1$model$parms.emission$mu,
#                                         sigma = mod.hmm.k2.v1$model$parms.emission$sigma),
#                       dens.emis = dmvnorm.hsmm)
#   train <- simulate(true.par, nsim = nrow(newDF), seed = b, rand.emis = rmvnorm.hsmm)
#   mod.boot <- hmmfit(train, true.par, mstep = mstep.mvnorm)
#   mu.boot[,b] <- unlist(mod.boot$model$parms.emission$mu)
#   sigma.boot[,b] <- unlist(mod.boot$model$parms.emission$sigma)
# }
# apply(mu.boot,1,mean)
# apply(mu.boot,1,sd)
# apply(sigma.boot,1,mean)
# apply(sigma.boot,1,sd)
# t(matrix(apply(sigma.boot,1,mean),ncol = 6))

#################################BOOTSTRAP######################################
bootstrap.N.HMM = function(data, nclust, R){
  master = list()
  cores <- 1
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  master<-append(master,foreach(num = 2:nclust) %dopar% {
    require(mhsmm)

    hsmmAIC <- function(model){
      if (class(model)=='hmm'){
        m = model$K
      } else if (class(model)=='hsmm'){
        m = model$J
      }
      k = length(model$model$parms.emission)
      p = m^2 + (k*m) - 1
      logL = max(model$loglik)
      AIC = (2 * logL) + (2*p)
      return(AIC)
    }

    hsmmBIC <- function(model){
      if (class(model)=='hmm'){
        m = model$K
      } else if (class(model)=='hsmm'){
        m = model$J
      }
      nObs = length(model$yhat)
      k = length(model$model$parms.emission)
      p = (m^2 + (k*m)) - 1
      logL = max(model$loglik)
      BIC = (2 * logL) + (p*log(nObs))
      return(BIC)
    }

    initialization <- kmeans(data,num)
    kmeans.means <- initialization$centers
    km.means.list = list()
    sigma = list()
    for (n in 1:num){
      sigma[[n]] = cov(data[initialization$cluster==n,])
      km.means.list[[n]] = kmeans.means[n,]
    }
    K <- num
    start.val <- hmmspec(init = rep(1/K, K),
                         trans = matrix(1/K, nrow = K, ncol = K),
                         parms.emis = list(mu=km.means.list,
                                           sigma=sigma),
                         dens.emis = dmvnorm.hsmm)
    mod.hmm <- hmmfit(matrix(unlist(data),ncol=2),
                      start.val, mstep = mstep.mvnorm)
    mu.boot <- matrix(NA,K*2,R)
    sigma.boot <- matrix(NA,K*4,R)
    for (b in 1:R)
    {
      true.par <- hmmspec(init = mod.hmm$model$init,
                          trans = mod.hmm$model$transition,
                          parms.emis = list(mu = mod.hmm$model$parms.emission$mu,
                                            sigma = mod.hmm$model$parms.emission$sigma),
                          dens.emis = dmvnorm.hsmm)

      train <- simulate(true.par, nsim = 100, seed = b, rand.emis = rmvnorm.hsmm)

      mod.boot <- hmmfit(train, true.par, mstep = mstep.mvnorm)
      mu.boot[,b] <- unlist(mod.boot$model$parms.emission$mu)
      sigma.boot[,b] <- unlist(mod.boot$model$parms.emission$sigma)
    }



    return(list(k = num,
                model = mod.hmm,
                means.mean = t(matrix(apply(mu.boot,1,mean),ncol = K*2)),
                means.sd = t(matrix(apply(mu.boot,1,sd),ncol = K*2)),
                sigma.mean = t(matrix(apply(sigma.boot,1,mean),ncol = K*2**2)),
                sigma.sd = t(matrix(apply(sigma.boot,1,sd),ncol = K*2**2)),
                AIC=hsmmAIC(mod.hmm),
                BIC=hsmmBIC(mod.hmm)))

  },length(master))
  stopCluster(cl)
  return(master)
}

HMM <- bootstrap.N.HMM(data = newDF, nclust = 4, R = 50)
save(HMM,file="N-HMM-2-4.RData",safe=T)

################################################################################

ggplot(StockReturns2, aes(x=NASDAQ)) +
  geom_histogram(aes(y =..density..),
                 color="black", fill="white", binwidth=.1) +
  stat_function(fun = dnorm, args = list(mean = result[[1]]$means.mean[1],
                                         sd = result[[1]]$sigma.mean[1]), colour="red")+
  stat_function(fun = dnorm, args = list(mean = result[[1]]$means.mean[3],
                                         sd = result[[1]]$sigma.mean[5]), colour="green")

ggplot(StockReturns2, aes(x=NASDAQ)) +
  geom_histogram(aes(y =..density..),
                 color="black", fill="white", binwidth=.1) +
  stat_function(fun = dnorm, args = list(mean = result[[2]]$means.mean[1],
                                         sd = result[[2]]$sigma.mean[1]), colour="red")+
  stat_function(fun = dnorm, args = list(mean = result[[2]]$means.mean[3],
                                         sd = result[[2]]$sigma.mean[5]), colour="green")+
  stat_function(fun = dnorm, args = list(mean = result[[2]]$means.mean[5],
                                         sd = result[[2]]$sigma.mean[9]), colour="yellow")

ggplot(StockReturns2, aes(x=NASDAQ)) +
  geom_histogram(aes(y =..density..),
                 color="black", fill="white", binwidth=.1) +
  stat_function(fun = dnorm, args = list(mean = result[[3]]$model$model$parms.emission$mu[[1]][1],
                                         sd = result[[3]]$model$model$parms.emission$sigma[[1]][1]), colour="red")+
  stat_function(fun = dnorm, args = list(mean = result[[3]]$model$model$parms.emission$mu[[2]][1],
                                         sd = result[[3]]$model$model$parms.emission$sigma[[2]][1]), colour="green")+
  stat_function(fun = dnorm, args = list(mean = result[[3]]$model$model$parms.emission$mu[[3]][1],
                                         sd = result[[3]]$model$model$parms.emission$sigma[[3]][1]), colour="yellow")+
  stat_function(fun = dnorm, args = list(mean = result[[3]]$model$model$parms.emission$mu[[4]][1],
                                         sd = result[[3]]$model$model$parms.emission$sigma[[4]][1]), colour="blue")

ggplot(StockReturns2, aes(x=NASDAQ)) +
  geom_histogram(aes(y =..density..),
                 color="black", fill="white", binwidth=.1) +
  stat_function(fun = dnorm, args = list(mean = result[[5]]$means.mean[1],
                                         sd = result[[5]]$sigma.mean[1]), colour="red")+
  stat_function(fun = dnorm, args = list(mean = result[[5]]$means.mean[3],
                                         sd = result[[5]]$sigma.mean[5]), colour="green")+
  stat_function(fun = dnorm, args = list(mean = result[[5]]$means.mean[5],
                                         sd = result[[5]]$sigma.mean[9]), colour="yellow")+
  stat_function(fun = dnorm, args = list(mean = result[[5]]$means.mean[7],
                                         sd = result[[5]]$sigma.mean[13]), colour="blue")+
  stat_function(fun = dnorm, args = list(mean = result[[5]]$means.mean[9],
                                         sd = result[[5]]$sigma.mean[17]), colour="purple")+
  stat_function(fun = dnorm, args = list(mean = result[[5]]$means.mean[11],
                                         sd = result[[5]]$sigma.mean[21]), colour="orange")
##########

ggplot(StockReturns2, aes(x=ESTX50)) +
  geom_histogram(aes(y =..density..),
                 color="black", fill="white", binwidth=.1) +
  stat_function(fun = dnorm, args = list(mean = result[[1]]$means.mean[2],
                                         sd = result[[1]]$sigma.mean[4]), colour="red")+
  stat_function(fun = dnorm, args = list(mean = result[[1]]$means.mean[4],
                                         sd = result[[1]]$sigma.mean[8]), colour="green")

ggplot(StockReturns2, aes(x=ESTX50)) +
  geom_histogram(aes(y =..density..),
                 color="black", fill="white", binwidth=.1) +
  stat_function(fun = dnorm, args = list(mean = result[[2]]$model$model$parms.emission$mu[[1]][2],
                                         sd = result[[2]]$model$model$parms.emission$sigma[[1]][4]), colour="red")+
  stat_function(fun = dnorm, args = list(mean = result[[2]]$model$model$parms.emission$mu[[2]][2],
                                         sd = result[[2]]$model$model$parms.emission$sigma[[2]][4]), colour="green")+
  stat_function(fun = dnorm, args = list(mean = result[[2]]$model$model$parms.emission$mu[[3]][2],
                                         sd = result[[2]]$model$model$parms.emission$sigma[[3]][4]), colour="yellow")

ggplot(StockReturns2, aes(x=ESTX50)) +
  geom_histogram(aes(y =..density..),
                 color="black", fill="white", binwidth=.1) +
  stat_function(fun = dnorm, args = list(mean = result[[3]]$model$model$parms.emission$mu[[1]][2],
                                         sd = result[[3]]$model$model$parms.emission$sigma[[1]][4]), colour="red")+
  stat_function(fun = dnorm, args = list(mean = result[[3]]$model$model$parms.emission$mu[[2]][2],
                                         sd = result[[3]]$model$model$parms.emission$sigma[[2]][4]), colour="green")+
  stat_function(fun = dnorm, args = list(mean = result[[3]]$model$model$parms.emission$mu[[3]][2],
                                         sd = result[[3]]$model$model$parms.emission$sigma[[3]][4]), colour="yellow")+
  stat_function(fun = dnorm, args = list(mean = result[[3]]$model$model$parms.emission$mu[[4]][2],
                                         sd = result[[3]]$model$model$parms.emission$sigma[[4]][4]), colour="blue")

################################################################################
initialization <- kmeans(newDF,3)
kmeans.means <- initialization$centers
sigma1 <- cov(newDF[initialization$cluster==1,])
sigma2 <- cov(newDF[initialization$cluster==2,])
sigma3 <- cov(newDF[initialization$cluster==3,])

K <- 3
init <- rep(1/K, K)
P <- matrix(c(0, .5, .5, .5, 0, .5, .5, .5, 0), nrow = K)
B <- list(
  mu = list(kmeans.means[1,],kmeans.means[2,],kmeans.means[3,]),
  sigma = list(sigma1,sigma2,sigma3))

d <- list(lambda = c(10, 30, 60), shift = c(1, 1, 1), type = "poisson")

true.par <- hsmmspec(init=init,
                  transition = P,
                  parms.emis = B,
                  sojourn = d,
                  dens.emis = dmvnorm.hsmm,
                  mstep = mstep.mvnorm)

train <- simulate(true.par, nsim = nrow(newDF), seed = 1234, rand.emis = rmvnorm.hsmm)
model <- hsmmfit(train, true.par, mstep = mstep.mvnorm)


plot(model,xlim=c(0,200))
plot(newDF,col=model$yhat)
hist(StockReturns2$NASDAQ, breaks = 1000, freq=F)
lines(y=dnorm(seq(-10,10,length=1000),
              model$model$parms.emission$mu[[1]][1],
              model$model$parms.emission$sigma[[1]][1]),x=seq(-10,10,length=1000),col='red',lwd=2)
lines(y=dnorm(seq(-10,10,length=1000),
              model$model$parms.emission$mu[[2]][1],
              model$model$parms.emission$sigma[[2]][1]),x=seq(-10,10,length=1000),col='green',lwd=2)
lines(y=dnorm(seq(-10,10,length=1000),
              model$model$parms.emission$mu[[3]][1],
              model$model$parms.emission$sigma[[3]][1]),x=seq(-10,10,length=1000),col='yellow',lwd=2)
lines(y=dnorm(seq(-10,10,length=1000),
              model$model$parms.emission$mu[[4]][1],
              model$model$parms.emission$sigma[[4]][1]),x=seq(-10,10,length=1000),col='yellow',lwd=2)

hist(StockReturns2$ESTX50, breaks = 1000, freq=F)
lines(y=dnorm(seq(-10,10,length=1000),
              model$model$parms.emission$mu[[1]][2],
              model$model$parms.emission$sigma[[1]][4]),x=seq(-10,10,length=1000),col='red',lwd=2)
lines(y=dnorm(seq(-10,10,length=1000),
              model$model$parms.emission$mu[[2]][2],
              model$model$parms.emission$sigma[[2]][4]),x=seq(-10,10,length=1000),col='green',lwd=2)
lines(y=dnorm(seq(-10,10,length=1000),
              model$model$parms.emission$mu[[3]][2],
              model$model$parms.emission$sigma[[3]][4]),x=seq(-10,10,length=1000),col='yellow',lwd=2)



####################################N-HSMM-GAMMA################################
initialization <- kmeans(newDF,3)
kmeans.means <- initialization$centers
sigma1 <- cov(newDF[initialization$cluster==1,])
sigma2 <- cov(newDF[initialization$cluster==2,])
sigma3 <- cov(newDF[initialization$cluster==3,])

K <- 3
init <- rep(1/K, K)
P <- matrix(1 / (K - 1), ncol = K, nrow = K)
diag(P) <- 0
B <- list(
  mu = list(kmeans.means[1,],kmeans.means[2,],kmeans.means[3,]),
  sigma = list(sigma1,sigma2,sigma3))

d <- list(shape = c(2, 5, 7), scale = c(2, 2, 2), type = "gamma")

true.par <- hsmmspec(init=init,
                  transition = P,
                  parms.emis = B,
                  sojourn = d,
                  dens.emis = dmvnorm.hsmm,
                  mstep = mstep.mvnorm)

train <- simulate(true.par, nsim = nrow(newDF), seed = 1234, rand.emis = rmvnorm.hsmm)
model <- hsmmfit(train, true.par, mstep = mstep.mvnorm)

plot(model$loglik, type = "b", ylab = "Log-likelihood", xlab = "Iteration")
plot(model,xlim=c(0,200))
plot(newDF,col=model$yhat)
hist(StockReturns2$NASDAQ, breaks = 1000, freq=F)
lines(y=dnorm(seq(-10,10,length=1000),
              model$model$parms.emission$mu[[1]][1],
              model$model$parms.emission$sigma[[1]][1]),x=seq(-10,10,length=1000),col='red',lwd=2)
abline(h = 0, v = model$model$parms.emission$mu[[1]][1], col = "blue")
lines(y=dnorm(seq(-10,10,length=1000),
              model$model$parms.emission$mu[[2]][1],
              model$model$parms.emission$sigma[[2]][1]),x=seq(-10,10,length=1000),col='green',lwd=2)
abline(h = 0, v = model$model$parms.emission$mu[[2]][1], col = "blue")
lines(y=dnorm(seq(-10,10,length=1000),
              model$model$parms.emission$mu[[3]][1],
              model$model$parms.emission$sigma[[3]][1]),x=seq(-10,10,length=1000),col='yellow',lwd=2)
abline(h = 0, v = model$model$parms.emission$mu[[3]][1], col = "blue")

hist(StockReturns2$ESTX50, breaks = 1000, freq=F)
lines(y=dnorm(seq(-10,10,length=1000),
              model$model$parms.emission$mu[[1]][2],
              model$model$parms.emission$sigma[[1]][4]),x=seq(-10,10,length=1000),col='red',lwd=2)
lines(y=dnorm(seq(-10,10,length=1000),
              model$model$parms.emission$mu[[2]][2],
              model$model$parms.emission$sigma[[2]][4]),x=seq(-10,10,length=1000),col='green',lwd=2)
lines(y=dnorm(seq(-10,10,length=1000),
              model$model$parms.emission$mu[[3]][2],
              model$model$parms.emission$sigma[[3]][4]),x=seq(-10,10,length=1000),col='yellow',lwd=2)


####################################N-HSMM-POISSON##############################
bootstrap.N.HSMM.POISSON <- function(data, nclust, R) {
  master = list()
  cores <- ifelse(detectCores() > 6,
                  cores <- 3,
                  cores <- 3)
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  master <- append(master, foreach(num = 2:nclust) %dopar% {
    require(mhsmm)

    hsmmAIC <- function(model) {
      if (class(model) == 'hmm') {
        m = model$K
      } else if (class(model) == 'hsmm') {
        m = model$J
      }
      k = length(model$model$parms.emission)
      p = m ^ 2 + (k * m) - 1
      logL = max(model$loglik)
      AIC = (2 * logL) + (2 * p)
      return(AIC)
    }

    hsmmBIC <- function(model) {
      if (class(model) == 'hmm') {
        m = model$K
      } else if (class(model) == 'hsmm') {
        m = model$J
      }
      nObs = length(model$yhat)
      k = length(model$model$parms.emission)
      p = (m ^ 2 + (k * m)) - 1
      logL = max(model$loglik)
      BIC = (2 * logL) + (p * log(nObs))
      return(BIC)
    }

    initialization <- kmeans(data, num)
    kmeans.means <- initialization$centers
    km.means.list = list()
    sigma.list = list()

    for (n in 1:num) {
      sigma.list[[n]] = cov(data[initialization$cluster == n, ])
      km.means.list[[n]] = kmeans.means[n, ]
    }

    K <- num
    init <- rep(1 / K, K)
    P <- matrix(1 / (K - 1), ncol = K, nrow = K)
    diag(P) <- 0
    B <- list(mu = km.means.list,
              sigma = sigma.list)

    d <-list(
      lambda = sample(1:500, K),
      shift = sample(10:100, K),
      type = "poisson")


    start.val <- hsmmspec(
      init = init,
      transition = P,
      parms.emis = B,
      sojourn = d,
      dens.emis = dmvnorm.hsmm,
      mstep = mstep.mvnorm)

    mod.hsmm <- hsmmfit(matrix(unlist(data),ncol = 2), start.val, mstep = mstep.mvnorm)
    mu.boot <- matrix(NA, K * 2, R)
    sigma.boot <- matrix(NA, K * 4, R)
    for (b in 1:R){
      true.par <- hsmmspec(
        init = mod.hsmm$model$init,
        transition = mod.hsmm$model$transition,
        parms.emis = list(
          mu = mod.hsmm$model$parms.emission$mu,
          sigma = mod.hsmm$model$parms.emission$sigma),
        sojourn = list(
          lambda = mod.hsmm$model$sojourn$lambda,
          shift = mod.hsmm$model$sojourn$shift,
          type = mod.hsmm$model$sojourn$type),
        dens.emis = dmvnorm.hsmm,
        mstep = mstep.mvnorm)

      train <- simulate(true.par, nsim = 100, seed = b, rand.emis = rmvnorm.hsmm)

      mod.boot <- hsmmfit(train, true.par, mstep = mstep.mvnorm)
      mu.boot[, b] <- unlist(mod.boot$model$parms.emission$mu)
      sigma.boot[, b] <- unlist(mod.boot$model$parms.emission$sigma)
    }

    return(list(
      k = num,
      model = mod.hsmm,
      means.mean = t(matrix(apply(mu.boot, 1, mean), ncol = K * 2)),
      means.sd = t(matrix(apply(mu.boot, 1, sd), ncol = K * 2)),
      sigma.mean = t(matrix(apply(sigma.boot, 1, mean), ncol = K * 2 ** 2)),
      sigma.sd = t(matrix(apply(sigma.boot, 1, sd), ncol = K * 2 ** 2)),
      AIC = hsmmAIC(mod.hsmm),
      BIC = hsmmBIC(mod.hsmm)))
  }, length(master))

  stopCluster(cl)
  return(master)
}

HSMM <- bootstrap.N.HSMM.POISSON(
    data = newDF,
    nclust = 4,
    R = 50)

save(HSMM, file = "N-HSMM-2-4-Poisson.RData", safe = T)

################################################################################
bootstrap.N.HSMM.GAMMA <- function(data, states, R) {
  master = list()
  cores <- ifelse(detectCores() > 6,
                  cores <- 3,
                  cores <- 3)
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  master <- append(master, foreach(num = 2:states) %dopar% {
    require(mhsmm)

    hsmmAIC <- function(model) {
      if (class(model) == 'hmm') {
        m = model$K
      } else if (class(model) == 'hsmm') {
        m = model$J
      }
      k = length(model$model$parms.emission)
      p = m ^ 2 + (k * m) - 1
      logL = max(model$loglik)
      AIC = (2 * logL) + (2 * p)
      return(AIC)
    }

    hsmmBIC <- function(model) {
      if (class(model) == 'hmm') {
        m = model$K
      } else if (class(model) == 'hsmm') {
        m = model$J
      }
      nObs = length(model$yhat)
      k = length(model$model$parms.emission)
      p = (m ^ 2 + (k * m)) - 1
      logL = max(model$loglik)
      BIC = (2 * logL) + (p * log(nObs))
      return(BIC)
    }

    initialization <- kmeans(data, num)
    kmeans.means <- initialization$centers
    km.means.list = list()
    sigma.list = list()

    for (n in 1:num) {
      sigma.list[[n]] = cov(data[initialization$cluster == n, ])
      km.means.list[[n]] = kmeans.means[n, ]
    }

    K <- num
    P <- matrix(1 / (K - 1), ncol = K, nrow = K)
    diag(P) <- 0
    start.val <- hsmmspec(
      init = rep(1 / K, K),

      transition = P,

      parms.emis = list(mu = km.means.list,
                        sigma = sigma.list),

      sojourn = list(shape = as.vector(sample(1:100,K)),
                     scale = as.vector(sample(1:10,K)),
                     type = "gamma"),

      dens.emis = dmvnorm.hsmm,

      mstep = mstep.mvnorm)

    mod.hsmm <- hsmmfit(matrix(unlist(data),ncol = 2),
                        start.val, mstep = mstep.mvnorm, maxit = 10)
    mu.boot <- matrix(NA, K * 2, R)
    sigma.boot <- matrix(NA, K * 4, R)
    for (b in 1:R){

      true.par <- hsmmspec(
        init = mod.hsmm$model$init,
        transition = mod.hsmm$model$transition,
        parms.emis = list(
          mu = mod.hsmm$model$parms.emission$mu,
          sigma = mod.hsmm$model$parms.emission$sigma),
        sojourn = list(
          shape = mod.hsmm$model$sojourn$shape,
          scale = mod.hsmm$model$sojourn$scale,
          type = mod.hsmm$model$sojourn$type),
        dens.emis = dmvnorm.hsmm,
        mstep = mstep.mvnorm)

      try({train <- simulate(true.par, nsim = 15, seed = b, rand.emis = rmvnorm.hsmm)},
          silent = T)


      mod.boot <- hsmmfit(train, true.par, mstep = mstep.mvnorm, maxit = 10)
      mu.boot[, b] <- unlist(mod.boot$model$parms.emission$mu)
      sigma.boot[, b] <- unlist(mod.boot$model$parms.emission$sigma)
    }

    return(list(
      k = num,
      model = mod.hsmm,
      means.mean = t(matrix(apply(mu.boot, 1, mean), ncol = K * 2)),
      means.sd = t(matrix(apply(mu.boot, 1, sd), ncol = K * 2)),
      sigma.mean = t(matrix(apply(sigma.boot, 1, mean), ncol = K * 2 ** 2)),
      sigma.sd = t(matrix(apply(sigma.boot, 1, sd), ncol = K * 2 ** 2)),
      AIC = hsmmAIC(mod.hsmm),
      BIC = hsmmBIC(mod.hsmm)))
  }, length(master))

  stopCluster(cl)
  return(master)
}

GAMMA.HSMM <-
  bootstrap.N.HSMM.GAMMA(
    data = newDF,
    states = 4,
    R = 20)

save(GAMMA.HSMM, file = "N-HSMM-2-4-Gamma.RData", safe = T)
