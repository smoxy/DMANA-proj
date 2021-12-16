###################################LOAD#########################################
load(url("https://github.com/smoxy/DMANA-proj/blob/main/HMM-HSMM/code/Examples_L31.RData?raw=true"))
source("https://raw.githubusercontent.com/MatteoFasulo/Rare-Earth/main/R/util/coreFunctions.R")
loadPackages(c('ggplot2','gamlss','tidyverse','tidyquant','magrittr','tseries','MVN','mhsmm','doParallel'))
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
    scale_x_date(date_labels = "%Y", date_breaks="2 year") +
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
initialization <- kmeans(newDF,3)
kmeans.means <- initialization$centers
sigma1 <- cov(newDF[initialization$cluster==1,])
sigma2 <- cov(newDF[initialization$cluster==2,])
sigma3 <- cov(newDF[initialization$cluster==3,])

K <- 3
start.val <- hmmspec(init = rep(1/K, K),
                     trans = matrix(1/K, nrow = K, ncol = K),
                     parms.emis = list(mu=list(kmeans.means[1,],kmeans.means[2,],kmeans.means[3,]),
                                       sigma=list(sigma1,sigma2,sigma3)),
                     dens.emis = dmvnorm.hsmm)
mod.hmm.k2.v1 <- hmmfit(matrix(unlist(newDF),ncol=2),
                        start.val, mstep = mstep.mvnorm)
plot(newDF,col=mod.hmm.k2.v1$yhat)

B <- 1000 # replicates
mu.boot <- matrix(NA,6,B)
sigma.boot <- matrix(NA,12,B)
for (b in 1:B)
{
  true.par <- hmmspec(init = mod.hmm.k2.v1$model$init,
                      trans = mod.hmm.k2.v1$model$transition,
                      parms.emis = list(mu = mod.hmm.k2.v1$model$parms.emission$mu,
                                        sigma = mod.hmm.k2.v1$model$parms.emission$sigma),
                      dens.emis = dmvnorm.hsmm)
  train <- simulate(true.par, nsim = nrow(newDF), seed = b, rand.emis = rmvnorm.hsmm)
  mod.boot <- hmmfit(train, true.par, mstep = mstep.mvnorm)
  mu.boot[,b] <- unlist(mod.boot$model$parms.emission$mu)
  sigma.boot[,b] <- unlist(mod.boot$model$parms.emission$sigma)
}
apply(mu.boot,1,mean)
apply(mu.boot,1,sd)
apply(sigma.boot,1,mean)
apply(sigma.boot,1,sd)
t(matrix(apply(sigma.boot,1,mean),ncol = 6))

#################################BOOTSTRAP######################################

bootstrap = function(data, nclust, R){
  master = list()
  registerDoParallel(5)
  master<-append(master,foreach(num = 2:nclust) %dopar% {
    require(mhsmm)

    hsmmAIC <- function(model){
      m = model$K
      k = length(model$model$parms.emission)
      p = m^2 + (k*m) - 1
      logL = max(model$loglik)
      AIC = (-2 * logL) + (2*p)
      return(AIC)
    }

    hsmmBIC <- function(model){
      nObs = length(model$yhat)
      m = model$K
      k = length(model$model$parms.emission)
      p = (m^2 + (k*m)) - 1
      logL = max(model$loglik)
      BIC = (-2 * logL) + (p*log(nObs))
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
      train <- simulate(true.par, nsim = nrow(data), seed = b, rand.emis = rmvnorm.hsmm)
      mod.boot <- hmmfit(train, true.par, mstep = mstep.mvnorm)
      mu.boot[,b] <- unlist(mod.boot$model$parms.emission$mu)
      sigma.boot[,b] <- unlist(mod.boot$model$parms.emission$sigma)
    }

    print(paste("Finished state: ",num,sep=""))

    return(list(k = num,
                model = mod.hmm,
                means.mean = t(matrix(apply(mu.boot,1,mean),ncol = K*2)),
                means.sd = t(matrix(apply(mu.boot,1,sd),ncol = K*2)),
                sigma.mean = t(matrix(apply(sigma.boot,1,mean),ncol = K*2**2)),
                sigma.sd = t(matrix(apply(sigma.boot,1,sd),ncol = K*2**2)),
                AIC=hsmmAIC(mod.hmm),
                BIC=hsmmBIC(mod.hmm)))

  },length(master))
}

result <- bootstrap(data = newDF, nclust = 3, R = 100)
stopImplicitCluster()
save(result,file="HMM-Result.RData",safe=T)

################################################################################

ggplot(StockReturns2, aes(x=NASDAQ)) +
  geom_histogram(aes(y =..density..),
                 color="black", fill="white", binwidth=.1) +
  stat_function(fun = dnorm, args = list(mean = -0.15234507, sd = 8.3006376), colour="red")+
  stat_function(fun = dnorm, args = list(mean = 0.03785593 , sd = 1.6453569), colour="green")+
  stat_function(fun = dnorm, args = list(mean = 0.14786969  , sd = 0.5317181), colour="blue")

ggplot(StockReturns2, aes(x=ESTX50)) +
  geom_histogram(aes(y =..density..),
                 color="black", fill="white", binwidth=.1) +
  stat_function(fun = dnorm, args = list(mean = -0.11874290, sd = 8.9762954), colour="red")+
  stat_function(fun = dnorm, args = list(mean = -0.06172499, sd = 2.1158288), colour="green")+
  stat_function(fun = dnorm, args = list(mean = 0.11453419, sd = 0.6186066), colour="blue")

################################################################################
initialization <- kmeans(newDF,3)
kmeans.means <- initialization$centers
sigma1 <- cov(newDF[initialization$cluster==1,])
sigma2 <- cov(newDF[initialization$cluster==2,])
sigma3 <- cov(newDF[initialization$cluster==3,])

J <- 3
init <- rep(1/J, J)
P <- matrix(c(0, .5, .5, .5, 0, .5, .5, .5, 0), nrow = J)
B <- list(
  mu = list(kmeans.means[1,],kmeans.means[2,],kmeans.means[3,]),
  sigma = list(sigma1,sigma2,sigma3))

d <- list(lambda = c(10, 30, 60), shift = c(10, 100, 30), type = "poisson")

model <- hsmmspec(init=init,
                  transition = P,
                  parms.emis = B,
                  sojourn = d,
                  dens.emis = dmvnorm.hsmm,
                  mstep = mstep.mvnorm)

model <- hsmmfit(matrix(unlist(newDF),ncol=2),
                model, mstep = mstep.mvnorm)
plot(model)
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
