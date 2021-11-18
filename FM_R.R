library(mclust)
data("faithful")
plot(faithful)
x <- faithful$waiting
hist(x,breaks=20,freq=FALSE)
clustering <- Mclust(x,G=2,modelNames = "V")
seq.values <- seq(min(x),max(x),length=2000)
den1 <- dnorm(seq.values,mean=clustering$parameters$mean[1],
              sd=sqrt(clustering$parameters$variance$sigmasq[1]))
den2 <- dnorm(seq.values,mean=clustering$parameters$mean[2],
              sd=sqrt(clustering$parameters$variance$sigmasq[2]))
pi1 <- clustering$parameters$pro[1]
pi2 <- clustering$parameters$pro[2]
den <- pi1*den1 + pi2*den2
lines(seq.values,den)
lines(seq.values,pi1*den1,col=2)
lines(seq.values,pi2*den2,col=3)
#
# Generate data from a finite mixture of Gaussians
#
gen.gauss.mixt <- function(N=1000,prior,mu,sigma2,start.seed=1234)
{
  set.seed(start.seed)
  K <- length(prior) # number of clusters
  cluster <- sample(1:K,prob = prior, size=N, replace=TRUE) # assignment variable, z
  x <- matrix(NA,nrow=N,ncol=1)
  for(i in 1:N)
  {
    x[i] <- rnorm(1,mean=mu[cluster[i]],sd=sqrt(sigma2[cluster[i]]))
  }
  return(list(x = x, cluster = cluster))
}
# Example
mu <- c(2,-2)
sigma2 <- c(1,1)
prior <- c(.6,.4)
data <- gen.gauss.mixt(N=2000,prior=prior,mu=mu,sigma2 = sigma2)
par(mfrow=c(1,2))
hist(data$x,breaks=20)
plot(data$x,col=data$cluster)
library(mclust)
mixt.mclust <- Mclust(data$x,G=2,model="V")
mixt.mclust$parameters$pro
mixt.mclust$parameters$mean
mixt.mclust$parameters$variance$sigmasq
#
# EM for univariate mixture of Gaussian
em.gauss.mixt <- function(data,prior,mu,sigma2,eps=10^-4)
{
  K <- length(prior)
  N <- length(data)
  dens <- matrix(NA,nrow=N,ncol=K)
  for(k in 1:K)
  {
    dens[,k] <- dnorm(data,mean = mu[k],sd=sqrt(sigma2[k]))
  }
  dif <- Inf
  w <- matrix(NA,nrow=N,ncol=K)
  log.lik.old <- -Inf
  iter <- 0
  while(dif>eps)
  {
    iter <- iter+1
    print(iter)
      # E-step
  for(k in 1:K)
  {
    w[,k] <- prior[k]*dens[,k]
  }
    w <- w/apply(w,1,sum) # Posterior weights
   #
   # M-step
   prior <- apply(w,2,sum)/N # prior = apply(w,2,mean)
   #
   for(k in 1:K)
   {
     mu[k] <- sum(w[,k]*data)/sum(w[,k])
     sigma2[k] <- sum(w[,k]*(data-mu[k])^2)/sum(w[,k])
     dens[,k] <- dnorm(data,mean = mu[k], sd = sqrt(sigma2[k]))
   }
   log.lik <- sum(log(dens%*%prior))
   dif <- log.lik-log.lik.old
   log.lik.old <- log.lik
   print(dif)
   print(log.lik)
  }
  npar <- (K-1) + K + K# prior probs + means + sigma2
  AIC <- -2*log.lik+2*npar
  BIC <- -2*log.lik+log(N)*npar
  return(list(posterior=w,mu=mu,sigma2=sigma2,prior=prior,llk = log.lik, AIC = AIC, BIC = BIC))
}
#
mixt2 <- em.gauss.mixt(data=data$x,prior=c(.5,.5),mu=c(7,-4),sigma2=c(1,1))
clustering2 <- apply(mixt2$posterior,1,which.max)
table(data$cluster,clustering2)
table(data$cluster,mixt.mclust$classification)
#
ini.par <- kmeans(data$x,centers = 3)
mixt3 <- em.gauss.mixt(data=data$x,prior=table(ini.par$cluster)/2000,mu=ini.par$centers,sigma2=c(1,1,1))
#
# Comparing models using the diabetes data
library(mclust)
data("diabetes")
x.dist=dist(diabetes[,-1],method="euclidean")
hc.single=hclust(x.dist,method="single",members=NULL)
hc.complete=hclust(x.dist,method="complete",members=NULL)
hc.average=hclust(x.dist,method="average",members=NULL)
hc.ward=hclust(x.dist,method="ward.D",members=NULL)
k.means <- kmeans(diabetes[,-1],centers=3,nstart = 20)
mixt.diab <- Mclust(diabetes[,-1],G=1:4)
plot(hc.single)
plot(hc.complete)
plot(hc.average)
plot(hc.ward)
hicluste1=cutree(hc.single,k=4)
hicluste2=cutree(hc.complete,k=4)
hicluste3=cutree(hc.average,k=3)
hicluste4=cutree(hc.ward,k=3)
table(diabetes[,1],hicluste1)
table(diabetes[,1],hicluste2)
table(diabetes[,1],hicluste3)
table(diabetes[,1],hicluste4)
adjustedRandIndex(diabetes[,1],hicluste1)
adjustedRandIndex(diabetes[,1],hicluste2)
adjustedRandIndex(diabetes[,1],hicluste3)
adjustedRandIndex(diabetes[,1],hicluste4)
table(diabetes[,1],k.means$cluster)
adjustedRandIndex(diabetes[,1],k.means$cluster)
table(diabetes[,1],mixt.diab$classification)
adjustedRandIndex(diabetes[,1],mixt.diab$classification)
#
# Model selection: the faithful data
mixt2.faith <- Mclust(faithful,G=2,model="VVV")
mixt3.faith <- Mclust(faithful,G=3,model="VVV")
mixt4.faith <- Mclust(faithful,G=4,model="VVV")
mixt.faith <- Mclust(faithful,G=1:4)

mixt2.cn.faith <- CNmixt(faithful, G = 1:3, parallel = TRUE, seed =1234, model = "VVV")

# EM algorithm for multivariate Gaussian mixtures
em.gauss.mixt.mult <- function(data,prior,mu,Sigma,eps=10^-4)
{
  # mu must be a matrix K \times P
  # Sigma must be an array P \times P \times K
  # with P = number of variables
  # and K = number of clusters
  library(mvtnorm)
  K <- length(prior)
  N <- dim(data)[1]
  P <- dim(data)[2]
  dens <- matrix(NA,nrow=N,ncol=K)
  for(k in 1:K)
  {
    dens[,k] <- dmvnorm(data,mean = mu[k,],sigma=Sigma[,,k])
  }
  dif <- Inf
  w <- matrix(NA,nrow=N,ncol=K)
  log.lik.old <- -Inf
  iter <- 0
  while(dif>eps)
  {
    iter <- iter+1
    print(iter)
    # E-step
    for(k in 1:K)
    {
      w[,k] <- prior[k]*dens[,k]
    }
    w <- w/apply(w,1,sum) # Posterior weights
    #
    # M-step
    prior <- apply(w,2,sum)/N # prior = apply(w,2,mean)
    #
    for(k in 1:K)
    {
      mu[k,] <- colSums(w[,k]/sum(w[,k])*data)
      Sigma[,,k]    <- crossprod(sqrt(w[,k]/sum(w[,k]))*(data-matrix(rep(mu[k,],each=N),ncol=P)))
      dens[,k] <- dmvnorm(data,mean = mu[k,],sigma=Sigma[,,k])
    }
    log.lik <- sum(log(dens%*%prior))
    dif <- log.lik-log.lik.old
    log.lik.old <- log.lik
    print(dif)
    print(log.lik)
  }
  npar <- (K-1) + K*P + K*(P*(P+1)/2)# prior probs + means + sigma2
  AIC <- -2*log.lik+2*npar
  BIC <- -2*log.lik+log(N)*npar
  return(list(posterior=w,mu=mu,Sigma=Sigma,prior=prior,llk = log.lik, AIC = AIC, BIC = BIC))
}
#
k.means.faith <- kmeans(faithful,centers=2,nstart = 20)
mu <- matrix(k.means.faith$centers,2,2)
prior <- table(k.means.faith$cluster)/272
Sigma <- array(NA,c(2,2,2))
Sigma[,,1] <- var(faithful[k.means.faith$cluster==1,])
Sigma[,,2] <- var(faithful[k.means.faith$cluster==2,])
# The faithful data are store as list(), check typeof(faithful)
# To run the "home-made" code, data must be a matrix object
faithful.matrix <- matrix(unlist(faithful),nrow=272,ncol=2)
model <- em.gauss.mixt.mult(faithful.matrix,prior=prior,mu=mu,Sigma=Sigma)
# Atypical data analysis
library("ContaminatedMixt")
library("mnormt")
p <- 2
set.seed(16)
n1 <- n2 <- 200
X1 <- rmnorm(n = n1, mean = rep(2, p), varcov = diag(c(5, 0.5)))
X2 <- rmnorm(n = n2, mean = rep(-2, p), varcov = diag(c(5, 0.5)))
bad <- matrix(runif(n = 20, min = -20, max = 20), nrow = 10, ncol = 2)
X <- rbind(X1, X2, bad)
res1 <- CNmixt(X, G = 1:3, parallel = TRUE, seed = 2)
library(teigen)
res2 <- teigen(X, G=1:3, models = "all", verbose = FALSE)
res2$allbic
#
#
# Variable selection
library(clustvarsel) # mclust
library(sparcl) # sparse k-means
library(vscc) # teigen and mclust
library(SelvarMix) # Lasso as in Zhou et al. (2009)
require(Rmixmod)
require(glasso)
#
#x.data <- banknote[,-1]
x.data <- wines[,-1]
clust.vscc.teigen <- vscc(x.data, G=1:5, automate = "teigen", initial = NULL, train = NULL, forcereduction = FALSE)
head(clust.vscc.mclust$topselected) #Show preview of selected variables
table(wines[,1], clust.vscc.mclust$initialrun$classification) #Clustering results on full data set
table(wines[,1], clust.vscc.mclust$bestmodel$classification) #Clustering results on reduced data set
#
clust <- clustvarsel(x.data,G=1:5)
clust$subset
table(clust$model$classification,wine[,1])
#
clust.b <- clustvarsel(x.data,G=1:5,direction = "backward")
clust.b$subset
table(clust.b$model$classification,wine[,1])
#
clust.mod <- Mclust(x.data[,clust$subset],G=1:5)
table(clust.mod$classification,banknote[,1])
# choose tuning parameter
km.perm <- KMeansSparseCluster.permute(x.data,K=2,wbounds=seq(3,7,len=15),nperms=5)
print(km.perm)
plot(km.perm)
# run sparse k-means
km.out <- KMeansSparseCluster(x.data,K=2,wbounds=km.perm$bestw)
print(km.out)
plot(km.out)
table(km.out[[1]]$Cs,banknote[,1])
#
lambda <- seq(0.1, 100, length = 25)
rho <- seq(1, 2, length=2)
lasso <- SelvarClustLasso(x=x.data, nbcluster=1:3,criterio="ICL",lambda=lambda,rho=rho)
summary(lasso)
lasso$S
lasso$R
lasso$U
lasso$nbcluster
table(lasso$partition,banknote[,1])
#
y.data <- wines[,-1]
wine.mclust <- Mclust(y.data,G=1:5)
summary(wine.mclust)
table(wine[,1],wine.mclust$classification)
adjustedRandIndex(wine[,1],wine.mclust$classification)
wine.contaminated <- CNmixt(y.data, G = 1:5, parallel = TRUE, seed =1234)
wine.teigen <- teigen(y.data, G=1:5, models = "all", verbose = TRUE)
summary(wine.teigen)
table(wine[,1],wine.teigen$iclresults$classification)
adjustedRandIndex(wine[,1],wine.teigen$iclresults$classification)
wine.kmeans <- kmeans(y.data,centers = 3,nstart = 20)
table(wine[,1],wine.kmeans$cluster)
adjustedRandIndex(wine[,1],wine.kmeans$cluster)
miss.clus <- wine[,1]-wine.kmeans$cluster
miss.clus.color <- rep(1,length(miss.clus))
miss.clus.color[miss.clus!=0] <- 2
wine.pgmm <- pgmmEM(y.data,rG=3,rq=2,icl=TRUE)
table(wines[,1],wine.pgmm$map)
adjustedRandIndex(wines[,1],wine.pgmm$map)
clust.vscc <- vscc(y.data, G=1:5, automate = "mclust", initial = NULL, train = NULL, forcereduction = FALSE)
table(wine[,1], clust.vscc$initialrun$classification) #Clustering results on full data set
table(wine[,1], clust.vscc$bestmodel$classification) #Clustering results on reduced data set
wine.clustvarsel <- clustvarsel(y.data,G=1:5)
wine.clustvarsel$subset
table(wine[,1],wine.clustvarsel$model$classification)
adjustedRandIndex(wine[,1],wine.clustvarsel$model$classification)
km.perm <- KMeansSparseCluster.permute(y.data,K=3,wbounds=seq(3,7,len=15),nperms=5)
km.out <- KMeansSparseCluster(y.data,K=3,wbounds=km.perm$bestw)
table(wines[,1],km.out[[1]]$Cs)
adjustedRandIndex(wines[,1],km.out[[1]]$Cs)
#
#
n <- 200      # sample size
pro <- 0.5    # mixing proportion
mu1 <- c(0,0) # mean vector for the first cluster
mu2 <- c(3,3) # mean vector for the second cluster
sigma1 <- matrix(c(1,0.5,0.5,1),2,2)       # covar matrix for the first cluster
sigma2 <- matrix(c(1.5,-0.7,-0.7,1.5),2,2) # covar matrix for the second cluster
X <- matrix(0, n, 5, dimnames = list(NULL, paste0("X", 1:5)))
set.seed(1234) # for replication
u <- runif(n)
Class <- ifelse(u < pro, 1, 2)
X[u < pro, 1:2]  <- MASS::mvrnorm(sum(u < pro), mu = mu1, Sigma = sigma1)
X[u >= pro, 1:2] <- MASS::mvrnorm(sum(u >= pro), mu = mu2, Sigma = sigma2)
X[, 3] <- X[, 1] + rnorm(n)
X[, 4] <- rnorm(n, mean = 1.5, sd = 2)
X[, 5] <- rnorm(n, mean = 2, sd = 1)
