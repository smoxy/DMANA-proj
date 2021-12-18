bootstrap.N.HMM = function(data, nclust, R){
  master = list()
  cores <- ifelse(detectCores()>6,
                  cores <- 3,
                  cores <- 3)
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

HMM <- bootstrap.N.HMM(data = newDF, nclust = 4, R = 500)
save(HMM,file="N-HMM-2.RData",safe=T)
