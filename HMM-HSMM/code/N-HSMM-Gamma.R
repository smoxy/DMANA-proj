bootstrap.N.HSMM.GAMMA <- function(data, nclust, R) {
  master = list()
  cores <- ifelse(detectCores() > 6,
                  cores <- 6,
                  cores <- 4)
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
                        start.val, mstep = mstep.mvnorm)
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

GAMMA.HSMM <-
  bootstrap.N.HSMM.GAMMA(
    data = newDF,
    nclust = 4,
    R = 500)

save(GAMMA.HSMM, file = "N-HSMM-2-4-Gamma.RData", safe = T)
