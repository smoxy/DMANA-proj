if(!"pacman"%in%installed.packages()[,"Package"]) install.packages("pacman")
pacman::p_load('biotools', 'NbClust', "resemble","Rcpp")
#############
preProcessing <- function(df) {
  df[3:10] <- scale(df[3:10])
  df <- df[-c(1, 2)]
  return(df)
}

catReg <- function(df) {
  regionNames <- c("Sud", "Isola", "Nord")
  for (i in unique(df$Region)) {
    df$Region[df$Region == i] <- regionNames[i]
  }
  return(df)
}

catArea <- function(df) {
  areaNames <-
    c("NApu",
      "Cala",
      "SApu",
      "Sici",
      "ISar",
      "CSar",
      "ELig",
      "WLig",
      "Umbr")
  for (i in unique(df$Area)) {
    df$Area[df$Area == i] <- areaNames[i]
  }
  return(df)
}
load("C:/Users/media/OneDrive/Desktop/University/Data Mining&Analytics/3Presentazione/ClusterData_L31 (1).RData")
df <- catReg(olive)
df <- catArea(df)
df <- preProcessing(df)
rm(coffee)
rm(f.voles)
rm(x2)
################################summary######################################
summary(olive[, -c(1, 2)])
plot(as.data.frame(scale(olive[, -c(1, 2)])))
plot(as.data.frame(scale(olive[, -c(1, 2)])), col = as.factor(olive$Region))
##############################distance########################
combDist <- function(distance, methods, data) {
  k <- 0
  results <- list()
  for (i in 1:length(distance)){
    ifelse(distance[i] == "minkowski",
           dist <- dist(df, method = distance[i], p = 3),
           ifelse(distance[i] == "mahalanobis",
                  dist <- D2.dist(data = df, cov = cov(df)),
                  dist <- dist(df, method = distance[i])))
  for (j in 1:length(methods)){
        k <- k + 1
        results[[k]] <- hclust(dist, method = methods[j])
  }
  }
  return(results)
}
results <- combDist(c("euclidean", "manhattan", "minkowski","mahalanobis"),c("single", "complete", "average", "ward"),df)
##############################NbClust################################


bestCluster <- function(data, cols, nclust){
  clustering <- matrix(NA, nrow = nrow(data), ncol = cols)
  for (j in 1:cols){
    clustering[, j] <- cutree(results[[j]], k = nclust)
  }
  return(clustering)
}

optimalCluster <- function(data, distance, method, indexes){
  o <- 0
  optimal <- list()
  for (i in 1:length(distance)){
    for (j in 1:length(method)){
      for (k in 1:length(indexes)){
        o <- o + 1
        optimal[[o]]<- NbClust::NbClust(data,method=method[j],distance=distance[i],index=indexes[k])
      }
    }
  }
  return(optimal)
}
mahala <- f_diss(df, diss_method = "mahalanobis")
NbClust::NbClust(df, diss=as.dist(f_diss(df, diss_method = "mahalanobis", 
                                         center = TRUE)), distance = NULL, method="ward.D")
