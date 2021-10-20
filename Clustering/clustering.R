stdVar <- function(vars, df){
  df <- na.omit(df)
  attach(df)
  for (i in vars){
    df[i] <- scale(df[i])
  }
  return(df)
}
stdOlive <- stdVar(c("Palmitic","Palmitoleic","Stearic","Oleic","Linoleic","Linolenic","Arachidic","Eicosenoic"),olive)
attach(stdOlive)

# Determine number of clusters
wss <- (nrow(stdOlive)-1)*sum(apply(stdOlive,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(stdOlive,
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

# K-Means Cluster Analysis
fit <- kmeans(stdOlive, 5) # 5 cluster solution
# get cluster means
aggregate(stdOlive,by=list(fit$cluster),FUN=mean)
# append cluster assignment
mydata <- data.frame(stdOlive, fit$cluster)

# Ward Hierarchical Clustering
d <- dist(stdOlive, method = "euclidean") # distance matrix
fit <- hclust(d, method="ward.D")
plot(fit) # display dendogram
groups <- cutree(fit, k=5) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters
rect.hclust(fit, k=5, border="red")

# Ward Hierarchical Clustering with Bootstrapped p values
library(pvclust)
fit <- pvclust(stdOlive, method.hclust="ward.D",
               method.dist="euclidean")
plot(fit) # dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit, alpha=.95)

# Model Based Clustering
library(mclust)
fit <- Mclust(stdOlive)
plot(fit) # plot results
summary(fit) # display the best model
