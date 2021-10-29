if(!"pacman"%in%installed.packages()[,"Package"]) install.packages("pacman")
pacman::p_load('biotools', 'NbClust', "resemble","Rcpp","stargazer","ggplot2","dplyr",
               "factoextra","cluster","dendextend","corrplot","igraph","fpc","clValid")
################################################################################
load("Z:/DesktopC/LUMSA/2/Data Mining/Clustering/ClusterData_L31.RData")

preProcessing <- function(df){
  df[3:10] <- scale(df[3:10])
  df <- df[-c(1,2)]
  return(df)
}

catReg <- function(df){
  regionNames <- c("Sud","Isola","Nord")
  for (i in unique(df$Region)){
    df$Region[df$Region == i] <- regionNames[i]
  }
  return(df)
}

catArea <- function(df){
  areaNames <- c("NApu","Cala","SApu","Sici","ISar","CSar","ELig","WLig","Umbr")
  for (i in unique(df$Area)){
    df$Area[df$Area == i] <- areaNames[i]
  }
  return(df)
}

df <- catReg(olive)
df <- catArea(df)
df <- preProcessing(df)
rm(coffee)
rm(f.voles)
rm(x2)
olive <- catReg(olive)
olive <- catArea(olive)

table(olive$Area)
table(olive$Region)

olive %>%
  count(region = factor(Region), area = factor(Area)) %>%
  mutate(pct = prop.table(n)) %>% 
  ggplot(aes(x = area, y = pct, fill = area, label = scales::percent(pct))) + 
  geom_col(position = 'dodge') + 
  geom_text(position = position_dodge(width = .9),
            vjust = -0.5, 
            size = 3) + 
  scale_y_continuous(name = "Percentage")+
  scale_x_discrete(name = "Area",labels=c("NApu","Cala","SApu","Sici","ISar","CSar","ELig","WLig","Umbr"))+
  scale_fill_discrete(name="Area",labels=c("NApu","Cala","SApu","Sici","ISar","CSar","ELig","WLig","Umbr"))+
  scale_fill_hue(l=40, c=35)+
  theme(legend.position = "none")

olive %>%
  count(region = factor(Region)) %>%
  mutate(pct = prop.table(n)) %>% 
  ggplot(aes(x = region, y = pct, fill = region, label = scales::percent(pct))) + 
  geom_col(position = 'dodge') + 
  geom_text(position = position_dodge(width = .9),
            vjust = -0.5, 
            size = 3) + 
  scale_y_continuous(name = "Percentage")+
  scale_x_discrete(name = "Region",labels=c("Sud","Isola","Nord"))+
  scale_fill_discrete(name="Region",labels=c("Sud","Isola","Nord"))+
  scale_fill_hue(l=40, c=35)+
  theme(legend.position = "none")
  
################################summary#########################################
summary(olive[, -c(1, 2)])
plot(as.data.frame(scale(olive[, -c(1, 2)])))
plot(as.data.frame(scale(olive[, -c(1, 2)])), col = as.factor(olive$Region))

# Between variations
olive %>% group_by(Area) %>%
  select(Palmitic,Palmitoleic,Stearic,Oleic,Linoleic,Linolenic,Arachidic,Eicosenoic) %>% 
  summarize_all(mean) %>% 
  as.data.frame %>% 
  select(-Area) %>%
  stargazer(type = "text")

# Within variations
olive %>% group_by(Area) %>%
  select(Palmitic,Palmitoleic,Stearic,Oleic,Linoleic,Linolenic,Arachidic,Eicosenoic) %>% 
  mutate_all(function(x) {x - mean(x)}) %>%
  as.data.frame %>% 
  select(-Area) %>%
  stargazer(type = "text", omit.summary.stat = "mean")
##############################distance##########################################
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
results <- combDist(c("euclidean", "manhattan", "minkowski","mahalanobis"),
                    c("single", "complete", "average", "ward.D"),df)

par(mfrow=c(4,2))
plot(results[[2]])
plot(results[[4]])
plot(results[[6]])
plot(results[[8]])
plot(results[[10]])
plot(results[[12]])
plot(results[[14]])
plot(results[[16]])

fviz_dend(results[[16]],k=3,cex=.5, k_colors = c("#2E9FDF","#00AFBB","#E7B800"),color_labels_by_k = T, rect=T)
##############################hierarchical######################################

bestCluster <- function(data, cols, nclust){
  clustering <- matrix(NA, nrow = nrow(data), ncol = cols)
  for (j in 1:cols){
    clustering[, j] <- cutree(results[[j]], k = nclust)
  }
  return(clustering)
}

clustering <- bestCluster(df, length(results), 9)
fviz_cluster(list(data=df,cluster=clustering[,4]),
             palette=c(1:9),
             ellipse.type = "convex",
             repel=T,
             show.clust.cent = F,
             ggtheme=theme_minimal())
res.agnes <- agnes(x=olive[,-c(1,2)], stand=T,metric="euclidean",method = "ward")
fviz_dend(res.agnes, ces=.6, k=3)
fviz_dend(res.agnes, ces=.6, k=9)
dendList <- dendlist(as.dendrogram(results[[2]]),as.dendrogram(results[[4]]))
tanglegram(as.dendrogram(results[[2]]),as.dendrogram(results[[4]])) #!!!!!
dend1 <- df %>% D2.dist(data = ., cov = cov(.)) %>% hclust("single") %>% as.dendrogram
dend2 <- df %>% D2.dist(data = ., cov = cov(.)) %>% hclust("complete") %>% as.dendrogram
dend3 <- df %>% D2.dist(data = ., cov = cov(.)) %>% hclust("average") %>% as.dendrogram
dend4 <- df %>% D2.dist(data = ., cov = cov(.)) %>% hclust("ward.D") %>% as.dendrogram
dendList <- dendlist("Single"=dend1,"Complete"=dend2,"Average"=dend3,"Ward"=dend4)
cors <- cor.dendlist(dendList)
corrplot(cors, "pie","lower")
fviz_dend(dend4, k=3,
          k_colors = "jco",
          type="phylogenic",
          repel=T,
          phylo_layout = "layout.gem")



optimalCluster <- function(data, distance, method, indexes, dissMatrix){
  o <- 0
  optimal <- list()
  for (i in 1:length(distance)){
    for (j in 1:length(method)){
      for (k in 1:length(indexes)){
        o <- o + 1
        ifelse(distance[i]=="mahalanobis",
               optimal[[o]] <- NbClust(data,diss = dissMatrix, distance = NULL, method=method[j],index=indexes[k],max.nc=9),
               optimal[[o]] <- NbClust(data,method=method[j],distance=distance[i],index=indexes[k],max.nc=9))
      }
    }
  }
  return(optimal)
}

optimal <- optimalCluster(df, distance=c("euclidean", "manhattan", "minkowski","mahalanobis"),
                              method=c("single", "complete", "average", "ward.D"),
                              indexes = c("alllong"),
                              dissMatrix = as.dist(f_diss(df, diss_method ="mahalanobis",center =T)))

clusterDF <- function(optimal){
  result <- data.frame()
  for (i in 1:length(optimal)){
    result[i,1] <- optimal[[i]]$Best.nc[1]
    result[i,2] <- optimal[[i]]$Best.nc[2]
  }
  rownames(result) <- c("EU-S","EU-C","EU-A","EU-W","MAN-S","MAN-C","MAN-A","MAN-W","MIN-S","MIN-C","MIN-A","MIN-W","MAH-S","MAH-C","MAH-A","MAH-W")
  return(result)
}
clusterStats5 <- clusterDF(optimal)

fviz_nbclust(df, kmeans, nstart=50, method = "wss")+
  geom_vline(xintercept = 5, linetype=2)+
  theme_classic()

fviz_nbclust(df, kmeans, nstart=50, method = "silhouette")+
  theme_classic()

fviz_nbclust(df, kmeans, nstart=50, method = "gap_stat", nboot = 500)+
  theme_classic()

set.seed(123)
km.res <- kmeans(df, 5, nstart=50)
print(km.res)
fviz_cluster(km.res, data=df,
             palette= c(1:5),
             ellipse.type="euclid",
             star.plot=T,
             repel=T,
             ggtheme=theme_minimal())

km.res <- kmeans(df, 7, nstart=50)
print(km.res)
fviz_cluster(km.res, data=df,
             palette= c(1:7),
             ellipse.type="euclid",
             star.plot=T,
             repel=T,
             ggtheme=theme_minimal())

pam(x=D2.dist(data = df, cov = cov(df)), k = 5, stand = T)

fviz_nbclust(df, pam, nstart=50, method = "wss")+
  theme_classic()

fviz_nbclust(df, pam, nstart=50, method = "silhouette")+
  theme_classic()

fviz_nbclust(df, pam, method = "gap_stat", nboot = 500)+
  theme_classic()

pam.res <- pam(x=df,metric ='euclidean', k = 5, stand = T)
fviz_cluster(pam.res,
             palette=c(1:5),
             ellipse.type = "t",
             repel=T,
             ggtheme=theme_classic())

pam.res <- pam(x=df,metric ='manhattan', k = 5, stand = T)
fviz_cluster(pam.res,
             palette=c(1:5),
             ellipse.type = "t",
             repel=T,
             ggtheme=theme_classic())

hc.res <- eclust(df, "hclust", k=5, hc_metric = "euclidean", hc_method = "ward.D", graph = F)
fviz_dend(hc.res, show_labels = F, palette="jco",as.ggplot=T)

km.res <- eclust(df, "kmeans", k=5, nstart=50, graph=F)
fviz_cluster(km.res, geom="point",ellipse.type = "norm", palette="jco",ggtheme=theme_minimal())

pam.res <- eclust(df, "pam", k=5, nstart=50, graph=F)
fviz_cluster(pam.res, geom="point",ellipse.type = "norm", palette="jco",ggtheme=theme_minimal())

fviz_silhouette(hc.res, palette="jco",ggtheme=theme_classic())
fviz_silhouette(km.res, palette="jco",ggtheme=theme_classic())
fviz_silhouette(pam.res, palette="jco",ggtheme=theme_classic())


table(olive$Area, hc.res$cluster)
cluster.stats(d=dist(df),as.numeric(olive$Area), hc.res$cluster)$vi

table(olive$Area, km.res$cluster)
cluster.stats(d=dist(df),as.numeric(olive$Area), km.res$cluster)$vi

table(olive$Area, pam.res$cluster)
cluster.stats(d=dist(df),as.numeric(olive$Area), pam.res$cluster)$vi

optimalScores(clValid(df, 2:9, clMethods = c("hierarchical","kmeans","pam"),
                      validation = "stability", 
                      metric="euclidean",
                      method = "ward"))

optimalScores(clValid(df, 2:9, clMethods = c("hierarchical","kmeans","pam"),
                      validation = "internal", 
                      metric="euclidean",
                      method = "ward"))

optimalScores(clValid(df, 2:9, clMethods = c("hierarchical","kmeans","pam"),
                      validation = "stability", 
                      metric="manhattan",
                      method = "ward"))

optimalScores(clValid(df, 2:9, clMethods = c("hierarchical","kmeans","pam"),
                      validation = "internal", 
                      metric="manhattan",
                      method = "ward"))

optimalScores(clValid(df, 2:9, clMethods = c("hierarchical","kmeans","pam"),
                      validation = "stability", 
                      metric="correlation",
                      method = "ward"))

optimalScores(clValid(df, 2:9, clMethods = c("hierarchical","kmeans","pam"),
                      validation = "internal", 
                      metric="correlation",
                      method = "ward"))

res.hk <- hkmeans(df, hc.metric = "euclidean", hc.method = "ward.D", 5)
res.hk <- hkmeans(df, hc.metric = "manhattan",  hc.method = "ward.D", 5)
fviz_dend(res.hk, cex=.6, palette="jco",rect=T, rect_border="jco",rect_fill=T)
fviz_cluster(res.hk, palette="jco",repel=T,ggtheme=theme_classic())

cluster.stats(d=dist(df),as.numeric(olive$Area), res.hk$cluster)$vi
