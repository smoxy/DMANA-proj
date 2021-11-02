if(!"pacman"%in%installed.packages()[,"Package"]) install.packages("pacman")
pacman::p_load('biotools', 'NbClust', "resemble","Rcpp","stargazer","ggplot2","dplyr",
               "factoextra","cluster","dendextend","corrplot","igraph","fpc","clValid","HDMD")
################################################################################
load("Z:/DesktopC/LUMSA/2/Data Mining/Clustering/ClusterData_L31.RData")

preProcessing <- function(df){
  df[3:10] <- scale(df[3:10])
  df <- df[-c(1,2)]
  return(df)
}

catReg <- function(df){
  regionNames <- c("Sud","Sardegna","Centro-Nord")
  for (i in unique(df$Region)){
    df$Region[df$Region == i] <- regionNames[i]
  }
  return(df)
}

catArea <- function(df){
  areaNames <- c("NPuglia","Calabria","SPuglia","Sicilia","InSardegna","CoSardegna","EastLiguria","WestLiguria","Umbria")
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
dissMatrix <- as.dist(pairwise.mahalanobis(olive[,-c(1,2)], grouping = c(1:nrow(olive)), cov = cov(olive[,-c(1,2)]))$distance)



combDist <- function(distance, methods, data) {
  k <- 0
  results <- list()
  for (i in 1:length(distance)){
    ifelse(distance[i] == "minkowski",
           dist <- dist(df, method = distance[i], p = 3),
           ifelse(distance[i] == "mahalanobis",
                  dist <- dissMatrix,
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

par(mfrow=c(2,2))
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
dend1 <- df %>% dist(method = "euclidean") %>% hclust("complete") %>% as.dendrogram
dend2 <- df %>% dist(method = "euclidean") %>% hclust("ward.D") %>% as.dendrogram
dend3 <- df %>% dist(method = "manhattan") %>% hclust("complete") %>% as.dendrogram
dend4 <- df %>% dist(method = "manhattan") %>% hclust("ward.D") %>% as.dendrogram
dend5 <- df %>% dist(method = "minkowski",p = 3) %>% hclust("complete") %>% as.dendrogram
dend6 <- df %>% dist(method = "minkowski",p = 3) %>% hclust("ward.D") %>% as.dendrogram
dend7 <- olive[,-c(1,2)] %>% D2.dist(data = ., cov = cov(.)) %>% hclust("complete") %>% as.dendrogram
dend8 <- olive[,-c(1,2)] %>% D2.dist(data = ., cov = cov(.)) %>% hclust("ward.D") %>% as.dendrogram

#Complete
dendList1 <- dendlist("EU-C"=dend1,
                     "MA-C"=dend3,
                     "MI-C"=dend5,
                     "MH-C"=dend7)
cors1 <- cor.dendlist(dendList1)
corrplot(cors1, "pie","lower")

#Ward
dendList2 <- dendlist("EU-W"=dend2,
                     "MA-W"=dend4,
                     "MI-W"=dend6,
                     "MH-W"=dend8)
cors2 <- cor.dendlist(dendList2)
corrplot(cors2, "pie","lower")

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
               optimal[[o]] <- NbClust(data,diss = dissMatrix, distance = NULL, method=method[j],index=indexes[k],min.nc = 3,max.nc=9),
               optimal[[o]] <- NbClust(data,method=method[j],distance=distance[i],index=indexes[k],min.nc = 3,max.nc=9))
      }
    }
  }
  return(optimal)
}

optimal <- optimalCluster(df, distance=c("euclidean", "manhattan", "minkowski","mahalanobis"),
                              method=c("single", "complete", "average", "ward.D"),
                              indexes = c("gap"),
                              dissMatrix = dissMatrix)

clusterDF <- function(optimal){
  result <- data.frame()
  for (i in 1:length(optimal)){
    result[i,1] <- optimal[[i]]$Best.nc[1]
    result[i,2] <- optimal[[i]]$Best.nc[2]
  }
  rownames(result) <- c("EU-S","EU-C","EU-A","EU-W","MAN-S","MAN-C","MAN-A","MAN-W","MIN-S","MIN-C","MIN-A","MIN-W","MAH-S","MAH-C","MAH-A","MAH-W")
  return(result)
}
clusterStats4 <- clusterDF(optimal)
################################ OPTIMAL CLUST N ###############################
5
##############################Partitioning######################################
fviz_nbclust(df, kmeans, nstart=50, method = "wss", k.max=9)+
  geom_vline(xintercept = 5, linetype=2)+
  theme_classic() # 5

fviz_nbclust(df, kmeans, nstart=50, method = "silhouette", k.max=9)+
  theme_classic() # 5

fviz_nbclust(df, kmeans, nstart=50, method = "gap_stat", nboot = 500, k.max=9)+
  theme_classic() # 5

set.seed(123)
km.res <- kmeans(df, 5, nstart=50)
print(km.res)
fviz_cluster(km.res, data=df,
             palette="jco",
             ellipse.type="euclid",
             star.plot=T,
             repel=T,
             ggtheme=theme_minimal())

set.seed(123)
km.res <- kmeans(df, 6, nstart=50)
print(km.res)
fviz_cluster(km.res, data=df,
             palette="jco",
             ellipse.type="euclid",
             star.plot=T,
             repel=T,
             ggtheme=theme_minimal())

set.seed(123)
km.res <- kmeans(df, 7, nstart=50)
print(km.res)
fviz_cluster(km.res, data=df,
             palette="jco",
             ellipse.type="euclid",
             star.plot=T,
             repel=T,
             ggtheme=theme_minimal())


################################## PAM #########################################
#Using Mahalanobis
fviz_nbclust(df, pam, diss = dissMatrix, nstart=50, method = "wss", k.max=9)+
  theme_classic() # 7

fviz_nbclust(df, pam, diss = dissMatrix, nstart=50, method = "silhouette", k.max=9)+
  theme_classic() # 7

fviz_nbclust(df, pam, diss = dissMatrix, method = "gap_stat", nboot = 500, k.max=9)+
  theme_classic() # 5 o 7

pam.res.5 <- pam(x=dissMatrix,diss=T,nstart=50, k = 5)
pam.res.5$data <- olive[,-c(1,2)]
fviz_cluster(pam.res.5,
             palette="jco",
             ellipse.type = "t",
             repel=T,
             ggtheme=theme_classic())

pam.res.7 <- pam(x=dissMatrix,diss=T,nstart=50, k = 7)
pam.res.7$data <- olive[,-c(1,2)]
fviz_cluster(pam.res.7,
             palette="jco",
             ellipse.type = "t",
             repel=T,
             ggtheme=theme_classic())

hc.res <- eclust(df, "hclust", k=5, hc_metric = "euclidean", hc_method = "ward.D", graph = T, nboot=500)
fviz_dend(hc.res, show_labels = F, palette="jco",as.ggplot=T)

km.res <- eclust(df, "kmeans", k=5, nstart=50, graph=F)
fviz_cluster(km.res, geom="point",ellipse.type = "norm", palette="jco",ggtheme=theme_minimal(), graph = T, nboot=500)

pam.res <- eclust(df, "pam", k=5, nstart=50, graph=F)
fviz_cluster(pam.res, geom="point",ellipse.type = "norm", palette="jco",ggtheme=theme_minimal(), graph = T, nboot=500)

fviz_silhouette(hc.res, palette="jco",ggtheme=theme_classic())+labs(caption="hclust")
fviz_silhouette(km.res, palette="jco",ggtheme=theme_classic())+labs(caption="kmeans")
fviz_silhouette(pam.res, palette="jco",ggtheme=theme_classic())+labs(caption="PAM")

################################################################################
load("Z:/DesktopC/LUMSA/2/Data Mining/Clustering/ClusterData_L31.RData")

table(olive$Area, hc.res$cluster)
cluster.stats(d=dissMatrix,as.numeric(olive$Area), hc.res$cluster)$vi

table(olive$Area, km.res$cluster)
cluster.stats(d=dissMatrix,as.numeric(olive$Area), km.res$cluster)$vi

table(olive$Area, pam.res$cluster)
cluster.stats(d=dissMatrix,as.numeric(olive$Area), pam.res$cluster)$vi

optimalScores(clValid(df, 3:9, clMethods = c("hierarchical","kmeans","pam"),
                      validation = "stability", 
                      metric="euclidean",
                      method = "ward"))

optimalScores(clValid(df, 3:9, clMethods = c("hierarchical","kmeans","pam"),
                      validation = "internal", 
                      metric="euclidean",
                      method = "ward"))

optimalScores(clValid(df, 3:9, clMethods = c("hierarchical","kmeans","pam"),
                      validation = "stability", 
                      metric="correlation",
                      method = "ward"))

optimalScores(clValid(df, 3:9, clMethods = c("hierarchical","kmeans","pam"),
                      validation = "internal", 
                      metric="correlation",
                      method = "ward"))

res.hk <- hkmeans(df, hc.metric = "euclidean", hc.method = "ward.D", 5)
fviz_dend(res.hk, cex=.6, palette="jco",rect=T, rect_border="jco",rect_fill=T)
fviz_cluster(res.hk, palette="jco",repel=T,ggtheme=theme_classic())

table(olive$Area, res.hk$cluster)
cluster.stats(d=dissMatrix,as.numeric(olive$Area), res.hk$cluster)$vi
