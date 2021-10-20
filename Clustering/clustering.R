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

x.dist1=dist(stdOlive,method="euclidean")
x.dist2=dist(stdOlive,method="manhattan")
x.dist3=dist(stdOlive,method="minkowski")
hc.single1=hclust(x.dist1,method="single",members=NULL)
hc.single2=hclust(x.dist2,method="single",members=NULL)
hc.single3=hclust(x.dist3,method="single",members=NULL)
hc.ward=hclust(x.dist1,method="ward.D",members=NULL)
plot(hc.single1)
plot(hc.single2)
plot(hc.single3)
plot(hc.ward)
hicluste1=cutree(hc.single1,k=3)
hicluste2=cutree(hc.single2,k=3)
hicluste3=cutree(hc.single3,k=3)
hiward=cutree(hc.ward,k=3)
plot(stdOlive, col=hicluste1,main="Single Linkage - Euclidean distance")
plot(stdOlive, col=hicluste2,main="Single Linkage - Manhattan distance")
plot(stdOlive, col=hicluste3,main="Single Linkage - Minkowski distance")
plot(stdOlive, col=hiward,main="Ward's minimum variance - Euclidean distance")
#
table(hicluste1,olive[,1])
table(hicluste2,olive[,1])
table(hicluste3,olive[,1])
table(hiward,olive[,1])
#
k.olive <- kmeans(stdOlive,centers=3,nstart = 20)
plot(stdOlive,col=k.olive$cluster)
t(table(olive[,1],k.olive$cluster))

table(olive[,1]) #true

###############################INTERNET#########################################
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization

distance <- get_dist(stdOlive)
fviz_dist(distance, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))
kOlive <- kmeans(stdOlive, centers = 3, nstart = 25)
kOlive
fviz_cluster(kOlive, data = stdOlive)
stdOlive %>%
  as_tibble() %>%
  mutate(cluster = kOlive$cluster,
         state = row.names(stdOlive)) %>%
  ggplot(aes(Palmitic, Palmitoleic, color = factor(cluster), label = olive$Region)) +
  geom_text()
kOlive1 <- kmeans(stdOlive, centers = 1, nstart = 25)
kOlive2 <- kmeans(stdOlive, centers = 2, nstart = 25)
kOlive3 <- kmeans(stdOlive, centers = 3, nstart = 25)
# plots to compare
p1 <- fviz_cluster(kOlive1, geom = "point", data = stdOlive) + ggtitle("k = 1")
p2 <- fviz_cluster(kOlive2, geom = "point",  data = stdOlive) + ggtitle("k = 2")
p3 <- fviz_cluster(kOlive3, geom = "point",  data = stdOlive) + ggtitle("k = 3")
library(gridExtra)
grid.arrange(p1, p2, p3, nrow = 2)

set.seed(123)
fviz_nbclust(stdOlive, kmeans, method = "wss")
fviz_nbclust(stdOlive, kmeans, method = "silhouette")
set.seed(123)
gap_stat <- clusGap(stdOlive, FUN = kmeans, nstart = 25,
                    K.max = 10, B = 50)
# Print the result
print(gap_stat, method = "firstmax")
fviz_gap_stat(gap_stat)
set.seed(123)
wss <- kmeans(stdOlive, 10, nstart = 25)
sil <- kmeans(stdOlive, 5, nstart = 25)
final <- kmeans(stdOlive, 9, nstart = 25)
print(final)
fviz_cluster(final, data = stdOlive)
stdOlive %>%
  mutate(Cluster = final$cluster) %>%
  group_by(Cluster) %>%
  summarise_all("mean")

p1 <- fviz_cluster(wss, geom = "point", data = stdOlive) + ggtitle("k = 5")
p2 <- fviz_cluster(sil, geom = "point",  data = stdOlive) + ggtitle("k = 9")
p3 <- fviz_cluster(final, geom = "point",  data = stdOlive) + ggtitle("k = 10")
library(gridExtra)
grid.arrange(p1, p2, p3, nrow = 2)
wss
sil
final
####################################BOOK########################################
install.packages(c("cluster","factorextra"))
dist.eucl <- dist(df, method = "euclidean")
round(as.matrix(dist.eucl)[1:9, 1:9], 1)
dist.cor <- get_dist(df, method = "pearson")
round(as.matrix(dist.cor)[1:9, 1:9], 1)
fviz_dist(dist.eucl)
fviz_nbclust(df, kmeans, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2)
set.seed(123)
km.res <- kmeans(df, 4, nstart = 25)
print(km.res)
aggregate(df, by=list(cluster=km.res$cluster),mean)
dd <- cbind(olive, cluster=km.res$cluster)
head(dd)
km.res$size
km.res$centers
fviz_cluster(km.res, data=df,
             palette=c("#2E9FDF","#00AFBB","#E7B800","#FC4E07"),
             ellipse.type = "euclid",
             star.plot=TRUE,
             repel=TRUE,
             ggtheme = theme_minimal())

fviz_nbclust(df, pam, method = "silhouette")+
  theme_classic()
pam.res <- pam(df, 5)
print(pam.res)
dd <- cbind(df, cluster=pam.res$cluster)
pam.res$medoids
head(pam.res$cluster)
fviz_cluster(pam.res,
             palette=c("#2E9FDF","#00AFBB","#E7B800","#FC4E07","purple"),
             ellipse.type = "t",
             repel=T,
             ggtheme = theme_classic())

fviz_nbclust(df, clara, method = "silhouette")+
  theme_classic()
clara.res <- clara(df, 5, samples=50, pamLike=TRUE)
print(clara.res)
dd <- cbind(df, cluster=clara.res$cluster)
head(clara.res$clustering,10)
clara.res$medoids
fviz_cluster(clara.res,
             palette=c("#2E9FDF","#00AFBB","#E7B800","#FC4E07","purple"),
             ellipse.type = "t",
             geom="point", pointsize = 1,
             ggtheme = theme_classic())
res.hc <- hclust(d = dist.eucl, method = "ward.D2")
fviz_dend(res.hc, cex=.5)
res.hc2 <- hclust(dist.eucl, method = "average")
grp <- cutree(res.hc, k=5)
head(grp, n=5)
table(grp)
fviz_dend(res.hc, k=5,
          cex = 0.5,
          k_colors = c("#2E9FDF","#00AFBB","#E7B800","#FC4E07","purple"),
          color_labels_by_k = TRUE,
          rect = TRUE)
fviz_cluster(list(data=df, cluster=grp),
             palette=c("#2E9FDF","#00AFBB","#E7B800","#FC4E07","purple"),
             ellipse.type = "convex",
             repel=T,
             show.clust.cent = F,
             ggtheme = theme_minimal())
res.agnes <- agnes(x=df,
                   stand = T,
                   metric = "euclidean",
                   method = "ward")
res.diana <- diana(x=df, 
                   stand = T,
                   metric = "euclidean")
fviz_dend(res.agnes, cex=0.6, k=5)
install.packages("dendextend")
library(dendextend)
res.dist <- dist(df, method = "euclidean")
hc1 <- hclust(res.dist,method = "average")
hc2 <- hclust(res.dist,method = "ward.D2")
dend1 <- as.dendrogram(hc1)
dend2 <- as.dendrogram(hc2)
dend_list <- dendlist(dend1, dend2)
tanglegram(dend1, dend2)
cor.dendlist(dend_list, method = "cophenetic")
cor.dendlist(dend_list, method = "baker")
cor_cophenetic(dend1, dend2)
cor_bakers_gamma(dend1,dend2)
dend1 <- df %>% dist %>% hclust("complete") %>% as.dendrogram()
dend2 <- df %>% dist %>% hclust("single") %>% as.dendrogram()
dend3 <- df %>% dist %>% hclust("average") %>% as.dendrogram()
dend4 <- df %>% dist %>% hclust("centroid") %>% as.dendrogram()
dend_list <- dendlist("Complete"=dend1,"Single"=dend2,"Average"=dend3,"Centroid"=dend4)
cors <- cor.dendlist(dend_list)
round(cors, 2)
corrplot::corrplot(cors, "pie", "lower")
dd <- dist(df, method = "euclidean")
hc <- hclust(dd, method = "ward.D2")
require("igraph")
fviz_dend(hc, k=4, k_colors = "jco",
          type="phylogenic",repel = T)
fviz_dend(hc, k=4, k_colors = "jco",
          type="phylogenic",repel = T,
          phylo_layout = "layout.gem")
install.packages("clustertend")
fviz_pca_ind(prcomp(df),habillage = olive$Area, palette = "jco",
             geom = "point",ggtheme=theme_classic(),legend="bottom")
fviz_pca_ind(prcomp(df),geom = "point",ggtheme=theme_classic())
km.res1 <- kmeans(df, 5)
fviz_cluster(list(data=df, cluster=km.res1$cluster),
             ellipse.type = "norm",
             geom = "point",
             stand=F,
             palette="jco",
             ggtheme = theme_minimal())
fviz_dend(hclust(dist(df)),k=5, k_colors = "jco",
          as.ggplot=T, show_labels = F)
set.seed(123)
hopkins(df, n=nrow(df)-1)
fviz_dist(dist(df),show_labels = F)
#install.packages("NbClust")
library(NbClust)
fviz_nbclust(df, kmeans, method = "wss")+
  geom_vline(xintercept = 4, linetype=2)
fviz_nbclust(df, kmeans, method = "silhouette")
set.seed(123)
fviz_nbclust(df, kmeans, nstart=25, method = "gap_stat", nboot = 500)
nb <- NbClust(df, distance = "euclidean", min.nc = 2, max.nc = 10, method = "kmeans")
fviz_nbclust(nb)
km.res <- eclust(df, "kmeans",k=5, nstart=25, graph = F)
fviz_cluster(km.res, geom="point",ellipse.type = "norm",palette="jco",ggtheme=theme_minimal())
hc.res <- eclust(df, "hclust",k=5, hc_metric = "euclidean",hc_method = "ward.D2",graph = F)
fviz_dend(hc.res,show_labels = F, palette = "jco",as.ggplot=T)
fviz_silhouette(km.res, palette="jco",ggtheme=theme_minimal())
silinfo <- km.res$silinfo
names(silinfo)
head(silinfo$widths[, 1:3],10)
silinfo$clus.avg.widths
silinfo$avg.width
km.res$size
library(fpc)
km_stats <- cluster.stats(dist(df), km.res$cluster)
km_stats$dunn
km_stats
table(olive$Area, km.res$cluster)
areas <- as.numeric(olive$Area)
clust_stats <- cluster.stats(d = dist(df),
                             areas,
                             km.res$cluster)
clust_stats$corrected.rand
clust_stats$vi
pam.res <- eclust(df, "pam",k=5, graph = F)
table(olive$Area, pam.res$cluster)
cluster.stats(d = dist(df),
                  areas,
                  instkm.res$cluster)$vi
res.hc <- eclust(df, "hclust",k=5,graph = F)
table(olive$Area, pam.res$cluster)
cluster.stats(d = dist(df),
                  areas,
                  res.hc$cluster)$vi
install.packages("clValid")
library(clValid)
clmethods <- c("hierarchical","kmeans","pam")
intern <- clValid(df, nClust = 2:6,
                  clMethods = clmethods, validation = "internal")
summary(intern)
stab <- clValid(df, nClust = 2:6,
                  clMethods = clmethods, validation = "stability")
summary(intern)
optimalScores(stab)
install.packages("pvclust")
library(pvclust)
set.seed(123)
res.pv <- pvclust(df, method.dist = "cor",
                  method.hclust = "average",nboot = 10)
plot(res.pv, hang=-1, cex=0.5)
pvrect(res.pv)
clusters <- pvpick(res.pv)
clusters
install.packages("parallel")
library(parallel)
c1 <- makeCluster(2, type="PSOCK")
res.pv <- parPvclust(c1, df, nboot = 1000)
stopCluster(c1)
res.hk<-hkmeans(df, 5)
fviz_cluster(res.hk, palette="jco",repel = T, ggtheme=theme_classic())
fviz_dend(res.hk, cex=0.6, palette="jco",
          rect=T, rect_border = "jco",rect_fill = T)
res.fanny <- fanny(df, 5)
head(res.fanny$membership, 5)
fviz_silhouette(res.fanny, palette="jco",ggtheme=theme_minimal())
library(mclust)
mc <- Mclust(df)
summary(mc)
mc$modelName
mc$G
head(mc$z, 30)
head(mc$classification, 30)
fviz_mclust(mc, "BIC",palette="jco")
fviz_mclust(mc, "classification",geom="point",palette="jco",poitsize=1.5)
fviz_mclust(mc, "uncertainty",palette="jco")
table(olive[,2],mc$classification)
