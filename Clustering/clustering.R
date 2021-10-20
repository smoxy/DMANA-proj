load("Z:/DesktopC/LUMSA/2/Data Mining/Clustering/ClusterData_L31.RData")

preProcessing <- function(df){
  df <- df[, c(-1,-2)]
  df = data.frame(scale(df))
  return(df)
}

catReg <- function(df){
  #areaNames <- c("NApu","Cala","SApu","Sici","ISar","CSar","ELig","WLig","Umbr")
  regionNames <- c("Sud","Isola","Nord")
  for (i in unique(df$Region)){
    df$Region[df$Region == i] <- regionNames[i]
  }
  return(df)
}
stdOlive <- preProcessing(olive)
olive <- catReg(olive)

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

################################################################################
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
