###################################LOAD#########################################
source("https://raw.githubusercontent.com/MatteoFasulo/Rare-Earth/main/R/util/coreFunctions.R")
loadPackages(c('pacman','ggplot2','sn','corrplot','FactoMineR','FactoInvestigate','Factoshiny','factoextra','caret'))
pacman::p_loaded()
################################################################################
data(wines)
################################################################################
corMatrix <- cor(wines[,-1])
plotCorr <- function(corMatrix){
  col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
  return(corrplot(corMatrix,
                  type="upper",
                  order="FPC",
                  sig.level = 0.01,
                  insig = "blank",
                  title = 'Correlation Plot',
                  tl.cex = .8,
                  mar=c(0,0,1,0)))
}
plotCorr(corMatrix)
################################################################################
year <- as.numeric(substr(rownames(wines), 6, 7))
wines[,'year'] <- year

res <- PCA(wines, quanti.sup=NULL, quali.sup=c(1,29))

fviz_screeplot(res, addlabels = TRUE, ylim = c(0, 26))

summary(res, nbelements=Inf)

fviz_pca_ind(res,
             col.ind="contrib",
             geom = c("point","text"),
             gradient.cols = c("#94bdff","#1b7fcc",'#00062e'),
             repel=T,
             select.ind= list(cos2 = .6),
             ggtheme=theme_minimal())

fviz_pca_var(res, col.var = "contrib",
             gradient.cols = c("#94bdff","#1b7fcc",'#00062e'),
             select.var= list(cos2 = 10),
             repel = T,
             ggtheme = theme_minimal())

fviz_pca_biplot(res, label = "var", habillage=wines$wine,
                addEllipses=TRUE, ellipse.level=0.90,
                ggtheme = theme_minimal())

plotellipses(model=res, autoLab = "auto")

fviz_pca_biplot(res, axes=c(3,4), label = "var", habillage=wines$wine,
                addEllipses=TRUE, ellipse.level=0.90,
                ggtheme = theme_minimal())

plotellipses(model=res, axes=c(3,4), autoLab = "auto")
################################################################################
set.seed(3033)
intrain <- createDataPartition(y = wines$wine, p=0.7, list = FALSE)
training <- wines[intrain,]
testing <- wines[-intrain,]
dim(training); dim(testing);
training[["wine"]] = factor(training[["wine"]])
trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
set.seed(3333)
knn_fit <- train(wine ~., data = training, method = "knn",
                 trControl=trctrl,
                 preProcess = c("center", "scale"),
                 tuneLength = 10)
knn_fit
plot(knn_fit)
test_pred <- predict(knn_fit, newdata = testing)
confusionMatrix(test_pred, testing$wine)



#dimdesc(res)
#dimdesc(res, proba=0.2)
#Investigate(res)
#PCAshiny(res)
