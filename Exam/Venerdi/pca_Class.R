# PCA
res<-PCA(pbcseq2, quanti.sup=c(5,7), quali.sup=c(1,2,3,4,6,8,17,18,19))
fviz_screeplot(res, addlabels = TRUE, ylim = c(0, 32))
summary(res, nbelements=10)
fviz_pca_ind(res,
             col.ind="cos2",
             geom = c("point","text"),
             gradient.cols = c("#94bdff","#1b7fcc",'#00062e'),
             repel=T,
             select.ind= list(cos2 = .8),
             ggtheme=theme_minimal())

fviz_pca_var(res, col.var = "cos2",
             gradient.cols = c("#94bdff","#1b7fcc",'#00062e'),
             select.var= list(cos2 = 10),
             repel = T,
             ggtheme = theme_minimal())

# DIM 1 = Funzionalità epatica
# DIM 2 = Capacità di coagulazione del sangue

fviz_pca_biplot(res, label = "var", habillage=pbcseq2$status,
                addEllipses=TRUE, ellipse.level=0.90,
                ggtheme = theme_minimal(),
                select.var= list(cos2 = 10))

fviz_pca_var(res, col.var = "cos2",
             axes = c(3,4),
             gradient.cols = c("#94bdff","#1b7fcc",'#00062e'),
             select.var= list(cos2 = 10),
             repel = T,
             ggtheme = theme_minimal())

pc.comp <- res$ind$coord
pc.comp1 <- pc.comp[,1]
pc.comp2 <- pc.comp[,2]
X = cbind(pc.comp1, pc.comp2)
X = as.data.frame(cbind(status=pbcseq2$status, X))
X[which(X$status==0),1] <- "Censored.Transplant"
X[which(X$status==1),1] <- "Dead"
X$status <- as.factor(X$status)

intrain <- createDataPartition(y = pbcseq2$status, p=0.7, list = FALSE)

trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3, classProbs=T, summaryFunction=twoClassSummary)
set.seed(1234)
knn_fit.PC12 <- caret::train(status ~ ., data = X[intrain,], method = "rf",
                             trControl=trctrl,
                             preProcess = c("center", "scale"))
print(knn_fit.PC12)
plot(knn_fit.PC12)
test_pred <- predict(k.means, newdata = X[-intrain,])
confusionMatrix(test_pred, as.factor(X[-intrain,'status']))
