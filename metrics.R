confusion <- matrix(sample.int(10, size = 2*2, replace = TRUE), nrow = 2, ncol = 2)

metrics <- function(matrix){
  diag <- c()
  rowsum <- c()
  colsum <- c()
  for (i in 1:nrow(matrix)){
    diag <- append(diag,matrix[i,i])
    rowsum <- append(rowsum,sum(matrix[i,]))
    colsum <- append(colsum,sum(matrix[,i]))
  }
  TP <- diag(matrix)
  FP <- colSums(matrix) - diag(matrix)
  FN <- rowSums(matrix) - diag(matrix)
  TN <- sum(matrix) - (FP + FN + TP)
  accuracy <- sum(TP) / sum(TP + FN)
  sensitivity <- TP / (TP + FN)
  specificity <- TN / (TN+FP) 
  precision <- TP / (TP + FP)
  recall <- TP / (TP + FN)
  beta <- as.numeric(readline(prompt="Beta: "))
  f1 <- ifelse(precision + recall == 0, 0, 2 / (beta * (1/precision) + (beta * 1/recall)))
  f1[is.na(f1)] <- 0
  
  cat("Acc:", accuracy)
  cat("\n")
  cat("Sensitivity:", sensitivity)
  cat("\n")
  cat("Specificity:", specificity)
  cat("\n")
  cat("Precision:", precision)
  cat("\n")
  cat("Recall:", recall)
  cat("\n")
  cat("F1:", f1)
}
