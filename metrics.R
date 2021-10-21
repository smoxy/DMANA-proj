confusion <- matrix(sample.int(10, size = 4*4, replace = TRUE), nrow = 4, ncol = 4)

metrics <- function(matrix){
  diag <- c()
  rowsum <- c()
  colsum <- c()
  for (i in 1:nrow(matrix)){
    diag <- append(diag,matrix[i,i])
    rowsum <- append(rowsum,sum(matrix[i,]))
    colsum <- append(colsum,sum(matrix[,i]))
  #cat(paste( ,sep="\n"))
  #return(df)
  #cat(paste("Accuracy:",diag/rowsum))
  }
  #print(sum(diag))
  print("Acc")
  print(sum(diag)/sum(rowsum))
}

