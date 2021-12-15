##' 2021-04-25
##' Removing duplicated gene symbols

rm_Dup_ID <- function(matrix){
  expreSet <- matrix
  rowMean <- apply(expreSet[,-1], 1, function(x){mean(as.numeric(x), na.rm = TRUE)})
  expreSet <- expreSet[order(rowMean,decreasing = TRUE),]
  expreSet <- expreSet[!duplicated(expreSet[,1]),]
  expreSet <- expreSet[!expreSet[,1]=="",]
  rownames(expreSet) <- expreSet[,1]
  expreSet <- expreSet[,-1]
  return(expreSet)
}
