##' 2021-04-24
##' Preparing Gene Expression Data: 
##' Use affy package to do RMA of the CEL files:
require("affy")
require("dplyr")
require("stringi")
require("GEOquery")

do_RMA_Symbol <- function(cancer_type, accession_number, platforms){
  untar(paste0("./01_data/_raw_", cancer_type, "/", accession_number, "_RAW.tar"), 
        exdir = paste0("./01_data/_raw_", cancer_type, "_extracted_files/", accession_number, "/"))

  data <- justRMA(celfile.path = paste0("./01_data/_raw_", cancer_type, "_extracted_files/", accession_number, "/"))
  
  #Extracting gene expression matrix:
  exprSet <- exprs(data)
  colnames(exprSet) <- strsplit(colnames(exprSet), "[_]") %>% sapply(., "[",1)

  gse <- getGEO(accession_number, GSEMatrix=TRUE)
  pdata <- pData(gse[[1]])
  pdata <- pdata[pdata$geo_accession%in%colnames(exprSet),]
  write.csv(pdata, paste0("./01_data/", accession_number, "_clinical.csv"))
  
  gpl <- getGEO(platforms, destdir=paste0("./01_data/_raw_", cancer_type, "_extracted_files/", accession_number, "/"))
  featuredata <- Table(gpl)
  probe2Symbol <- featuredata[,c("ID","Gene Symbol")]
  rownames(probe2Symbol) <- probe2Symbol[,"ID"]
  probe2Symbol$"Gene Symbol 2" <- probe2Symbol$`Gene Symbol` %>% 
    strsplit(., " /// ") %>% 
    sapply(., "[",1)
  exprSet <- cbind("Gene Symbol"=probe2Symbol[match(rownames(probe2Symbol), rownames(exprSet)),3],exprSet)
  exprSet <- exprSet[!is.na(exprSet[,1]),] 
  
  mydata <- rm_Dup_ID(exprSet)
  write.table(mydata, paste0("./01_data/", accession_number, "_rma_symbol.txt"), sep = "\t", col.names = NA)
}




