
DEG_anlysis <- function(expr, pdata, type){
  
  if(type=="CC"){
    cat("Identification of DE genes between cervical cancer and normal tissues\n")
  }else{
    cat("Identification of DE genes between endometrial cancer and normal tissues\n")
  }
  cat("Sample size: ", ncol(expr), "\n")
  cat("The number of genes: ", nrow(expr), "\n")
  
  colnames(expr) <- sapply(strsplit(colnames(expr), "_"), "[", 1)
  
  # Zscore transformation
  mydata <- t(scale(t(expr), center = TRUE, #for cols (genes)
                    scale = TRUE))
  mydata <- mydata[,match(colnames(mydata), pdata$geo_accession)]
  
  # Using Limma package to identify differential expressed genes
  
  require("limma")
  condition <- factor(pdata$tissues)
  design <- model.matrix(~0+condition)
  colnames(design) <- levels(condition)
  rownames(design) <- pdata$geo_accession
  contrast.matrix <- makeContrasts(Cancer-Normal, levels = colnames(design))
  
  # lmFit() method
  # This method will fit a linear model (defined in design) to the data 
  # to calculate the mean expression level in the subgroup.
  fit <- lmFit(mydata, design)
  
  # contrasts.fit() method
  # Now you have to tell limma which groups you want to compare.
  fit1 <- contrasts.fit(fit, contrast.matrix)
  
  # eBayes() method (it has performed a moderated t-test on each gene.)
  # Now limma is ready to perform the statistical test to compare the groups.
  fit2 <- eBayes(fit1)
  # fit2$coefficients # retrieve the log fold changes of each gene 
  list <- topTable(fit2, n=nrow(fit2), adjust="BH")
  write.csv(list, paste0("./02_table/DEG_", type, "_Limma.csv"))
  cat("The DEG analysis is done.")
}