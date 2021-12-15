##' 2021-04-24
##' Downloading of gene expression profiling (microarray) for cervical cancer

rm(list = ls())
Sys.setenv(LANGUAGE="en")
setwd("C:/Users/flora/Desktop/gynecological cancers") 

source("./_function/doRmaSymbol.R")
source("./_function/rmDupID.R")

listSamples <- list(Tpyes = rep(c("cervical_cancer", "endometrial_cancer"), c(2,3)),
                    GEO = c("GSE9750", "GSE7803", "GSE63514",
                            "GSE17025", "GSE115810", "GSE36389"),
                    Platforms = c( "GPL96", "GPL96", "GPL570",
                                  "GPL570", "GPL96", "GPL96"))

for(i in 1:6){
  print(listSamples$Tpyes[i])
  print(listSamples$GEO[i])
  print(listSamples$Platforms[i])
  do_RMA_Symbol(listSamples$Tpyes[i], listSamples$GEO[i], listSamples$Platforms[i])
}

do_RMA_Symbol("endometrial_cancer", "GSE36389", "GPL96")
