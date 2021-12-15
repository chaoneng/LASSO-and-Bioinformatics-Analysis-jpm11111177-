

listSamples <- c("GSE9750", "GSE7803", "GSE63514",
                 "GSE17025", "GSE115810", "GSE36389")
filenames <- NULL
for(i in listSAmples){
  filename <- paste0("./01_data/", i, "_clinical.csv")
  filenames <- c(filenames, filename)
}

All <- lapply(filenames, function(i){
  read.csv(i, header = TRUE, sep = ",", row.names = 1)
})
names(All) <- listSAmples