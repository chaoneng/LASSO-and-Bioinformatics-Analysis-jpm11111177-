## 2021-04-25
## Creating a vennDiagram of the DE genes

rm(list = ls())
Sys.setenv(LANGUAGE = "en")
setwd("C:/Users/flora/Desktop/gynecological cancers") 

#-----------------------------------------------------------------#
# Data input
#-----------------------------------------------------------------#

library("dplyr")

filenames <- list.files(path = "./02_table/", pattern = "_Limma.csv", full.names=TRUE)
All <- lapply(filenames, function(i){
  read.csv(i, header = TRUE, sep = ",", row.names = 1)
})
names(All) <- c("CC", "EC")
CC <- All[["CC"]]
EC <- All[["EC"]]

library("VennDiagram")
x <- list(
  CC %>% filter(abs(logFC)>1&adj.P.Val<0.05) %>% rownames() %>% unlist(),
  EC %>% filter(abs(logFC)>1&adj.P.Val<0.05) %>% rownames() %>% unlist())
# intersect(x[[1]],x[[2]])

grid.newpage()
draw.pairwise.venn(920, 843, 78, 
                   category = c("CC vs. normal", "EC vs. normal"), 
                   lty = rep("blank", 2), 
                   fill = c("red", "green"), 
                   alpha = rep(0.6, 2),
                   cex = 1.2, # Vector giving the size for each area label
                   cat.cex = 1, # Vector giving the size for each category name
                   cat.fontface = 2, # Vector giving the fontface for each category name
                   cat.pos = c(0, 0),
                   cat.dist = rep(0.025, 2), 
                   scaled = TRUE)


png("./03_figure/venn.png",        
    width = 5*800,        # 5 x 300 pixels
    height = 5*800,
    res = 900)             # 300 pixels per inch 
draw.pairwise.venn(920, 843, 78, 
                   category = c("CC vs. normal", "EC vs. normal"), 
                   lty = rep("blank", 2), 
                   fill = c("red", "green"), 
                   alpha = rep(0.6, 2),
                   cex = 1.2, # Vector giving the size for each area label
                   cat.cex = 1, # Vector giving the size for each category name
                   cat.fontface = 2, # Vector giving the fontface for each category name
                   cat.pos = c(0, 0),
                   cat.dist = rep(0.025, 2), 
                   scaled = TRUE)
dev.off()
