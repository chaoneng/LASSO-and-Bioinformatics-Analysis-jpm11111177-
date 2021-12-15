##' 2021-04-25
##' Generating a volcano plot

rm(list = ls())
Sys.setenv(LANGUAGE = "en")
setwd("C:/Users/flora/Desktop/gynecological cancers") 

library("cowplot")
library("dplyr")
source("./_function/plotVolcano.R")
#-----------------------------------------------------------------#
# Data input
#-----------------------------------------------------------------#

filenames <- list.files(path = "./02_table/", pattern = "_Limma.csv", full.names=TRUE)
All <- lapply(filenames, function(i){
  read.csv(i, header = TRUE, sep = ",", row.names = 1)
})
names(All) <- c("CC", "EC")

CC <- All[["CC"]]
p1 <- volcano_plot(CC, 1, 0.05, "CC")
EC <- All[["EC"]]
p2 <- volcano_plot(EC, 1, 0.05, "EC")

#-----------------------------------------------------------------#
# Combining plots
#-----------------------------------------------------------------#

legend <- get_legend(
  # create some space to the left of the legend
  p1 + theme(legend.box.margin = margin(0, 0, 0, 10)))
prow <- plot_grid(p1 + theme(legend.position="none"), 
                  p2 + theme(legend.position="none"),
                  nrow=1, labels = "AUTO")
pt <- plot_grid(prow, legend, rel_widths = c(.8,.15))
pt

png("./03_figure/volcano.png",    # create PNG for the heat map        
    width = 12*800,        # 5 x 300 pixels
    height = 5*800,
    res = 800)             # 300 pixels per inch 
pt
dev.off()
