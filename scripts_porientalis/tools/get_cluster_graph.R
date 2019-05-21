#!/usr/bin/Rscript

list.of.packages<-c("tidyverse")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
file1 = as.character(args[1])
file2 = as.character(args[2])

cluster_table=read.delim(args[1], header=FALSE) 
levels_plot=scan(args[2], character(), quote = '')

pdf( file="cluster_heatmap.pdf", width = 12, height = 10 )

ggplot(data =cluster_table, aes(x=factor(V1, level=levels_plot), y=factor(V2, level=levels_plot), fill=V3)) + geom_tile(colour="white",size=0.25) + theme_minimal() + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + theme(axis.title.x=element_blank(), axis.title.y=element_blank())

