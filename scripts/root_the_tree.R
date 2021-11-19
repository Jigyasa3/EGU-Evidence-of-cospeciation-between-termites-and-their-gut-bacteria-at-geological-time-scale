#!/usr/bin/env Rscript

#USAGE- Rscript 
args <- commandArgs(TRUE)
tree <- args[1] #original markergene tree
outgroup <- args[2] #eg-"Fibrobacter"
output1 <- args[3] #rooted tree

##on markergene trees with sister groups-
library(ape)
library(dplyr)
library(treeio)
library(stringr)

#1. root the trees on the first tip corresponding to the outgroup phyla-
cogtree<-read.tree(file=tree)
d<-data.frame(label=cogtree$tip.label)

outgroupname<-d %>%dplyr::filter(str_detect(label, as.vector(outgroup)))
outgroupname[1,]

cogtree2<-root(cogtree,as.vector(outgroupname[1,]))
is.rooted(cogtree2)
write.tree(cogtree2,file=output1)
#----------------------------------------------------------------------------------------------
