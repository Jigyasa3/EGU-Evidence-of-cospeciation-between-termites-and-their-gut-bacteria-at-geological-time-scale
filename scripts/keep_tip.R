#!/usr/bin/env Rscript

#USAGE- Rscript [tree] [clade] [outputfile]
args <- commandArgs(TRUE)
tree <- args[1]
clade <- args[2]
outputfile <- args[3]

library(ape)
tree<-read.tree(file=tree)

clade<-read.csv(file=clade)
colnames(clade)<-c("X","V1")
#clade$V2<-gsub("$","_1",clade$V1) #if the headers are based on nt.sequence file. They lack "_1" at the end.

library(dplyr)
library(stringr)
#remove GTDB ids from termite specific clades-
clade<-clade%>%mutate(toremove=ifelse(str_detect(V1,"^RS_GCF|^GB_GCA"),"remove",as.character(V1)))
clade2<-clade%>%filter(toremove!="remove")

tree2<-keep.tip(tree,as.vector(clade2$V1))

write.tree(tree2,file=outputfile)

