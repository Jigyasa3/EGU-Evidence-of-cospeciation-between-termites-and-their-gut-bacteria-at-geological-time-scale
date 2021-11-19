#!/usr/bin/env Rscript

#USAGE- Rscript [tree1] [output_file1] [output_file2]
args <- commandArgs(TRUE)
symbionttree <- args[1]
outputtree<-args[2] #symbionttree-${i}
outputfile<-args[3] #symbiontheader-${i}

##step1-rename symbiont tree file with hostname_1
library(ape)
symbionttree<-read.tree(file=symbionttree)

d<-data.frame(label=symbionttree$tip.label)
#d$label2<-gsub("d__.*$","",d$label)
#d$label2 <- sub("--[^--]+$", "", d$label2)
d$runnumber<-gsub("--.*$","",d$label)

hosttree<-read.tree("/flash/BourguignonU/Jigs/markergene_paper_scripts/scripts/rna12-16S_202samples_newnames_feb2020.nwk")
e<-data.frame(label=hosttree$tip.label)
e$runnumber<-gsub("--.*$","",e$label)

d2<-merge(d,e,by="runnumber")

library(dplyr)
d2<-d2 %>%group_by(label.y) %>%mutate(onemore = paste(label.y,row_number(),sep="_"))

##output- newtree, header
library(ggtree)
library(treeio)
symbionttree2<-rename_taxa(symbionttree, d2, label.x, onemore)
write.tree(symbionttree2,file=outputtree)

d2<-as.data.frame(d2)
d3<-d2%>%select(onemore)
write.csv(d3,file=outputfile)

##step2-get sequences ids out-

#awk -F"," '{print $2}' symbiontheader-d__Bacteria_p__Spirochaetota__p__UBP6__p__Deferribacterota__COG0552_tips_4.txt | sed 's/"//g' > 2-symbiontheader-d__Bacteria_p__Spirochaetota__p__UBP6__p__Deferribacterota__COG0552_tips_4.txt
#or
#for i in symbiontheader-d__*;do awk -F"," '{print $2}' ${i} | sed 's/"//g' > 2-${i};done

##step3-run "newhosttree.R" script to get hosttree with "0" branch lengths equal to symbionttree-
#IN_DIR="/flash/BourguignonU/Jigs/markers/cophylogeny"
#Rscript ${IN_DIR}/newhosttree.R 2-symbiontheader-d__Bacteria_p__Spirochaetota__p__UBP6__p__Deferribacterota__COG0552_tips_4.txt ${IN_DIR}/rna12-16S_202samples_newnames_feb2020.nwk hosttree-d__Bacteria_p__Spirochaetota__p__UBP6__p__Deferribacterota__COG0552_tips_4.nwk
#or
#for i in 2-symbiontheader-d__*.txt;do filename=`echo ${i}| sed 's/2-symbiontheader-//g'`; filename2=`echo ${filename}| sed 's/.txt//g'`; Rscript ${IN_DIR}/newhosttree.R ${i} ${IN_DIR}/rna12-16S_202samples_newnames_feb2020.nwk hosttree-${filename2};done
