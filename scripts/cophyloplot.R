#!/usr/bin/env Rscript
#USAGE- Rscript
args <- commandArgs(TRUE)
tree <- args[1] #tips_tree
output <- args[2] #pdf format

##Generate cophyloplot in pdf format-

library(ape)
library(dplyr)
library(stringr)

termite_tree0<-read.tree("/flash/BourguignonU/Jigs/markergene_paper_scripts/scripts/rna12-16S_202samples_newnames_feb2020.nwk") #load the old tree
d<-data.frame(label=termite_tree0$tip.label)
d$runnumber<-unlist(lapply(strsplit(as.character(d$label),split="--"),"[",1))

#read the sample info-
sample_names<-read.csv("/flash/BourguignonU/Jigs/markergene_paper_scripts/scripts/sample_names_ids_withotherfactors_tomedit.202samples.csv",header=TRUE)
sample_names2<-sample_names%>%dplyr::select(runnumber,termitespecies2)
#rename termite tree-
library(treeio)
d<-merge(d,sample_names2,by="runnumber")

termite_tree<-rename_taxa(termite_tree0,d,label,termitespecies2)
d<-data.frame(label=termite_tree$tip.label)
d$runnumber<-unlist(lapply(strsplit(as.character(d$label),split="--"),"[",1))

#read the lineage2 colors used in chapter1-
lineagecolors<-read.csv("/flash/BourguignonU/Jigs/markergene_paper_scripts/scripts/lineage_itolcolors_202samples.csv")
lineagecolors<-lineagecolors%>%dplyr::select(runnumber,V2)
colnames(lineagecolors)<-c("runnumber","color")
#----------------------------------------
#tree1-
#read the symbiont tree-
symbiont_tree_1<-read.tree(file=tree) ##this is original file!symbiont_tree_1<-read.tree(file="ntsequences/nt_sequences_alltrees_june2021/rep_cogs_trees/d__Archaea_p__Methanobacteriota__p__Thermoplasmatota__COG0552_tips_1.nwk")
e<-data.frame(label=symbiont_tree_1$tip.label)
e$runnumber<-gsub("--.*$","",e$label)
e$proteins<-gsub("^.*--","",e$label)
e$proteins<-gsub("_d__.*$","",e$proteins)
e$taxa<-gsub("^.*d__","",e$label)
e$phyla<-gsub("_c__.*$","",e$taxa)
e$phyla<-gsub("^.*p__","p_",e$phyla)
e$class<-gsub("_o__.*$","",e$taxa)
e$class<-gsub("^.*c__","c_",e$class)
e$order<-gsub("_f__.*$","",e$taxa)
e$order<-gsub("^.*o__","o_",e$order)
e$family<-gsub("_g__.*$","",e$taxa)
e$family<-gsub("^.*f__","f_",e$family)
e$genus<-gsub("_s__.*$","",e$taxa)
e$genus<-gsub("^.*g__","g_",e$genus)

e$genus<-gsub("^.*f__","f_",e$genus)
e$genus<-gsub("^.*o__","o_",e$genus)
e$genus<-gsub("^.*c__","c_",e$genus)
e$genus<-gsub("_1$","",e$genus)

#rename symbiont tree-
e<-merge(e,sample_names2,by="runnumber")
e$newlabel<-paste(e$termitespecies2,e$proteins,e$genus,sep="--")

symbiont_tree_2<-rename_taxa(symbiont_tree_1,e,label,newlabel)

#read the PACO output file-
#res_1<-read.csv(file="ntsequences/nt_sequences_alltrees_june2021/rep_cogs_trees/significant_tips-d__Bacteria_p__Spirochaetota__p__UBP6__p__Deferribacterota__COG0552_tips_4.nwk",header=TRUE)
#colnames(res_1)<-c("samplenames","paco_residues")

#merge "res" df with d2-
#res_1$samplenames<-as.character(res_1$samplenames)
#res_1$runnumber<-unlist(lapply(strsplit(res_1$samplenames,split="--"),"[",1))
#res2_1<-merge(res_1,sample_names,by.x="runnumber",by.y="runnumber")  #get the host names "label" column.

#get bacterial names out in a "vector"
#res2_1$binname<-gsub("^.*?_","",res2_1$samplenames) #get everything after first underscore
#res2_2<-merge(res2_1,e,by.x="binname",by.y="label")

#extract termite species from host tree with representatives in symbiont tree-
d2<-merge(d,e,by="runnumber")
termitetree_1<-keep.tip(termite_tree,as.vector(unique(d2$label.x)))

#add colored lines corresponding to termite lineages
#res2_2<-merge(res2_2,lineagecolors,by.x="runnumber.x",by.y="runnumber")
d2<-merge(d2,lineagecolors,by.x="runnumber",by.y = "runnumber")
d2$color<-as.character(d2$color)

assoc1<-d2%>%select(label.x,newlabel)
##plot1-phytools-
#library(phytools)
source("/flash/BourguignonU/Jigs/markergene_paper_scripts/scripts/cophylo.R")
obj<-cophylo(termitetree_1,symbiont_tree_2,assoc = assoc1)
#obj<-cophylo(termitetree_1,symbiont_tree_2,assoc = assoc1,rotate=FALSE)
#obj<-cophylo(termitetree_1,symbiont_tree_2,assoc = assoc1,rotate.multi=TRUE,methods=c("pre","post"))
#plot(obj,part=0.4,fsize=c(0.5,0.5),link.lwd=2,link.lty="solid",link.col=d2$color)

#plot2-ape-https://www.rdocumentation.org/packages/ape/versions/4.1/topics/cophyloplot
#ape::cophyloplot(termitetree_1, symbiont_tree_2, assoc = assoc1, use.edge.length = TRUE, space = 50,
#                 length.line = 1, gap = 2, type = "phylogram", rotate = FALSE,
#                 col = res2_2$color, lwd = 2, lty = par("lty"),
#                 show.tip.label = TRUE, font = 0.6)

pdf(file=output,height=7,width=14)
plot(obj,part=0.4,fsize=c(0.6,0.6),link.lwd=2,link.lty="solid",link.col=d2$color)
dev.off()

#svg(file=output,height=9,width=13)
#plot(obj,part=0.4,fsize=c(0.6,0.6),link.lwd=2,link.lty="solid",link.col=d2$color)
#dev.off()
