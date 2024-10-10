#!/usr/bin/env Rscript

#USAGE- Rscript [tree] [termitetree] [output_file1]
args <- commandArgs(TRUE)
tree <- args[1]
termitetree<-args[2]
output_file1 <- args[3]


library(dplyr)
library(tidyr)
library(plyr)
library(ape)

#read the symbiont tree-
symbiont_tree<-read.tree(file=tree)
symbiont_tree2<-symbiont_tree


e<-data.frame(label=symbiont_tree2$tip.label)
e$label<-as.character(e$label)

#read the termite tree
termite_tree_new<-read.tree(file=termitetree)
#is.binary(termite_tree_new) # true
#is.ultrametric(termite_tree_new) #true
#is.rooted(termite_tree_new) #true

d<-data.frame(label=termite_tree_new$tip.label)
d$label<-as.character(d$label)

library(vegan)

H.matrix<-cophenetic(termite_tree_new)
E.matrix<-cophenetic(symbiont_tree2)

#generate the HE.matrix-
HE.matrix <- data.frame(matrix(ncol = nrow(e), nrow = nrow(d))) #columns=nrow(e), rows=nrow(d)
rownames(HE.matrix) <- d$label
colnames(HE.matrix)<-e$label

HE.matrix2<-HE.matrix
#remove the gene names from column names (everything after "--")
#rownames(HE.matrix2) <- gsub(x = rownames(HE.matrix2), pattern = "_.*$", replacement = "")
#names(HE.matrix2) <- gsub(x = names(HE.matrix2), pattern = "_.*$", replacement = "")
#names(HE.matrix2) <- gsub(x = names(HE.matrix2), pattern = "_", replacement = "-")


#match rownames and column names-
out <- outer(row.names(HE.matrix2), colnames(HE.matrix2), `==`) #find cells that are common
dimnames(out) <- dimnames(HE.matrix2)
out<-as.data.frame(out)
colnames(out)<-colnames(HE.matrix) #change the column names to original tip.labels
rownames(out)<-rownames(HE.matrix)
out <- as.matrix(1*out) #convert logical to numeric and into a matrix


#-------------------------------------------------------
#run paco-

library(paco)
#ver2-https://www.uv.es/cophylpaco/index.html (if ver1 gives an error)

PACo.dV <- function (H.dist, P.dist, HP.bin) {  
  HP.bin <- which(HP.bin > 0, arr.in=TRUE) 
  H.PCo <- pcoa(sqrt(H.dist), correction="none")$vectors  
  P.PCo <- pcoa(sqrt(P.dist), correction="none")$vectors  
  H.PCo <- H.PCo[HP.bin[,1],]  
  P.PCo <- P.PCo[HP.bin[,2],]  
  list (H.PCo = H.PCo, P.PCo = P.PCo) } 

CP.D=E.matrix
N.D=H.matrix
NCP=out
PACo.fit <- PACo.dV(N.D, CP.D, NCP)
NCP.proc <- procrustes(PACo.fit$H.PCo, PACo.fit$P.PCo) 
m2.obs <- NCP.proc$ss  
N.perm = 1000
P.value = 0 
set.seed(2)  
NLinks = sum(NCP) 
HP <- diag(NLinks) 

for (n in c(1:N.perm)) 
  { if (NLinks <= nrow(NCP) | NLinks <= ncol(NCP))     
    {  flag2 <- TRUE  
    while (flag2 == TRUE)  
      {  NCP.perm <- t(apply(NCP,1,sample)) 
      if(any(colSums(NCP.perm) == NLinks)) flag2 <- TRUE else flag2 <- FALSE 
      }   
    } else { NCP.perm <- t(apply(NCP,1,sample))}
    PACo.perm <- PACo.dV(N.D, CP.D, NCP.perm) 
    m2.perm <- procrustes(PACo.perm$H.PCo, PACo.perm$P.PCo)$ss  
    if (m2.perm <= m2.obs) {P.value = P.value + 1}  
  } 
P.value <- P.value/N.perm 


write.csv(P.value,file=output_file1)
