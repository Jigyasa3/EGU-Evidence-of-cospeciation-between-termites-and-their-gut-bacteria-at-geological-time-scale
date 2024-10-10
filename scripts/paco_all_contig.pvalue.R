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
rownames(HE.matrix2) <- gsub(x = rownames(HE.matrix2), pattern = "--.*", replacement = "")
names(HE.matrix2) <- gsub(x = names(HE.matrix2), pattern = "--.*", replacement = "")
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
D<-prepare_paco_data(H=H.matrix,P=E.matrix,HP=out)
D<-add_pcoord(D,correction="cailliez")
D<-PACo(D,nperm=10000,seed=12,method="backtrack",symmetric=FALSE) #"r0" algorithm is used when host maintains the symbionts evolution. If not known use "backtracking" or "swaps"
#symmetric=FALSE: one group is not assumed to track the evolution of the other.
D<-paco_links(D) #a jackknife procedure to estimate the degree of individual interactions
res<-as.data.frame(residuals_paco(D$proc)) #residuals of each interactions
links<-as.data.frame(D$jackknife)
D.pvalue<-D$gof #output D$gof$ss ->m^2xy=56.92 #gives the p-value of overall phylogenetic coevolution.

write.csv(D.pvalue,file=output_file1)

#----------------------------------------------------------
#NOTE-OR run the script "paco2.R" if the above script gives the following error-"Error in if (nrow(X) != nrow(Y)) stop(gettextf("matrices have different number of rows: %d and %d",  : 
  argument is of length zero"
