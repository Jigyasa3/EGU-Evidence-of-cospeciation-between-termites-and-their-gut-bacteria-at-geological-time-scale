#!/usr/bin/env Rscript

#USAGE- Rscript [tree1] [tree2] [output_file1] [output_file2]
args <- commandArgs(TRUE)
tree1 <- args[1]
tree2<-args[2]
output_file1 <- args[3]
output_file2<- args [4]

## run "before_treedist.R" to generate the tree1 and tree2 files
library(ape)
tree1<-read.tree(file=tree1)
tree2<-read.tree(file=tree2)

#install.packages("rlang")
library(rlang) #0.4.10 version
library(devtools)
library(rJava) #0.9-13
library(htmltools) #0.5.1.1
#install.packages("TreeDist")
library(TreeDist)
library(TreeTools)
#install_github("ms609/TreeDistData")
library(TreeDistData)
library(TreeSearch)

#method1-generate random trees for the host tree-<USING GENERALIZED RF method>
tree1<-unroot(tree1) #based on https://github.com/ms609/TreeDist/issues/58
nRep <- 100000 # Use more replicates for more accurate estimate of expected value
randomTrees <- lapply(logical(nRep), function (x) RandomTree(tree1$tip.label))
randomDists <- ClusteringInfoDistance(tree1, randomTrees, normalize = TRUE)
expectedCID <- mean(randomDists)


dist12 <- ClusteringInfoDistance(tree1, tree2, normalize = TRUE)
# Now count the number of random trees that are this similar to tree1
nThisSimilar <- sum(randomDists < dist12)
pValue <- nThisSimilar / nRep

write.csv(pValue,file=output_file1)

#----------------------------------------------------------------------------------
#method2-generate random trees for the host tree-<USING NYE method>
nRep <- 100000 # Use more replicates for more accurate estimate of expected value
randomTrees <- lapply(logical(nRep), function (x) RandomTree(tree1$tip.label))
randomDists <- NyeSimilarity(tree1, randomTrees, normalize = TRUE,similarity = FALSE)
expectedCID <- mean(randomDists)


dist12 <- NyeSimilarity(tree1, tree2, normalize = TRUE,similarity = FALSE)
# Now count the number of random trees that are this similar to tree1
nThisSimilar <- sum(randomDists < dist12)
pValue2 <- nThisSimilar / nRep

write.csv(pValue2,file=output_file2)
