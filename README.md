# EGU-Evidence-of-cospeciation-between-termites-and-their-gut-bacteria-at-geological-time-scale

This repository contains the scripts to run the analysis performed for the manuscript *Evidence-of-cospeciation-between-termites-and-their-gut-bacteria-at-geological-time-scale*

The data analysis occurs in the following steps-

##step1. Root the tree with outgroup bacterial sequences-
```module load R/3.6.1```

#<Rscript root_the_tree.R treefile outgroupname outputfile>

```
 Rscript root_the_tree.R ../data/ry_GTR_ry_protein-nucleo-uniq-d__Bacteria_p__Firmicutes_A__p__Fusobacteriota__COG0012.fasta.treefile Fusobacteriota ../data/output/rooted-ry_GTR_ry_protein-nucleo-uniq-d__Bacteria_p__Firmicutes_A__p__Fusobacteriota__COG0012.fasta.treefile

 Rscript root_the_tree.R ../data/ry_GTR_ry_protein-nucleo-uniq-d__Bacteria_p__Firmicutes_A__p__Fusobacteriota__COG0552.fasta.treefile Fusobacteriota ../data/output/rooted-ry_GTR_ry_protein-nucleo-uniq-d__Bacteria_p__Firmicutes_A__p__Fusobacteriota__COG0552.fasta.treefile

  Rscript root_the_tree.R ../data/ry_GTR_ry_protein-nucleo-uniq-d__Bacteria_p__Spirochaetota__p__UBP6__p__Deferribacterota__COG0012.fasta.treefile UBP6 ../data/output/rooted-ry_GTR_ry_protein-nucleo-uniq-d__Bacteria_p__Spirochaetota__p__UBP6__p__Deferribacterota__COG0012.fasta.treefile

  Rscript root_the_tree.R ../data/ry_GTR_ry_protein-nucleo-uniq-d__Bacteria_p__Spirochaetota__p__UBP6__p__Deferribacterota__COG0552.fasta.treefile UBP6 ../data/output/rooted-ry_GTR_ry_protein-nucleo-uniq-d__Bacteria_p__Spirochaetota__p__UBP6__p__Deferribacterota__COG0552.fasta.treefile
```
  
##step2. Get termite-specific clusters with atleast 10 termite species-
###A-If there is one outgroup sequences-
```
  Rscript termitespecificclades_10species.R ../data/output/rooted-ry_GTR_ry_protein-nucleo-uniq-d__Bacteria_p__Firmicutes_A__p__Fusobacteriota__COG0012.fasta.treefile p__Fusobacteriota ../data/output/clade-Firmicutes_A-COG0012.txt ../data/output/Firmicutes_A-COG0012_tips

  Rscript termitespecificclades_10species.R ../data/output/rooted-ry_GTR_ry_protein-nucleo-uniq-d__Bacteria_p__Firmicutes_A__p__Fusobacteriota__COG0552.fasta.treefile p__Fusobacteriota ../data/output/clade-Firmicutes_A-COG0552.txt ../data/output/Firmicutes_A-COG0552_tips
```
  
###B-If there are two outgroup sequences-
```
  Rscript termitespecificclades_10species_2sisters.R ../data/output/rooted-ry_GTR_ry_protein-nucleo-uniq-d__Bacteria_p__Spirochaetota__p__UBP6__p__Deferribacterota__COG0012.fasta.treefile p__UBP6 p__Deferribacterota ../data/output/clade-Spirochaetota-COG0012.txt ../data/output/Spirochaetota-COG0012_tips

  Rscript termitespecificclades_10species_2sisters.R ../data/output/rooted-ry_GTR_ry_protein-nucleo-uniq-d__Bacteria_p__Spirochaetota__p__UBP6__p__Deferribacterota__COG0552.fasta.treefile p__UBP6 p__Deferribacterota ../data/output/clade-Spirochaetota-COG0552.txt ../data/output/Spirochaetota-COG0552_tips
```
  
##step3. Get the subtree corresponding to extracted tips
```
  for i in Spirochaetota-COG0012_tips*csv;do filename=`echo ${i}| sed 's/_.csv//g'` ; Rscript ../../scripts/keep_tip.R rooted-ry_GTR_ry_protein-nucleo-uniq-d__Bacteria_p__Spirochaetota__p__UBP6__p__Deferribacterota__COG0012.fasta.treefile ${i} ${filename}.nwk;done
```
  
