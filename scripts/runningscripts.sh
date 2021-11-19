
##step1. Root the tree with outgroup bacterial sequences-
module load R/3.6.1
#<Rscript root_the_tree.R treefile outgroupname outputfile>

Rscript root_the_tree.R ../data/ry_GTR_ry_protein-nucleo-uniq-d__Bacteria_p__Firmicutes_A__p__Fusobacteriota__COG0012.fasta.treefile Fusobacteriota ../data/output/rooted-ry_GTR_ry_protein-nucleo-uniq-d__Bacteria_p__Firmicutes_A__p__Fusobacteriota__COG0012.fasta.treefile
Rscript root_the_tree.R ../data/ry_GTR_ry_protein-nucleo-uniq-d__Bacteria_p__Firmicutes_A__p__Fusobacteriota__COG0552.fasta.treefile Fusobacteriota ../data/output/rooted-ry_GTR_ry_protein-nucleo-uniq-d__Bacteria_p__Firmicutes_A__p__Fusobacteriota__COG0552.fasta.treefile
Rscript root_the_tree.R ../data/ry_GTR_ry_protein-nucleo-uniq-d__Bacteria_p__Spirochaetota__p__UBP6__p__Deferribacterota__COG0012.fasta.treefile UBP6 ../data/output/rooted-ry_GTR_ry_protein-nucleo-uniq-d__Bacteria_p__Spirochaetota__p__UBP6__p__Deferribacterota__COG0012.fasta.treefile
Rscript root_the_tree.R ../data/ry_GTR_ry_protein-nucleo-uniq-d__Bacteria_p__Spirochaetota__p__UBP6__p__Deferribacterota__COG0552.fasta.treefile UBP6 ../data/output/rooted-ry_GTR_ry_protein-nucleo-uniq-d__Bacteria_p__Spirochaetota__p__UBP6__p__Deferribacterota__COG0552.fasta.treefile

##step2. Get termite-specific clusters with atleast 10 termite species-
#A-If there is one outgroup sequences-
Rscript termitespecificclades_10species.R ../data/output/rooted-ry_GTR_ry_protein-nucleo-uniq-d__Bacteria_p__Firmicutes_A__p__Fusobacteriota__COG0012.fasta.treefile p__Fusobacteriota ../data/output/clade-Firmicutes_A-COG0012.txt ../data/output/Firmicutes_A-COG0012_tips
Rscript termitespecificclades_10species.R ../data/output/rooted-ry_GTR_ry_protein-nucleo-uniq-d__Bacteria_p__Firmicutes_A__p__Fusobacteriota__COG0552.fasta.treefile p__Fusobacteriota ../data/output/clade-Firmicutes_A-COG0552.txt ../data/output/Firmicutes_A-COG0552_tips

#B-If there are two outgroup sequences-
Rscript termitespecificclades_10species_2sisters.R ../data/output/rooted-ry_GTR_ry_protein-nucleo-uniq-d__Bacteria_p__Spirochaetota__p__UBP6__p__Deferribacterota__COG0012.fasta.treefile p__UBP6 p__Deferribacterota ../data/output/clade-Spirochaetota-COG0012.txt ../data/output/Spirochaetota-COG0012_tips
Rscript termitespecificclades_10species_2sisters.R ../data/output/rooted-ry_GTR_ry_protein-nucleo-uniq-d__Bacteria_p__Spirochaetota__p__UBP6__p__Deferribacterota__COG0552.fasta.treefile p__UBP6 p__Deferribacterota ../data/output/clade-Spirochaetota-COG0552.txt ../data/output/Spirochaetota-COG0552_tips

##step3. Get the subtree corresponding to extracted tips
for i in Spirochaetota-COG0012_tips*csv;do filename=`echo ${i}| sed 's/_.csv//g'` ; Rscript ../../scripts/keep_tip.R rooted-ry_GTR_ry_protein-nucleo-uniq-d__Bacteria_p__Spirochaetota__p__UBP6__p__Deferribacterota__COG0012.fasta.treefile ${i} ${filename}.nwk;done

##step4. cophylo plot
Rscript cophyloplot.R ../data/output/Spirochaetota-COG0012_tips_1.nwk ../data/output/Spirochaetota-COG0012_tips_1.pdf

##step5. cophylo stats 
#A. PACo-
Rscript paco_all_contig.pvalue.R ../data/output/Spirochaetota-COG0012_tips_1.nwk rna12-16S_202samples_newnames_feb2020.nwk ../data/output/overall-pvalue-Spirochaetota-COG0012_tips_1.nwk.txt

#B. Treedist-
#convert the symbiont tipnames to hostname_1, hostname_2 etc.-(or see "before_treedist.R")
Rscript before_treedist.R ../data/output/Spirochaetota-COG0012_tips_1.nwk ../data/output/symbionttree-Spirochaetota-COG0012_tips_1.nwk ../data/output/symbiontheader-Spirochaetota-COG0012_tips_1.txt

#get column 2 out from the header file-(or see "before_treedist.R")
awk -F"," '{print $2}' symbiontheader-Spirochaetota-COG0012_tips_1.txt | sed 's/"//g' > 2-symbiontheader-Spirochaetota-COG0012_tips_1.txt

#generate a new hosttree with same no.of tips as symbionts of zero branch lengths-(or see "before_treedist.R")
Rscript newhosttree.R ../data/output/2-symbiontheader-Spirochaetota-COG0012_tips_1.txt rna12-16S_202samples_newnames_feb2020.nwk ../data/output/hosttree-Spirochaetota-COG0012_tips_1.nwk

#run treedist-
Rscript treedist.R hosttree-Spirochaetota-COG0012_tips_1.nwk symbionttree-Spirochaetota-COG0012_tips_1.nwk generalized_rf_pvalue-Spirochaetota-COG0012_tips_1.nwk.txt nye-pvalue-Spirochaetota-COG0012_tips_1.nwk.txt
