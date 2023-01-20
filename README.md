# EGU-Evidence-of-cospeciation-between-termites-and-their-gut-bacteria-at-geological-time-scale

This repository contains the scripts to run the analysis performed for the manuscript *Evidence-of-cospeciation-between-termites-and-their-gut-bacteria-at-geological-time-scale*

The data analysis occurs in the following steps- The folder names and filenames are same as in the repository. The scripts are in ```script/``` folder and data in ```data/``` folder.

### mtDNA beast2 tree-

## Run mitofinder on contigs extracted by BLAST search-
```
module load MitoFinder/1.4
module load python/2.7.18
module load ncbi-blast/2.10.0+
module load mitfi/0.1


DB="/bucket/BourguignonU/SimonH/databases/ANNOTATION_MT_GENOMES/MITOFINDER_DB"

mitofinder --seqid ${file2} --assembly ${file1} --refseq ${DB}/sequence.gb --organism 5 --processors 14  --blast-eval 0.00001 --blast-identity-nucl 50 --blast-identity-prot 40 --blast-size 30 --min-contig-size 300 --max-contig-size 30000

```

## Extract the outputs from mitofinder-
```
#Get all the gff files to the same directory-
for i in mtdna*;do cp ${i}/${i}_MitoFinder_mitfi_Final_Results/*.gff mitofinder_output/${i}.gff; done

## only print out lines that contain "gene"
for i in *gff;do awk '!/gene/' $i > ${i}_new1;done
##replace column3 with column9
for i in *_new1 ;do awk '{$3=$9;$9=$1}1'  $i > ${i}_new2;done
##change file to Tab-deliminated format
for i in *_new2 ;do awk '{ for(i=1;i<=NF;i++){if(i==NF){printf("%s\n",$NF);}else {printf("%s\t",$i)}}}' ${i} > ${i}_new3;done
##Only print out first 9 columns
for i in *_new3 ;do awk '{print$1,$2,$3,$4,$5,$6,$7,$8,$9}' ${i} > ${i}_new4;done
##change file to Tab-deliminated format
for i in *_new4;do awk '{ for(i=1;i<=NF;i++){if(i==NF){printf("%s\n",$NF);}else {printf("%s\t",$i)}}}' ${i} > ${i}_final.gff;done

##Extract Genes separately
##Make sure that fasta header and col1 of bed file match!
for i in *gff;do awk '{print $1; exit}' ${i}> ${i}-header.txt;done
cat *header.txt > allheadersgff.txt

for i in mtdna*fasta;do seqkit replace -p "(.+)$" -r "{kv}" -k ../2-original_mitofinder_headers.txt ${i} > mitofinder_output/2-${i};done
for i in 2-mtdna272-*.fasta;do filename=`echo ${i}|awk -F".fasta" '{print $1}'`; filename2=`echo ${filename}| awk -F"2-mtdna" '{print "mtdna"$2}'`; bedtools getfasta -fi ${i} -bed ${filename2}.gff_new1_new2_new3_new4_final.gff -name > ${filename2}_separate_genes.fasta; done

<on one file> 
bedtools getfasta -fi 2-mtdna301-59.fasta -bed mtdna301-59.gff_new1_new2_new3_new4_final.gff -name > mtdna301-59_se
parate_genes.fasta

##Generate a gene specific file for all mtdna genes
#Get the 38 mtdna gene ids out-
grep ">" mtdna230-24_separate_genes.fasta > individual_genes/all_mtdna_genenames.txt #230-24 has all 38 genes!
sed -i 's/mtdna.*$//g' all_mtdna_genenames.txt

#Generate separate files with all sampleids for each mtdna geneid-
<on one file> grep ">Name=tRNA-Asn::" *separate_genes.fasta | sed 's/^.*.fasta:>//g' > individual_genes/ids-tRNA-Asn.txt
while read line;do newname=`echo ${line}|awk -F">Name=" '{print $2}'`; newname2=`echo ${newname}| awk -F "::" '{print $1}'` ; grep "${line}" *separate_genes.fasta | sed 's/^.*.fasta:>//g' > individual_genes/ids-${newname2}.txt  ;done < individual_genes/all_mtdna_genenames.txt

#Extract the fasta sequences out for each mtdna geneid-
cat *_separate_genes.fasta > all_separate_genes.fasta
for i in ids-*.txt;do filename=`echo ${i}| awk -F"ids-" '{print $2}'`; filename2=`echo ${filename}| awk -F".txt" '{print $1}'` ; seqtk subseq ../all_separate_genes.fasta ${i} > ${filename2}.fasta ;done

##Rename the fasta headers with original filename-
#Remove unnecessary info from fasta headers-
for i in *fasta;do sed 's/>.*::mtdna/>/g' ${i} > newnames/2-${i};done
for i in 2-*.fasta;do sed 's/\.1.*$//g' ${i} > 3-${i};done

#Remove the 3 samples from all files- <301-30, 301-39, 272-20>
nano idstoremove.txt
for i in 3-2-*.fasta;do awk 'BEGIN{while((getline<"idstoremove.txt")>0)l[">"$1]=1}/^>/{f=!l[$1]}f' ${i} > 4-${i};done

#Get the old names and new names in tab separated file-
awk -F"," '{print $1"\t" $16}' sample_names_ids_withotherfactors_tomedit.202samples_nov2022.csv > 2-sample_names_ids_withotherfactors_tomedit.202samples_nov2022.csv

#rename the fasta files-
for i in 4-3-2-*fasta;do seqkit replace -p "(.+)$" -r "{kv}" -k ../../../../../2-sample_names_ids_withotherfactors_tomedit.202samples_nov2022.csv ${i} > named-${i};done
```

## Generate smatrix to run BEAST-
```
##run mafft
module load mafft/7.305

#mafft --maxiterate 1000 --localpair --amino ${IN_DIR}/${file1} > ${OUT_DIR}/aligned-${file1} #for short no. of sequences
mafft --auto ${IN_DIR}/${file1} > ${OUT_DIR}/aligned-${file1} #for large no. of sequences

##Make sure that aligned files have extension ".fas"
/home/j/jigyasa-arora/bin/rename/rename "s/.fasta/.fas/" aligned*.fasta

##Separate out the mtDNA genes into three folders- protein/ rrna/ trna/
##Generate a separate smatrix for all three folders-

perl /home/j/jigyasa-arora/local/FASconCAT/FASconCAT_v1.11.pl -s

cp protein/FcC_smatrix_prot.fas .
cp rrna/FcC_smatrix_rrna.fas .
cp trna/FcC_smatrix_trna.fas .

perl /home/j/jigyasa-arora/local/FASconCAT/FASconCAT_v1.11.pl -s #smatrix of all rrna, trna, proteins combined!
```

## Run BEAST2.4.8
```
1. Generate a NEXUS file from FcC_smatrix.fas alignment file in GENEIOUS PRIME
2. In the NEXUS file add the following partitions at the end-

begin assumptions;
        charset prot_1st=1-25156\2;
        charset prot_2nd=2-25156\2;
        charset rrna=25157-30349;
        charset trna=30350-34553;
end;

3. Script to run BEAST2-
#SBATCH --partition=largemem
#SBATCH --time=35-0
#SBATCH --mem=140G

module load beast2/v2.4.8
java -jar /apps/free72/beast2/v2.4.8/lib/beast.jar -seed 13579 -overwrite 18Aug2019_termites_Serriconstrained.xml

4. Script to run BEAST2 Treeannotator-
#SBATCH --partition=largemem
#SBATCH --time=19-0
#SBATCH --mem=700G

module load beast2/v2.4.8
##NOTE-check TRACER to get how much percent of generations need to be removed for burn-in. Sometime it varies from the one defined in .xml file
java -Xms550m -Xmx550g -Djava.library.path="$BEAST/lib" -cp "$BEAST/lib/launcher.jar" beast.app.treeannotator.TreeAnnotatorLauncher -burnin 10 rna12-16S.trees rna12-16S_202samples_concensus.trees


```
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

### Marker gene trees-

## step0. Check extractmarkergenes.md and GTDBver95_outgroups.md followed by generate_phyla_specific_trees.md in the scripts/ folder

## step1. Root the tree with outgroup bacterial sequences-

```module load R/3.6.1```

```
Rscript root_the_tree.R ../data/ry_GTR_ry_protein-nucleo-uniq-d__Bacteria_p__Firmicutes_A__p__Fusobacteriota__COG0012.fasta.treefile Fusobacteriota ../data/output/rooted-ry_GTR_ry_protein-nucleo-uniq-d__Bacteria_p__Firmicutes_A__p__Fusobacteriota__COG0012.fasta.treefile

Rscript root_the_tree.R ../data/ry_GTR_ry_protein-nucleo-uniq-d__Bacteria_p__Firmicutes_A__p__Fusobacteriota__COG0552.fasta.treefile Fusobacteriota ../data/output/rooted-ry_GTR_ry_protein-nucleo-uniq-d__Bacteria_p__Firmicutes_A__p__Fusobacteriota__COG0552.fasta.treefile

Rscript root_the_tree.R ../data/ry_GTR_ry_protein-nucleo-uniq-d__Bacteria_p__Spirochaetota__p__UBP6__p__Deferribacterota__COG0012.fasta.treefile UBP6 ../data/output/rooted-ry_GTR_ry_protein-nucleo-uniq-d__Bacteria_p__Spirochaetota__p__UBP6__p__Deferribacterota__COG0012.fasta.treefile

Rscript root_the_tree.R ../data/ry_GTR_ry_protein-nucleo-uniq-d__Bacteria_p__Spirochaetota__p__UBP6__p__Deferribacterota__COG0552.fasta.treefile UBP6 ../data/output/rooted-ry_GTR_ry_protein-nucleo-uniq-d__Bacteria_p__Spirochaetota__p__UBP6__p__Deferribacterota__COG0552.fasta.treefile
```
  
## step2. Get termite-specific clusters with atleast 10 termite species-
### A-If there is one outgroup sequences-
```Rscript termitespecificclades_10species.R ../data/output/rooted-ry_GTR_ry_protein-nucleo-uniq-d__Bacteria_p__Firmicutes_A__p__Fusobacteriota__COG0012.fasta.treefile p__Fusobacteriota ../data/output/clade-Firmicutes_A-COG0012.txt ../data/output/Firmicutes_A-COG0012_tips

Rscript termitespecificclades_10species.R ../data/output/rooted-ry_GTR_ry_protein-nucleo-uniq-d__Bacteria_p__Firmicutes_A__p__Fusobacteriota__COG0552.fasta.treefile p__Fusobacteriota ../data/output/clade-Firmicutes_A-COG0552.txt ../data/output/Firmicutes_A-COG0552_tips
```
  
### B-If there are two outgroup sequences-
```
Rscript termitespecificclades_10species_2sisters.R ../data/output/rooted-ry_GTR_ry_protein-nucleo-uniq-d__Bacteria_p__Spirochaetota__p__UBP6__p__Deferribacterota__COG0012.fasta.treefile p__UBP6 p__Deferribacterota ../data/output/clade-Spirochaetota-COG0012.txt ../data/output/Spirochaetota-COG0012_tips

Rscript termitespecificclades_10species_2sisters.R ../data/output/rooted-ry_GTR_ry_protein-nucleo-uniq-d__Bacteria_p__Spirochaetota__p__UBP6__p__Deferribacterota__COG0552.fasta.treefile p__UBP6 p__Deferribacterota ../data/output/clade-Spirochaetota-COG0552.txt ../data/output/Spirochaetota-COG0552_tips
```
  
## step3. Get the subtree corresponding to extracted tips
```
for i in Spirochaetota-COG0012_tips*csv;do filename=`echo ${i}| sed 's/_.csv//g'` ; Rscript ../../scripts/keep_tip.R rooted-ry_GTR_ry_protein-nucleo-uniq-d__Bacteria_p__Spirochaetota__p__UBP6__p__Deferribacterota__COG0012.fasta.treefile ${i} ${filename}.nwk;done
```

## step4. cophylo plot
```
Rscript cophyloplot.R ../data/output/Spirochaetota-COG0012_tips_1.nwk ../data/output/Spirochaetota-COG0012_tips_1.pdf
```

## step5. cophylo stats
### A. PACo-
```
Rscript paco_all_contig.pvalue.R ../data/output/Spirochaetota-COG0012_tips_1.nwk rna12-16S_202samples_newnames_feb2020.nwk ../data/output/overall-pvalue-Spirochaetota-COG0012_tips_1.nwk.txt
```
### B. Treedist-
#### convert the symbiont tipnames to hostname_1, hostname_2 etc.-(or see "before_treedist.R")
```
Rscript before_treedist.R ../data/output/Spirochaetota-COG0012_tips_1.nwk ../data/output/symbionttree-Spirochaetota-COG0012_tips_1.nwk ../data/output/symbiontheader-Spirochaetota-COG0012_tips_1.txt
```
#### get column 2 out from the header file-(or see "before_treedist.R")
```
awk -F"," '{print $2}' symbiontheader-Spirochaetota-COG0012_tips_1.txt | sed 's/"//g' > 2-symbiontheader-Spirochaetota-COG0012_tips_1.txt
```
#### generate a new hosttree with same no.of tips as symbionts of zero branch lengths-(or see "before_treedist.R")
```
Rscript newhosttree.R ../data/output/2-symbiontheader-Spirochaetota-COG0012_tips_1.txt rna12-16S_202samples_newnames_feb2020.nwk ../data/output/hosttree-Spirochaetota-COG0012_tips_1.nwk
```
#### run treedist-
```
Rscript treedist.R hosttree-Spirochaetota-COG0012_tips_1.nwk symbionttree-Spirochaetota-COG0012_tips_1.nwk generalized_rf_pvalue-Spirochaetota-COG0012_tips_1.nwk.txt nye-pvalue-Spirochaetota-COG0012_tips_1.nwk.txt
``` 

## step6. Get TSCs common across marker genes-
```
#1. Run script "getcommonclades_termitegtdbfunction.R" on RSTUDIO as the library HOME doesn't install on DEIGO
#2. In Excel-   Convert rows to columns
                Manually combine columns if any TSC is common between them
                Remove duplicates
                Save file
```
