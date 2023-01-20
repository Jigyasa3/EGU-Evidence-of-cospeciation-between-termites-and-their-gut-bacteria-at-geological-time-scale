##Revisions2-Paper figures
setwd("C:/Users/81704/Dropbox (OIST)/OIST/Tom'sUnitProject/TomUnit_project/manuscripts/paper2_markergenes/manuscript_submittedProceedingB/Revisions/figures")

library(dplyr)
library(ape)
library(stringr)
library(tidyr)

#1. Figure1

##pvalue data-
data<-read.csv("allcophylo_beasttree_pvalue_analysis.csv")
data<-data%>%filter(str_detect(filename,"COG0172"))

data$tsc<-paste0("TSC",row.names(data)) #add TSC names

#add complete taxonomy-
data$taxonomy<-gsub("^.*_p_","",data$taxonomy)
data$taxonomy<-gsub("_c_"," ",data$taxonomy)
data$taxonomy<-gsub("_o_"," ",data$taxonomy)
data$taxonomy<-gsub("_f_"," ",data$taxonomy)
data$taxonomy<-gsub("_g_"," ",data$taxonomy)

data$final_name<-paste0(data$tsc," ", data$taxonomy)
data$taxonomy<-NULL
data$tsc<-NULL
data$filename<-NULL

data_long<-gather(data, metric, coevolution_sig, -c(final_name), factor_key=TRUE)
data_long<-data_long%>%mutate(coevolution_sig2=ifelse(coevolution_sig<0.05 & coevolution_sig>0.01,"*",ifelse(coevolution_sig<0.01 & coevolution_sig>0.001,"**",ifelse(coevolution_sig<0.001,"***",ifelse(metric=="Transfer.rate","Transfer.rate","notsig")))))

#orginal values data-
data2<-read.csv("allcophylo_beasttree_orginalvalues_analysis.csv")
data2<-data2%>%filter(str_detect(filename,"COG0172"))

data2$tsc<-paste0("TSC",row.names(data2)) #add TSC names

#add complete taxonomy-
data2$taxonomy<-gsub("^.*_p_","",data2$taxonomy)
data2$taxonomy<-gsub("_c_"," ",data2$taxonomy)
data2$taxonomy<-gsub("_o_"," ",data2$taxonomy)
data2$taxonomy<-gsub("_f_"," ",data2$taxonomy)
data2$taxonomy<-gsub("_g_"," ",data2$taxonomy)

data2$final_name<-paste0(data2$tsc," ", data2$taxonomy)
data2$taxonomy<-NULL
data2$tsc<-NULL
data2$filename<-NULL
colnames(data2)<-c("Nye","GRF","PACO","final_name")

data2$Nye<-as.numeric(data2$Nye)
data2$GRF<-as.numeric(data2$GRF)
data2_long<-gather(data2, metric, value, -c(final_name), factor_key=TRUE)

#combine two dfs-
alldata<-merge(data2_long,data_long,by=c("final_name","metric"),all.y=TRUE)
alldata<-alldata%>%mutate(value=ifelse(metric=="Transfer.rate",coevolution_sig,value))

#plot-https://stackoverflow.com/questions/67746044/r-heatmap-with-circles
library(ggplot2)
library(forcats)

#color_values<-c("0"="grey","0.25"="#CDAA7F","0.5"="#03C04a","1"="dark green","NA"="black")
color_values2<-c("notsig"="grey","*"="#CDAA7F","**"="#03C04a","***"="dark green","Transfer.rate"="black")
#joined_long3.2<-joined_long3.2%>%mutate(coevolution_sig2_colors=ifelse(coevolution_sig=="not-significant","grey",ifelse(coevolution_sig=="*","#CDAA7F",ifelse(coevolution_sig=="**","#03C04a",ifelse(coevolution_sig=="***","dark green","black")))))

#for COG0552-
alldata$final_name <- factor(alldata$final_name, levels=c("TSC1 Thermoplasmatota Thermoplasmata Methanomassiliicoccales Methanomethylophilaceae","TSC2 Thermoproteota Bathyarchaeia B26-1 UBA233","TSC3 Actinobacteriota Coriobacteriia Coriobacteriales","TSC4 Actinobacteriota Actinomycetia","TSC5 Bacteroidota Bacteroidia","TSC6 Bacteroidota Bacteroidia Bacteroidales","TSC7 Bacteroidota Bacteroidia","TSC8 Campylobacterota Campylobacteria Campylobacterales Sulfurospirillaceae","TSC9 Desulfobacterota Desulfovibrionia Desulfovibrionales Desulfovibrionaceae","TSC10 Desulfobacterota Desulfarculia Adiutricales Adiutricaceae Adiutrix","TSC11 Fibrobacterota Fibrobacteria Fibrobacterales Fibrobacteraceae Fibromonas","TSC12 Fibrobacterota Chitinivibrionia Chitinivibrionales","TSC13 Fibrobacterota Chitinivibrionia Chitinivibrionales","TSC14 Firmicutes_A Clostridia Oscillospirales Ruminococcaceae","TSC15 Firmicutes_A Clostridia Peptostreptococcales Anaerovoracaceae","TSC16 Firmicutes_A Clostridia Oscillospirales Oscillospiraceae","TSC17 Firmicutes_A Clostridia Oscillospirales","TSC18 Firmicutes_A Clostridia Oscillospirales Ruminococcaceae","TSC19 Planctomycetota Planctomycetes Pirellulales","TSC20 Proteobacteria Gammaproteobacteria Burkholderiales","TSC21 Proteobacteria Gammaproteobacteria Burkholderiales Rhodocyclaceae","TSC22 Proteobacteria Gammaproteobacteria Burkholderiales Rhodocyclaceae","TSC23 Proteobacteria Gammaproteobacteria Burkholderiales Burkholderiaceae","TSC24 Spirochaetota Spirochaetia Treponematales Treponemataceae_B Treponema_E","TSC25 Spirochaetota Spirochaetia Treponematales Treponemataceae_B Treponema_E","TSC26 Spirochaetota Spirochaetia Treponematales Treponemataceae_B","TSC27 Spirochaetota Spirochaetia Treponematales Treponemataceae_B"))

#for COG0172-
#alldata$final_name<-factor(alldata$final_name,levels=c("TSC1 Thermoplasmatota Thermoplasmata Methanomassiliicoccales Methanomethylophilaceae","TSC2 Thermoproteota Bathyarchaeia B26-1 UBA233","TSC3 Acidobacteriota Holophagae Holophagales Holophagaceae","TSC4 Actinobacteriota Coriobacteriia","TSC5 Actinobacteriota Actinomycetia Propionibacteriales Propionibacteriaceae","TSC6 Bacteroidota Bacteroidia","TSC7 Bacteroidota Bacteroidia","TSC8 Bacteroidota Bacteroidia","TSC9 Bacteroidota Bacteroidia","TSC10 Bacteroidota Bacteroidia","TSC11 Bacteroidota Bacteroidia Bacteroidales","TSC12 Bacteroidota Bacteroidia Bacteroidales","TSC13 Bacteroidota Bacteroidia Bacteroidales","TSC14 Bacteroidota Bacteroidia Bacteroidales","TSC15 Desulfobacterota Desulfarculia Adiutricales Adiutricaceae Adiutrix","TSC16 Elusimicrobiota Endomicrobia Endomicrobiales Endomicrobiaceae","TSC17 Fibrobacterota Fibrobacteria Fibrobacterales Fibrobacteraceae Fibromonas","TSC18 Fibrobacterota Chitinivibrionia Chitinivibrionales Chitinispirillaceae","TSC19 Firmicutes Bacilli Lactobacillales","TSC20 Planctomycetota Planctomycetes Pirellulales","TSC21 Proteobacteria Gammaproteobacteria Burkholderiales","TSC22 Proteobacteria Alphaproteobacteria Rickettsiales Anaplasmataceae Wolbachia","TSC23 Proteobacteria Gammaproteobacteria Burkholderiales","TSC24 Proteobacteria Gammaproteobacteria Burkholderiales Burkholderiaceae","TSC25 Spirochaetota Spirochaetia Treponematales","TSC26 Spirochaetota Spirochaetia Treponematales","TSC27 Spirochaetota Spirochaetia Treponematales","TSC28 Spirochaetota Spirochaetia Treponematales","TSC29 Spirochaetota Spirochaetia Treponematales","TSC30 Spirochaetota Spirochaetia Treponematales"))

p1<-alldata%>%ggplot(aes(x=metric, y=reorder(as.factor(final_name), desc(as.factor(final_name))), fill = coevolution_sig2))+#, size = coevolution_sig)) +
  geom_point(shape = 21, stroke = 0,aes(size=value)) +
  #scale_x_discrete(position = "top", labels=c("Generalized Robinson-Foulds","Nye et al.","PACo","Rate of transfer")) +
  #scale_y_discrete(labels=(levelnames2))+
  #scale_fill_gradientn(colors=c("grey","#CDAA7F","#03C04a","dark green"),breaks=c(0,0.25,0.5,1),labels = c("Not Significant", "* p<0.05","** p<0.01","*** p<0.001"))+
  scale_fill_manual(values=color_values2,labels = c("Not Significant", "* p<0.05","** p<0.01","*** p<0.001","Transfer values"))+
  theme_minimal()+scale_size(range = c(1.5,8)) +
  theme(legend.position = "bottom", 
        panel.grid.major = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),axis.text.y.left = element_text(hjust = 0, size=8)) +
  #guides(size = guide_legend(override.aes = list(fill = NA, color = "white", stroke = .25), 
  #                           label.position = "bottom",
  #                           title.position = "top", 
  #                           order = 1),
  #       fill = guide_colorbar(ticks.colour = NA, title.position = "top", order = 2)) +
  labs(size = "Cophylogenetic metrics", fill = "Significance:", x = NULL, y = NULL)

#p2<-ggplot(joined_long3.2,aes(x=totaltips,y=as.factor(final_name)))+   # this makes the plot horizontal
#  geom_bar(stat = "identity",position = "dodge")+theme_classic()+scale_x_continuous(breaks=c(0,50,100,150,200,250,1200))+scale_y_discrete(labels=rev(levelnames2))
#p2.2<-p2+theme(axis.text.y=element_blank())  

#library(ggpubr)
#library(grid)
#figure <- ggarrange(p1, NULL, p2.2,
#                    nrow = 1, align="hv",widths = c(1, 0.05, 1),
#                    common.legend = TRUE,legend="bottom")

#library(svglite)
#ggsave(file="COG0552trees_figure1.svg", plot=p1, width=8, height=8)


pdf(file="COG0172trees_figure1.pdf",height=7,width=11)
p1
dev.off()

#-------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------
#2.Figure2,3 cophylo plots

#-------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
#3. Figure3A

#BEASTTREE-
allmarkergene_dtl<-read.csv("allTSCs_beastmttree_generaxoutput.txt",header=TRUE,sep="\t")
#allmarkergene_dtl$V1<-gsub(".log","",allmarkergene_dtl$V1)
#allmarkergene_dtl$V1<-gsub("2-folder-","",allmarkergene_dtl$V1)

#colnames(allmarkergene_dtl)<-c("filename","D","L","T")

#allmarkergene_dtl$D<-gsub("D=","",allmarkergene_dtl$D)
#allmarkergene_dtl$L<-gsub("L=","",allmarkergene_dtl$L)
#allmarkergene_dtl$T<-gsub("T=","",allmarkergene_dtl$T)
allmarkergene_dtl$D<-as.numeric(allmarkergene_dtl$D)
allmarkergene_dtl$L<-as.numeric(allmarkergene_dtl$L)
allmarkergene_dtl$T<-as.numeric(allmarkergene_dtl$T)

#thirtytrees<-c("Spirochaetota_COG0552_tips_1","Fibrobacterota_COG0552_tips_1","Firmicutes_A_COG0552_tips_5","Firmicutes_A_COG0552_tips_1","Bacteroidota_COG0552_tips_1","Bacteroidota_COG0552_tips_3","Proteobacteria_COG0552_tips_1","Fibrobacterota_COG0552_tips_2","Thermoplasmatota_COG0552_tips_1","Proteobacteria_COG0552_tips_4","Methanobacteriota_COG0552_tips_1","Thermoproteota_COG0552_tips_1","Acidobacteriota_COG0552_tips_1","Actinobacteriota_COG0552_tips_1","Actinobacteriota_COG0552_tips_2","Bacteroidota_COG0552_tips_2","Bacteroidota_COG0552_tips_4","Campylobacterota_COG0552_tips_1","Desulfobacterota_COG0552_tips_1","Desulfobacterota_COG0552_tips_2","Fibrobacterota_COG0552_tips_3","Firmicutes_A_COG0552_tips_2","Firmicutes_A_COG0552_tips_3","Firmicutes_A_COG0552_tips_4","Planctomycetota_COG0552_tips_1","Proteobacteria_COG0552_tips_2","Proteobacteria_COG0552_tips_3","Spirochaetota_COG0552_tips_2","Spirochaetota_COG0552_tips_3","Spirochaetota_COG0552_tips_4")
thirtytrees<-allmarkergene_dtl%>%filter(str_detect(filename,"COG0172"))

allmarkergene_30<-allmarkergene_dtl%>%filter(filename %in% thirtytrees$filename)
allmarkergene_30$geneinfo<-rep("bacterial",times=nrow(allmarkergene_30))


mtdna_dtl<-read.csv("allmtdna_beastmttree_generaxoutput.txt",header=FALSE)
mtdna_dtl$V1<-gsub("2-folder-rooted-mttree-aligned-named-4-3-2-","",mtdna_dtl$V1)
mtdna_dtl$V1<-gsub(".fas.log","",mtdna_dtl$V1)
colnames(mtdna_dtl)<-c("filename","D","L","T")
mtdna_dtl$D<-gsub("D=","",mtdna_dtl$D)
mtdna_dtl$L<-gsub("L=","",mtdna_dtl$L)
mtdna_dtl$T<-gsub("T=","",mtdna_dtl$T)
mtdna_dtl$D<-as.numeric(mtdna_dtl$D)
mtdna_dtl$L<-as.numeric(mtdna_dtl$L)
mtdna_dtl$T<-as.numeric(mtdna_dtl$T)
mtdna_dtl$geneinfo<-rep("mtdna",times=nrow(mtdna_dtl))
#mtdna_dtl$cophylo<-rep("nonsignificant",times=nrow(mtdna_dtl))

#combine both dfs together-
alldtl<-rbind(allmarkergene_30,mtdna_dtl)%>%as.data.frame()
alldtl2<-alldtl%>%dplyr::filter(!str_detect(filename,"^tRNA")) #removing trna from the file

#add significant or not info-
data<-read.csv("allcophylo_beasttree_pvalue_analysis.csv")
data<-data%>%filter(str_detect(filename,"COG0172"))

data<-data%>%mutate(allsig=ifelse(GRF<0.05 & Nye <0.05 & PACO <0.05,"significant","notsignificant"))
data<-data%>%select(filename,allsig)

alldtl2<-merge(alldtl2,data,by="filename",all.x=TRUE)

write.csv(alldtl2,file="DTL_COG0172microbial_mtdna_beasttree.csv")

##plot1-Transfer rate
library(ggplot2)
p<-ggplot(alldtl2, aes(x=geneinfo, y=T,col=allsig)) + 
  geom_jitter(position=position_jitter(0.2),size=2)#geom_dotplot(binaxis='y', stackdir='center',stackratio=1.5, dotsize=1.2)
p1<-p+ylim(0,1)+theme_bw()+ggtitle("Transfer Rate against Beast tree") +
  xlab("Genes") + ylab("Rate of Transfer")+scale_x_discrete(labels= c("microbial markergenes","termite mtDNA genes"))

#--------------------------------------------------------------------------------------
#IQTREE-
allmarkergene_dtl<-read.csv("allTSCs_iqtree_mttree_generaxoutput.txt",header=TRUE,sep="\t")
#allmarkergene_dtl$V1<-gsub(".log","",allmarkergene_dtl$V1)
#allmarkergene_dtl$V1<-gsub("2-folder-","",allmarkergene_dtl$V1)

#colnames(allmarkergene_dtl)<-c("filename","D","L","T")

#allmarkergene_dtl$D<-gsub("D=","",allmarkergene_dtl$D)
#allmarkergene_dtl$L<-gsub("L=","",allmarkergene_dtl$L)
#allmarkergene_dtl$T<-gsub("T=","",allmarkergene_dtl$T)
allmarkergene_dtl$D<-as.numeric(allmarkergene_dtl$D)
allmarkergene_dtl$L<-as.numeric(allmarkergene_dtl$L)
allmarkergene_dtl$T<-as.numeric(allmarkergene_dtl$T)

#thirtytrees<-c("Spirochaetota_COG0552_tips_1","Fibrobacterota_COG0552_tips_1","Firmicutes_A_COG0552_tips_5","Firmicutes_A_COG0552_tips_1","Bacteroidota_COG0552_tips_1","Bacteroidota_COG0552_tips_3","Proteobacteria_COG0552_tips_1","Fibrobacterota_COG0552_tips_2","Thermoplasmatota_COG0552_tips_1","Proteobacteria_COG0552_tips_4","Methanobacteriota_COG0552_tips_1","Thermoproteota_COG0552_tips_1","Acidobacteriota_COG0552_tips_1","Actinobacteriota_COG0552_tips_1","Actinobacteriota_COG0552_tips_2","Bacteroidota_COG0552_tips_2","Bacteroidota_COG0552_tips_4","Campylobacterota_COG0552_tips_1","Desulfobacterota_COG0552_tips_1","Desulfobacterota_COG0552_tips_2","Fibrobacterota_COG0552_tips_3","Firmicutes_A_COG0552_tips_2","Firmicutes_A_COG0552_tips_3","Firmicutes_A_COG0552_tips_4","Planctomycetota_COG0552_tips_1","Proteobacteria_COG0552_tips_2","Proteobacteria_COG0552_tips_3","Spirochaetota_COG0552_tips_2","Spirochaetota_COG0552_tips_3","Spirochaetota_COG0552_tips_4")
thirtytrees<-allmarkergene_dtl%>%filter(str_detect(filename,"COG0172"))

allmarkergene_30<-allmarkergene_dtl%>%filter(filename %in% thirtytrees$filename)
allmarkergene_30$geneinfo<-rep("bacterial",times=nrow(allmarkergene_30))


mtdna_dtl<-read.csv("allmtdna_iqtree_mttree_generaxoutput.txt",header=FALSE)
mtdna_dtl$V1<-gsub("2-folder-rooted-mttree-aligned-named-4-3-2-","",mtdna_dtl$V1)
mtdna_dtl$V1<-gsub(".fas.log","",mtdna_dtl$V1)
colnames(mtdna_dtl)<-c("filename","D","L","T")
mtdna_dtl$D<-gsub("D=","",mtdna_dtl$D)
mtdna_dtl$L<-gsub("L=","",mtdna_dtl$L)
mtdna_dtl$T<-gsub("T=","",mtdna_dtl$T)
mtdna_dtl$D<-as.numeric(mtdna_dtl$D)
mtdna_dtl$L<-as.numeric(mtdna_dtl$L)
mtdna_dtl$T<-as.numeric(mtdna_dtl$T)
mtdna_dtl$geneinfo<-rep("mtdna",times=nrow(mtdna_dtl))
#mtdna_dtl$cophylo<-rep("nonsignificant",times=nrow(mtdna_dtl))

#combine both dfs together-
alldtl_iqtree<-rbind(allmarkergene_30,mtdna_dtl)%>%as.data.frame()
alldtl_iqtree2<-alldtl_iqtree%>%dplyr::filter(!str_detect(filename,"^tRNA")) #removing trna from the file


#add significant or not info-
data<-read.csv("allcophylo_uceiqtree_pvalue_analysis.csv")
data<-data%>%filter(str_detect(filename,"COG0172"))

data<-data%>%mutate(allsig=ifelse(GRF<0.05 & Nye <0.05 & PACO <0.05,"significant","notsignificant"))
data<-data%>%select(filename,allsig)

alldtl_iqtree2<-merge(alldtl_iqtree2,data,by="filename",all.x=TRUE)

write.csv(alldtl_iqtree2,file="DTL_COG0172microbial_mtdna_iqtree.csv")

##plot1-Transfer rate
library(ggplot2)
q<-ggplot(alldtl_iqtree2, aes(x=geneinfo, y=T,col=allsig)) + 
  geom_jitter(position=position_jitter(0.2),size=2)#geom_dotplot(binaxis='y', stackdir='center',stackratio=1.5, dotsize=1.2)
q1<-q+ylim(0,1)+theme_bw()+ggtitle("Transfer Rate against UCE tree") +
  xlab("Genes") + ylab("Rate of Transfer")+scale_x_discrete(labels= c("microbial markergenes","termite mtDNA genes"))


##All plots together-
library(ggpubr)
figure <- ggarrange(p1, q1,
                    #labels = c("A", "B"),
                    ncol = 2, nrow = 1,
                    common.legend = TRUE, legend = "right")

pdf(file="TransferRate_COG0172_withcolors_microbial_mtdna_genes.pdf",height=6,width=10)
figure
dev.off()




