#!/usr/bin/env Rscript
#USAGE- Rscript
args <- commandArgs(TRUE)
tree <- args[1] #rootedtree
sister <- args[2] #p__Fibrobacterota
output2 <- args[3] #no.of termitespecific clades >=10 species
suffix1 <- args[4] #suffix for termitespecific clades tips files

##Get no.of sequences per markergene
#for one markergene-
library(ape)
library(dplyr)
library(treeio)
library(stringr)

#get user defined sistergroup info-
cogtree<-read.newick(file=tree)
d<-data.frame(label=cogtree$tip.label)
d$phyla<-gsub("^.*_d__","d__",d$label)

sistergroup<-as.vector(as.character(sister)) #user defined
print(sistergroup)
#d<-d%>%mutate(groupinfo=ifelse(str_detect(label,"^RS_GCF|^GB_GCA"),"notphylaofinterest",ifelse(str_detect(phyla,sistergroup),"notphylaofinterest","phylaofinterest")))
#-------------------------------------------------------------------------------------------
#2.Get no.of termite specific clades >=10 termite species per markergene
##get termite species info for each tip-
d<-data.frame(label=cogtree$tip.label)
#get "termitegut" info from labels-
d<-d%>%mutate(info2=ifelse(str_detect(label,sistergroup),as.character(label),ifelse(str_detect(label,"termitegut"),"termitegut",as.character(label))))
#d$info<-gsub(".*termitegut_1","termitegut",d$label) #get termite or not info
d<-d%>%mutate(info2=ifelse(str_detect(label,"^RS_GCF_900167335_1|^RS_GCF_001639295_1|^RS_GCF_001639285_1"),"termitegut",as.character(d$info2))) #add the termitegut GTBD ids

d<-d%>%dplyr::mutate(outgroupinfo=ifelse(str_detect(label,"^RS_GCF|^GB_GCA"),"gtdb","markers")) #get gtdb outgroup or markergene info
d$runnumber<-gsub("--.*$","",d$label)
d$runnumber<-gsub("_d.*$","",d$runnumber)

d<-d%>%dplyr::mutate(info=ifelse(outgroupinfo=="markers","termitegut",ifelse(outgroupinfo=="gtdb" & info2 !="termitegut","others",as.character(info2))))


sampleinfo<-read.csv("/flash/BourguignonU/Jigs/markergene_paper_scripts/scripts/sample_names_ids_withotherfactors_tomedit.202samples.csv")
sampleinfo<-sampleinfo%>%dplyr::select(runnumber,termitespecies)
d2<-merge(d,sampleinfo,by="runnumber",all.x=TRUE)

##get node info for termite specific clades-
cogtree2<-cogtree
x <- tidytree::as_tibble(cogtree2) #%>% mutate(isTip = ! .data$node %in% .data$parent)

gtdbtermite<-function(parentnode){
  library(caper)
  library(phytools)
  parentnode_tips<-clade.members(parentnode, cogtree2, tip.labels=TRUE)%>%as.data.frame()
  colnames(parentnode_tips)<-c("tips")
  parentnode_tips<-merge(parentnode_tips,d2,by.x="tips",by.y="label")
  parentnode_tips_info<-parentnode_tips%>%dplyr::select(info)%>%dplyr::group_by(info)%>%dplyr::summarize(n())
  parentnode_tips_info<-as.data.frame(parentnode_tips_info)
  colnames(parentnode_tips_info)<-c("info","n")
  
  parentnode_tips_outgroupinfo<-parentnode_tips%>%dplyr::select(outgroupinfo)%>%dplyr::group_by(outgroupinfo)%>%dplyr::summarize(n())
  parentnode_tips_outgroupinfo<-as.data.frame(parentnode_tips_outgroupinfo)
  colnames(parentnode_tips_outgroupinfo)<-c("outgroupinfo","n")
  
  parentnode_tips_species<-parentnode_tips%>%dplyr::select(termitespecies)%>%dplyr::group_by(termitespecies)%>%dplyr::summarize(n())
  parentnode_tips_species<-na.omit(parentnode_tips_species)
  parentnode_tips_species<-as.data.frame(parentnode_tips_species)
  colnames(parentnode_tips_species)<-c("species","n")
  
  if(nrow(na.omit(parentnode_tips_species))>=10 && parentnode_tips_info$info!="others"){
    print(paste0("This parent node= ", parentnode, "has more than equal to 10 termite species"))
    return(parentnode)
    
  }else{
    print(paste0("This parent node= ", parentnode, "doesn't contain 10 or more termite species. Ignore it."))
    parentnode <- c() #else statement is TRUE, then return an empty vector
    return(parentnode)  
    
  }
}

h_lista<-as.list(x$parent)
gtdbtermite_list<-lapply(h_lista,gtdbtermite)

if (class(gtdbtermite_list)=="list") {
  gtdbtermite_df2<-do.call("rbind", gtdbtermite_list)%>%as.data.frame()%>%unique() #these are nodeids to remove!
} else {
  gtdbtermite_df2<-(gtdbtermite_list)%>%as.data.frame()%>%unique()
  colnames(gtdbtermite_df2)<-c("V1")
}

#----------------------------------------------------------------------------------------------
#compare every pair of nodes for checkpoint1-
gtdbtermite_df2_matrix_combinations <- expand.grid(gtdbtermite_df2$V1,gtdbtermite_df2$V1)
gtdbtermite_df2_matrix_combinations$Var1<-as.numeric(gtdbtermite_df2_matrix_combinations$Var1)
gtdbtermite_df2_matrix_combinations$Var2<-as.numeric(gtdbtermite_df2_matrix_combinations$Var2)

gtdbtermite_df2_matrix_combinations<-gtdbtermite_df2_matrix_combinations%>%dplyr::mutate(repeats=ifelse(Var1==Var2,"repeat","unique"))
gtdbtermite_df2_matrix_combinations<-gtdbtermite_df2_matrix_combinations%>%dplyr::filter(repeats!="repeat")  
gtdbtermite_df2_matrix_combinations$repeats<-NULL

gtdbtermite_checkpoint1<-function(node1,node2){
  node1_tips<-clade.members(node1, cogtree2, tip.labels=TRUE)%>%as.data.frame()
  colnames(node1_tips)<-c("V1")
  
  node2_tips<-clade.members(node2, cogtree2, tip.labels=TRUE)%>%as.data.frame()
  colnames(node2_tips)<-c("V1")
  
  node12_tips<-merge(node1_tips,node2_tips,by="V1") #get common tips between two nodes
  
  if(nrow(node12_tips)==nrow(node1_tips)){
    print(paste0("Remove the node ", node1))
    return(node1) #get the nodeid to remove from the list
  }else if (nrow(node12_tips)==nrow(node2_tips)){
    print(paste0("Remove the node ", node2))
    return(node2) #get the nodeid to remove from the list
  } else{
    #print(paste0("Keep both nodes ", node1 ," and", node2))
  }
  
  
}

h_list1<-as.list(gtdbtermite_df2_matrix_combinations$Var1)
h_list2<-as.list(gtdbtermite_df2_matrix_combinations$Var2)
gtdbtermite_checkpoint1_list<-mapply(gtdbtermite_checkpoint1,h_list1,h_list2)

if (class(gtdbtermite_checkpoint1_list)=="list") {
  gtdbtermite_checkpoint1_df2<-do.call("rbind", gtdbtermite_checkpoint1_list)%>%as.data.frame()%>%unique() #these are nodeids to remove!
} else {
  gtdbtermite_checkpoint1_df2<-(gtdbtermite_checkpoint1_list)%>%as.data.frame()%>%unique()
  colnames(gtdbtermite_checkpoint1_df2)<-c("V1")
}


#Remove the node labels in "gtdbtermite_checkpoint1_df2" from "gtdbtermite_df2" as they are represented by a larger clade-
gtdbtermite_df2_final<-gtdbtermite_df2%>%dplyr::filter(!V1 %in% gtdbtermite_checkpoint1_df2$V1)

#--------------------------
#checkpoint2- remove clades corresponding to the sistergroup-

sistergroupfind<-function(nodenumber,sistergroup){
  a<-clade.members(nodenumber, cogtree2, tip.labels=TRUE)%>%as.data.frame()
  colnames(a)<-c("tips")
  a$taxa<-gsub("_c__.*$","",a$tips)
  a$taxa<-gsub("^.*_d","d",a$taxa)
  a<-a%>%dplyr::mutate(phyla=ifelse(str_detect(taxa, sistergroup),"sistergroup","maingroup"))
  b<-a%>%dplyr::select(phyla)%>%dplyr::group_by(phyla)%>%dplyr::count()
  b<-as.data.frame(b)
  if (b$phyla=="sistergroup") {
    print(paste0("the node ", nodenumber, "needs to be removed. It consists of sister phyla."))
    return(nodenumber)
  } else{
    print(paste0("the node", nodenumber, "has phyla of interest"))
    nodenumber<-c() #empty
    return(nodenumber)
  }
} #get the nodes to remove

nodenumber_list<-list() #initialize a list
nodenumber_list<-as.list(gtdbtermite_df2_final$V1) #add values to empty list
sistergroupfind_list<-mapply(sistergroupfind,nodenumber_list,as.vector(as.character(sistergroup)))


if (class(sistergroupfind_list)=="list") {
  sistergroupfind_list_df2<-do.call("rbind", sistergroupfind_list)%>%as.data.frame()%>%unique() #these are nodeids to remove!
} else {
  sistergroupfind_list_df2<-(sistergroupfind_list)%>%as.data.frame()%>%unique()
  colnames(sistergroupfind_list_df2)<-c("V1")
}


gtdbtermite_df2_final<-gtdbtermite_df2_final%>%dplyr::filter(!V1 %in% sistergroupfind_list_df2$V1)

outputfile2<-nrow(gtdbtermite_df2_final) #no.of termitespecific clades in the markergene
write.csv(outputfile2,file=output2)

#------------------------------------------------------------------------------------------------------
##get the sequences-
finalnodes_df2<-gtdbtermite_df2_final #load the output file "gtdbtermite_df2_final" to get dataframe format
finalnodes_df3<-as.data.frame(t(finalnodes_df2))

clades_tips<-function(nodeinfo){
  a<-clade.members(nodeinfo, cogtree2, tip.labels=TRUE)%>%as.data.frame()
  return(a)
}

node_list<-list() #initialize a list
node_list<-as.list(finalnodes_df2$V1) #add values to empty list

clademembers_list<-lapply(node_list,clades_tips)


#write to csv-
rowvector<-c(1 : length(finalnodes_df3))
write2csv<-function(data,file) write.csv(data,file)
library(plyr)
plyr::m_ply(cbind(data=clademembers_list,file=paste(as.vector(suffix1),rowvector,".csv",sep="_")),write2csv)

