
setwd("C:/Users/81704/Documents")

##prepare the data-
library(ape)
library(dplyr)
library(stringr)
library(plyr)

filename<-list.files(pattern= "rooted-ry_GTR_ry_protein_nucleo-uniq-d__Bacteria_p__Spirochaetota") #use ".+" to use multiple patterns together. for script1

newcolumns_function<-function(filename){
  repname_vector<-as.data.frame(filename) #get clade names from filename
  repname_vector$filename<-gsub("^.*COG","COG",repname_vector$filename)
  #repname_vector$filename<-gsub("_gtdbtermiteclades_april2021_","_",repname_vector$filename)
  repname_vector$filename<-gsub("_.csv","",repname_vector$filename)
  
  dat<-read.csv(file=filename)
  colnames(dat)<-c("X","tips")
  dat$clade<-rep(repname_vector$filename,times=nrow(dat))
  return(dat)
}


newcolumns_list<-lapply(filename,newcolumns_function)
allcogs_Bacteroidota<-do.call("rbind", newcolumns_list)%>%as.data.frame()%>%unique()

allcogs_Bacteroidota<-allcogs_Bacteroidota%>%dplyr::mutate(info=ifelse(str_detect(tips, "(^[0-9]+)"),"termite",ifelse(str_detect(tips, "^GB_GCA|^RS_GCF"),"gtdb","others" )))
gtdbcounts<-allcogs_Bacteroidota%>%dplyr::filter(info=="gtdb")%>%dplyr::select(tips)%>%dplyr::group_by(tips)%>%count() #which GTDB ids are repeated across COGs
colnames(gtdbcounts)<-c("tips","n")
gtdbcounts<-gtdbcounts%>%filter(n>1)


#get a final list/df with all combinations-

complist<-c() #create an empty list
for (i in 1:nrow(gtdbcounts)){
  tipstoextract<-as.data.frame(gtdbcounts[i,1])
  colnames(tipstoextract)<-c("column")
  tipstoextract<-gsub("__d.*$","",tipstoextract$column)
  comp<-allcogs_Bacteroidota%>%dplyr::filter(str_detect(tips, tipstoextract))%>%dplyr::select(clade)
  comp2<-t(comp)%>%as.data.frame()
  complist[[i]]<-comp2

}

complistdf<-plyr::rbind.fill(complist)%>%as.data.frame()%>%unique()
write.csv(complistdf,file="Spirochaetota_commonclades.csv")
