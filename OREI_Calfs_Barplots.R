library(tidyverse)
library(ggfortify)
library(compositions)
library(stats)
library(DESeq2)
#install.packages("BiocManager")
library(BiocManager)
#BiocManager::install("DESeq2")
library(ggrepel)
#BiocManager::install("metagenomeSeq")
library(metagenomeSeq)
#BiocManager::install("ALDEx2")
library(ALDEx2)

#Load in data
amr<-read.csv("~/Desktop/NoyesLab/OREI_Calfs/Noyes_Project_019_AMR++_Results/SNPconfirmed_AMR_analytic_matrix.csv")
metadata<-read.csv("~/Desktop/NoyesLab/OREI_Calfs/sample_data_06_08_2023.csv")
rownames(metadata)<-metadata$SampleId
#dim(amr)
#separate amr classifications into own columns
amr_classification<-amr %>% separate_wider_delim(cols = gene_accession,delim='|',names = c("accession","type" ,"class", "mechanism","group"),too_many="merge")

#Manually fix bugged Rifampin
amr_classification$group[amr_classification$accession=="MEG_8674"]<-
  paste0(amr_classification$accession[amr_classification$accession=="MEG_8674"],"|RequiresSNPConfirmation")

#Exclude low sensitivity SNP confirmed AMR
amr_classification<-amr_classification[!grepl("RequiresSNPConfirmation",amr_classification$group),]

#View(amr_classification)
#Check on any 0 counts
#colnames(amr)[colSums(as.matrix(amr[,2:dim(amr)[2]]))==0]
#amr_classification[rowSums(amr_classification[,6:dim(amr_classification)[2]])==0,5])



#####
#####BEGIN COUNT PROCESSING FUNCTION
#####
class_process<- function(x) {
  
  ###Group by and sum according to classification of AMR
  #Agg ***sum*** accross classifications is correct
  mechsums<-aggregate(.~mechanism,data=x[,c(4,6:dim(x)[2])],FUN=sum)
  classsums<-aggregate(.~class,data=x[,c(3,6:dim(x)[2])],FUN=sum)
  groupsums<-aggregate(.~group,data=x[,c(5,6:dim(x)[2])],FUN=sum)
  rownames(mechsums)<-mechsums$mechanism
  rownames(classsums)<-classsums$class
  rownames(groupsums)<-groupsums$group
  
  #mechsums<-rowsum(amr_classification[,c(6:dim(amr_classification)[2])],group=amr_classification$mechanism);dim(mechsums) #Sum of counts by mech
  #classsums<-rowsum(amr_classification[,c(6:dim(amr_classification)[2])],group=amr_classification$class);dim(classsums) #Sum of counts by class
  #groupsums<-rowsum(amr_classification[,c(6:dim(amr_classification)[2])],group=amr_classification$group);dim(groupsums) #Sum of counts by group
  
  
  ###Data cleaning and exclusion
  #Keep only real USDA samples, no controls
  mechsums<-mechsums[,grepl("USDA",colnames(mechsums))] #Keep only real samples
  classsums<-classsums[,grepl("USDA",colnames(classsums))]     
  groupsums<-groupsums[,grepl("USDA",colnames(groupsums))]
  
  #Reformat names to match metadata
  colnames(mechsums)<-gsub("(USDA\\d+)_.*",'\\1',colnames(mechsums))
  colnames(classsums)<-gsub("(USDA\\d+)_.*",'\\1',colnames(classsums))
  colnames(groupsums)<-gsub("(USDA\\d+)_.*",'\\1',colnames(groupsums))
  
  ###Remove 0 counts
  mechsums<-mechsums[rowSums(mechsums)!=0,]
  #sum(rowSums(mechsums)==0) #Fixed it
  classsums<-classsums[rowSums(classsums)!=0,]
  #sum(rowSums(classsums)==0) #Fixed it
  groupsums<-groupsums[rowSums(groupsums)!=0,]
  #sum(rowSums(groupsums)==0) #Fixed it

  
  #transpose for analysis
  mech_t<-as.data.frame(t(mechsums) )
  class_t<-as.data.frame(t(classsums) )
  group_t<-as.data.frame(t(groupsums) )
  
  mat_list<-list(mech_t,class_t,group_t)
  names(mat_list)<-c("mech","class","group")
  
  return(mat_list)
  
}

amr_list<-class_process(amr_classification)



#View(amr_list$class)



####Maybe start function here, input is class_process object

###Merge transposed dataset with metadata
amr_list$class$SampleId<-rownames(amr_list$class)
class_t_meta<-amr_list$class %>% left_join(metadata[c("SampleId","CaseOrControl","FarmId")],by=c("SampleId") ) #Keep certain columns
#sum(is.na(class_t_farm$WeeksToInfection)) #Confirmed no unmatched

colnames(class_t_meta)
class_t_meta<-class_t_meta[,!( colnames(class_t_meta) %in% c("SampleId") )]  #Remove all unwanted metadata fields
#class_t_farm_clr<-class_t_farm_clr[,c(1:48,52,60)]  #Remove all unwanted metadata fields

class_agg<-aggregate(.~CaseOrControl + FarmId,data=class_t_meta,FUN = mean)

#Convert feature counts to proportion of reads per sample
class_agg[,!(colnames(class_agg) %in% c("CaseOrControl","FarmId"))]<-t( apply(   #Needs transpose to fit
  class_agg[,!(colnames(class_agg) %in% c("CaseOrControl","FarmId"))],1,
  function(x) x/sum(x,na.rm=TRUE)
) )

#Long format convert all features
class_agg_long<-gather(class_agg,amr,count,Acetate_resistance:Zinc_resistance) 



#####Filter down to only the Top 10
#Group dataframe into Farm and CaseOrControl, take Top 10 counts of each
top_10<-class_agg_long %>% group_by(FarmId,CaseOrControl) %>% 
  arrange(desc(count),.by_group = T) %>% top_n(10,count)
#View(top_10)

#Add an Other category, proportion of all other classes
others<-cbind( data.frame(amr=rep("Other",8)),
       aggregate(count~CaseOrControl + FarmId,data=top_10,function(x) 1-sum(x)) )
top_10<-rbind(top_10,others)



###Make classic barplot
ggplot(data=top_10,aes(x=interaction(FarmId,CaseOrControl),y=count,fill=amr)) + 
  geom_bar(stat="identity") +
  theme_minimal() +
  scale_fill_manual(name="AMR Class", values=colorRampPalette(brewer.pal(11,"Spectral"))(15) ) +
  theme(legend.title = element_text(size=10)) +
  theme(axis.text.x = element_text(angle = -45,hjust=0.15,vjust=-0.4)) +
  theme(legend.text = element_text(size=7)) +
  theme(legend.key.height= unit(3, 'mm'),
        legend.key.width= unit(5, 'mm')) +
  xlab("Case or Control by Farm") +
  ylab("AMR Class Proportion")













#Group by on farm and case/control
class_agg<-aggregate(.~CaseOrControl + FarmId,data=class_t_farm_clr,FUN = mean)
#Group by on farm or other variables  -- group by sum or group by mean??

#remove and store metadata
rownames(class_agg)<-paste(class_agg[,1],class_agg[,2])
meta<-class_agg[,1:2] #set it aside for after clr
class_agg<-class_agg[,-c(1,2)]
#class_clr<-clr(class_agg) #perform clr normalization

heatmap(as.matrix(t(class_agg)),cexRow = 0.5,Rowv=NA) #This is correct
heatmap(as.matrix(class_agg),cexRow = 0.5)#,Colv=NA)


###Check which iron_resistance mechanisms/groups are most common

#legend(x="right", legend=c("min", "med", "max"),fill=heat.colors(3))
#heatmap(as.matrix(t(class_clr)),cexRow = 0.5)
#dim(class_agg)
#View(as.matrix(class_agg))

###Create GGplot class heatmap by converting to long format
class_agg_meta<-cbind(meta,as.data.frame(class_agg)) #reattach metadata
class_agg_long<-gather(class_agg_meta,amr,count,Acetate_resistance:Zinc_resistance) #Long format
#View(class_agg_long)
#View(class_agg)
#View(class_agg_meta)
#GGplot heatmap
ggplot(data=class_agg_long,aes(x=CaseOrControl,y=amr,fill=count)) + geom_tile() +
  facet_wrap(~FarmId,ncol=4) + scale_fill_viridis_c(option="magma") 





