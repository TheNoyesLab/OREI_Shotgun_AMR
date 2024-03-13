library(ggfortify)
#library(compositions)
library(stats)
library(compositions)
library(DESeq2)
#install.packages("BiocManager")
library(BiocManager)
#BiocManager::install("DESeq2")
library(ggrepel)
#BiocManager::install("metagenomeSeq")
library(metagenomeSeq)
#BiocManager::install("ALDEx2")
library(tidyverse)
#devtools::install_github("ggloor/ALDEx_bioc")
library(ALDEx2)
library(ggpubr)
#install.packages("ggplotify")
library(ggplotify)
library(vegan)






snps<-read.csv("~/Desktop/NoyesLab/OREI_Calfs/Final_Meta_SNPs.csv")
metadata<-read.csv("~/Desktop/NoyesLab/OREI_Calfs/sample_data_06_08_2023.csv")
rownames(metadata)<-metadata$SampleId
dim(snps)
#separate amr classifications into own columns
amr_classification<-snps %>% separate_wider_delim(cols = SNP_accession,delim='|',
                                                  names = c("accession","type" ,"class", "mechanism","group"),too_many="merge")
amr_classification<-as.data.frame(amr_classification)
#colSums(amr_classification[,6:854])
#View(amr_classification)
#####
#####BEGIN DATA CLEANING AND EXCLUSION FUNCTION
#####
clean_and_exclude<- function(amr_class,DIMfilter){
  
  class_cols<-c("accession","type","class","mechanism","group") #classification columns to keep
  
  ###SNP Confirmation issues
  #Manually fix bugged Rifampin
  #amr_class$group[amr_class$accession=="MEG_8674"]<-
  #  paste0(amr_class$accession[amr_class$accession=="MEG_8674"],"|RequiresSNPConfirmation")
  
  #Exclude low sensitivity SNP confirmed AMR
  #amr_class<-amr_class[!grepl("RequiresSNPConfirmation",amr_class$group),]
  
  
  ###Naming schemes
  #Keep only real USDA samples, no controls
  amr_class<-amr_class[,colnames(amr_class) %in% class_cols | grepl("USDA",colnames(amr_class))]
  
  #Reformat names to match metadata
  colnames(amr_class)<-gsub("(USDA\\d+)_.*",'\\1',colnames(amr_class))
  #sum(colnames(amr_class) %in% metadata$SampleId) == dim(amr_class)[2] #Confirm matching sample_ids
  
  
  ###Filter on DIM 
  # if (DIMfilter == "after"){  #Days in Milk after birth
  #   new_DIM<-metadata$SampleId[metadata$DIM>0]               #List of IDs
  #   amr_class<-amr_class[,colnames(amr_class) %in% class_cols | colnames(amr_class) %in% new_DIM]  #Keep columns if in filtered list
  #   
  # } else if(DIMfilter == "before"){ #Days in milk before birth
  #   new_DIM<-metadata$SampleId[metadata$DIM<=0]              #List of IDs
  #   amr_class<-amr_class[,colnames(amr_class) %in% class_cols | colnames(amr_class) %in% new_DIM]  #Keep columns if in filtered list
  # } else if(DIMfilter =="all"){}
  
  ###Exclude 0-read samples
  #zerocol<-which(colSums(amr_class[,6:dim(amr_class)[2]])==0) #vector of 0 columns 
  #amr_class<-amr_class[, -(zerocol+5) ] #Remove 0 columns, shift 5 to match columns
  #Find a way to operationalize this
  
  return(amr_class)
  
}


amr_class_clean_all<-clean_and_exclude(amr_classification,"all")
dim(amr_class_clean_all)




#####
#####Remove SNP-Confirmation Only Features
#####

annots<-read.csv("~/Desktop/NoyesLab/OREI_Calfs/megares_annotations_v3.00.csv")
#View(annots)

snponly<-annots[grepl("RequiresSNPConfirmation",annots$header),] #Only snp confirmed
trimsnp<-gsub("(MEG_\\d+)\\|.*","\\1",snponly$header) #Extract accession #

amr_class_clean_all<-amr_class_clean_all[!(amr_class_clean_all$accession %in% trimsnp),]
dim(amr_class_clean_all)


#####
#####BEGIN COUNT PROCESSING FUNCTION
#####
class_process<- function(x) {
  
  ###Group by and sum according to classification of AMR
  #Agg ***sum*** accross classifications is correct
  mechsums<-aggregate(.~mechanism,data=x[,c(4,6:dim(x)[2])],FUN=sum)
  classsums<-aggregate(.~class,data=x[,c(3,6:dim(x)[2])],FUN=sum)
  #groupsums<-aggregate(.~group,data=x[,c(5,6:dim(x)[2])],FUN=sum)
  
  #Convert column to rowname, delete columns
  rownames(mechsums)<-mechsums$mechanism
  rownames(classsums)<-classsums$class
  #rownames(groupsums)<-groupsums$group
  mechsums<- mechsums %>% dplyr::select(-c(mechanism))
  classsums<- classsums %>% dplyr::select(-c(class))
  #groupsums<- groupsums %>% dplyr::select(-c(group))
  
  #mechsums<-rowsum(amr_classification[,c(6:dim(amr_classification)[2])],group=amr_classification$mechanism);dim(mechsums) #Sum of counts by mech
  #classsums<-rowsum(amr_classification[,c(6:dim(amr_classification)[2])],group=amr_classification$class);dim(classsums) #Sum of counts by class
  #groupsums<-rowsum(amr_classification[,c(6:dim(amr_classification)[2])],group=amr_classification$group);dim(groupsums) #Sum of counts by group
  
  
  
  ###Remove 0 counts
  mechsums<-mechsums[rowSums(mechsums)!=0,]
  #sum(rowSums(mechsums)==0) #Fixed it
  classsums<-classsums[rowSums(classsums)!=0,]
  #sum(rowSums(classsums)==0) #Fixed it
  #groupsums<-groupsums[rowSums(groupsums)!=0,]
  #sum(rowSums(groupsums)==0) #Fixed it
  
  
  #transpose for analysis
  mech_t<-as.data.frame(t(mechsums) )
  class_t<-as.data.frame(t(classsums) )
  #group_t<-as.data.frame(t(groupsums) )
  
  #Put into object
  mat_list<-list(mech_t,class_t)#,group_t)
  names(mat_list)<-c("mechanism","class")#,"group")
  
  return(mat_list)
  
}


amr_list<-class_process(amr_class_clean_all)






#####
#####GRAB FULL METADATA
#####
amr_meta<-amr_list$mechanism #Load in amr_list for join
amr_meta$SampleId<-rownames(amr_meta)

#Left join to grab metadata fields
amr_meta<-left_join(amr_meta,metadata[c("SampleId","CaseOrControl","FarmId","DIM","Batch")],by="SampleId") %>%
  dplyr::select("SampleId","CaseOrControl","FarmId","DIM","Batch") %>%  #keep only metadata
  column_to_rownames("SampleId") %>%
  mutate("Calving" = ifelse(DIM>=0,"Postnatal","Prenatal")) #Divide DIM into case-control

#Confirm metadata is in same order, identical rownames regardless of type
#identical(rownames(amr_list$class),rownames(amr_meta))



#####
#####SPLIT DATA INTO FARMS
#####
amr_farm<-lapply(amr_list,function(x) split(x,amr_meta$FarmId) )    #Split counts
amr_farm_meta<-split(amr_meta,amr_meta$FarmId)   #Split metadata

#Confirm metadata is in same order, identical rownames regardless of type, Farm
#identical(rownames(amr_farm$mechanism$`Farm A`),rownames(amr_farm_meta$'Farm A'))
#identical(rownames(amr_farm$class$`Farm B`),rownames(amr_farm_meta$'Farm B'))






################################START BARPLOT SECTION -- WHOLE DATASET
####Maybe start function here, input is FULL class_process object

###Merge transposed dataset with metadata
amr_list$class$SampleId<-rownames(amr_list$class)
class_t_meta<-amr_list$class %>% left_join(metadata[c("SampleId","CaseOrControl","FarmId")],by=c("SampleId") ) #Keep certain columns
#sum(is.na(class_t_farm$WeeksToInfection)) #Confirmed no unmatched

#colnames(class_t_meta)
class_t_meta<-class_t_meta[,!( colnames(class_t_meta) %in% c("SampleId") )]  #Remove all unwanted metadata fields
#class_t_farm_clr<-class_t_farm_clr[,c(1:48,52,60)]  #Remove all unwanted metadata fields

###Aggregate counts by Case-Control & Farm, *Average
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
#Create & reorder interaction variable for x-axis
top_10$FarmCase<-factor(interaction(top_10$FarmId,top_10$CaseOrControl),
                        levels= c("Farm A.Control","Farm A.Case","Farm B.Control","Farm B.Case",
                                  "Farm C.Control","Farm C.Case","Farm D.Control","Farm D.Case")
)

#barplot
ggplot(data=top_10,aes(x=FarmCase,y=count,fill=amr)) + 
  geom_bar(stat="identity") +
  theme_minimal() +
  scale_fill_manual(name="AMR Class", values=colorRampPalette(brewer.pal(11,"Spectral"))(16) ) +
  theme(legend.title = element_text(size=10)) +
  theme(axis.text.x = element_text(angle = -45,hjust=0.15,vjust=-0.4)) +
  theme(legend.text = element_text(size=7)) +
  theme(legend.key.height= unit(3, 'mm'),
        legend.key.width= unit(5, 'mm')) +
  ggtitle("Top 10 AMR Classes") +
  theme(plot.title = element_text(hjust=0.5)) +
  xlab("Case or Control by Farm") +
  ylab("AMR Class Proportion")






#####
#####START PCA SECTION
#####

###Apply CLR normalization to SPLIT amr object
amr_clr<-lapply(amr_farm$mechanism,function(x) as.data.frame(clr(x)) )

#Apply PCA to clr-transformed data
amr_pca<-lapply(amr_clr,function(x) prcomp(x)) #use only numeric values

A<-autoplot(amr_pca$`Farm A`,data = amr_farm_meta$`Farm A`,colour = "DIM") + scale_color_gradient2(midpoint=0,low="royalblue4",mid="brown2",high="yellow") + theme_minimal() + ggtitle("Farm A") + theme(plot.title = element_text(hjust=0.5))
#loadings=T,loadings.label=T,loadings.size=0.1,loadings.label.size=2,loadings.color="blue",loadings.label.color="black")
B<-autoplot(amr_pca$`Farm B`,data = amr_farm_meta$`Farm B`,colour = "DIM") + scale_color_gradient2(midpoint=0,low="royalblue4",mid="brown2",high="yellow") + theme_minimal() + ggtitle("Farm B") + theme(plot.title = element_text(hjust=0.5))
C<-autoplot(amr_pca$`Farm C`,data = amr_farm_meta$`Farm C`,colour = "DIM") + scale_color_gradient2(midpoint=0,low="royalblue4",mid="brown2",high="yellow") + theme_minimal() + ggtitle("Farm C") + theme(plot.title = element_text(hjust=0.5))
D<-autoplot(amr_pca$`Farm D`,data = amr_farm_meta$`Farm D`,colour = "DIM") + scale_color_gradient2(midpoint=0,low="royalblue4",mid="brown2",high="yellow") + theme_minimal() + ggtitle("Farm D") + theme(plot.title = element_text(hjust=0.5))


ggarrange(A,B,C,D,common.legend = T,legend="right")

















#######################PERMANOVA
#prcomp used for running PCA
#Aitchison Distance matrix was computed using the vegdist function in ‘vegan’
#Aitchison distance matrix (euclidean distance of a centered-log ratio scaled abundance matrix) was computed and used as input to the adonis function in the ‘vegan’ package.
#test if the differences between time period distance were greater than those within-time distances using the adonis function in ‘vegan’


#SPLIT DATASET
amr_dist<-lapply(amr_farm$mechanism, function(x) vegdist(as.data.frame(clr(x)),method="euclidean") )


amr_test_list<-list()
for(i in names(amr_dist)){
  amr_test<-adonis2(amr_dist[[i]] ~ DIM + CaseOrControl + Batch,data=amr_farm_meta[[i]],method = "aitchison",parallel=4)
  amr_test_list[[i]]<-amr_test
}

amr_test_list

###Conclusion
#Days in Milk is significant for all 4 farms
#Days in Milk explains 1.8-4.1% of the variance in the farms
#CaseOrControl is significant for Farm B & D ONLY
#No notable batch effects
#Clear trends in Farm A and B, but seems to offer minimal explanatory value
#PCA may not pick up on nuanced associations, Random Forest?




#FULL DATASET
amr_dist<-lapply(amr_list, function(x) vegdist(as.data.frame(clr(x)),method="euclidean") )

amr_test<-lapply(amr_dist,function(x) adonis2(x ~ DIM + CaseOrControl + FarmId + Batch,data=amr_meta,method = "aitchison",parallel=4))

#Adonis PERMANOVA test results
amr_test$mechanism
#amr_test$class
#amr_test$group

#All independent variables are significant
#FarmId explains roughly 6% of variance
#DIM explains roughly 1.5% of variance
#CaseOrControl explains roughly 0.2% of variance
#Batch effects explain 0.2%

