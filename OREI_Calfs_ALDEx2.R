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



#Load in data
amr<-read.csv("~/Desktop/NoyesLab/OREI_Calfs/Noyes_Project_019_AMR++_Results/SNPconfirmed_AMR_analytic_matrix.csv")
metadata<-read.csv("~/Desktop/NoyesLab/OREI_Calfs/sample_data_06_08_2023.csv")
rownames(metadata)<-metadata$SampleId
#dim(amr)

#separate amr classifications into own columns
amr_classification<-amr %>% separate_wider_delim(cols = gene_accession,delim='|',
                                                 names = c("accession","type" ,"class", "mechanism","group"),too_many="merge")
amr_classification<-as.data.frame(amr_classification)



View(amr_classification)

#####
#####BEGIN DATA CLEANING AND EXCLUSION FUNCTION
#####
clean_and_exclude<- function(amr_class,DIMfilter){
  
  class_cols<-c("accession","type","class","mechanism","group") #classification columns to keep
  
  ###SNP Confirmation issues
  #Manually fix bugged Rifampin
  amr_class$group[amr_class$accession=="MEG_8674"]<-
    paste0(amr_class$accession[amr_class$accession=="MEG_8674"],"|RequiresSNPConfirmation")
  
  #Exclude low sensitivity SNP confirmed AMR
  amr_class<-amr_class[!grepl("RequiresSNPConfirmation",amr_class$group),]

  
  ###Naming schemes
  #Keep only real USDA samples, no controls
  amr_class<-amr_class[,colnames(amr_class) %in% class_cols | grepl("USDA",colnames(amr_class))]

  #Reformat names to match metadata
  colnames(amr_class)<-gsub("(USDA\\d+)_.*",'\\1',colnames(amr_class))
  #sum(colnames(amr_class) %in% metadata$SampleId) == dim(amr_class)[2] #Confirm matching sample_ids

  
  ###Filter on DIM 
  if (DIMfilter == "after"){  #Days in Milk after birth
    new_DIM<-metadata$SampleId[metadata$DIM>0]               #List of IDs
    amr_class<-amr_class[,colnames(amr_class) %in% class_cols | colnames(amr_class) %in% new_DIM]  #Keep columns if in filtered list
    
  } else if(DIMfilter == "before"){ #Days in milk before birth
    new_DIM<-metadata$SampleId[metadata$DIM<=0]              #List of IDs
    amr_class<-amr_class[,colnames(amr_class) %in% class_cols | colnames(amr_class) %in% new_DIM]  #Keep columns if in filtered list
  } else if(DIMfilter =="all"){}


  ###Exclude 0-read samples
  zerocol<-which(colSums(amr_class[,6:dim(amr_class)[2]])==0) #vector of 0 columns 
  amr_class<-amr_class[, -(zerocol+5) ] #Remove 0 columns, shift 5 to match columns
  #Find a way to operationalize this

  return(amr_class)
  #sum(colSums(amr_classification[,6:dim(amr_classification)[2]])==0)
  #View(amr_classification)
  #Check on any 0 counts
  #colnames(amr)[colSums(as.matrix(amr[,2:dim(amr)[2]]))==0]
  #amr_classification[rowSums(amr_classification[,6:dim(amr_classification)[2]])==0,5])

}



#####
#####BEGIN COUNT PROCESSING FUNCTION
#####
class_process<- function(x) {

  ###Group by and sum according to classification of AMR
  #Agg ***sum*** accross classifications is correct
  mechsums<-aggregate(.~mechanism,data=x[,c(4,6:dim(x)[2])],FUN=sum)
  classsums<-aggregate(.~class,data=x[,c(3,6:dim(x)[2])],FUN=sum)
  groupsums<-aggregate(.~group,data=x[,c(5,6:dim(x)[2])],FUN=sum)
  
  #Convert column to rowname, delete columns
  rownames(mechsums)<-mechsums$mechanism
  rownames(classsums)<-classsums$class
  rownames(groupsums)<-groupsums$group
  mechsums<- mechsums %>% dplyr::select(-c(mechanism))
  classsums<- classsums %>% dplyr::select(-c(class))
  groupsums<- groupsums %>% dplyr::select(-c(group))
  
  #mechsums<-rowsum(amr_classification[,c(6:dim(amr_classification)[2])],group=amr_classification$mechanism);dim(mechsums) #Sum of counts by mech
  #classsums<-rowsum(amr_classification[,c(6:dim(amr_classification)[2])],group=amr_classification$class);dim(classsums) #Sum of counts by class
  #groupsums<-rowsum(amr_classification[,c(6:dim(amr_classification)[2])],group=amr_classification$group);dim(groupsums) #Sum of counts by group
  
  
  
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
  names(mat_list)<-c("mechanism","class","group")
  
  return(mat_list)
  
}






amr_class_clean_pre<-clean_and_exclude(amr_classification,"before")
amr_class_clean_post<-clean_and_exclude(amr_classification,"after")
amr_class_clean_all<-clean_and_exclude(amr_classification,"all")

amr_list<-class_process(amr_class_clean_all)

################################END OBJECT CREATION SECTION






################################START BARPLOT SECTION
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
#Create & reorder interaction variable for x-axis
top_10$FarmCase<-factor(interaction(top_10$FarmId,top_10$CaseOrControl),
                        levels= c("Farm A.Control","Farm A.Case","Farm B.Control","Farm B.Case",
                                  "Farm C.Control","Farm C.Case","Farm D.Control","Farm D.Case")
)

#barplot
ggplot(data=top_10,aes(x=FarmCase,y=count,fill=amr)) + 
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




#####
#####GRAB METADATA
#####
#Note that you really only need one copy of the metadata... But I did 3 anyways

amr_meta<-amr_list #Rename use in loop

#Add sample_id to each dataset
for(i in names(amr_meta)) {
  amr_meta[[i]]$SampleId<-rownames(amr_meta[[i]])
}

#Merge to grab metadata
amr_meta<-lapply(amr_meta,function(x){ left_join(x,metadata[c("SampleId","CaseOrControl","FarmId","DIM","Batch")],by="SampleId") %>%
      dplyr::select("SampleId","CaseOrControl","FarmId","DIM","Batch") %>%  #keep only metadata
      column_to_rownames("SampleId") }  #Add rownames back
                   )








#####
#####START PCA SECTION
#####

###Apply CLR normalization to amr object
amr_clr<-lapply(amr_list,function(x) as.data.frame(clr(x)) )

#Apply PCA to clr-transformed data
amr_pca<-lapply(amr_clr,function(x) prcomp(x)) #use only numeric values


autoplot(amr_pca$mechanism,data = amr_meta$mechanism,colour = "DIM") +
  scale_color_gradient2(midpoint=0,low="royalblue4",mid="brown2",high="yellow") +
  facet_wrap(~FarmId)
autoplot(amr_pca$class,data = amr_meta$class,colour = "DIM") +
  scale_color_gradient2(midpoint=0,low="royalblue4",mid="brown2",high="yellow") +
  facet_wrap(~FarmId)
autoplot(amr_pca$group,data = amr_meta$group,colour = "DIM") +
  scale_color_gradient2(midpoint=0,low="royalblue4",mid="brown2",high="yellow") +
  facet_wrap(~FarmId)

#Editors note:
#Chris used negative binomial model to associate S. aureus vs other taxa
#Also used the model w/emmeans to calculate *adjusted OTU abundances* How? What does this mean?
#Came up with adjusted richness/diversity metrics


#######################PERMANOVA
#prcomp used for running PCA
#Aitchison Distance matrix was computed using the vegdist function in ‘vegan’
#Aitchison distance matrix (euclidean distance of a centered-log ratio scaled abundance matrix) was computed and used as input to the adonis function in the ‘vegan’ package.
#test if the differences between time period distance were greater than those within-time distances using the adonis function in ‘vegan’

amr_dist<-lapply(amr_list, function(x) vegdist(as.data.frame(clr(x)),method="euclidean") )

amr_test<-lapply(amr_dist,function(x) adonis2(x ~ DIM + CaseOrControl + FarmId + Batch,data=amr_meta$mechanism,method = "aitchison",parallel=4))

#Adonis PERMANOVA test results
amr_test$mechanism
amr_test$class
amr_test$group

#All independent variables are significant
#FarmId explains roughly 8% of variance
#DIM explains roughly 0.5-0.6% of variance
#CaseOrControl explains roughly 0.3% of variance
#Minimal batch effects




################################START ALDEX2 SECTION
aldex_amr<- function(amr_list,amr) {
  
  amrdf<-amr_list[[amr]] #Call data for the type of amr
  
  ###ALDEX2 Preprocessing/Formatting
  #Make SampleId column for join
  amrdf$SampleId<-rownames(amrdf)
 
  #Join on metadata to get Case/Control status
  CorC<-amrdf %>% left_join(metadata[c("SampleId","CaseOrControl")],by=c("SampleId") ) %>% 
    dplyr::select('SampleId','CaseOrControl')  #Extract CaseOrControl in correct order for ALDEX2
  CorC

  
  #Transpose class dataframe and exclude SampleId column
  aldex_classdf<-t( amrdf[-c(which(colnames(amrdf) == "SampleId"))] )
  
  class.aldex <- aldex(aldex_classdf, CorC$CaseOrControl, mc.samples=128, test="t", effect=TRUE, 
                       include.sample.summary=FALSE, denom="all", verbose=FALSE, paired.test=T, gamma=NULL)
  
  
  #Plot 3 significance testing graphs
  par(mfrow=c(1,3))
  aldex.plot(class.aldex, type="MA", test="welch", xlab="Log-ratio abundance",
             ylab="Difference", main='Bland-Altman plot',called.cex = 3 )
  aldex.plot(class.aldex, type="MW", test="welch", xlab="Dispersion",
             ylab="Difference", main='Effect plot',called.cex = 3)
  aldex.plot(class.aldex, type="volcano", test="welch", xlab="Difference",
             ylab="-1(log10(q))", main= 'Volcano plot',ylim = c(0,1.5),xlim = c(-2,2),called.cex = 3)  #volcano plot
  #mtext(paste(toupper(amr),"Signficance Plots"), side = 3,line=2, outer = TRUE)     #Add title with amr type
  
  
}

#Make plots and convert to ggplot
p<-as.ggplot( function() aldex_amr(amr_list,"class"),scale = 0.9 )
q<-as.ggplot( function() aldex_amr(amr_list,"mechanism"),scale = 0.9 )
r<-as.ggplot( function() aldex_amr(amr_list,"group"),scale = 0.9 )

#Add borders to plots
p<-p+theme(panel.border = element_rect(colour = "black", fill=NA, size=1) )
q<-q+theme(panel.border = element_rect(colour = "black", fill=NA, size=1) )
r<-r+theme(panel.border = element_rect(colour = "black", fill=NA, size=1) )
p
q
r
ggarrange(p,q,r,labels=c("Class","Mech","Group"),nrow = 3)



#aldex_amr(amr_list$class)
#aldex_amr(amr_list$mechansim)
#aldex_amr(amr_list$group)



#######################% of Total Reads
#View(class_agg_long)










