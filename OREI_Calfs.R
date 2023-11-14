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


#Load in data
amr<-read.csv("~/Desktop/NoyesLab/OREI_Calfs/Noyes_Project_019_AMR++_Results/SNPconfirmed_AMR_analytic_matrix.csv")
metadata<-read.csv("~/Desktop/NoyesLab/OREI_Calfs/sample_data_06_08_2023.csv")
rownames(metadata)<-metadata$SampleId
#dim(amr)
#separate amr classifications into own columns
amr_classification<-amr %>% separate_wider_delim(cols = gene_accession,delim='|',names = c("accession","type" ,"class", "mechanism","group"),too_many="merge")
#View(amr_classification)

#Manually fix bugged Rifampin
amr_classification$group[amr_classification$accession=="MEG_8674"]<-
  paste0(amr_classification$accession[amr_classification$accession=="MEG_8674"],"|RequiresSNPConfirmation")

#Exclude low sensitivity SNP confirmed AMR
#amr_classification<-amr_classification[!grepl("RequiresSNPConfirmation",amr_classification$group),]

#Check on any 0 counts
#colnames(amr)[colSums(as.matrix(amr[,2:dim(amr)[2]]))==0]
#amr_classification[rowSums(amr_classification[,6:dim(amr_classification)[2]])==0,5])



#####
#####BEGIN DIFFABS FUNCTION
#####
man_diffabs<- function(x) {

###Group by and sum according to classification of AMR
#Should I aggregate with a sum, median, or mean?
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
#View(mech_t)

######Try some MANUAL differential abundance
#####Diff abundance w/o normalization

#mech_t_clr<-as.data.frame(clr(mech_t)) #CLR normalization
#class(mech_t_clr)
###Merge transposed datasets with metadata
mech_t$SampleId<-rownames(mech_t)
mech_t_farm<-mech_t %>% left_join(metadata,by=c("SampleId") )
#colnames(mech_t_farm)
#Keep only desired metadata columns
mechname<-names(mech_t_farm[,1:(dim(mech_t)[2]-1) ] )
##names(mech_t_farm[,!(names(mech_t_farm) %in% c("CaseOrControl","FarmId"))])
mechname
wilcoxon_p<-c()
logfc<-c()
for (i in mechname){
  result<-wilcox.test(mech_t_farm[,i] ~ CaseOrControl, data=mech_t_farm)
  logfc[[i]]<-mean( log2(mech_t_farm[mech_t_farm$CaseOrControl=="Case",i]+0.1) ) -
    mean( log2(mech_t_farm[mech_t_farm$CaseOrControl=="Control",i]+0.1) )
  wilcoxon_p[[i]]<-result$p.value
}
###Validated model of diff abundance
#FitZig: 0 inflated Gaussian model
#metagenomeSeq R package
#Includes diff abundance with CLR normallization
#ANCOM, ALDEX2 Normalization procedures
wilcoxon_p<-data.frame(class=names(wilcoxon_p),p_raw=unlist(wilcoxon_p),logfc=unlist(logfc))
wilcoxon_p$p_adj<-p.adjust(wilcoxon_p$p_raw,method = "fdr")

#View(wilcoxon_p)

wilcoxon_p$differential_abs<-ifelse(wilcoxon_p$p_adj<0.05 & wilcoxon_p$logfc>1,"UP",
                                    ifelse(wilcoxon_p$p_adj<0.05 & wilcoxon_p$logfc< -1,"DOWN","NO") )
wilcoxon_p$label<-ifelse(wilcoxon_p$differential_abs!="NO",wilcoxon_p$class,NA)

ssDA<-wilcoxon_p$label[!is.na(wilcoxon_p$label)]

#plot(wilcoxon_p$logfc,-log10(wilcoxon_p$p_adj))
plotty<-ggplot(data=wilcoxon_p,aes(x=logfc,y=-log10(p_adj),colour=differential_abs,label=label )) + geom_point() +
    geom_hline(yintercept=-log10(0.05),col="red") + geom_vline(xintercept=c(-1,1),col="red") + 
    geom_label_repel(size=3) + scale_color_manual(values=c("blue", "black", "red")) + 
    xlim(-2,2)  + theme_minimal()

return(list(ssDA,plotty) )

}

#####
#####END DIFFABS FUNCTION
#####

#Include SNP Confirmed
j<-man_diffabs(amr_classification) 
j[[1]]
j[[2]]
#Exclude SNP Confirmed
k<-man_diffabs(amr_classification[!grepl("RequiresSNPConfirmation",amr_classification$group),])
k[[1]]
k[[2]]
View(amr_classification)
k[[1]][ !( k[[1]] %in% j[[1]] ) ]  #These are added when SNPs removed
j[[1]][ !( j[[1]] %in% k[[1]] ) ]  #These are dropped when SNPs removed

#Next steps:
#Try out different tools (specifically ANCOM, ALDEX2)
#Turn above manual DA into a function w/graph output











#######Diff abundance w/normalization

mech_t_clr<-as.data.frame(clr(mech_t)) #CLR normalization
class(mech_t_clr)
###Merge transposed datasets with metadata
mech_t_clr$SampleId<-rownames(mech_t_clr)
mech_t_farm_clr<-mech_t_clr %>% left_join(metadata,by=c("SampleId") )

mechname<-names(mech_t_farm_clr[,c(1:145)] )
  ##names(mech_t_farm[,!(names(mech_t_farm) %in% c("CaseOrControl","FarmId"))])
mechname
wilcoxon_p<-c()
logfc<-c()
for (i in mechname){
  result<-wilcox.test(mech_t_farm_clr[,i] ~ CaseOrControl, data=mech_t_farm_clr)
  logfc[[i]]<-mean(mech_t_farm_clr[mech_t_farm_clr$CaseOrControl=="Case",i]) -
    mean(mech_t_farm_clr[mech_t_farm_clr$CaseOrControl=="Control",i])
  wilcoxon_p[[i]]<-result$p.value
}
wilcoxon_p<-data.frame(class=names(wilcoxon_p),p_raw=unlist(wilcoxon_p),logfc=unlist(logfc))
wilcoxon_p$p_adj<-p.adjust(wilcoxon_p$p_raw,method = "fdr")

View(wilcoxon_p)

plot(wilcoxon_p$logfc,-log10(wilcoxon_p$p_adj))





#######Heatmap creation below
#mech_t_clr$SampleId<-rownames(mech_t_clr)
#mech_t$SampleId<-rownames(mech_t)
#mech_t_farm<-mech_t %>% left_join(metadata,by=c("SampleId") )

# #Remove all unwanted metadata fields
# #colnames(mech_t_farm)
# #sum(is.na(mech_t_farm$WeeksToInfection)) #Confirmed no unmatched
# mech_t_farm_clr<-mech_t_farm_clr[,c(1:145,149,157)]  
# head(mech_t_farm_clr)
# ###Group by on farm and case/control
# #Should I aggregate with sum median or mean?
# mech_agg<-aggregate(.~CaseOrControl + FarmId,data=mech_t_farm,FUN = mean)
# #Group by on farm or other variables  -- group by sum or group by mean??
# 
# ###Make Heatmap with wide data
# rownames(mech_agg)<-paste(mech_agg[,1],mech_agg[,2])
# mech_agg<-mech_agg[,-c(1,2)]
# #Apply centered log ratio normalization
# #mech_clr<-as.matrix(clr(mech_agg)) #clr applies to rows, within sample
# 
# heatmap(as.matrix(t(mech_agg)),cexRow = 0.3)
# #heatmap(as.matrix(t(mech_clr)),cexRow = 0.3) #Normalized heatmap



########Try the heatmap with class
class_t_clr<-as.data.frame(clr(class_t)) #CLR normalization

###Merge transposed dataset with metadata
class_t_clr$SampleId<-rownames(class_t_clr)
class_t_farm_clr<-class_t_clr %>% left_join(metadata,by=c("SampleId") )
#sum(is.na(class_t_farm$WeeksToInfection)) #Confirmed no unmatched
colnames(class_t_farm_clr)
class_t_farm_clr<-class_t_farm_clr[,c(1:44,48,56)]  #Remove all unwanted metadata fields
#class_t_farm_clr<-class_t_farm_clr[,c(1:48,52,60)]  #Remove all unwanted metadata fields


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



#####End of heatmaps and diff abs
















###Microbiome count normalization
#Centered log-ratio (divide taxa counts by geometric mean)
mech_clr<-clr(mechsums)
class_clr<-clr(classsums)
group_clr<-clr(groupsums)
#View(mech_clr)
#Note, might break some assumptions of PCA

#NB Editors note:
#Chris used negative binomial model to associate S. aureus vs other taxa
#Also used the model w/emmeans to calculate *adjusted OTU abundances* How? What does this mean?
#Came up with adjusted richness/diversity metrics

#CLR Editor's note:
#Chris used CLR to normalize abundances for heatmaps





metadata$SampleId[!(metadata$SampleId %in% mech_t$SampleId)] %in% colnames(mechsums)

colnames(amr)[grepl("USDA2186",colnames(amr))]

inner_meta<-inner_join(mech_t,metadata,by="SampleId")

#PCA by amr mech
pcamech<-prcomp(mech_t[,-c(dim(mech_t)[2])])
autoplot(pcamech,data=inner_meta,colour="CaseOrControl") + facet_wrap(~FarmId,ncol=2)

#PCA by amr class
pcaclass<-prcomp(clr(class_t[,-c(dim(class_t)[2])]))
autoplot(pcaclass,data=inner_meta,colour="CaseOrControl") + facet_wrap(~FarmId,ncol=2)

pcaclass<-prcomp(clr(class_t[,-c(dim(class_t)[2])]))
autoplot(pcaclass,data=inner_meta,colour="FarmId") + facet_wrap(~CaseOrControl,ncol=2)

#PCA by amr group
pcagroup<-prcomp(clr(group_t[,-c(dim(group_t)[2])]))
autoplot(pcagroup,data=inner_meta,colour="CaseOrControl") + facet_wrap(~FarmId,ncol=2)

min(mech_clr)
max(mech_clr)










--------------------------------



###Group by and sum according to classification of AMR
#Should I aggregate with a sum, median, or mean?
#Agg ***sum*** accross classifications is correct
mechsums<-aggregate(.~mechanism,data=amr_classification[,c(4,6:dim(amr_classification)[2])],FUN=sum)
classsums<-aggregate(.~class,data=amr_classification[,c(3,6:dim(amr_classification)[2])],FUN=sum)
groupsums<-aggregate(.~group,data=amr_classification[,c(5,6:dim(amr_classification)[2])],FUN=sum)
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
#View(mech_t)

######Try some MANUAL differential abundance
#####Diff abundance w/o normalization

#mech_t_clr<-as.data.frame(clr(mech_t)) #CLR normalization
#class(mech_t_clr)
###Merge transposed datasets with metadata
mech_t$SampleId<-rownames(mech_t)
mech_t_farm<-mech_t %>% left_join(metadata,by=c("SampleId") )
#colnames(mech_t_farm)
#Keep only desired metadata columns
mechname<-names(mech_t_farm[,1:(dim(mech_t)[2]-1) ] )
##names(mech_t_farm[,!(names(mech_t_farm) %in% c("CaseOrControl","FarmId"))])
mechname
wilcoxon_p<-c()
logfc<-c()
for (i in mechname){
  result<-wilcox.test(mech_t_farm[,i] ~ CaseOrControl, data=mech_t_farm)
  logfc[[i]]<-mean( log2(mech_t_farm[mech_t_farm$CaseOrControl=="Case",i]+0.1) ) -
    mean( log2(mech_t_farm[mech_t_farm$CaseOrControl=="Control",i]+0.1) )
  wilcoxon_p[[i]]<-result$p.value
}
###Validated model of diff abundance
#FitZig: 0 inflated Gaussian model
#metagenomeSeq R package
#Includes diff abundance with CLR normallization
#ANCOM, ALDEX2 Normalization procedures
wilcoxon_p<-data.frame(class=names(wilcoxon_p),p_raw=unlist(wilcoxon_p),logfc=unlist(logfc))
wilcoxon_p$p_adj<-p.adjust(wilcoxon_p$p_raw,method = "fdr")

#View(wilcoxon_p)

wilcoxon_p$differential_abs<-ifelse(wilcoxon_p$p_adj<0.05 & wilcoxon_p$logfc>1,"UP",
                                    ifelse(wilcoxon_p$p_adj<0.05 & wilcoxon_p$logfc< -1,"DOWN","NO") )
wilcoxon_p$label<-ifelse(wilcoxon_p$differential_abs!="NO",wilcoxon_p$class,NA)



plot(wilcoxon_p$logfc,-log10(wilcoxon_p$p_adj))
ggplot(data=wilcoxon_p,aes(x=logfc,y=-log10(p_adj),colour=differential_abs,label=label )) + geom_point() +
  geom_hline(yintercept=-log10(0.05),col="red") + geom_vline(xintercept=c(-1,1),col="red") + 
  geom_label_repel(size=3) + scale_color_manual(values=c("blue", "black", "red")) + 
  xlim(-2,2)  + theme_minimal()



