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
dim(amr)
#separate amr classifications into own columns
amr_classification<-amr %>% separate_wider_delim(cols = gene_accession,delim='|',names = c("accession","type" ,"class", "mechanism","group"),too_many="merge")
#View(amr_classification)

#Check on any 0 counts
#colnames(amr)[colSums(as.matrix(amr[,2:dim(amr)[2]]))==0]
#amr_classification[rowSums(amr_classification[,6:dim(amr_classification)[2]])==0,5])


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

###Remove 0 counts
mechsums<-mechsums[,colSums(mechsums)!=0]
#sum(rowSums(mechsums)==0) #Fixed it
classsums<-classsums[,colSums(classsums)!=0]
#sum(rowSums(classsums)==0) #Fixed it
groupsums<-groupsums[,colSums(groupsums)!=0]
#sum(rowSums(groupsums)==0) #Fixed it



sum(colSums(groupsums)==0)
head(groupsums)





##########################METAGENOMESEQ DIFFABS


#Remove samples not present in both datasets
sum(!(colnames(classsums) %in% metadata$SampleId))
sum(!(metadata$SampleId %in% colnames(classsums)))
metadata<-metadata[metadata$SampleId %in% colnames(classsums),]

#order the samples properly bc touchy MRexperiment object
metadata_order<-metadata[order(rownames(metadata)),]
classsums_order<-classsums[,order(colnames(classsums))]

#Add a small number to prevent normalization issues
#classsums_order[classsums_order[1:dim(classsums_order)[1],]==0]<-0.001
#Pump in some changes to verify significance testing
#cases<-metadata_order$SampleId[metadata_order$CaseOrControl=="Case"]
#classsums_order[,colnames(classsums_order) %in% cases] <-classsums_order[,colnames(classsums_order) %in% cases] +0.01
#View(classsums_order)


###Convert our metadata into metagenomeseq-specific object
phen_metadata <- AnnotatedDataFrame(metadata_order)
#Combine new metadata with count matrix into new object
#Count matrix should be AMR in rows, samples in cols
#no transpose needed
experiment<-newMRexperiment(counts=classsums_order,phenoData=phen_metadata)

classsums_order[5,]

###TIME SERIES ON ONE FEATURE FOR NOW
###
res1<-fitTimeSeries(experiment,feature=5,
                    class='CaseOrControl',time='DIM',id='SampleId',
                    B=10,norm=T)
#Normalization is Cumulative Sum Normalization
#Must account for 0 values, add small number to each

res1[1] #p-value is always = 1/(B+1) -- Something is busted
#Confirmed that p-value is in face 1/(B+1) in the source code

#Create factor of Yes/No Case vs Control
classInfo <- factor(res1$data$class)

res2<- fitTimeSeries(experiment,feature=6,
                     class='CaseOrControl',time='DIM',id='SampleId',
                     B=10,norm=T,log=F)

plotClassTimeSeries(res1,pch=21,bg=classInfo)
plotTimeSeries(res1)

par(mfrow=c(1,1))
plotClassTimeSeries(res1,pch=21,bg=classInfo)
plotTimeSeries(res2)
plotClassTimeSeries(res2,pch=21,bg=classInfo)
###TIME SERIES ON ONE FEATURE FOR NOW
###




###LOOP THROUGH ALL FEATURES
###
feats<-1:(dim(experiment)[[1]]-35)
#experiment[1]
#Basically a for loop pumping into a list
timeSeriesFits <- lapply(feats,function(i){
  fitTimeSeries(obj=experiment,
                feature=i,
                class="CaseOrControl",
                id="SampleId",
                time="DIM",
                #lvl='class',
                #C=.3,# a cutoff for 'interesting' 
                B=10, # B is the number of permutations and should clearly not be 1
                norm=T)
})
#names(timeSeriesFits) <- classes


# Removing classes of bacteria without a potentially
# interesting time interval difference.
#looping through this: timeSeriesFits[[i]][[1]]
#Which is: timeSeriesFits[[classification]][[significance]]
alltimeSeriesFits <- lapply(timeSeriesFits,function(i){i[[1]]})
sstimeSeriesFits<-alltimeSeriesFits[-grep("No",alltimeSeriesFits)] #Extract only statistically significant entries

#rowbind all significant periods
sstimeSeriesFits_df <- do.call(rbind,sstimeSeriesFits) 
sstimeSeriesFits_df

#Bonferroni testing adjustment
pvalues = sstimeSeriesFits_df[,"p.value"]
adjPvalues = p.adjust(pvalues,"bonferroni")
sstimeSeriesFits_df = cbind(sstimeSeriesFits_df,adjPvalues)

head(sstimeSeriesFits_df)


###There are no statistically significant time points between Case and Control for 
###ANY of the AMR features
timeSeriesFits[[1]]

#Plot a couple of the fits
par(mfrow=c(2,1))
plotClassTimeSeries(timeSeriesFits[[7]],pch=21,bg=classInfo)
plotTimeSeries(timeSeriesFits[[7]])


##########################METAGENOMESEQ DIFFABS