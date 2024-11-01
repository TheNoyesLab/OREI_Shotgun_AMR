#!/bin/bash

#Establish file paths


#straindir='/home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/strain_analysis'
#test_reads='/home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/strain_analysis/test_reads'
#assembly="/home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/strain_analysis/spades_output/K55/scaffolds.fasta"
#reads="$straindir/test_reads" 
#filenames="$straindir/final_read_list.txt" 


#database='/home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/staphA_reference/ncbi_dataset/data/GCF_000013425.1/StaphA_reference.fasta'
#reads='/scratch.global/fermx014/data/elder099/Noyes_Project_019/NonHostFastq'
straindir='/home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/strain_analysis'
rundir="$straindir/coassembly_run_v4"
assemblies="$rundir/assemblies"
reads="$straindir/test_reads"
coReads="$reads/coReads"
indiv_reads="$rundir/indiv_read_list.txt"
coRead_list="$rundir/coRead_list.txt"
binning="$rundir/binning"
work="$binning/work_files"


# Activate conda
. /home/noyes046/elder099/anaconda3/etc/profile.d/conda.sh
. /home/noyes046/elder099/anaconda3/etc/profile.d/mamba.sh

#Activate conda environment
mamba activate metawrap #Change this once conda is fixed


#metawrap binning -o $straindir/binning_round1 -t 50 -a $assembly --metabat2 $test_reads/USDA871_S124.non.host.R_1.fastq $test_reads/USDA871_S124.non.host.R_2.fastq
mkdir -p "$rundir/binning"


#####
#####BINNING INDIVIDUAL ASSEMBLIES
#####
#for ((i=1;i<=10;i++))
cat $indiv_reads | while read file
#tail $indiv_reads -n +2 | while read file
do
        echo "$file"
	mkdir "$binning/work_files" #Create work directory
	

	###Align original reads to generated scaffold
	echo "alignment"
	#Index scaffold file
	bwa index $assemblies/spades_${file}_output/scaffolds.fasta
	
	#Make SAM of dehosted reads vs Scaffold
	bwa mem -t 80 $assemblies/spades_${file}_output/scaffolds.fasta $reads/$file.R1.fastq.gz $reads/$file.R2.fastq.gz > $work/$file.sam
	
	#Make BAM from SAM
	samtools view -@ 50 -bS $work/$file.sam > $work/$file.bam
	
	echo "sorting"
	#Sort BAM file
	samtools sort -@ 50 $work/$file.bam -o $work/${file}_sorted.bam

	echo "depth file"
	#Grab depth file from BAM
	jgi_summarize_bam_contig_depths -t 50 --outputDepth $work/depth.txt $work/${file}_sorted.bam




	###Run actual binning using metabat2
	echo "Binning"
	echo "Bins directory: $binning/metabat2/${file}_bins"
        mkdir -p $binning/metabat2/${file}_bins
	
	metabat2 -i $assemblies/spades_${file}_output/scaffolds.fasta -a $work/depth.txt -o $binning/metabat2/${file}_bins/${file}_bin -t 80 -m 1500
	

	###Complete binning script by deleting work files
	cd $binning
	rm -fr work_files

done



#####
#####BINNING COASSEMBLIES
#####

cat $coRead_list | while read file
#tail $coRead_list -n +2 | while read file
do
        echo "$file"
	mkdir "$binning/work_files" #Create work directory
	

	###Align original reads to generated scaffold
	echo "alignment"
	#Index scaffold file
	bwa index $assemblies/spades_${file}_output/scaffolds.fasta
	
	#Make SAM of dehosted reads vs Scaffold
	bwa mem -t 80 $assemblies/spades_${file}_output/scaffolds.fasta $coReads/$file.R1.fastq.gz $coReads/$file.R2.fastq.gz > $work/$file.sam
	
	#Make BAM from SAM
	samtools view -@ 50 -bS $work/$file.sam > $work/$file.bam
	
	echo "sorting"
	#Sort BAM file
	samtools sort -@ 50 $work/$file.bam -o $work/${file}_sorted.bam

	echo "depth file"
	#Grab depth file from BAM
	jgi_summarize_bam_contig_depths -t 50 --outputDepth $work/depth.txt $work/${file}_sorted.bam




	###Run actual binning using metabat2
	echo "Binning"
	echo "Bins directory: $binning/metabat2/${file}_bins"
	mkdir -p $binning/metabat2/${file}_bins
	
	metabat2 -i $assemblies/spades_${file}_output/scaffolds.fasta -a $work/depth.txt -o $binning/metabat2/${file}_bins/${file}_bin -t 80 -m 1500
	

	###Complete binning script by deleting work files
	cd $binning
	rm -fr work_files

done

