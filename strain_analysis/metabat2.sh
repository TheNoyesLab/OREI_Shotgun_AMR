#!/bin/bash

#Establish file paths


straindir='/home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/strain_analysis'
#test_reads='/home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/strain_analysis/test_reads'
#assembly="/home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/strain_analysis/spades_output/K55/scaffolds.fasta"
reads="$straindir/test_reads" 
filenames="$straindir/final_read_list.txt" 

# Activate conda
. /home/noyes046/elder099/anaconda3/etc/profile.d/conda.sh
. /home/noyes046/elder099/anaconda3/etc/profile.d/mamba.sh

#Activate conda environment
mamba activate metawrap

#Note: Input reads must be uncompressed .fastq and have "R_1" "R_2" structure
#Future: switch to a more flexible tool (just metabat2?)

#metawrap binning -o $straindir/binning_round1 -t 50 -a $assembly --metabat2 $test_reads/USDA871_S124.non.host.R_1.fastq $test_reads/USDA871_S124.non.host.R_2.fastq
mkdir -p "$straindir/binning_v3"

#for ((i=1;i<=10;i++))
#cat $filenames | while read file
tail $filenames -n +2 | while read file
do
        echo "$file"
	
	work="$straindir/binning_v3/work_files"
	mkdir "$straindir/binning_v3/work_files"
	###Align original reads to generated scaffold
	
	echo "alignment"
	#Index scaffold file
	bwa index $straindir/spades_${file}_output/scaffolds.fasta
	
	#Make SAM of dehosted reads vs Scaffold
	bwa mem -t 50 $straindir/spades_${file}_output/scaffolds.fasta $reads/$file.R1.fastq.gz $reads/$file.R2.fastq.gz > $straindir/binning_v3/work_files/$file.sam
	
	#Make BAM from SAM
	samtools view -bS $straindir/binning_v3/work_files/$file.sam > $straindir/binning_v3/work_files/$file.bam
	
	echo "sorting"
	#Sort BAM file
	samtools sort -@ 50 $straindir/binning_v3/work_files/$file.bam -o $straindir/binning_v3/work_files/${file}_sorted.bam

	echo "depth file"
	#Grab depth file from BAM
	jgi_summarize_bam_contig_depths --outputDepth $straindir/binning_v3/work_files/depth.txt $straindir/binning_v3/work_files/${file}_sorted.bam
	
	echo "Bins directory: $straindir/binning_v3/${file}_bins"


	###Run actual binning using metabat2
	echo "Binning"
	mkdir $straindir/binning_v3/${file}_bins
	metabat2 -i $straindir/spades_${file}_output/scaffolds.fasta -a $straindir/binning_v3/work_files/depth.txt -o $straindir/binning_v3/${file}_bins/${file}_bin -t 50 -m 1500
	
	
	#Complete binning script by deleting work files
	cd $straindir/binning_v3
	rm -fr work_files

	#head -n 1 $straindir/spades_${file}_output/scaffolds.fasta
	#head -n 1 $reads/$file.R_1.fastq
	#head -n 1 $reads/$file.R_2.fastq

      # echo $i
done
