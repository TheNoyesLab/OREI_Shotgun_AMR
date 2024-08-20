#!/bin/bash

#Establish file paths


straindir='/home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/strain_analysis'
#test_reads='/home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/strain_analysis/test_reads'
assembly="/home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/strain_analysis/spades_output/K55/scaffolds.fasta"
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


#for ((i=1;i<=10;i++))
#cat $filenames | while read file
tail $filenames -n +2 | while read file
do
        echo "$file"
	
	echo "Bins directory: $straindir/binning_898_996/${file}_bins"
	#Note: Input reads must be uncompressed .fastq and have "R_1" "R_2" structure
	#Future: switch to a more flexible tool (just metabat2?)
	
	#Run metabat2 using scaffolds/assemblies and original reads
	metawrap binning -o $straindir/binning_898_996/${file}_bins -t 50 -a $straindir/spades_${file}_output/scaffolds.fasta  --metabat2 $reads/$file.R_1.fastq $reads/$file.R_2.fastq
	

	#head -n 1 $straindir/spades_${file}_output/scaffolds.fasta
	#head -n 1 $reads/$file.R_1.fastq
	#head -n 1 $reads/$file.R_2.fastq

      # echo $i
done
