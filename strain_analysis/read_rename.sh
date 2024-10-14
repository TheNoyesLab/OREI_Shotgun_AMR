#!/bin/bash

#reads='/scratch.global/fermx014/data/elder099/Noyes_Project_019/NonHostFastq'
straindir='/home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/strain_analysis'
reads="$straindir/test_reads"
assembly="/home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/strain_analysis/spades_output/K55/scaffolds.fasta"
test_reads='/home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/strain_analysis/test_reads'
filenames="$straindir/final_read_list.txt"


# Activate conda
#. /home/noyes046/elder099/anaconda3/etc/profile.d/conda.sh
#. /home/noyes046/elder099/anaconda3/etc/profile.d/mamba.sh


###Fix read names to match Sample_ID
cd $test_reads

tail $filenames -n +2 | while read file
do
        #Remove _S127 etc suffix
        echo ${file}.R1.fastq.gz | sed "s/_S[0-9]*//"
        mv "${file}.R1.fastq.gz" "`echo ${file}.R1.fastq.gz | sed "s/_S[0-9]*//"`"

        echo ${file}.R2.fastq.gz | sed "s/_S[0-9]*//"
        mv "${file}.R2.fastq.gz" "`echo ${file}.R2.fastq.gz | sed "s/_S[0-9]*//"`"
done 
