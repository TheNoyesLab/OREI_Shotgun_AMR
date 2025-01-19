#!/bin/bash

#Establish file paths
rundir="$straindir/coassembly_run_v4"
assemblies="$rundir/assemblies"
coassemblies="$assemblies/coassemblies"
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

#tail $indiv_reads -n +2 | while read file
#do

run_MaxBin.pl -contig $assemblies/spades_${file}_output/scaffolds.fasta -reads $reads/$file.R1.fastq.gz -reads2 $reads/$file.R2.fastq.gz -threads 32 -min_contig_length 1500 -out $binning/maxbin2/${file}_bins/${file}_bin

#done



#cat $coRead_list | while read file
#do

run_MaxBin.pl -contig $assemblies/spades_${file}_output/scaffolds.fasta -reads $coReads/$file.R1.fastq.gz -reads2 $coReads/$file.R2.fastq.gz -threads 32 -min_contig_length 1500 -out $binning/maxbin2/${file}_bins/${file}_bin

#done

#/home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/strain_analysis/spades_USDA898_S127.non.host_output
