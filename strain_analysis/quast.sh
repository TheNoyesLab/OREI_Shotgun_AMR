#!/bin/bash



#database='/home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/staphA_reference/ncbi_dataset/data/GCF_000013425.1/StaphA_reference.fasta'
#reads='/scratch.global/fermx014/data/elder099/Noyes_Project_019/NonHostFastq'
straindir='/home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/strain_analysis'
reads="$straindir/test_reads"
filenames="$straindir/final_read_list.txt"

# Activate conda
. /home/noyes046/elder099/anaconda3/etc/profile.d/conda.sh
. /home/noyes046/elder099/anaconda3/etc/profile.d/mamba.sh

mamba activate quast

metaquast.py $straindir/spades_USDA898_S127.non.host_output/scaffolds.fasta $straindir/spades_USDA996_S67.non.host_output/scaffolds.fasta -o $straindir/quast_898_996_output -m 1000 -t 32
#Run metaquast on 898 and 996 metagenomic contigs
#minimum contig length of 1000, 32 threads

