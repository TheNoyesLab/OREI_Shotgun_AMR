#!/bin/bash

#Establish file paths

assembly="/home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/strain_analysis/spades_output/K55/scaffolds.fasta"
straindir='/home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/strain_analysis'
test_reads='/home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/strain_analysis/test_reads'

# Activate conda
. /home/noyes046/elder099/anaconda3/etc/profile.d/conda.sh
. /home/noyes046/elder099/anaconda3/etc/profile.d/mamba.sh

#Activate conda environment
mamba activate metawrap


metawrap binning -o $straindir/binning_round1 -t 50 -a $assembly --metabat2 $test_reads/USDA871_S124.non.host.R_1.fastq $test_reads/USDA871_S124.non.host.R_2.fastq


