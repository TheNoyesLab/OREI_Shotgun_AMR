#!/bin/bash

#Establish file paths

assembly="/home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/strain_analysis/spades_output/K55/scaffolds.fasta"
straindir='/home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/strain_analysis'
test_reads='/home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/strain_analysis/test_reads'
final_bins='/home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/strain_analysis/binning_898_996/dRep_output/dereplicated_genomes/'
gtdbtk='/home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/strain_analysis/gtdb-tk_output'

# Activate conda
. /home/noyes046/elder099/anaconda3/etc/profile.d/conda.sh
. /home/noyes046/elder099/anaconda3/etc/profile.d/mamba.sh

#Activate conda environment
mamba activate mag_taxon #fasttree is here

###Build phylogenetic tree of *user input genomes only*
zcat $gtdbtk/gtdbtk.bac120.user_msa.fasta.gz | fasttree  > $straindir/output_tree2.nwk
