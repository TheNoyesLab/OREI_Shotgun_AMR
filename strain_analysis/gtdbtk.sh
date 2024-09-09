#!/bin/bash

#Establish file paths

assembly="/home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/strain_analysis/spades_output/K55/scaffolds.fasta"
straindir='/home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/strain_analysis'
test_reads='/home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/strain_analysis/test_reads'
final_bins='/home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/strain_analysis/binning_898_996/dRep_output/dereplicated_genomes/'


# Activate conda
. /home/noyes046/elder099/anaconda3/etc/profile.d/conda.sh
. /home/noyes046/elder099/anaconda3/etc/profile.d/mamba.sh

#Activate conda environment
mamba activate mag_taxon #gtdb-tk is here

gtdbtk classify_wf --genome_dir $final_bins --out_dir $straindir/gtdb-tk_output --extension fa --cpus 32 --mash_db $straindir/gtdb-tk_output/mash_db
