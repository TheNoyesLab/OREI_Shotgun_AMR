#!/bin/bash

#Establish file paths

assembly="/home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/strain_analysis/spades_output/K55/scaffolds.fasta"
straindir='/home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/strain_analysis'
test_reads='/home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/strain_analysis/test_reads'

# Activate conda
. /home/noyes046/elder099/anaconda3/etc/profile.d/conda.sh
. /home/noyes046/elder099/anaconda3/etc/profile.d/mamba.sh

#Activate conda environment
mamba activate strains #dRep is in here

#Create list of bins for input to dRep
#Future: Edit this so it's a for loop concatting all of the different bin paths

#Enter directory of bins
cd $straindir/binning_898_996/USDA898_S127.non.host_bins/metabat2_bins
#Change bin names to add Sample ID so they are unique
for f in * ; do mv -- "$f" "USDA898_$f" ; done
#Export list of bins paths to file
ls -d -1 "$PWD/"*[!unbinned].fa* > $straindir/binning_898_996/bin_list.txt

#Enter directory of bins
cd $straindir/binning_898_996/USDA996_S67.non.host_bins/metabat2_bins
#Chagne bin names to add Sample ID so they are unique
for f in * ; do mv -- "$f" "USDA996_$f" ; done
#Export list of bins paths to file, overwriting
ls -d -1 "$PWD/"*[!unbinned].fa* >> $straindir/binning_898_996/bin_list.txt

dRep dereplicate -p 50 $straindir/binning_898_996/dRep_output -g "$straindir/binning_898_996/bin_list.txt" -pa 0.95 -sa 0.99


