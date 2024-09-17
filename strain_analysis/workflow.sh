#!/bin/bash

#reads='/scratch.global/fermx014/data/elder099/Noyes_Project_019/NonHostFastq'
straindir='/home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/strain_analysis'
reads="$straindir/test_reads"
assembly="/home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/strain_analysis/spades_output/K55/scaffolds.fasta"
test_reads='/home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/strain_analysis/test_reads'
final_bins='/home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/strain_analysis/binning_898_996/dRep_output/dereplicated_genomes/'


# Activate conda
. /home/noyes046/elder099/anaconda3/etc/profile.d/conda.sh
. /home/noyes046/elder099/anaconda3/etc/profile.d/mamba.sh


###Establish workspace/filepaths
#Move to strain workspace
cd $straindir

#Grab full file names
ls $reads > file_nm.txt
#Trim file names, keeping only Sample IDs
#sed -n 's/.non.host.R[12].fastq.gz//p' file_nm.txt > file_name_trim.txt
sed -n 's/.R[12].fastq.gz//p' file_nm.txt > file_name_trim.txt
#Sort and deduplicate list
sort -u file_name_trim.txt > final_read_list.txt

rm file_nm.txt
rm file_name_trim.txt


###Run individual assemblies
mamba activate strains

#for ((i=1;i<=10;i++))
#cat $filenames | while read file
tail $filenames -n +2 | while read file
do
       # ngless map.ngl -j 32 /scratch.global/fermx014/data/elder099/Noyes_Project_019/NonHostFastq /scratch.global/elder099/StaphA/final_file_names.txt
        echo "$file"

        spades.py --meta -t 64 -1 $reads/$file.R1.fastq.gz -2 $reads/$file.R2.fastq.gz -o $straindir/spades_${file}_output

        #Read online that want to set memory limit high (e.g. 500GB) and fewer threads (e.g. 16)
      # echo $i
done


