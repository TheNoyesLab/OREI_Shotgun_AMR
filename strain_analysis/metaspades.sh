#!/bin/bash



#database='/home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/staphA_reference/ncbi_dataset/data/GCF_000013425.1/StaphA_reference.fasta'
#reads='/scratch.global/fermx014/data/elder099/Noyes_Project_019/NonHostFastq'
straindir='/home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/strain_analysis'
reads="$straindir/test_reads"
filenames="$straindir/final_read_list.txt"

# Activate conda
. /home/noyes046/elder099/anaconda3/etc/profile.d/conda.sh
conda activate strains

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
