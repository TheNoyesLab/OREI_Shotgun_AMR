#!/bin/bash

#reads='/scratch.global/fermx014/data/elder099/Noyes_Project_019/NonHostFastq'
straindir='/home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/strain_analysis'
rundir='/home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/strain_analysis/coassembly_run_v4'
reads="$straindir/test_reads"

#Move to strain workspace
cd $straindir
mkdir -p coassembly_run_v4
cd $rundir

#Grab full file names
ls $reads > file_nm.txt
#Trim file names, keeping only Sample IDs
#sed -n 's/.non.host.R[12].fastq.gz//p' file_nm.txt > file_name_trim.txt
sed -n 's/.R[12].fastq.gz//p' file_nm.txt > file_name_trim.txt
#Sort and deduplicate list
sort -u file_name_trim.txt > indiv_read_list.txt #full list of individual reads

rm file_nm.txt
rm file_name_trim.txt


