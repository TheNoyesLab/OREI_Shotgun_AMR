#!/bin/bash

#####
#####START METASNV
#####
echo "Start MetaSNV"

###Start conda
. /home/noyes046/elder099/anaconda3/etc/profile.d/conda.sh

#Start metasnv environment
conda activate metasnv

#Create list of BAMs using full path
cd outputs
ls -d -1 "$PWD/"** > ../meta_sample_list.txt
#Create list of BAMs using file name
ls > ../sample_headers.txt

cd ..

#sub_db="$Bench/SNP_Injector/DB_noSNP_subset${j}.fasta"


#/usr/bin/time -v -o $Bench/S$j/M$i/Meta_Time.txt metaSNV.py $Bench/S$j/M$i/Meta_Out \
#        		$Bench/S$j/M$i/meta_sample_list.txt \
#        		$sub_db --threads 32

OREI="/scratch.global/elder099/Noyes_Project_019/all_pools"
megares="/home/noyes046/shared/databases/megares_v2.0/megares_full_database_v2.00.fasta"

metaSNV.py $OREI/Meta_Out $OREI/meta_sample_list.txt $megares --threads 48

#mv $Bench/S$j/M$i/Meta_Time.txt $Bench/S$j/M$i/Meta_Out/Meta_Time.txt #Move from outside in because Meta_Out didn't exist
#cat $Bench/S$j/M$i/Meta_Out/snpCaller/called_SNPs.best_split_* | sed -n -e 's/[0-9][0-9]*|//g; /[ACTG]/ s/|.*//p' > $Bench/S$j/M$i/Meta_Out/Meta_Out_Fix.csv #Fix to look like normal csv

cat $OREI/Meta_Out/snpCaller/called_SNPs.best_split_*  > $OREI/Meta_Out/Meta_Out_Fix.txt


