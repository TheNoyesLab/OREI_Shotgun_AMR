#!/bin/bash

#Grab full file names
ls ./NonHostFastq/ > file_nm.txt
#Trim file names, keeping only Sample IDs
sed -n 's/.non.host.R[12].fastq.gz//p' file_nm.txt > file_name_trim.txt
#Sort and deduplicate list
sort -u file_name_trim.txt > final_file_names.txt

rm file_nm.txt
rm file_name_trim.txt

