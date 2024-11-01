#!/bin/bash

#reads='/scratch.global/fermx014/data/elder099/Noyes_Project_019/NonHostFastq'
straindir='/home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/strain_analysis'
rundir='/home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/strain_analysis/coassembly_run_v4'
test_reads='/home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/strain_analysis/test_reads'
coReads="$test_reads/coReads"

# Activate conda
#. /home/noyes046/elder099/anaconda3/etc/profile.d/conda.sh
#. /home/noyes046/elder099/anaconda3/etc/profile.d/mamba.sh

#Make coReads directory
cd $test_reads
mkdir -p coReads

###Read in csv file with groupings

if [ -f "$rundir/coassembly_groupR1.txt" ] ; then
    rm "$rundir/coassembly_group.txt"
    rm "$rundir/coassembly_groupR1.txt"
    rm "$rundir/coassebmly_groupR2.txt"
fi

for i in {1..10} #Loop through each cow group (317 groups)
do
	#Read through CSV file
	#cat $rundir/cow_coassembly_groups.csv | while IFS="," read -r col1 col2 col3 col4 col5
	cat $rundir/test_coassembly_groups.csv | while IFS="," read -r col1 col2 col3 col4 col5
	do
		#Output column 2 (Sample ID) and append to file
		if [ $col1 = $i ]; then
			#echo $test_reads/${col2}.non.host.R1.fastq.gz
			echo $col2 | tr -d '"' >> $rundir/coassembly_group.txt
			echo $test_reads/${col2}.non.host.R1.fastq.gz | tr -d '"' >> $rundir/coassembly_groupR1.txt
			echo $test_reads/${col2}.non.host.R2.fastq.gz | tr -d '"' >> $rundir/coassembly_groupR2.txt

    		fi
		

	done


	#Concatenate reads for each existing group
	if [ -f "$rundir/coassembly_groupR1.txt" ] ; then
		echo "next"

		cat $rundir/coassembly_group.txt
		cat $rundir/coassembly_groupR1.txt
		cat $rundir/coassembly_groupR2.txt
	        
		###Algo for concatenating reads
        	#cat (list of reads in each group here)
        	coread_name=$(cat $rundir/coassembly_group.txt | tr '\n' '_')
		echo $coread_name

		#Concatenate forward and reverse reads
		xargs cat < $rundir/coassembly_groupR1.txt > $coReads/${coread_name}.non.host.R1.fastq.gz
		xargs cat < $rundir/coassembly_groupR2.txt > $coReads/${coread_name}.non.host.R2.fastq.gz
		
		
		rm "$rundir/coassembly_group.txt"	
		rm "$rundir/coassembly_groupR1.txt"
		rm "$rundir/coassembly_groupR2.txt"
	fi
done


###Create final list of reads w/coReads

#Move to strain workspace
cd $rundir

#Grab full file names
ls $coReads > coRead_list.txt

#Trim file names, keeping only Sample IDs
#sed -n 's/.non.host.R[12].fastq.gz//p' file_nm.txt > file_name_trim.txt
sed -n 's/.R[12].fastq.gz//p' coRead_list.txt > coRead_list_trim.txt

#Sort and deduplicate list
sort -u coRead_list_trim.txt > coRead_list.txt #full list of coReads

rm coRead_list_trim.txt

###Final list of reads + coReads
cat indiv_read_list.txt coRead_list.txt > final_read_list.txt