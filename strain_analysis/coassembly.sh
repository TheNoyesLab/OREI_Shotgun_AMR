#!/bin/bash

#reads='/scratch.global/fermx014/data/elder099/Noyes_Project_019/NonHostFastq'
straindir='/home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/strain_analysis'
reads="$straindir/test_reads"
assembly="/home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/strain_analysis/spades_output/K55/scaffolds.fasta"
test_reads='/home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/strain_analysis/test_reads'
final_bins='/home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/strain_analysis/binning_898_996/dRep_output/dereplicated_genomes/'


# Activate conda
#. /home/noyes046/elder099/anaconda3/etc/profile.d/conda.sh
#. /home/noyes046/elder099/anaconda3/etc/profile.d/mamba.sh


###Read in csv file with groupings

if [ -f "$straindir/coassembly_group.txt" ] ; then
    rm "$straindir/coassembly_group.txt"
fi

for i in {1..4} #Loop through each cow group (317 groups)
do
	#Read through CSV file
	cat $straindir/cow_coassembly_groups.csv | while IFS="," read -r col1 col2 col3 col4 col5
	do
		#Output column 2 (Sample ID) and append to file
		if [ $col1 = $i ]; then
			echo $col2 | tr -d '"' >> $straindir/coassembly_group.txt
    		fi
		#echo "I got:$col1|$col2"

		#####
		#####EITHER: Create 317 text files 
		#####OR: Run assemblies as files are created, then delete them (cleaner server, uglier code)


	done
	if [ -f "$straindir/coassembly_group.txt" ] ; then
		echo "beep"
		cat $straindir/coassembly_group.txt
		
		rm "$straindir/coassembly_group.txt"
	fi
	#rm $straindir/coassembly_group.txt
done
