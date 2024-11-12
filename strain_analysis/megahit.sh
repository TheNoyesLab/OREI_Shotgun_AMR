#!/bin/bash



#database='/home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/staphA_reference/ncbi_dataset/data/GCF_000013425.1/StaphA_reference.fasta'
#reads='/scratch.global/fermx014/data/elder099/Noyes_Project_019/NonHostFastq'
straindir='/scratch.global/elder099'
rundir='/scratch.global/elder099/strains_run_v5'
#straindir='/home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/strain_analysis'
#rundir="$straindir/coassembly_run_v4"
assemblies="$rundir/assemblies"
#reads="$straindir/test_reads"
reads='/scratch.global/fermx014/help/elder099/Noyes_Project_019/NonHostFastq'
coReads="$reads/coReads"
indiv_reads="$rundir/indiv_read_list.txt"
#indiv_reads="$rundir/tmp_read_list.txt"
coRead_list="$rundir/coRead_list.txt"

# Activate conda
. /home/noyes046/elder099/anaconda3/etc/profile.d/conda.sh
. /home/noyes046/elder099/anaconda3/etc/profile.d/mamba.sh

#Spades is in strains environment
mamba activate metawrap #megahit is in here

#Make assemblies directory
cd $rundir
mkdir -p assemblies


###MEGAHIT individual reads assembly
#for ((i=1;i<=10;i++))
#iterator
n=1

cat $indiv_reads | while read file
#tail $indiv_reads -n +839 | while read file
do
        echo "$file"
	
	#Only run assembly on samples that haven't been run
	assembly_path="$assemblies/spades_${file}_output"	
	if [ ! -d "$assembly_path" ]; then
	
	#Run assembly w/MEGAHIT
	megahit -m 0.3 -t 80 -1 $reads/$file.R1.fastq.gz -2 $reads/$file.R2.fastq.gz -o $assemblies/spades_${file}_output 
	#spades.py --meta -t 128 -1 $reads/$file.R1.fastq.gz -2 $reads/$file.R2.fastq.gz -o $assemblies/spades_${file}_output
	
	#Remove excess intermediates
	cd $assemblies/spades_${file}_output
	rm -fr intermediate_contigs
	cd $rundir

	#Read online that want to set memory limit high (e.g. 500GB) and fewer threads (e.g. 16)
      	#echo $i
      	
      	n=$(($n + 1))
	if [ $((n % 20)) == 0 ]; then
		echo $n
		echo "File $file is done assembling. This is file Number:$n" | mail -s "Assembly Progress" "elder099@umn.edu"

	fi
	

	fi
done


###METASPADES coReads coassembly
#for ((i=1;i<=10;i++))
#cat $coRead_list | while read file
#tail $coRead_list -n +2 | while read file
#do
#       # ngless map.ngl -j 32 /scratch.global/fermx014/data/elder099/Noyes_Project_019/NonHostFastq /scratch.global/elder099/StaphA/final_file_names.txt
#        echo "$file"
#	
#	spades.py --meta -t 32 -1 $coReads/$file.R1.fastq.gz -2 $coReads/$file.R2.fastq.gz -o $assemblies/spades_${file}_output
#	
#	#Read online that want to set memory limit high (e.g. 500GB) and fewer threads (e.g. 16)
#      # echo $i
#done

