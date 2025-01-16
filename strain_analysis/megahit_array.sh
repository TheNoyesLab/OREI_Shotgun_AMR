#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --mem=200g
#SBATCH -t 18:00:00
#SBATCH --mail-type=END 
#

#mins18:00:00
#database='/home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/staphA_reference/ncbi_dataset/data/GCF_000013425.1/StaphA_reference.fasta'
#reads='/scratch.global/fermx014/data/elder099/Noyes_Project_019/NonHostFastq'
straindir='/scratch.global/elder099'
rundir='/scratch.global/elder099/strains_run_v5'
#straindir='/home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/strain_analysis'
#rundir="$straindir/coassembly_run_v4"
assemblies="$rundir/assemblies"
#reads="$straindir/test_reads"
#reads='/scratch.global/fermx014/help/elder099/Noyes_Project_019/NonHostFastq'
reads='/scratch.global/fermx014/data/elder099/elder099_2024-12-10_Noyes_Project_019/NonHostFastq'
coReads="$reads/coReads"
indiv_reads="$rundir/indiv_read_list.txt"
#indiv_reads="$rundir/tmp_read_list.txt"
coRead_list="$rundir/coRead_list.txt"

# Activate conda
. /home/noyes046/elder099/anaconda3/etc/profile.d/conda.sh
. /home/noyes046/elder099/anaconda3/etc/profile.d/mamba.sh

#Spades is in strains environment
mamba activate metawrap #megahit is in here

cd $rundir
file=$(sed -n ${SLURM_ARRAY_TASK_ID}p $indiv_reads)

#du -sh $reads/$SAMPLE.R1.fastq.gz
#du -sh $reads/$SAMPLE.R2.fastq.gz


#Only run assembly on samples that haven't been run
assembly_path="$assemblies/spades_${file}_output"
if [ ! -d "$assembly_path" ]; then

	#Run assembly w/MEGAHIT
	megahit -m 0.9 -t 20 -1 $reads/$file.R1.fastq.gz -2 $reads/$file.R2.fastq.gz -o $assemblies/spades_${file}_output
	#spades.py --meta -t 128 -1 $reads/$file.R1.fastq.gz -2 $reads/$file.R2.fastq.gz -o $assemblies/spades_${file}_output

	#Remove excess intermediates
	cd $assemblies/spades_${file}_output
	rm -fr intermediate_contigs
	cd $rundir

else
	echo "File $file exists!"

fi


