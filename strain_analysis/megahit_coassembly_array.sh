#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --mem=200g
#SBATCH -t 18:00:00
#SBATCH --mail-type=END 
#

straindir='/scratch.global/elder099'
rundir='/scratch.global/elder099/strains_run_v5'
assemblies="$rundir/assemblies"
coassemblies="$assemblies/coassemblies"
#reads='/scratch.global/fermx014/help/elder099/Noyes_Project_019/NonHostFastq'
reads='/scratch.global/fermx014/data/elder099/elder099_2024-12-10_Noyes_Project_019/NonHostFastq'
#coReads="$reads/coReads"
coReads="$rundir/coReads"
indiv_reads="$rundir/indiv_read_list.txt"
#indiv_reads="$rundir/tmp_read_list.txt"
coRead_list="$rundir/coRead_list.txt"

# Activate conda
. /home/noyes046/elder099/anaconda3/etc/profile.d/conda.sh
. /home/noyes046/elder099/anaconda3/etc/profile.d/mamba.sh

#Spades is in strains environment
mamba activate metawrap #megahit is in here

cd $rundir
file=$(sed -n ${SLURM_ARRAY_TASK_ID}p $coRead_list)

cd $assemblies
mkdir -p coassemblies
#du -sh $coReads/${file}.R1.fastq.gz


#Only run assembly on samples that haven't been run
#assembly_path="$assemblies/spades_${file}_output"
#if [ ! -d "$assembly_path" ]; then

#	#Run assembly w/MEGAHIT
#	megahit -m 0.9 -t 20 -1 $reads/$file.R1.fastq.gz -2 $reads/$file.R2.fastq.gz -o $assemblies/spades_${file}_output
#	#spades.py --meta -t 128 -1 $reads/$file.R1.fastq.gz -2 $reads/$file.R2.fastq.gz -o $assemblies/spades_${file}_output
#
#	#Remove excess intermediates
#	cd $assemblies/spades_${file}_output
#	rm -fr intermediate_contigs
#	cd $rundir
#
#else
#	echo "File $file exists!"
#
#fi

###METASPADES coReads coassembly

#Only run assembly on samples that haven't been run
assembly_path="$coassemblies/spades_${file}_output"
if [ ! -d "$assembly_path" ]; then

       #Run assembly w/MEGAHIT
       megahit -m 0.9 -t 20 -1 $coReads/$file.R1.fastq.gz -2 $coReads/$file.R2.fastq.gz -o $coassemblies/spades_${file}_output

       #Remove excess intermediate files
	cd $coassemblies/spades_${file}_output
        rm -fr intermediate_contigs
else
       echo "File $file exists!"

fi
