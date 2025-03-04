#!/bin/bash

#Establish file paths
rundir='/scratch.global/elder099/strains_run_v5'
reads='/scratch.global/fermx014/help/elder099/Noyes_Project_019/NonHostFastq'
assemblies="$rundir/assemblies"
coassemblies="$assemblies/coassemblies"
coReads="$rundir/coReads"
indiv_reads="$rundir/indiv_read_list.txt"
coRead_list="$rundir/coRead_list.txt"
binning="$rundir/binning"
work="$binning/work_files"


# Activate conda
. /home/noyes046/elder099/anaconda3/etc/profile.d/conda.sh
. /home/noyes046/elder099/anaconda3/etc/profile.d/mamba.sh

#Activate conda environment
mamba activate metawrap #Change this once conda is fixed


mkdir -p "$rundir/binning"


#####
#####BINNING INDIVIDUAL ASSEMBLIES
#####
cat $indiv_reads | while read file
#tail $indiv_reads -n +2 | while read file
do
        echo "$file"
	mkdir "$binning/work_files" #Create work directory
	
	
	assembly_path="$assemblies/spades_${file}_output"
	
	#Only run binning on existing assemblies
	if [ -d "$assembly_path" ]; then
	
		###Align original reads to generated scaffold
		echo "alignment"
		#Index scaffold file
		bwa index $assemblies/spades_${file}_output/final.contigs.fa
	
		#Make SAM of dehosted reads vs Scaffold
		bwa mem -t 80 $assemblies/spades_${file}_output/final.contigs.fa $reads/$file.R1.fastq.gz $reads/$file.R2.fastq.gz > $work/$file.sam
	
		#Make BAM from SAM
		samtools view -@ 80 -bS $work/$file.sam > $work/$file.bam
	
		echo "sorting"
		#Sort BAM file
		samtools sort -@ 80 $work/$file.bam -o $work/${file}_sorted.bam

		echo "depth file"
		#Grab depth file from BAM
		jgi_summarize_bam_contig_depths --outputDepth $work/depth.txt $work/${file}_sorted.bam



	
		###Run actual binning using metabat2
		
		echo "Bins directory: $binning/metabat2/${file}_bins"
        	mkdir -p $binning/metabat2/${file}_bins
		
		metabat2 -i $assemblies/spades_${file}_output/final.contigs.fa -a $work/depth.txt -o $binning/metabat2/${file}_bins/${file}_bin -t 80 -m 1500
	

	else
		#If not present, do nothing
  		echo "File does not exist!"
	
	
	fi
	

	###Complete binning script by deleting work files
	cd $binning
	rm -fr work_files

done



#####
#####BINNING COASSEMBLIES
#####

#cat $coRead_list | while read file
#tail $coRead_list -n +2 | while read file
#do
#        echo "$file"
#	mkdir "$binning/work_files" #Create work directory
	

	###Align original reads to generated scaffold
#	echo "alignment"
	#Index scaffold file
#	bwa index $assemblies/spades_${file}_output/scaffolds.fasta
	
	#Make SAM of dehosted reads vs Scaffold
#	bwa mem -t 80 $assemblies/spades_${file}_output/scaffolds.fasta $coReads/$file.R1.fastq.gz $coReads/$file.R2.fastq.gz > $work/$file.sam
	
	#Make BAM from SAM
#	samtools view -@ 50 -bS $work/$file.sam > $work/$file.bam
	
#	echo "sorting"
	#Sort BAM file
#	samtools sort -@ 50 $work/$file.bam -o $work/${file}_sorted.bam

#	echo "depth file"
	#Grab depth file from BAM
#	jgi_summarize_bam_contig_depths -t 50 --outputDepth $work/depth.txt $work/${file}_sorted.bam




	###Run actual binning using metabat2
#	echo "Binning"
#	echo "Bins directory: $binning/metabat2/${file}_bins"
#	mkdir -p $binning/metabat2/${file}_bins
	
#	metabat2 -i $assemblies/spades_${file}_output/scaffolds.fasta -a $work/depth.txt -o $binning/metabat2/${file}_bins/${file}_bin -t 80 -m 1500
	

	###Complete binning script by deleting work files
#	cd $binning
#	rm -fr work_files

#done

