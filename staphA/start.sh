#!/bin/bash

database='/home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/staphA_reference/ncbi_dataset/data/GCF_000013425.1/StaphA_reference.fasta'
reads='/scratch.global/fermx014/data/elder099/Noyes_Project_019/NonHostFastq'
filenames='/scratch.global/elder099/StaphA/final_file_names.txt'
align='/scratch.global/elder099/StaphA/bwa_output'

#for ((i=1;i<=10;i++))
cat $filenames | while read file
do
       # ngless map.ngl -j 32 /scratch.global/fermx014/data/elder099/Noyes_Project_019/NonHostFastq /scratch.global/elder099/StaphA/final_file_names.txt
	echo "$file"	
       
	bwa mem -t 64 -R "@RG\tID:$file\tSM:2Strain\tPL:Illumina" $database \
                      $reads/$file.R1.fastq.gz $reads/$file.R2.fastq.gz |
		      samtools view --threads 16 -Sb > $align/${file}_aligned.bam

                #Alignment post-processing
                #samtools view -S -b $Align/Sim${i}_Ref.sam > $Align/Sim${i}_Ref.bam --threads 16
                samtools sort $align/${file}_aligned.bam -o $align/${file}_sorted.bam --threads 16
                samtools index $align/${file}_sorted.bam
                #bedtools bamtobed -i $Align/Sim${i}_Ref_sorted.bam > $Align/Sim${i}_Ref_sorted.bed

 
      # echo $i
done



