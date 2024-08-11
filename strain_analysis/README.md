# Strain Analysis Directory

*Steps for Strain Analysis*
- Assembly:
  - Use your favorite assembler on each sample individually
  - Optional: also perform coassembly by merging all reads from all samples (catches low-abundance microbes) N/A
- Genome Binning:
  - Bin each assembly separately
- Dereplication:
  - Use dRep on all bins together to extract similar strains across samples
- Downstream Analysis
  - Perform downstream analyses on de-repped list of genomes


# Captain's Log:


Finished assembly of USDA USDA871_S124.non.host reads. Next step is to test them out, see if they look good. 2024-07-19
 
 * Corrected reads are in /home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/strain_analysis/spades_output/corrected/
 * Assembled contigs are in /home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/strain_analysis/spades_output/contigs.fasta
 * Assembled scaffolds are in /home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/strain_analysis/spades_output/scaffolds.fasta
 * Paths in the assembly graph corresponding to the contigs are in /home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/strain_analysis/spades_output/contigs.paths
 * Paths in the assembly graph corresponding to the scaffolds are in /home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/strain_analysis/spades_output/scaffolds.paths
 * Assembly graph is in /home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/strain_analysis/spades_output/assembly_graph.fastg
 * Assembly graph in GFA format is in /home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/strain_analysis/spades_output/assembly_graph_with_scaffolds.gfa

======= SPAdes pipeline finished.

Completed QC of assemblies using QUAST. 2024-07-19

Finished binning of USDA USDA871_S124.non.host reads. Next step is to dereplicate all bins (each bin?). 2024-08-10
 * Bin output of MetaWrap (which ran Metabat2) is in: /home/noyes046/elder099/OREI_Shotgun_AMR_Analyses/strain_analysis/binning_round1/metabat2_bins/.fa
 * All bins are multi-fastas in one directory
 * MetaWrap doesn't work with zipped reads, may need to pick another tool later -- runtime around 25 mins for one samples at 50 cpus
 * dRep failed when running on one bin, and subsequently failed when using multiple
 * dRep in strains conda env

