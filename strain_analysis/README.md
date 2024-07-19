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
