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

