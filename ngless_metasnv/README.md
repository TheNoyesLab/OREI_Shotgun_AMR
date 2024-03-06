# Scripts to Run NGLess Alignment and MetaSNV Variant Calling

## Step 1: filepath_create.sh
- Create list of file names using sample_ids for input to NGLess
## Step 2: map.ngl 
- Script in pythonic NGLess language. Peforms BWA alignment
## Step 3: start.sh
- Run NGLess. Must be run in a loop to align all samples, despite parallelization
## Step 4: metasnv.sh
- Runs file setup then runs MetaSNV variant calling on all BAM files (based on mpileup)
## Step 5: meta_clean.py
- Processes MetaSNV outputs into a a wide-format count matrix of SNPs


*Not all samples run yet, working on subset*
