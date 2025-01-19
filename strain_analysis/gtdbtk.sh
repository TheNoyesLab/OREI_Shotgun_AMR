#!/bin/bash

#Establish file paths

straindir='/scratch.global/elder099'
rundir='/scratch.global/elder099/strains_run_v5'
reads='/scratch.global/fermx014/help/elder099/Noyes_Project_019/NonHostFastq'
coReads="$rundir/coReads"
assemblies="$rundir/assemblies"
coassemblies="$assemblies/coassemblies"
indiv_reads="$rundir/indiv_read_list.txt"
coRead_list="$rundir/coRead_list.txt"
all_read_list="$rundir/final_read_list.txt"
binning="$rundir/binning"
final_bins="$binning/dRep_output/dereplicated_genomes"


# Activate conda
. /home/noyes046/elder099/anaconda3/etc/profile.d/conda.sh
. /home/noyes046/elder099/anaconda3/etc/profile.d/mamba.sh

#Activate conda environment
mamba activate mag_taxon #gtdb-tk is here

gtdbtk classify_wf --genome_dir $final_bins --out_dir $rundir/gtdb-tk_output --extension fa --cpus 60 --mash_db $rundir/gtdb-tk_output/mash_db
