#!/bin/bash

#Establish file paths
rundir='/scratch.global/elder099/strains_run_v5'
reads='/scratch.global/fermx014/help/elder099/Noyes_Project_019/NonHostFastq'
coReads="$rundir/coReads"
assemblies="$rundir/assemblies"
coassemblies="$assemblies/coassemblies"
indiv_reads="$rundir/indiv_read_list.txt"
coRead_list="$rundir/coRead_list.txt"
all_read_list="$rundir/final_read_list.txt"
binning="$rundir/binning"
work="$binning/work_files"


# Activate conda
. /home/noyes046/elder099/anaconda3/etc/profile.d/conda.sh
. /home/noyes046/elder099/anaconda3/etc/profile.d/mamba.sh

#Activate conda environment
mamba activate strains #dRep is in here

#Create list of bins for input to dRep

cd $binning
rm -f bin_list.txt

#Loop through all of the read names, including coReads/coBins
cat $indiv_reads | while read file
#tail $all_read_list -n +2 | while read file
do
	echo $file

        bin_path="$binning/metabat2/${file}_bins"
	
        #Only run binning on existing assemblies
        if [ -d "$bin_path" ]; then
	
		#No longer need to update names -- fixed upstream
		cd $binning/metabat2/${file}_bins

		#Pull full file path of all successful bins
		ls -d -1 "$PWD/"*[!unbinned].fa* >> $binning/bin_list.txt
	
	else
                #If not present, do nothing
                echo "File does not exist!"
        fi


done

#Enter directory of bins
#cd $straindir/binning_898_996/USDA898_S127.non.host_bins/metabat2_bins
#Change bin names to add Sample ID so they are unique
#for f in * ; do mv -- "$f" "USDA898_$f" ; done
#Export list of bins paths to file
#ls -d -1 "$PWD/"*[!unbinned].fa* > $straindir/binning_898_996/bin_list.txt

#Enter directory of bins
#cd $straindir/binning_898_996/USDA996_S67.non.host_bins/metabat2_bins
#Chagne bin names to add Sample ID so they are unique
#for f in * ; do mv -- "$f" "USDA996_$f" ; done
#Export list of bins paths to file, overwriting
#ls -d -1 "$PWD/"*[!unbinned].fa* >> $straindir/binning_898_996/bin_list.txt

###Dereplicate Bins 
#primary clustering 95%, secondary 99%
dRep dereplicate -p 80 $binning/dRep_output -g "$binning/bin_list.txt" -pa 0.95 -sa 0.99



