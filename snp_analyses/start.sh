#!/bin/bash

for ((i=1;i<=598;i++))
do
	ngless map.ngl -j 32 /scratch.global/elder099/Noyes_Project_019/all_pools/NonHostFastq /scratch.global/elder099/Noyes_Project_019/all_pools/final_file_names.txt
	echo $i
done




