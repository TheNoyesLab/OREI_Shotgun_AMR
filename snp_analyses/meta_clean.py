
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 10:07:43 2024
@author: Jared Young
"""

# metaSNV cleaner
## a tool to clean the output of metaSNV snp calling into a usable count matrix 

# load test data
headers = open('/scratch.global/elder099/Noyes_Project_019/all_pools/sample_headers.txt','rt')
data = open('/scratch.global/elder099/Noyes_Project_019/all_pools/Meta_Out/Meta_Out_Fix.txt','rt')
output = open('/scratch.global/elder099/Noyes_Project_019/all_pools/Final_Meta_SNPs.csv','wt')

# make lists to save data
names = list()
ref_pos = list()
ref_nuc = list()
snp_counts = list()

# loop to save data to set lists 
for line in data: 
    line = line.rstrip()
    line = line.split(sep = '\t')
    names.append(line[0])
    ref_pos.append(line[2])
    ref_nuc.append(line[3])
    snp_counts.append(line[5])
    
# extract SNP data 
# split by ann "."
nonref_nuc = list()
counts = list()
for line in snp_counts:
    line = line.rstrip()
    line = line.split(sep=".")
    nuc = line[0].split(sep="|")
    nonref_nuc.append(nuc[1])
    counts.append(line[1])

# extract headers 
heads = list()
heads.append("SNP_accession")
for line in headers:
    line = line.rstrip()
    heads.append(line)

# clean count data 
## make new list
clean_counts = list()
for entry in counts:
    entry = entry.rstrip()
    index = 0 
    result_string = ""
    for i in range(len(entry)):
        if i != index:
            result_string += entry[i]
    result_string = result_string.replace("|",",")
    clean_counts.append(result_string)
    
## make new row names 
## use enumerate(maybe)
meg = list()
compound = list()
rclass = list()
mechanism = list()
rgroup = list()

# Split name into principal components 
for entry in names: 
    entry = entry.rstrip()
    entry = entry.split(sep="|")
    meg.append(entry[0])
    compound.append(entry[1])
    rclass.append(entry[2])
    mechanism.append(entry[3])
    rgroup.append(entry[4])

# make new naming format 
new_name = list()
for (count,entry) in enumerate(meg):
    new_name.append(entry+"|"+compound[count]+"|"+rclass[count]+"|"+mechanism[count]+"|"+rgroup[count]+"|"+rgroup[count]+"_"+ref_nuc[count]+"_"+ref_pos[count]+"_"+nonref_nuc[count])

# print headers to output
for i,entry in enumerate(heads): 
    entry = entry.rstrip()
    if(i<len(heads)-1):
        print(entry,end=",",file=output)
    else:
        print(entry,file=output)
#print(entry,sep=",",file=output)
#print("",file=output)
#print(heads[0:],sep=",",file=output)

# print new name and counts to final output 
for (count,entry) in enumerate(new_name):
    print(entry,clean_counts[count],sep=",",file=output)

# close files
headers.close()
data.close()
output.close()
