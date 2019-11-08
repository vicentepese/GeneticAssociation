#!bin/bash

# Go to data
cd /scratch/users/vipese/GeneticAssociation/AA_GzipOutputs/

# Filter based on maximum and minimum position of the referred study 
echo "Filtering"
awk -v FS=, '{ if($1 >= 25997458 && $1 <= 33188223) print $0 | "gzip > CHR_6_filt.csv.gz"}' <(gzip -dc CHR_6.csv.gz) 
echo "Filtered"
