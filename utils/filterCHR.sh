#!bin/bash

# Go to data
cd /scratch/users/vipese/GeneticAssociation/AA_GzipOutputs/

# Filter based on maximum and minimum position of the referred study 
echo "Filtering lower MAF"
awk -v FS=, '{ if($2 >= 25997458) print $0}' CHR_6.csv > CHR6_filt
awk -v FS=, '{ if($2 <= 33188223) print $0}' CHR6_filt > CHR6_filt_def
awk -v FS=, '{print$0}' CHR_6.csv
