#!/bin/bash

#SBATCH --job-name=gzipFiles
#SBATCH -p owners
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH -t 10:00:00

# Go to folder 
cd /scratch/users/vipese/GeneticAssociation/Outputs/

# Get files 
files="$(ls)"

# Gzip each file 
for file in $files ; do
  echo "Gziping $file"
  gzip $file
  echo "$file succesfully gziped"
  echo "Moving $file to AAO folder"
  mv "$file.gz" ../Association_Analysis_Outputs
  echo "$file succesfully moved"
done 