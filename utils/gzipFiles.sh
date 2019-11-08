#!/bin/bash

#SBATCH --job-name=Gzip%a
#SBATCH -p owners
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=20G
#SBATCH --nodes=1
#SBATCH -t 05:00:00
#SBATCH --array=1-22

cd /scratch/users/vipese/GeneticAssociation/Outputs/
pigz CHR_$SLURM_ARRAY_TASK_ID.csv CHR_$SLURM_ARRAY_TASK_ID.csv.gz
mv CHR_$SLURM_ARRAY_TASK_ID.csv.gz ..//AA_GzipOutputs/CHR_$SLURM_ARRAY_TASK_ID.csv.gz