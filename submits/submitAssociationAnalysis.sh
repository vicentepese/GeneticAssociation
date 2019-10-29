#!/bin/bash

#SBATCH --job-name=AssociationAnalysis
#SBATCH -p mignot
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=20G
#SBATCH --output=../SlurmOutputs/AssociationAnalysis/AssociationAnalysis__CHR%a.out
#SBATCH --error=../SlurmOutputs/AssociationAnalysis/AssociationAnalysis__CHR%a.err
#SBATCH --nodes=1
#SBATCH -t 05:00:00

# Load the module 
ml R/3.5.1
cd .. 
Rscript --vanilla associationAnalysis.R $SLURM_ARRAY_TASK_ID
