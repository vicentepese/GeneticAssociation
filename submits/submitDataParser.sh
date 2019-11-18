#!/bin/bash

#SBATCH --job-name=dataParser
#SBATCH -p mignot
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --output=../SlurmOutputs/DataParser/dataParser_CHR%a.out
#SBATCH --error=../SlurmOutputs/DataParser/dataParser_CHR%a.err
#SBATCH --nodes=1
#SBATCH -t 03:00:00


cd ..
python3 dataParser.py -ChrIndex $SLURM_ARRAY_TASK_ID

