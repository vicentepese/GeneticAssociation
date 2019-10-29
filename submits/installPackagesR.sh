#!/bin/bash

#SBATCH --job-name=installPackagesR
#SBATCH -p owners
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --output=../SlurmOutputs/installPackagesR.out
#SBATCH --error=../SlurmOutputs/installPackagesR.err
#SBATCH --nodes=1
#SBATCH -t 01:00:00

# Load the module 
ml R/3.5.1
R --no-save << EOF
check.packages <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
        install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
}

# Usage example
packages<-c("optparse")
check.packages(packages)
EOF
