# GeneticAssociation

## Introduction 

This repository was designed to perform an Association Analysis (AA) between Single Nucleotide Polymorphisms (SNP) and protein sequencing expression genes in peripheral blood cells. Firstly, files (SNP and protein files) are pre-processed in `dataParser.py` to be further analyzed in `associationAnalysis.R`.

This pipeline has been designed to function in cloud computer services &ndash; such as the one available at Stanford University: Sherlock. It allows the usage of the so-called Job Arrays, in which each job is designated to a single chromosome. In this way, tasks are efficiently distributed over _slurm_ jobs, both for pre-processing and analysis. If no argument is provided to `dataParser.py` or `associationAnalysis.R`, the scripts will be ran in the chromosomes specified in `options.json`

## Pipeline

### Options 

The pipeline is designed in a way in which no input must be given. For the sake of simplification, all the settings are stored in `options.json`. Such a file is structured in a tree fashion, containing the paths of the home directory, and paths to the folders where the protein, genomic, processed data and outputs are stored. Such paths can be modified as pleased and adapted to the needs. 

The option _chrArray_ allows to select the chromosomes that are to be analyzed, in a _slurm task array_ fashion -- \[_number_\] for individual analysis, \[_number1_-_number2_\] for a range of chromosomes, or _"all"_ for the analysis of all the chromosomes. 

Other complementary options **must** be included, such as the file indicating the ID for each patient *GWAS_ID_file*, the file of the covariates data *covariateFile* and the path to the log output *logPath*.

Multiple options have been included for both preprocessing (*genOptions*) and AA (*analysis*) where thresholds for Minimum Allele Frequency (MAF), HWE, number of chromosomes to be filtered out, and p-values are selected. 

### Pre-processing 

In order to run the data pre-processing or data parsing, `rundataParser.sh` must be executed, after the options/settings have been selected. 

Preprocessing consist on six major steps:
1. **Data removal (optional)**: in this step, unnecessary (unused) files are deleted from the data folder. 
2. **Protein data pre-processing**: protein data is loaded and pre-processed. Each patient ID from the protein sequencing file is matched with their corresponding ID from the GWAS, as indicated by the file in options *GWAS_ID_file*.
3. **Find patient IDs in the GWAS**: since the GWAS contains the information of more patients than in the protein data and in a different order, the indexes of the patients that appear in the protein data are taken as a variable. 
4. **Patient ID reordering**: to match the order of the patients of the GWAS ID file and the protein file, the rows (patients) of the protein file are reordered as in the GWAS (same procedure is applied to the covariate data). This step and the previous one has been performed for each chromosome to avoid inconsistencies and errors and according to the *.sample* files for each chromosome. 
5. **Process impute (SNP) data**: to avoid excessive computational power, the file is open in *write text* mode, and for each line, the probability of belonging to a specific allele is computed only for the indexes taken in the third step, and saved in a separate file. De novo, the computational power is reduced. 
6. **Save files**: files are seved in the designed folder in `options.json`. 

With the performance of the aforementioned steps, the pre-processing of the data is done, and the analysis can be subsequently performed. 

### Association Analysis 

The AA is performed in `associationAnalysis.R` in a similar fashion as the pre-processing --- each chromosome will be sent as an independent *slurm* job &ndash; unless no argument is given, in which case the AA will be computed in the chromosomes specified in `options.json`  

The protein data is normalazed using the R `scale` command, and further saved in a separate file. A linear model (using the library `MatrixEQTL`) is fit using the SNPS data (*impute*), the protein data and the covariates. The results are saved in the folder defined in `options.json`

## Reproducibility check 

In order to corroborate the outcome, the results are compared to a similar study &ndash; [Genetic variation in MHC proteins is associated with T cell receptor expression biases](https://www.nature.com/articles/ng.3625#supplementary-information) by E.Sharon et al.. After downloading the data from the study,the chromosome number is computed based on the rsID for each gene. The dataset to compute the chromosome from the rsID can be found in [this link](https://atgu.mgh.harvard.edu/plinkseq/resources.shtml) named *refdb.dbsnp.gz*. Then, for the chromosomes in which the aforementioned studies found significance, matches between SNP (chromosome number and location)-protein are found. 

In this case, the aforementioned study only found significance association in chromosome 6, thus the reproducibility is solely checked in such chromosome. Our results replicate ~7% of the significant associations in the aforementioned study. Their results are fully available online. The data used to check the repdoducility was taken from the Supplementary Table 4 of the study. 

