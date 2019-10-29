#!/usr/bin/env Rscript

library(MatrixEQTL)
library(rjson)
library(readr)
library(ggplot2)
library(MASS) 
library(car)
library(gridExtra)

# Set working directory 
setwd("~/Projects/GeneticAssociation")

############# Set arguments ##############

args = commandArgs(trailingOnly=TRUE)

if (length(args) == 0 || args[1] == 'all'){
    chrArray = c(1:22)
} else{
    chrArray = list()
    for (arg in args){
      chrArray = c(chrArray, arg)
    }
}
CHR <- chrArray[[1]]
print(paste("The Association analysis will be run in chromosome ", as.character(CHR), sep = ''))

############# Initialization ###############

# Load options
options <- fromJSON(file = 'options.json')

# Get list of files 
processedFiles = list.files(options$path$processedFolder)

# Load prot data
chrProtFile = processedFiles[grep(paste('CHR',CHR,'_', sep = ''),processedFiles)]
chrProtFile = chrProtFile[grep('protData', chrProtFile)]
print(paste(" Loading ", chrProtFile, sep = ''))
prot.data <- read_csv(file = paste(options$path$processedFolder, chrProtFile, sep = ''), )
print(head(prot.data))
rownames(prot.data) <- prot.data$Protein
prot.data$Protein <- NULL


# Normalize data 
print("Normalizing data")
norm.prot.data <- scale(prot.data)

# Write if in options
if (options$analysis$saveNorm == TRUE){
    print(CHR)
    write.csv(norm.prot.data, file = paste(options$path$processedFolder, "CHR",as.character(CHR), "_", "normProtData.csv", sep = ''))
}

################## Analysis ###############

# Create model 
useModel <- MatrixEQTL::modelLINEAR

# Create outputfile
output_filename <- tempfile()

# P-value
pvalthres <- options$analysis$pvalthres

# Set error covariance matrix
errorCovariance <- numeric()

# Load genotype data (cols = patients, rows = snps)
chrSnpFile = processedFiles[grep('impute', processedFiles)][grep(paste('CHR',CHR,'_',sep=''), processedFiles[grep('impute', processedFiles)])]
snps = SlicedData$new();
snps$fileDelimiter = ",";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 0;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
try(
    snps$LoadFile(paste(options$path$processedFolder,chrSnpFile , sep = ''), delimiter = '')
)


# Remove based on MAF
mafthreshold <- options$analysis$maf 
maf.list = vector('list', length(snps))
for(sl in 1:length(snps)) {
slice = snps[[sl]];
maf.list[[sl]] = rowMeans(slice,na.rm=TRUE)/2;
maf.list[[sl]] = pmin(maf.list[[sl]],1-maf.list[[sl]]);
}
maf = unlist(maf.list)
cat('SNPs before filtering:',nrow(snps), "\n")
snps$RowReorder(maf>=mafthreshold);
cat('SNPs after filtering:',nrow(snps), "\n")


# Load gene expression data (cols = patients, rows = genes)
chrProtNormFile =  processedFiles[grep('normProt', processedFiles)][grep(paste('CHR',CHR,'_',sep=''), processedFiles[grep('normProt', processedFiles)])]
prot = SlicedData$new();
prot$fileDelimiter = ",";      # the TAB character
prot$fileOmitCharacters = "NA"; # denote missing values;
prot$fileSkipRows = 1;          # one row of column labels
prot$fileSkipColumns = 1;       # one column of row labels
prot$fileSliceSize = 2000;      # read file in slices of 2,000 rows
prot$LoadFile(paste(options$path$processedFolder, chrProtNormFile, sep = ''));


# Load covariate (cols = patients, rows = covariates)
chrCovarFile = processedFiles[grep('covar', processedFiles)][grep(paste('CHR',CHR,'_',sep=''), processedFiles[grep('covar', processedFiles)])]
cvrt = SlicedData$new();
cvrt$fileDelimiter = ",";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 0;       # one column of row labels
cvrt$fileSliceSize = 2000;      # read file in slices of 2,000 rows
cvrt$LoadFile(paste(options$path$processedFolder, chrCovarFile, sep = ''));

# Run analysis 
me = Matrix_eQTL_engine(
    snps = snps,
    gene = prot,
    cvrt = cvrt,
    output_file_name = output_filename,
    pvOutputThreshold = pvalthres,
    useModel = useModel, 
    errorCovariance = errorCovariance, 
    verbose = TRUE,
    pvalue.hist = TRUE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE
    );
print("Analysis finished")
print("Unlinking output")
unlink(output_filename);

## Results:
cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected eQTLs:', '\n');
head(me$all$eqtls)

# Write 
print("Saving file")
output_filename <- paste(options$path$outputs, "CHR_", as.character(CHR),".csv", sep = '')
write.csv(me$all$eqtls, output_filename)
print(cat("File saved in", output_filename ))





