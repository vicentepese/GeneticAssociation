import json
import os
import csv
import pandas as pd
import numpy as np
import argparse
import gzip
import logging
import sys
from collections import OrderedDict
import sys

def argsParser():

	# Define argument parser
	parser = argparse.ArgumentParser(description='A script that cleans the  oxforf gen file')

	# Add arguments
	parser.add_argument('-ChrIndex', help= 'int: number of the chromosome to be parsed. Only for slurm', required=False)

	# Initialize variables
	args = parser.parse_args()

	# If no CHR specified, analyze as specified in options
	if not args.ChrIndex:
		# Load options
		with open("options.json", "r") as jsonFile:
			options = json.load(jsonFile)
		if isinstance(options['chrArray']['CHRpreprocess'], int):
			chrArray = ["CHR" + str(options)]
		else: 
			chrArray = ["CHR" + str(s) for s in options['chrArray']['CHRpreprocess']]

	# Else, analyze the inputed chromosome (slurm)
	elif len(args.ChrIndex) > 1:
		chrArray = ["CHR" + str(args.ChrIndex)]
		print(" Len > 1" + chrArray)
	elif len(args.ChrIndex) == 1:
		print("Len = 1 " + chrArray)
		chrArray = ["CHR" + str(args.ChrIndex)]

	return chrArray


def importSample(options, sampleFile):

	with open(options['path']['genomFolder'] + sampleFile[0], 'r') as csv_file:
		csv_reader = csv.reader(csv_file, delimiter = ' ')
		samples = list()
		for row in csv_reader:
			if  csv_reader.line_num is 1 or csv_reader.line_num is 2:
				continue
			samples.append(row[1].split('_')[1]+row[1].split('_')[2])

	return samples

def importProtData(options):

	protData = OrderedDict()
	with open(options['path']['protFolder'] + "MergedProtData.csv", 'r') as csv_file:
		csv_reader = csv.reader(csv_file, delimiter = ',')
		for row in csv_reader:
			if csv_reader.line_num is 1:
				proteins = row[2:]
			else:
				protData[row[1].split(' ')[0] + row[1].split(' ')[1]] = row[2:]

	return protData, proteins

def importGWASID(options):

	GWAS_IDs = OrderedDict()
	with open(options['path']['protFolder'] + options['GWAS_ID_file']) as csvFile:
		csv_reader = csv.reader(csvFile)
		for row in csv_reader:
			GWAS_IDs[row[0].split(' ')[-1]] = row[0].split(' ')[1]

	return GWAS_IDs

def importCovar(options):

	coverData = OrderedDict()
	with open(options['path']['genomFolder'] + options['covariateFile']) as csvFile:
		csv_reader = csv.reader(csvFile)
		for row in csv_reader:
			if csv_reader.line_num is 1:
				continue
			coverData[row[0].split()[0].split('_')[1] + row[0].split()[0].split('_')[2]] = row[0].split()[3:]

	return coverData

def removeData(options):

	# Get list of files
	genomFiles = os.listdir(options["path"]["genomFolder"])

	# Check if files have not been filtered -- removed
	print("Checking if files were processed")
	if not any("slurm" in genomFile for genomFile in genomFiles):
		print("Files were processed")
		return

	# If files were not filtered, remove files that do not start with "CHR"
	else:
		print("Files not processed, removing unnecessary files")
		notCHRfiles = list(filter(lambda s: not (s.startswith("CHR")) , genomFiles))
		notCHRfiles = list(filter(lambda s: not(s.endswith(".sample")), notCHRfiles))

		for notCHRfile in notCHRfiles:
			os.remove(options['path']['genomFolder'] + notCHRfile)

		# Get list of files
		genomFiles = os.listdir(options['path']["genomFolder"])

		# Remove files that contian Shapeit and do not end with sample
		genomFiles = list(filter(lambda s: "Shapeit" in s or "IGHREGION" in s, genomFiles))
		for gfile in genomFiles:
			os.remove(options['path']['genomFolder'] + gfile)

		# Remove IGHREGION GENOTIPES

		print("Files have been processed")

def writeOrderedDict(data, filename):

	# Get keys (columns)
	keys = list(data.keys())

	# Get length
	length = len(data[keys[0]])

	#Write file
	i = 0
	with open(filename, 'w') as outfile:
		csv_writer = csv.writer(outfile)
		csv_writer.writerow(keys)
		while i < length:
			row = list()
			for key in keys:
				row.append(data[key][i])
			csv_writer.writerow(row)
			i += 1


def processProtData(options, protData):

	# Load GWAS ID
	GWAS_ID = importGWASID(options)

	# Find ID match
	processedProtData = OrderedDict()
	for id in list(protData.keys()):
		if id in list(GWAS_ID.keys()):
			processedProtData[GWAS_ID[id]] = protData[id]

	return(processedProtData)

def getIndexes(options, processedProdData, samples):

	# Get create True-False vector
	indexes = list()
	for sample in samples:
		if sample in list(processedProdData.keys()):
			indexes.append(True)
		else:
			indexes.append(False)

	# Get indxs
	indexes = np.where(np.array(indexes) == True)[0]

	# Array of IDS
	IDs = np.array(samples)[(indexes)]

	# Return
	return indexes, IDs

def reorderData(indexes, IDs, data):

	# Reorder
	reorData = OrderedDict()
	for id in IDs:
		reorData[id] = data[id]

	return reorData

def GzipFileHandler(FileName, Read=True):

	# If file is in gz format open with gzip
	if '.gz' in FileName:
		if Read is True:
			OpenFile = gzip.open(FileName, 'rb')

		else:
			OpenFile = gzip.open(FileName, 'wb')
	# Else open in text mode
	else:
		if Read is True:
			OpenFile = open(FileName, 'r')

		else:
			OpenFile = open(FileName, 'w')

	return OpenFile

def ConvertGenDos(Genos, indexes):

	# Convert the allele into a probability
	SampleN = len(Genos)
	SampleCheck = SampleN % 3
	assert SampleCheck == 0
	Sample = 0
	Dosages = ''
	for idx in indexes:
		AA, AB, BB=map(float, Genos[idx*3:idx*3+3])
		DosageInd=(0*AA)+AB+(2*BB)
		Sample += 1
		Dosages += str(DosageInd)+' '

	return Dosages


def ProcessStat(Stat, mafthreshold, hwethreshold, infothreshold):

	# Verbose
	logging.info('STAT FILE {} WITH THRESHOLDS AS FOLLOWS MAF {} HWE {} INFOSCORE {}'.format(Stat, mafthreshold, hwethreshold, infothreshold))

	# Initialize
	StatFile = GzipFileHandler(Stat)
	GoodSnps =0
	GoodIndex =[]

	# For gene file in the impute data
	for n, line in enumerate(StatFile):
		if n > 10:
			LineParse = line.strip().split(' ')
			if len(LineParse) == 37:
				info, all_maf, index = map(float, [LineParse[8], LineParse[30], LineParse[6]])
				hwe = set([True if i > hwethreshold else False for i in map(float, LineParse[32:36])])
				if not False in hwe and info >= infothreshold and all_maf >= mafthreshold:
					GoodSnps += 1
					GoodIndex.append(index)
	StatFile.close()
	prop='{0:.1f}'.format((float(GoodSnps)/float(n))*100)
	logging.info('FOUND QCED SNPS {} OUT OF {} IN {} FILE {} % PERCENT '.format(GoodSnps, n,  Stat, prop))

	return GoodIndex


def ProcessGenFile(GenFile, Stat, OutFile, Dosage, Exclude, chr, indexes, options):

	# Handle input if in gzip format, and output
	GenIn = GzipFileHandler(GenFile)
	if OutFile is None:
		OutGen = GzipFileHandler(options['path']['processedFolder'] + chr + "impute", Read=False)
	else:
		OutGen = GzipFileHandler(OutFile, Read=False)
	logging.info('READING THE GENOTYPE FILE {}'.format(GenFile))

	# Eclusion
	if Exclude == 1:
		Indices=ProcessStat(Stat=Stat,mafthreshold=options['genOptions']['mafthreshold'], hwethreshold=options['genOptions']['hwethreshold'], infothreshold=options['genOptions']['infothreshold'])

	# Initialize
	Lines = 0
	IndexLine =0
	LineBuffer =0

	# For each gene (line)
	for line in GenIn:
		Lines += 1
		IndexLine += 1

		# Verbose
		if Lines == 100000:
			LineBuffer += 1
			Lines =0
			logging.info('PROCESSED {} MILLION VARIANTS'.format(LineBuffer))

		# If exclude critiria 
		if Exclude == 1:
			if IndexLine in Indices: # if the line is present in the indices then it is a good snp having passed all the thresholds
				GenLine = line.strip().split(' ')
				RsidHeader = GenLine[2]
				if Dosage == 1:
					Genos = ConvertGenDos(GenLine[6:], Indices)
					ParsedLine = RsidHeader +' '+Genos.strip()
					OutGen.write(ParsedLine+'\n')
				else:
					OutGen.write(line)
					
		# If no exclusion, take all SNPS
		else:
			GenLine = str(line).split(' ')
			posLine = GenLine[2]
			if Dosage == 1:
				Genos = ConvertGenDos(GenLine[5:], indexes)
				ParsedLine = posLine + ' ' + Genos.strip()
				OutGen.write(ParsedLine+'\n')
			else:
				OutGen.write(line)

	# Close
	OutGen.close()
	GenIn.close()
	logging.info('PROCESSED IN TOTAL {} FROM {} AND DUMPED THE VARIANTS INTO {}'.format(IndexLine, GenFile, OutFile))

def main(chrArray):

	# Load options
	with open("options.json", "r") as jsonFile:
		options = json.load(jsonFile)

	# TODO: recheck how to do this
	# Remove unnecesary data
	removeData(options)

	# List files
	genomFiles = os.listdir(options['path']['genomFolder'])

	# Import protein files
	protData, proteins = importProtData(options)
	print('Protein data imported')

	# Import covariates
	covarData = importCovar(options)
	print('Covariate data imported')

	# Process protein data
	processedProtData = processProtData(options, protData)
	print('Protein data processed')

	for chr in chrArray:

		# Import sample
		sampleFile =  [s for s in genomFiles if "sample" in s and chr+'.' in s]
		samples = importSample(options, sampleFile)
		print('Sample file of %s imported' % (chr))

		# Match samples
		indexes,IDs = getIndexes(options, processedProtData, samples)

		# Reorder protein count file and write
		protDataReord = reorderData(indexes, IDs, processedProtData)
		protDataReord.update({"Protein": proteins})
		protDataReord.move_to_end("Protein", last = False)
		writeOrderedDict(protDataReord, options['path']['processedFolder'] + chr + '_' + 'protData.csv')
		print('Protein data reordered and saved')

		# Reorder covariates data and write
		covarDataReord = reorderData(indexes,IDs, covarData)
		writeOrderedDict(covarDataReord, options['path']['processedFolder'] + chr + '_' + 'covarData.csv')
		print('Covariate data reordered and saved')

		# Get impute file name
		genFile = [s for s in genomFiles if "impute" in s and chr+'_' in s][0]
		genPath = options['path']['genomFolder'] + genFile
		statFile = [s for s in genomFiles if "summary" in s and chr+'_' in s][0]
		statPath = options['path']['genomFolder'] + statFile

		# Initialize
		outFile = options['path']['processedFolder'] + chr + '_' + 'imputeData'
		Dosage = options['genOptions']['Dosage']
		Exclude = options['genOptions']['Exclude']
		Log = options['logPath'] + str(chr)
		logging.basicConfig(filename=Log,level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

		# Process genome file
		ProcessGenFile(genPath, statPath, outFile, Dosage, Exclude, chr, indexes, options)

if __name__ == "__main__":
	chrArray = argsParser()
	print(chrArray)
	main(chrArray)
