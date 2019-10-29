"""
scripts to clean an oxford gen file using thresholds from snptest summary file
@author: Aditya Ambati ambati@stanford.edu, Mignot Lab, Stanford University
"""
import gzip
import logging



def ProcessStat(Stat, mafthreshold, hwethreshold, infothreshold):
	logging.info('STAT FILE {} WITH THRESHOLDS AS FOLLOWS MAF {} HWE {} INFOSCORE {}'.format(Stat, mafthreshold, hwethreshold, infothreshold))
	StatFile = GzipFileHandler(Stat)
	GoodSnps =0
	GoodIndex =[]
	for n, line in enumerate(StatFile):
		if n > 10:
			LineParse = line.strip().split(' ')
			if len(LineParse) == 37:
				info, all_maf, index = map(float, [LineParse[8], LineParse[30], LineParse[6]])
				hwe = set([True if i > hwethreshold else False for i in map(float, LineParse[32:36])])
				if not False in hwe and info >= infothreshold and all_maf >= mafthreshold:
					GoodSnps += 1
					GoodIndex.append(index)
		# elif n == 10:
		# 	Header = line.strip()
		# 	Header = {i:j for i,j in enumerate(Header.split(' '))}
	StatFile.close()
	prop='{0:.1f}'.format((float(GoodSnps)/float(n))*100)
	logging.info('FOUND QCED SNPS {} OUT OF {} IN {} FILE {} % PERCENT '.format(GoodSnps, n,  Stat, prop))
	return GoodIndex


def ConvertGenDos(Genos):
	SampleN = len(Genos)
	SampleCheck = SampleN % 3
	assert SampleCheck == 0
	Sample = 0
	Dosages = ''
	for i in xrange(0, len(Genos), 3):
		AA, AB, BB=map(float, Genos[i:i+3])
		DosageInd=(0*AA)+AB+(2*BB)
		#assert DosageInd <= 2.0
		Sample += 1
		Dosages += str(DosageInd)+' '
	return Dosages


def ProcessGenFile(GenFile, Stat, OutFile,Dosage, Exclude):
	GenIn = GzipFileHandler(GenFile)
	OutGen = GzipFileHandler(OutFile, Read=False)
	logging.info('READING THE GENOTYPE FILE {}'.format(GenFile))
	if Exclude == 1:
		Indices=ProcessStat(Stat=Stat,mafthreshold=0.01,hwethreshold=5e-6,infothreshold=0.7)
	Lines = 0
	IndexLine =0
	LineBuffer =0
	for line in GenIn:
		Lines += 1
		IndexLine += 1
		if Lines == 100000:
			LineBuffer += 1
			Lines =0
			logging.info('PROCESSED {} MILLION VARIANTS'.format(LineBuffer))
		if Exclude == 1:
			if IndexLine in Indices: # if the line is present in the indices then it is a good snp having passed all the thresholds
				GenLine = line.strip().split(' ')
				RsidHeader = GenLine[2]#' '.join(GenLine[:6])
				if Dosage == 1:
					Genos = ConvertGenDos(GenLine[6:])
					ParsedLine = RsidHeader +' '+Genos.strip()
					OutGen.write(ParsedLine+'\n')
				else:
					OutGen.write(line)
		else:
			GenLine = line.strip().split(' ')
			RsidHeader = GenLine[2]#' '.join(GenLine[:6])
			if Dosage == 1:
				Genos = ConvertGenDos(GenLine[6:])
				ParsedLine = RsidHeader +' '+Genos.strip()
				OutGen.write(ParsedLine+'\n')
			else:
				OutGen.write(line)


	OutGen.close()
	GenIn.close()
	logging.info('PROCESSED IN TOTAL {} FROM {} AND DUMPED THE VARIANTS INTO {}'.format(IndexLine, GenFile, OutFile))



def main():
	import argparse
	parser = argparse.ArgumentParser(description='A script that cleans the  oxforf gen file')
	parser.add_argument('-StatFile', help='A summary stat file output from snptest summary stat', required=True)
	parser.add_argument('-Genotype', help='The Gen file to be cleaned', required=True)
	parser.add_argument('-Out', help='Name with path to which the cleaned Genos should be written', required=True)
	parser.add_argument('-Log', help='Name with path to which the log should be written', required=False)
	parser.add_argument('-Convert', choices=['1', '0'], help='Should the imputed best guess genotypes be converted to probabilities ranging from 0 to 2', required=False)
	parser.add_argument('-Exclude', choices=['1', '0'], help='Should the variants be excluded', required=False)
	args=parser.parse_args()
	GenFile = args.Genotype
	Stat = args.StatFile
	OutFile = args.Out
	Log = args.Log
	Dosage = args.Convert
	Exclude = args.Exclude
	if Log is not None:
		logging.basicConfig(filename=Log,level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
		logging.info('ARGUMENTS GIVEN {}'.format(args))
		logging.info('STARTED PROCESSING {} '.format(Stat))
	else:
		print 'LOGGING HAS BEEN TURNED OFF'
	if Dosage is not None:
		logging.info('CONVERT TO DOSAGES IS {} '.format(Dosage))
		Dosage = int(Dosage)
	else:
		Dosage =1
	
	if Exclude is not None:
		logging.info('EXCLUSION OF VARIANTS IS TURNED {} '.format(Exclude))
		Exclude = int(Exclude)
	else:
		Exclude = 1
		
	ProcessGenFile(GenFile=GenFile, Stat=Stat,OutFile=OutFile,Dosage=Dosage, Exclude=Exclude)


if __name__ == '__main__':main()

# import logging
# Stat='CHR1_QTL.stat.gz'
# logging.basicConfig(filename='test.log',level=logging.DEBUG, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
# logging.info('STARTED PROCESSING {} '.format(Stat))
# ProcessGenFile(GenFile='CHR1_test.gen.gz', Stat='CHR1_QTL.stat.gz',Dosage=True,Subset=False, OutFile='testcleaned.gen.gz')


