import json
import os
import csv
import pandas as pd
import numpy as np
import gzip
import sys
from collections import OrderedDict
import sys


def loadReproData(filename):

    # Create xl file
    xlFile = pd.ExcelFile(filename)

    # Create Ordered dictionnary
    pre_repro_data = {
        sheet_name: xlFile.parse(sheet_name)
        for sheet_name in xlFile.sheet_names
    }

    pre_repro_data = pre_repro_data['MHC long range associations']
    repro_data = OrderedDict()
    for val in pre_repro_data.values:
        repro_data[val[0]] = val[1:]

    return repro_data

def GzipFileHandler(FileName, Read=True):
    if '.gz' in FileName:
        if Read is True:
            OpenFile = gzip.open(FileName, 'rt', newline = '')

        else:
            OpenFile = gzip.open(FileName, 'wb')
    else:
        if Read is True:
            OpenFile = open(FileName, 'r')

        else:
            OpenFile = open(FileName, 'w')

    return OpenFile

def loadCHRData(filename):

    # Open file in text mode
    OpenFile = GzipFileHandler(filename, Read=True)
    print("Reading association analysis outputs. This might take several minutes.")
    reader = csv.reader(OpenFile)

    # Load file 
    CHR_data = OrderedDict()
    for row in reader:
        if row.index('') == 0:
            colnames = row[1:]
        else:
            



    return reader


def getMatches(CHR_data, repro_data):

    # Get genes
    snps_pos = CHR_data['snps']

    # Check if snps is among the significant results
    for snp in snps_pos:
        if snp in list(repro_data.keys()):

            # Get index 
            idx = list(repro_data.keys()).index(snp)

            # Check if 


def main():

    # Load options
    with open("optionsCheck.json", 'r') as jsonFile:
        options = json.load(jsonFile)

    # Load "reproducibility data"
    repro_data = loadReproData(options['file']["checkRepro"])

    # Load association analysis data
    CHR_data = loadCHRData('CHR_14.csv.gz')

    # Get matches
    getMatches(CHR_data, repro_data)

    # Load chromosome data
    print('stop')


main()