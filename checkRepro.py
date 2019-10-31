import json
import os
import csv
import pandas as pd
import numpy as np
import gzip
import sys
from collections import OrderedDict
import sys
import sqlite3


def loadReproData(filename):

    # Create xl file
    xlFile = pd.ExcelFile(filename)

    # Create Ordered dictionnary
    pre_repro_data = {
        sheet_name: xlFile.parse(sheet_name)
        for sheet_name in xlFile.sheet_names
    }

    pre_repro_data = pre_repro_data['MHC long range associations']
    repro_data = dict()
    for val in pre_repro_data.values:
        app = [val[0]]
        app.extend(val[2:])
        repro_data[val[1]] = app

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
    colnames = ["position", "gene", "statistic", "pvalue", "FDR", "beta"]
    next(reader)
    i = 0
    for row in reader:
        CHR_data[row[1].split(':')[0]] =  [row[1].split(':')[1]].extend(row[2:])
        i += 1

    return CHR_data, colnames


def rsID2CHRnum(repro_data):

    # Get Ids

    conn = sqlite3.connect('/scratch/users/vipese/GeneticAssociation/Resources/refdb.dbsnp')
    cur = conn.cursor()
    cur.execute("SELECT name FROM sqlite_master WHERE type='table';")

    cursor = conn.execute('select * from refvariants')
    names = list(map(lambda x: x[0], cursor.description))

    rsID_CHR_dict = OrderedDict()
    for row in cursor:
        if row[1] in list(repro_data.keys()):
            ins = repro_data[row[1]]
            ins = ins.insert(0, row[2])
            repro_data[row[1]] = ins

    return None


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

    # Convert rsID to chromosume number
    rsID2CHRnum(repro_data)


    # Load association analysis data
    CHR_data = loadCHRData('CHR_14.csv.gz')

    # Get matches
    getMatches(CHR_data, repro_data)

    # Load chromosome data
    print('stop')


main()