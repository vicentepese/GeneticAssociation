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

    # Reorganize dictionnary to create output
    pre_repro_data = pre_repro_data['MHC long range associations']
    repro_data = list()
    for val in pre_repro_data.values:
        repro_data.append(np.ndarray.tolist(val))

    return repro_data

def GzipFileHandler(FileName, Read=True):

    # If file is in gz format
    if '.gz' in FileName:
        if Read is True:
            OpenFile = gzip.open(FileName, 'rt', newline='')

        else:
            OpenFile = gzip.open(FileName, 'wb')

    # Else open in text format
    else:
        if Read is True:
            OpenFile = open(FileName, 'r')

        else:
            OpenFile = open(FileName, 'w')

    return OpenFile

def rsID2CHRnum(repro_data, options):

    # Connect to SQL dataset containing rsID and CHR number
    conn = sqlite3.connect(options["file"]["refdb"])
    cur = conn.cursor()
    cur.execute("SELECT name FROM sqlite_master WHERE type='table';")
    cursor = conn.execute('select * from refvariants')

    # Print: verbose
    print("Converting rsID to chromosome number. This might take several minutes.")

    # For each rsID in dataset, find if in reproData and take CHR number
    rsIds = np.array([item[1] for item in repro_data])
    for row in cursor:
        if row[1] in rsIds:

            # Find idx and append chromosome number at the beginning 
            idxs = np.where(rsIds == row[1])[0][0]
            for idx in np.nditer(idxs):
                repro_data[idx].insert(0, row[2])

    # Save data
    with open("Data/checkReproDataMod.csv" , 'w') as outFile:
        writer = csv.writer(outFile)
        writer.writerows(repro_data)
    print("rsID successfully converted to chromosome number")

    return repro_data
    

def getMatches(repro_data, options):

    # Open data from chromosome
    CHR_openFile = GzipFileHandler(options['file']['chrData'])

    # Reformat repro_data

    # Get positions 
    snp_pos = np.asarray([item[3] for item in repro_data])

    # For each SNP of the AA data, check if it matches the position of the chromosome
    # Run first with dicts, change to list of lists and see if there is any difference
    matches = list()
    print("Checking matches reproducibility. This might take several minutes")
    pos_match_count = 0

    for snp in CHR_openFile:

        # Find if there is a match in the SNP position
        if snp.split(",")[0] in snp_pos:

            # Find the idxs
            idxs = np.where(snp_pos == snp.split(",")[0])  

            pos_match_count += 1
            # print("%i position match found \n" % pos_match_count)

            # For each position match
            for idx in np.nditer(idxs):

                # Check if the gene matches as well and append
                if repro_data[idx][1] in snp.split(',')[1]:
                    print("Position %s matches with gene %s as %s" % (repro_data[idx][3], repro_data[idx][1], snp.split(',')[1]))
                    matches.append([repro_data[idx][3], repro_data[idx][1]] + snp.split(",")[1:] )
            
    print("Reproducibility checked")
    print("Number of matches found: %i" % len(matches))

    # Save matches 
    with open('reproOutput/matches.csv', 'w' ) as matches_out:
        writer = csv.writer(matches_out)
        writer.writerows(matches)
    
    repro_per =  len(matches) * 100 / len(repro_data) 
    print("The percentage of reproducibility is " + "{:.2%}".format(repro_per/100))
    

def main():

    # Load options
    with open("optionsCheck.json", 'r') as jsonFile:
        options = json.load(jsonFile)

    # If reproducibility data not processed, then process. Else load
    if "checkReproDataMod.csv" not in os.listdir(options['folder']['Data']):

        # Load original "reproducibility data"
        repro_data = loadReproData(options['file']["checkRepro"])

        # Convert rsID to chromosume number
        repro_data = rsID2CHRnum(repro_data, options)
    else:
        print("Reproducibility data loaded")
        with open(options['file']['reproDatamod'], 'r') as infile:
            csv_reader = csv.reader(infile)
            repro_data = list()
            for row in csv_reader:
                repro_data.append(row)


    # Get matches
    getMatches(repro_data, options)


main()
