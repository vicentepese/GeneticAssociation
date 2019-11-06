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
            OpenFile = gzip.open(FileName, 'rt', newline='')

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
        CHR_data[row[1]] = row[2:]

    return CHR_data, colnames


def rsID2CHRnum(repro_data, options):

    # Get Ids

    conn = sqlite3.connect(options["file"]["refdb"])
    cur = conn.cursor()
    cur.execute("SELECT name FROM sqlite_master WHERE type='table';")

    cursor = conn.execute('select * from refvariants')
    names = list(map(lambda x: x[0], cursor.description))

    print("Converting rsID to chromosome number. This might take several minutes.")
    for row in cursor:
        if row[1] in list(repro_data.keys()):
            ins = repro_data[row[1]]
            ins.insert(0, row[2])
            repro_data[row[1]] = ins

    # Save data
    with open("Data/checkReproDataMod.csv" , 'w') as outFile:
        for key in repro_data.keys():
            repro_data[key].insert(0, key)
            repro_data[key] = [str(itm) for itm in repro_data[key]]
            line = ", ".join(repro_data[key])
            outFile.write(line + "\n")

    print("rsID successfully converted to chromosome number")
    return repro_data
    

def getMatches(repro_data_old, options):


    # Open data from chromosome
    CHR_openFile = GzipFileHandler(options['file']['chrData'])

    # Reformat repro_data
    repro_data = dict()
    for key in list(repro_data_old.keys()):
        rp = [repro_data_old[key][1][1:]]
        for item in repro_data_old[key][3:]:
            rp.append(item)
        repro_data[repro_data_old[key][2][1:]] = rp

    # Get positions 
    snp_pos = list(repro_data.keys())
    # prot_genes = list(repro_data[key][1] for key in list(repro_data.keys()))


    # For each SNP of the AA data, check if it matches the position of the chromosome
    matches = dict()
    matches_nonGene = dict()
    pos_matches = {'Position': ['Gene', 'Repro_Gene']}
    print("Checking matches reproducibility. This might take several minutes")
    count = 1
    pos_match_count = 0
    gene_match = 0

    for snp in CHR_openFile:
        # print("Current line: %i" % count, end='\r')

        # Find if there is a match in the SNP position
        if any(snp.split(" ")[2] in pos for pos in snp_pos):
            pos_match_count += 1
            repro_data_gene = []
            repro_data_nongene = []
            print("%i position match found \n" % pos_match_count)

            # If the position and the gene matches
            repro_data_gene = []
            repro_data_nongene = []
            if repro_data[snp.split(" ")[2]][0] in snp.split("\"")[3]:
                repro_data_gene.append(repro_data[snp.split(" ")[2]])
            else:
                repro_data_nongene.append(repro_data[snp.split(" ")[2]])
            
            if len(repro_data_gene) > 0:
                rp = []
                for item in repro_data_gene[0]:
                    rp.append(item)
                matches[snp.split(" ")[2]] = [snp.split("\"")[3], rp]
            if len(repro_data_nongene) > 0:
                rp = [snp.split("\"")[3]]
                for item in repro_data_nongene[0]:
                    rp.append(item)
                matches_nonGene[snp.split(" ")[2]] = rp

            # # Check by rsID
            # for repro_rsID in list(repro_data.keys()):

            #     # Find rsID that matches the position
            #     if snp.split(" ")[2] in repro_data[repro_rsID[2]:

            #         # If there is a match in postion and in the gene
            #         if repro_data[repro_rsID][1] in snp.split("\"")[3]:
            #             repro_data_gene.append(repro_data[key][1][1:])

            #         # If there is a macth in the position but not in the gene 
            #         else: 
            #             repro_data_nongene.append(repro_data[key][1][1:])

            # # Save the match 
            # if len(repro_data_gene > 0):
            #     matches[snp.split(" ")[2]] = [snp.split("\"")[3], repro_data_gene]
            # if len(repro_data_nongene > 0):
            #     matches[snp.split(" ")[2]] = [snp.split("\"")[3], repro_data_nongene]

        count += 1
    print("Reproducibility checked")
    print("Number of matches found: %i" % len(matches))

    # Save matches 
    with open('matches.csv', 'w' ) as matches_out:
        for key, val in matches.items():
            matches_out.write(key + ", " + ', '.join(val) + ' \n')

    with open('matches_nonGene.csv', 'w' ) as matches_out:
        for key, val in matches_nonGene.items():
            matches_out.write(key + ", " + ', '.join(val) + ' \n')
    
    for key in list(matches.keys()):
        valn = [matches[key][0]]
        for val in matches[key][1]:
            valn.append(val)
        matches[key] = valn


    return matches        


def main():

    # Load options
    with open("optionsCheck.json", 'r') as jsonFile:
        options = json.load(jsonFile)

    # If reproducibility data not processed, then process. Else load
    if "checkReproDataMod.csv" not in os.listdir("./Data/"):
        # Load "reproducibility data"
        repro_data = loadReproData(options['file']["checkRepro"])

        # Convert rsID to chromosume number
        repro_data = rsID2CHRnum(repro_data, options)
    else:
        print("Reproducibility data loaded")
        with open(options['file']['reproDatamod'], 'r') as infile:
            csv_reader = csv.reader(infile)
            repro_data = dict()
            for row in csv_reader:
                repro_data[row[0]] = row[1:]


    # Get matches
    matches = getMatches(repro_data, options)

    # Load chromosome data
    print('stop')


main()

  # if any(repro_data[key][1] in snp.split("\"")[3] for key in list(repro_data.keys())):

            #     # If so, save the position
            #     gene_match += 1
            #     print("Position %s match with %s. Gene matches: %i \n" % (snp.split(" ")[2], snp.split("\"")[3], gene_match))
            #     repro_data_gene = []
            #     for key in list(repro_data.keys()):
            #         if snp.split(" ")[2] in repro_data[key][2]:
            #             repro_data_gene.append(repro_data[key][1][1:])

            #     # Save the match 
            #     matches[snp.split(" ")[2]] = [snp.split("\"")[3], repro_data_gene]
            # else:
            #     print("No match found for position %s \n" % snp.split(" ")[2] )
            #     repro_data_gene = []
            #     for key in list(repro_data.keys()):
            #         if snp.split(" ")[2] in repro_data[key][2] :
            #             repro_data_gene.append(repro_data[key][1][1:])
            #     pos_matches[snp.split(" ")[2]] = [snp.split("\"")[3], repro_data_gene]