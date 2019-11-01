import sqlite3
conn = sqlite3.connect('/scratch/users/vipese/GeneticAssociation/Data/refdb.dbsnp')
cur = conn.cursor()
cur.execute("SELECT name FROM sqlite_master WHERE type='table';")

cursor = conn.execute('select * from refvariants')
names = list(map(lambda x: x[0], cursor.description))

# missingPos = {line.strip().split('\t')[1]:tuple(line.strip().split('\t')) for line in open('MissingPositions.txt')}


missingPosdict = {}
for row in cursor:
    if row[1]:
        rsid = row[1]
        if rsid in missingPos:
            missingPosdict[rsid] = row[3]
            print('matched {} to position {}'.format(rsid, row[3]))


bimFile = open('GenRED.II.autosomal.Final.tped')
outBimfile = open('GenRED.II.autosomalSNPupdate.Final.tped', 'w')
excludeFile = open('GenRED.II.autosomalExclude.txt', 'w')
for row in bimFile:
    rowParse = row.strip().split(' ')
    if rowParse[1] in missingPos:
        if rowParse[1] in missingPosdict:
            rowParse[3] = str(missingPosdict.get(rowParse[1]))
            outBimfile.write(' '.join(rowParse)+'\n')
        else:
            excludeFile.write(row)
    else:
        outBimfile.write(row)
        #excludeFile.write(row)
excludeFile.close()
outBimfile.close()
bimFile.close()