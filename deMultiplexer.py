import counter, tools, sys, os
import pandas as pd
import re
import distance

def hitFinder(target, r1, r2):
    """
    Barcode added by Rea's protocol is found on the 5' end of the read. That barcode can be found in r1 if the read
    is short enough, but you'll definitely see it in read 2 since that starts from the 5' end. Using the illumina
    platform, r1 will need to be reverse complemented to be correct and r2 will need to not be. This function finds all
    reads with the barcode in r1 and then returns the corresponding r2 reads.
    :param target: Barcode to be found
    :param r1: list of reads produced by fastqparser
    :param r2: list of reads produced by fastqparser
    :return: list of reads
    """
    hits = []

    for i in range(0, len(r2)):
        if target in r2[i]:
            hits.append(r1[i])

    return hits

def hitFinder2 (target, r1, r2, ham=0):
    hits = []

    for i in range(0, len(r2)):

        if len(target) > len(r2[i]): continue

        if distance.hamming(target, r2[i][:len(target)]) <= ham:
            hits.append(r1[i])

    return hits

def deMultiplexer(name, barcode, ranMerLen, r1, r2, mismatch=0, ham=0):
    """
    Takes the name of the sample, the barcode sequence, and the length of the randomMer. Finds the barcode in read1,
    returns the corresponding reads in read2, performs the read count operation and outputs it to a folder with the
    samplename.csv
    :param name: experiment name
    :param barcode:
    :param ranMerLen:
    :param r1:
    :param r2:
    :param mismatch:
    :return:
    """
    print(name, barcode, ranMerLen)
    hits = hitFinder2(barcode, r1, r2, ham)

    counts = counter.uniqueFilter(readList=hits, barcodeLength=ranMerLen,barcodeMismatch=mismatch)
    return counts
def csv_from_excel(inLoc):
    xls = pd.ExcelFile(inLoc)
    df = xls.parse(sheetname="Sheet1", index_col=None, na_values=['NA'])
    df.to_csv('temp.csv', index=False)


def multipleLigationTrim(read, ranmer):
    patt = re.compile("AG[ATGC]{" + str(ranmer) + "}AGATCGGAAGAG")

    match = patt.search(read)

    if not match:
        return read
    else:
        return read[:match.start()]

########################################################################################################################

if __name__=="__main__":

    help = """
    #Takes 2 read files and a manifest file. Finds barcoded reads and outputs them as a counts file.
    #Mandatory arguments:
    #-r1: Read1 Location. Fastq format, gz compressed or not.
    #-r2: Read2 Location. Fastq format, gz compressed or not.
    #-m: Manifest location
    #-o: Folder to put output files into
    #-x(default 0): hamming distance
"""

    ham=0

    for x in range(0, len(sys.argv)):
        if sys.argv[x] == '-r1': r1Loc = sys.argv[x+1]
        if sys.argv[x] == '-r2': r2Loc = sys.argv[x+1]
        if sys.argv[x] == '-m' : manifestLoc = sys.argv[x+1]
        if sys.argv[x] == '-o' : outFolder = sys.argv[x+1]
        if sys.argv[x] == '-x' : ham = int(sys.argv[x+1])
        if sys.argv[x] == '-h': print(help);sys.exit()

    header = True #Does the manifest have a header

    flag = False
    if manifestLoc.endswith(".xlsx") or manifestLoc.endswith(".xls") or manifestLoc.endswith(".xlm"):
        flag = True
        csv_from_excel(manifestLoc)
    if flag == True: manifestLoc = "temp.csv"

    r1 = counter.fastqParser(r1Loc)
    r2 = counter.fastqParser(r2Loc, revComp=False)
    assert len(r1) == len (r2)

    manifest = open(manifestLoc, 'r')
    if header: next(manifest) #skips header
    for line in manifest:
        line = line.split(",")
        name = line[0]
        barcode = line[2].rstrip() + line[3].rstrip()
        ranMerLen = int(line[4]) + 2
        counts = deMultiplexer(name, barcode, ranMerLen, r1, r2, mismatch=1, ham=ham)
        counts2 = []
        for count in counts:
            temp_count = count[0][len(line[2].rstrip()):]
            temp_count = multipleLigationTrim(temp_count, ranMerLen-2)
            temp_count = multipleLigationTrim(temp_count, ranMerLen - 2)
            temp_count = multipleLigationTrim(temp_count, ranMerLen - 2)
            counts2.append([temp_count,count[1],count[2]])
        tools.CSVWriter(counts2,outLoc=outFolder+name+"_counts.csv",header="Sequence,TotalReads,UniqueReads")

    if flag: os.remove("temp.csv")



    





