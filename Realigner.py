import subprocess, os, sys, time
import counter, tools

def queryMaker(readList, outLoc="queryTemp.txt"):
    """
    Takes a counts object (nested list in the form [[read, total, unique],...]) and creates a file formatted for input
    into BLAST.
    :param readList: nested list in the form [[read, total, unique],...]
    :param outLoc: file location. temporary place to store the query file
    :return: None
    """
    """Takes a list of reads in the form of [[read, total, unique],...] and creates a file formatted for input into
    BLAST. """
    f= open(outLoc, 'w')
    for i in range(len(readList)):
        f.write(">"+str(i)+"\n")
        f.write(readList[i][0]+"\n")
    f.close()
    return True
def checkBLASTdb(dbName = "superset_withtRNA.fa"):
    """
    Checks if BLAST indices have been made. Creates the indices if they're not found. Works for mac and pc.
    :param dbName: title of fasta formatted database in databases folder
    :return: None
    """
    databaseloc = os.path.join(os.path.dirname(__file__),"databases/"+dbName)
    if not os.path.isfile(databaseloc+".nhr"):
        if sys.platform == "darwin": makeblastdbloc = os.path.join(os.path.dirname(__file__), "blast/makeblastdb")
        elif sys.platform.startswith("win"): makeblastdbloc = os.path.join(os.path.dirname(__file__), "blast/makeblastdb.exe")
        else:
            raise Exception ("System is not supported")
        subprocess.call([makeblastdbloc, "-in", databaseloc, "-dbtype", "nucl"])
    return True
def blaster(counts, database, outname = "output"):
    """
    Takes a list of counts in the format [[sequence, total reads, unique],...] then blasts each sequence against the
    RNASeqList (fasta formatted). Top match is appended to the end of each member of the count.
    :param counts: counts object nested list in the format of [[sequence, total reads, unique],...]
    :param dbLoc: database of reads to blast against
    :return: nested list of aligned counts. [[sequence, total reads, unique reads, blast alignment outfmt 6],...]
    """
    queryLoc = outname+"queryTemp.txt"
    blastOut = outname+"blastTemp.txt"

    #Set blast location
    if sys.platform == "darwin":
        blastLoc = os.path.join(os.path.dirname(__file__), "blast/blastn")
    elif sys.platform.startswith("win"):
        blastLoc = os.path.join(os.path.dirname(__file__), "blast/blastn.exe")
    else:
        raise Exception("System is not supported")

    checkBLASTdb(database)
    dbLoc = os.path.join(os.path.dirname(__file__), "databases/"+database)

    print ("Creating BLAST query")
    queryMaker(counts, outLoc =queryLoc)
    print(str(len(counts))+" queries.")
    open(blastOut,'w').close() #blast won't be happy if the file doesn't already exist
    #makes a call to blast everything in the query file against the database. Outputs to a temporary blast file.

    print("BLASTing...")
    start_time = time.time()
    subprocess.call([blastLoc, '-db', dbLoc , '-query', queryLoc, '-out', blastOut, '-outfmt', '6'])
    os.remove(queryLoc)
    print("BLAST Successful")
    print("BLAST took "+str(time.time()-start_time)+" seconds")

    #reads temporary blast file into hits
    f = open(blastOut, 'r')
    hits = f.readlines()
    f.close()
    #os.remove(blastOut)

    flag = -1
    for line in hits:
        #Each query was given a number by querymaker which corresponds to the index of each count
        #Each hit also has this index so a hit can be assigned to a count by that index
        templine = line.split('\t')
        i = int(templine[0])
        if i<=flag: continue
        counts[i].append(line.rstrip())
        flag=i
    namehits = [""]*len(counts)
    scores = [0]*len(counts)
    for line in hits:
        templine = line.split('\t')
        i = int(templine[0])
        score = float(templine[11].rstrip())
        name = templine[1]

        if score>=scores[i]:
            namehits[i]+=" | "+name
            scores[i]=score
    for i in range(len(namehits)):
        if namehits[i]: counts[i].append(namehits[i][3:])


    #each count now has its BLAST hit appended to it
    return counts
def seqFinder(name, seqList ):
    """
    Takes a gene name and matches it to the master RNA sequence list. Returns the sequence of that gene.
    :param name: str. gene name, must be exact match
    :param seqList: nested list. [[gene, sequence],...]. Can be produced from fasta formatted file using tools.seqParser
    :return: str. sequence of gene specified in name.
    """
    for i in range(len(seqList)):
        if name in seqList[i][0] == name:
            return seqList[i][1]
    print (name)
    raise Exception ("Couldn't find the sequence. No idea what happened in seqFinder. Shouldn't have happened")
def blastTailer(alignedCount, fastaList):
    """
    takes an aligned count produced by blaster function. Compares its 3' end to its target sequence and judges
    differences in length and sequence.
    :param alignedCount: list. [sequence, total reads, unique reads, blast alignment out fmt 6]
    :param fastaList: nested list. [[gene, sequence],...]. Can be produced from fasta formatted file using tools.seqParser
    :return: list. tail object. [sequence, unique reads, gene name, 3' end location, tail length, tail sequence]
    """
    sequence = alignedCount[0]
    uniqReads = alignedCount [2]
    if len(alignedCount)!=5: return [sequence, uniqReads, "no db match", "N/A", "N/A"]
    alignment = alignedCount[3]

    blastLine = alignment.split('\t')

    targetSeq = seqFinder(blastLine[1], fastaList)
    threeEnd = int(blastLine[9]) - len(targetSeq)
    tailLen = len(sequence) - int(blastLine[7])
    tailSeq = sequence[int(blastLine[7]):]

    return [sequence, uniqReads, alignedCount[4], threeEnd, tailLen, tailSeq,blastLine[4]]
def tailCalc (alignedCounts, dbName, outFolder="",outName="output"):
    """
    takes aligned counts and uses blastTailer algorithm to judge tail lengths. Outputs a list of tail objects.
    :param alignedCounts: nested list. [[sequence, total reads, unique reads, blast alignment out fmt 6],...]
    :param dbName: atabase of reads to blast against
    :param outFolder: str. folder to output files to
    :param outName: str. name to give otuput
    :return: nested list. tail objects. [[sequence, unique reads, gene name, 3' end location, tail
    length, tail sequence],...]
    """
    print('Generating tail files...')
    databaseloc = os.path.join(os.path.dirname(__file__), "databases/" + dbName)
    fastaList = tools.seqParser(databaseloc)
    tails = []
    i = 0
    for count in alignedCounts:
        if i%10000 == 0: print ("{:.0%}".format(i/len(alignedCounts)))
        tails.append(blastTailer(count,fastaList))
        i+=1

    tools.CSVWriter(tails, outFolder+outName+".tails", header="Sequence,UniqueReads,Gene,3Loc,TailLength,TailSequence,Mismatches")
    return tails


if __name__ == "__main__":

    base_dir = "/Users/tlshaw/Desktop/Rea/"

    database = "superset_withtRNA.fa"
    name = "output"

    counts = tools.countsin2(base_dir + "U1-v15+16-01_counts.csv_Human snRNA sequences-190111_Taildata.csv")

    print(counts[0])
    alignedCounts = blaster(counts, database, outname="U1-v15")
    tailCalc(alignedCounts, database, outFolder=base_dir, outName="U1-v15")