import subprocess, os

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
    return 1
def getBLASTdb(dbName = "superset_withtRNA.fa"):
    #if os = mac.....
    makeblastdbloc = os.path.join(os.path.dirname(__file__),"blast/makeblastdb.exe")
    databaseloc = os.path.join(os.path.dirname(__file__),"databases/"+dbName)
    subprocess.call([makeblastdbloc])








def blaster(counts, dbLoc, blastLoc = "blastn.exe", blastOut="blastTemp.txt"):
    """
    Takes a list of counts in the format [[sequence, total reads, unique],...] then blasts each sequence against the
    RNASeqList (fasta formatted). Matches are appended to the end of each member of the count. Returns this new appended count list.
    :param counts: counts object nested list in the format of [[sequence, total reads, unique],...]
    :param dbLoc: database of reads to blast against
    :param blastLoc:
    :param blastOut:
    :return:
    """
    queryMaker(counts, outLoc ="queryTemp.txt")
    open(blastOut,'w').close() #blast won't be happy if the file doesn't already exist
    #makes a call to blast everything in the query file against the database. Outputs to a temporary blast file.
    subprocess.call([blastLoc, '-db', dbLoc , '-query', "queryTemp.txt", '-out', 'blastTemp.txt', '-outfmt', '6', '-max_target_seqs', '1'])
    os.remove("queryTemp.txt")

    #reads temporary blast file into hits
    f = open("blastTemp.txt", 'r')
    hits = f.readlines()
    f.close()
    #os.remove(blastOut)

    for count in counts:
        count.append([]) #Adds a fourth entry to each count for add blast hits to.
    for line in hits:
        #Each query was given a number by querymaker which corresponds to the index of each count
        #Each hit also has this index so a hit can be assigned to a count by that index
        line = line.split('\t')
        i = int(line[0])
        counts[i][3].append(line[1])

    #each count now has the blast hits appended to it
    return counts


getBLASTdb()
