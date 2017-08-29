def CSVWriter (iterable, outLoc, header="", ):
    """
    Writes an iterable to a CSV file.
    :param iterable: List of list
    :param outLoc: file location. Where to place it.
    :param header: header of the CSV file
    :return: 1
    """
    if not iterable:
        print ("nothing to write")
        return 0

    out = open(outLoc, 'w')

    if header:
        out.write(header+'\n')

    #Only works if iterable is a nested list
    for member in iterable:
        for item in member:
            out.write(str(item)+',')
        out.write('\n')

    print("write to "+outLoc+" successful.")
    return 1
def seqParser(seqLoc):
    """
    Takes a FASTA formatted list of sequences and returns a properly formatted nested list
    :param seqLoc: fasta formatted file
    :return: nested list in the form [[name, seq],...]
    """
    f=open(seqLoc,'r')
    RNAseqs = f.readlines()
    f.close()
    RNAseqlist = []
    for i in range(0, len(RNAseqs)):
        if RNAseqs[i][0] == ">":
            RNAseqlist.append([RNAseqs[i].rstrip()[1:],RNAseqs[i+1].rstrip()])
    return RNAseqlist
def countsin(inLoc):
    """
    Takes saved count file and reads it into a counts nested list.
    :param inLoc: counts file
    :return: nested list. counts nested list. [[read, total number, unique number],...]
    """
    countFile = open(inLoc, "r").readlines()
    counts=[]
    for i in range(1, len(countFile)):
        temp = countFile[i].rstrip().split(",")
        counts.append([temp[0][8:], temp[1], temp[2]])
    return counts
