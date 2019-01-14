import pandas as pd
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
def tailParser(inLoc):
    """
    parses .tail file into a nested list usable by other modules
    :param inLoc: CSV input .tail file as produced by aligner.tailcalc
    :return: nested list. [[sequence, #reads, gene, 3'end, tail len, tail seq],...]
    """
    f = open(inLoc, 'r')
    tails = f.readlines()
    f.close()

    tailList = []

    for i in range(len(tails)):
        if i==0: continue #skips the header
        line = tails[i].rstrip().split(',')
        tailList.append(line)
    return tailList
def repeater(item, list, reps):
    '''Takes an item and a list and then adds a copy of the item to the list reps number of times.'''
    for i in range(reps):
        list.append(item)
    return

def pdTailMaker(inLoc):
    """
    Takes standard tail file and returns a pandas dataframe
    """
    tails = tailParser(inLoc)
    pdTails = []
    for tail in tails:
        type = tail[2][tail[2].find("|")+1:]
        name = tail[2][:tail[2].find("|")]
        repeater([name,tail[3],tail[4],tail[5],type],pdTails,int(tail[1]))
    df = pd.DataFrame(pdTails,columns=['Gene','3Loc','TailLength','TailSeq', 'Type'])
    df[['3Loc','TailLength']] = df[['3Loc','TailLength']].apply(pd.to_numeric,errors='coerce')
    return df

def Ten_vs_Eleven_Finder(read):
    print(read[::-1])

def countsin2(inLoc):
    """
    Takes saved count file and reads it into a counts nested list.
    :param inLoc: counts file
    :return: nested list. counts nested list. [[read, total number, unique number],...]
    """
    countFile = open(inLoc, "r").readlines()
    counts = []
    for i in range(1, len(countFile)):
        temp = countFile[i].rstrip().split(",")
        counts.append([temp[0], temp[1], temp[2]])
    return counts






if __name__=="__main__":
    Ten_vs_Eleven_Finder("CTCTCTAGAGGGCCGGCC")
