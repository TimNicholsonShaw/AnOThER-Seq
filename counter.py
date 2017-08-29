import gzip, distance, time
from tools import CSVWriter


def reverseComplement(read):
    """
    :param read: a read in string format
    :return: the reverse complement of that read in string format
    """
    comp=[]

    for j in read:
        if j == 'A': comp.append('T')
        elif j == 'T': comp.append('A')
        elif j == 'C': comp.append('G')
        elif j == 'G': comp.append('C')
        else: comp.append('N')
    comp.reverse()
    comp=''.join(comp)
    return(comp)
def fastqParser(fileLoc, revComp = True):
    """
    Reads a fastq formatted file and returns a list of reads in string format
    :param fileLoc: fastq file, reads are assumed to be represented as 4 lines, the second of which is the basecall
    :param revComp: default True, should the function return the reverse complement of the reads
    :return: list of reads in string format
    """
    if fileLoc[-3:] == ".gz":
        fastq = gzip.open(fileLoc, 'rt').readlines()
    else:
        fastq = open(fileLoc, 'r').readlines()
    reads =[]


    for i in range(len(fastq)):
            if revComp:
                if (i - 1) % 4 == 0:
                    reads.append(reverseComplement(fastq[i].rstrip()))
            else:
                if (i - 1) % 4 == 0:
                    reads.append(fastq[i].rstrip())
    if not reads[-1]: reads = reads[:-1]
    return reads
def uniqueFilter(readList, barcodeLength, barcodeMismatch =0):
    """
    Filters reads, collapses identical reads into a more compact object. Identifies duplicate reads by 3' barcode
    :param readList: list of reads
    :param barcodeLength: length of barcode at 3' end of read used to find unique reads. set to 0 if no barcode
    is used.
    :param barcodeMismatch: how similar can barcodes be and still call them duplicate reads. 0 means that they
    need to be exactly the same. 1 means if there is only one nt difference they will be called duplicates.
    Until the alogrithm is sped up, anything other than 0 will be intensely slow. Barcodes are used to judge
    how many of those reads are unique
    :return: List. returns counts list of the structure [[read sequence, total reads, unique reads],...]
    """
    barcoded = []
    counts = []
    if barcodeLength:
        for read in readList:
            barcoded.append([read[:-barcodeLength],read[-barcodeLength:]])
        barcoded = sorted(barcoded, key = lambda x:x[0])
        tempBarcodes = []
        currRead = barcoded[0][0]
        for i in range(len(barcoded)):
            if i %100000 == 0: print (i/len(barcoded))
            if barcoded[i][0] == currRead:
                tempBarcodes.append(barcoded[i][1])
            else:
                counts.append([currRead, len(tempBarcodes),uniqueNumber(tempBarcodes,mismatch=barcodeMismatch)])
                currRead = barcoded[i][0]
                tempBarcodes = [barcoded[i][1],]
        counts.append([currRead, len(tempBarcodes),uniqueNumber(tempBarcodes, mismatch=barcodeMismatch)])
        counts = sorted(counts, key = lambda x:x[2], reverse=True)
        return counts
    else:
        readList = sorted(readList)
        currRead = readList[0]
        n=0
        for i in range(len(readList)):
            if i%100000 == 0: print (i/len(readList))
            if readList[i] ==  currRead:
                n+=1
            else:
                counts.append([currRead,n, n])
                currRead = readList[i]
                n = 0
        counts.append([currRead,n,n])
        counts = sorted(counts, key = lambda x:x[2], reverse = True)
        return counts
def uniqueNumber(list, mismatch = 0):
    """
    Used to judge how many barcodes are unique. If mismatch is greater than 0, it uses hamming distance to
    judge how different items are. The algorithm is slower if the mismatch is greater than 0
    :param list: list of string items
    :param mismatch: how many differences do two items need to have before they are not called duplicates
    :return: int. How many of the items in that list are unique
    """
    if not list:
        return 0
    if not mismatch:
        return len(set(list))
    #This algorithm is slightly slower than above
    n=1
    list = sorted(list)
    for i in range(len(list)-1):
        if distance.hamming(list[i], list[i+1]) > mismatch: n+=1
    return n
def countReads (inLoc, barcodeLength, outFolder="",mismatch = 0, name = "output"):
    """
    Takes fastq file, removes duplicate reads, condenses sequencing data into a counts file. Returns a list
    object. Can also write into a counts file for faster processing later.
    :param inLoc: str. fastq formatted file
    :param outFolder: str. optional. folder to output into. Be sure to include delimiter at the end
    :param barcodeLength: int. length of the barcode at the 3' end of the read. Used to judge duplicate reads
    :param mismatch: int. how different do two barcodes need to be for them to be judged as separate reads
    :param name: str. what name to give all files. defaults to "output"
    :return: List. returns counts list of the structure [[read sequence, total reads, unique reads],...]. Also
    prints to outFolder location.
    """
    print('Reading in from ' + inLoc)
    reads = fastqParser(inLoc)
    print('Counting...')
    reads = uniqueFilter(reads, barcodeLength, barcodeMismatch = mismatch)
    print('Count successful')
    if outFolder:
        CSVWriter(reads, outFolder+name+".counts", header="Sequence,Total Reads,Unique Reads")
    return reads




start_time =time.time()
countReads ("siLuc2_S5_L001_R1_001.fastq.gz", outFolder="/Users/Lykke-AndersenLab/Desktop/",barcodeLength=13, mismatch = 1)
print(time.time()-start_time)

