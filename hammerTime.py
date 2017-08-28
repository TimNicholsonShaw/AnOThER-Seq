import gzip, distance, time


def reverseComplement(read):
    """Returns the reverse complement of a single read"""

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
    """Reads in a fastq formatted file and returns a list of list of reads
    Reads are assumed to be represented as 4 lines, the second line of which is the basecalls"""
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
    barcoded = []
    for read in readList:
        barcoded.append([read[:-barcodeLength],read[-barcodeLength:]])
    barcoded = sorted(barcoded, key = lambda x:x[0])
    counts = []
    tempBarcodes = []
    currRead = barcoded[0][0]
    for i in range(len(barcoded)):
        if barcoded[i][0] == currRead:
            tempBarcodes.append(barcoded[i][1])
        else:
            if not barcodeMismatch:
                counts.append([currRead, len(set(tempBarcodes))])
            else:
                counts.append([currRead, uniqueNumber(tempBarcodes,mismatch=barcodeMismatch)])
                currRead = barcoded[i][0]
                tempBarcodes = [barcoded[i][1],]
    counts.append([currRead, uniqueNumber(tempBarcodes, mismatch=barcodeMismatch)])
    counts = sorted(counts, key = lambda x:x[1], reverse=True)
    return counts
def uniqueNumber(list, mismatch = 0):
    unique = [list[0]]
    for barcode in list:
        flag = False
        for item in unique:
            if distance.hamming(barcode, item) <= mismatch:
                flag = True
                break
        if not flag:
            unique.append(barcode)
    return len(unique)









start_time =time.time()
reads = fastqParser("siLuc1_S1_L001_R1_001.fastq")
print('readin successful')
reads = uniqueFilter(reads, 12)
print(reads[:20])
print(time.time()-start_time)

