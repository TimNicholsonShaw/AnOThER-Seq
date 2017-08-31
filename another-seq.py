import counter, aligner,sys



########################################################################################################################



if __name__ == "__main__":
    database = "superset_withtRNA.fa"
    name = "output"
    allowedmismatch = 1

    help = """
    Performs the computationally taxing portions of this sequencing pipeline. Uses read 1, finds the random nucleotide
    sequence appended to the 3' end of the read during library prep and removes duplicate sequences. Creates a CSV type
    counts file of unique sequences in the form [Sequence, Total Reads, Unique Reads]. Then aligns it using BLAST
    to a database in the databases folder that can be specified. Using the alignment, judges the composition of the 3'
    end and outputs a tail file with that information. CSV formatted file with headers [sequence, unique reads,
    gene name, 3' end location, tail length, tail sequence]. Tail file can be used for downstream analysis.
    Mandatory arguments:
    -i: fastq file location. Use read 1.
    -o: Folder to put counts file and tail file into. Make sure it ends with the file delimiter
    -b: Length of random mer. Total length of added nucleotides. AG10 = 12, AG11 = 13

    Optional arguments
    -db: Name of database to use. Must be in databases folder. superset_withtRNA.fa is the default
    -n: Name to give files. Default is "output"
    -m: Allowed mismatch in the random barcode. Default is 1. If random barcode has 1 difference between eachother,
        they'll be called duplicates

    """

    for x in range(0, len(sys.argv)):
        if sys.argv[x] == '-i': inLoc = sys.argv[x+1]
        if sys.argv[x] == '-o': outFolder = sys.argv[x+1]
        if sys.argv[x] == '-b' : barcodeLength = int(sys.argv[x+1])
        if sys.argv[x] == '-db' : database = sys.argv[x+1]
        if sys.argv[x] == '-n': name = sys.argv[x+1]
        if sys.argv[x] == '-m': allowedmismatch = int(sys.argv[x+1])

        if sys.argv[x] == '-h': print(help);sys.exit()




    counts = counter.countReads(inLoc, barcodeLength, outFolder=outFolder, mismatch=allowedmismatch, name=name)
    alignedCounts=aligner.blaster(counts, database, outname = name)
    aligner.tailCalc(alignedCounts, database, outFolder=outFolder, outName=name)