import counter, aligner




if __name__ == "__main__":
    inLoc = "NoctWT_S3_L001_R1_001.fastq.gz"
    barcodeLength = 13
    outFolder = "/Users/Lykke-AndersenLab/Desktop/"
    database = "superset_withtRNA.fa"
    name = "NoctWT"
    allowedmismatch = 1


    counts = counter.countReads(inLoc, barcodeLength, outFolder=outFolder, mismatch=allowedmismatch, name=name)
    alignedCounts=aligner.blaster(counts, database, outname = name)
    aligner.tailCalc(alignedCounts, database, outFolder=outFolder, outName=name)