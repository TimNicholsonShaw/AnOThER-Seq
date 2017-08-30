import counter, aligner




if __name__ == "__main__":
    inLoc = "siNoct_S2_L001_R1_001.fastq.gz"
    barcodeLength = 12
    outFolder = "/Users/Lykke-AndersenLab/Dropbox/NoctData/"
    database = "superset_withtRNA.fa"
    name = "siNoct"
    allowedmismatch = 1


    counts = counter.countReads(inLoc, barcodeLength, outFolder=outFolder, mismatch=allowedmismatch, name=name)
    alignedCounts=aligner.blaster(counts, database, outname = name)
    aligner.tailCalc(alignedCounts, database, outFolder=outFolder, outName=name)