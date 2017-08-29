import tools
from scipy import stats
import numpy as np

def getNames(tails):
    """
    Returns a list of the names of all genes found in a tail file
    :param tails: nested list. Produced by tail parser from a .tails file
    :return: list of gene names
    """
    names = []
    for tail in tails:
        names.append (tail[2])
    return list(set(names))
def tailFilter(tails, max3End=-10, maxTail=10):
    """
    Filters a nested tail list produced from a .tails file by tailparser. Tail must have mapped to something. Tail
    can't be trimmed past the max 3 End. Tail can't be longer than the maxTail
    :param tails: nested list. Produced by tail parser from a .tails file
    :param max3End: Tail is filtered out if the read is trimmed past the 3' end
    :param maxTail: Tail is filtered out if its extension is beyond this length
    :return:
    """
    """Filters tails based on current criteria:
        1) Must have mapped to something
        2) Can't be chewed back farther than max3end (default is -10)"""

    tail1=[]
    tail2=[]
    for tail in tails:
        try:
            int(tail[3])
            tail1.append(tail)
        except:
            pass
    for tail in tail1: #must have mapped something closer to 3' end than -10
        if int(tail[3])+int(tail[4])>=max3End and int(tail[3])+int(tail[4])<=maxTail:
            tail2.append(tail)
    return tail2

def getTailComparison(tail1, tail2, transcript):
    """
    Takes two nested list tail data structures as provided by tailParser. Looks for transcript in each tail file,
    counts the total number of transcripts, gets the total tail (3loc + tail len), finds the average, standard
    deviation, and then finds a p-value between the two sets based
    :param tail1: nested list. Produced by tail parser from a .tails file
    :param tail2: nested list. Produced by tail parser from a .tails file
    :param transcript: gene name
    :return:
        Output tuple:
        0 - Transcript Name
        1 - Tail1: Number of transcripts
        2 - Tail1: Average of total tails
        3 - Tail1: Standard Deviation of total tails
        4 - Tail2: Number of transcripts
        5 - Tail2: Average of total tails
        6 - Tail2: Standard deviation of total tails
        7 - P-value between total tails
        """
    totalTails1 = []
    totalTails2 = []
    for tail in tail1:
        if transcript in tail[2]:
            tools.repeater(int(tail[3])+int(tail[4]), totalTails1, int(tail[1]))
    for tail in tail2:
        if transcript in tail[2]:
            tools.repeater(int(tail[3])+int(tail[4]), totalTails2, int(tail[1]))
    if not totalTails1 or not totalTails2:
        pCombo = "nan"
    else:
        pCombo = stats.ks_2samp(totalTails1, totalTails2)[1]
    return (transcript, len(totalTails1), np.average(totalTails1),\
           np.std(totalTails1),len(totalTails2), np.average(totalTails2),\
           np.std(totalTails2), pCombo)
def getCandidates(tail1, tail2, pThresh=0.05, minTranscripts=10, minDiff = 0.10, outFolder="", outName="output"):
    """
    Takes two nested list tail data structures. Gets the names of all the transcripts mapped in tail1. Looks for each
    individually in tail2 and gets the tail comparison. Then looks for hits based on the following criteria:
    1) p-value must be less than or equal to pthresh default is 0.05
    2) each tail data structure must contain at least mintranscripts number of transcripts (default is 10)
    3) The difference between the average tail lengths must avgDiff fraction different. default is 0.10
    :param tail1: nested list. Produced by tail parser from a .tails file
    :param tail2: nested list. Produced by tail parser from a .tails file
    :param pThresh: p-value must be less than or equal to pthresh default is 0.05
    :param minTranscripts: each tail data structure must contain at least mintranscripts number of
    transcripts (default is 10)
    :param minDiff: The difference between the average tail lengths must avgDiff fraction different. default is 0.10
    :return: nested list. Ordered by p-value. Each member of the nested list contains:
        0 - Transcript Name
        1 - Tail1: Number of transcripts
        2 - Tail1: Average of total tails
        3 - Tail1: Standard Deviation of total tails
        4 - Tail2: Number of transcripts
        5 - Tail2: Average of total tails
        6 - Tail2: Standard deviation of total tails
        7 - P-value between total tails
        """
    candidates = []
    transcripts = getNames(tail1)
    for transcript in transcripts:
        stat = getTailComparison(tail1, tail2, transcript)
        try: #fillers for when stat tests fails is "nan", probably a better way to do this, but this catches all the
            #exceptions when this pops up
            if stat[7] > pThresh: continue #not a candidate if the p-value is greater than the threshold
            if stat[1] < minTranscripts or stat[4] < minTranscripts: continue

            perDiff = abs((stat[2]-stat[5])/stat[2]) #fractional difference between averages
            if perDiff < minDiff: continue #need a min diff between the two averages, doesn't care longer or shorter

            candidates.append(stat)
        except:
            continue
    candidates = sorted(candidates, key=lambda x: x[7])
    if outFolder:
        tools.CSVWriter(candidates, outFolder+outName+"_candidates.csv", header="Transcript,1-Reads,1-Avg,\
        1-STDev,2-reads,2-Avg,2-STDev,P-value")

    return candidates
def typeCounter(tails, typeList=[]):
    """
    Takes nest tail list data structure and compares it to a list of types for transcripts. Your own type list can
    be provided or it can default to below. It then looks through the type list and counts how many of each type
    of transcript it sees. Returns a nested list of [[type, number of transcripts],...].
    :param tails: nested list. Produced by tail parser from a .tails file
    :param typeList: list. optional. can specify types to look for. automatically generated if not specified
    :return: nested list. [[type, number of transcripts],...]
    """
    typeCounts=[]
    if not typeList: #default types list
        typeList = ['Mt_rRNA', 'Mt_tRNA','gtRNA', 'miRNA', 'misc_RNA',\
                'rRNA', 'scRNA', 'snRNA', 'snoRNA',\
                'ribozyme', 'sRNA', 'scaRNA', 'lncRNA', \
                'no db match', 'Tail Analyzer Failed']

    for item in typeList:
        count = 0
        for tail in tails:
            if item in tail[2]:
                count += int(tail[1])
        typeCounts.append((item, count))
    typeCounts = sorted(typeCounts, key=lambda x:x[1], reverse=True) #sorts by total number of reads
    return typeCounts
def getCumulativeTail(tails, transcript, min=-10, max=10):
    """
    Takes a tail formatted list. Looks for the transcript and finds the cumulative tail plot of total tail length for
    that specific transcript.
    :param tails: nested list. Produced by tail parser from a .tails file
    :param transcript: str. gene name of transcript to plot
    :param min: minimum x coordinate to plot
    :param max: maximum x coordinate to plot
    :return: tuple of the form (X coords, y coords)
    """
    sepTails = []
    for tail in tails:
        if transcript in tail[2]:  #finds all occurences of the transcript and puts it into a new list.
            sepTails.append(tail)

    total = 0
    for tail in sepTails: total+=int(tail[1])  #counts the total number of reads that match the transcript

    x = list(range(min, max))  # creates the x values as a range between the max and min
    y=[]
    for num in x:
        count = 0
        for tail in sepTails:  #gets the cumulative tail percentage and adds it to the list y
            if (int(tail[3])+int(tail[4])) <= num: count += int(tail[1])
        try:
            y.append(count/total)
        except:
            y.append(0)

    return (x,y)