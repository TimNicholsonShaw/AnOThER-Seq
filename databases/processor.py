file = open("JamesmiRNA.csv",'r')
out = open("JamesmiRNA_Tim.csv",'w')
for line in file:
    line=line.rstrip()
    gene = line.split(",")[6]
    period=gene.find(".")
    flag = -1
    for i in range(period + 1, len(gene)):
        if gene[i].isalpha() or gene[i] == ".":
            flag = i
            break



    if flag==-1: line+=","+gene+"\n"
    else: line+=","+gene[:flag]+"\n"

    out.write(line)

