f = open("HFE_refseq_transcript.fasta")
contents = f.read()
reads = contents.rsplit(">")
reads.pop(0)
for animal in reads:
    splitFile = animal.split("\n", 1)
    name=""
    splitName=splitFile[0].split(":")
    for i in splitName:
        name=name+i
    splitName=name.split(".")
    name=""
    for i in splitName:
        name=name+i
    splitName=name.split()
    notFoundName=True
    latinName=""
    for i in splitName:
        if(notFoundName):
            if(i!="PREDICTED"):
                if(i.isalpha()):
                    if (latinName==""):
                        latinName=i+" "+splitName[splitName.index(i)+1]
                        notFoundName=False


    print(latinName)
    name = latinName + ".txt"
    unstrippedList = splitFile[1].split()
    cleanSeq=""
    for i in unstrippedList:
        cleanSeq=cleanSeq+i.strip()
    n = open(name,"w")
    n.write(cleanSeq)
