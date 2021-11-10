import pandas as pd

def getSpeciesDivergence(species_1, species_2):

    # Get abbreviation
    sci_to_abbr = {}
    abbr = pd.read_csv("data/species_abbr.csv")
    abbr_lst = abbr["Abbr."]
    sci_lst = abbr["Scientific Name"]
    i = 0
    for i in range(len(abbr_lst)):
        sci_to_abbr[sci_lst[i]] = [abbr_lst[i], i]
        i += 1

    # Get divergence matrix
    diverg_mat = pd.read_csv("data/species_diverg.csv")
    return diverg_mat[sci_to_abbr[species_1][0]][sci_to_abbr[species_2][1]]

def parseFasta(file_name, loc="./data"):

    # Read in data file
    f = open(loc+file_name, 'r')
    contents = f.read()
    f.close()

    # Divide up reads into individual contents
    reads = contents.rsplit(">")
    reads.pop(0)

    # Iterate through each gene read and create a labelled txt file for that species's gene
    for animal in reads:
        splitFile = animal.split("\n", 1)

        # Break down file
        name = ""
        splitName = splitFile[0].split(":")
        for i in splitName:
            name = name + i
        splitName = name.split(".")
        name = ""
        for i in splitName:
            name = name + i
        splitName = name.split()

        # Try to find Latin name in file
        notFoundName = True
        latinName = ""
        for i in splitName:
            if (notFoundName):
                if (i != "PREDICTED"):
                    if (i.isalpha()):
                        if (latinName == ""):
                            latinName = i + " " + splitName[splitName.index(i) + 1]
                            notFoundName = False

        # Write the file name as a combo of Latin and .txt suffix
        name = latinName + ".txt"

        # Strip away all whitespace from sequence, leaving just the sequence
        unstrippedList = splitFile[1].split()
        cleanSeq = ""
        for i in unstrippedList:
            cleanSeq = cleanSeq + i.strip()

        # Write the file with sequence only
        n = open(name, "w")
        n.write(cleanSeq)
        n.close()