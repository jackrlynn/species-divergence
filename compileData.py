import dataNavigation as nav
import tools.hashbox as hb
import tools.sequenceFeatures as sf
import tools.smithWaterman as sw
import numpy as np
import os
import csv


def compileDataNuc(seq1):
    # Nucleotide frequencies
    nfreqs = []
    nfreqs.append(sf.nucleoFreq(seq1))
    # nfreqs.append(sf.dinucleoFreq(seq1))

    # Complexity
    complex = []
    complex.append(sf.k1Complexity(seq1))
    complex.append(sf.k2Complexity(seq1))

    return flatten([nfreqs, complex])


def compileDataAlign(seq1, seq2):
    print("Landmark #1: DATA LOADED")
    align1, align2, matrix, hang, mismatch, skip = sw.smithWaterman(seq1, seq2, 2, -1, -2, True, True, True, True)
    align1, align2 = sw.fixFront(seq1, seq2, align1, align2)
    align1, align2 = sw.fixBack(seq1, seq2, align1, align2)

    mismatch = np.nan_to_num(mismatch)
    mismatch = mismatch.tolist()

    print("LANDMARK #2: ALIGNMENT COMPLETE")

    # k-mers 11-15
    kmers = []
    for k in range(11, 16):
        kmers.append(hb.getSimilarity(seq1, seq2, k))
        kmers.append(hb.getSmilarityExcludeAlignment(seq1, seq2, k, align1, align2))
        print(f'LANDMARK #3: {k}-MER SIMILARITY COMPLETE')

    print(hang)
    print(mismatch)
    print(skip)
    print(kmers)
    return flatten([hang, mismatch, skip, kmers])


# SOURCE: https://stackabuse.com/python-how-to-flatten-list-of-lists/
def flatten(arr):
    isFlat = False
    while (not isFlat):
        isFlat = True
        flat_arr = []
        for line in arr:
            if type(line) == list:
                for i in line:
                    flat_arr.append(i)
                isFlat = False
            else:
                flat_arr.append(line)
        arr = flat_arr
    return arr


def createTable(directory):
    genes = []
    table = []
    names = []
    holdRow = []
    i = 0
    os.chdir(directory)

    table.append(['Species Name', 'C freq', 'A freq', 'T freq', 'G freq', 'K1 Complexity', 'K2 Complexity'])

    cat = ['Overhang1: ', 'Overhang2: ', 'Fronthang1: ', 'Fronthang2: ', 'Backhang1: ', 'Backhang2: ',
           'AA Match: ', 'AC Mismatch: ', 'AG Mismatch: ', 'AT Mismatch: ', 'CA Mismatch: ', 'CC Match: ',
           'CG Mismatch: ',
           'CT Mismatch: ', 'GA Mismatch: ', 'GC Mismatch: ', 'GG Match: ', 'GT Mismatch: ', 'TA Mismatch: ',
           'TC Mismatch: ', 'TG Mismatch: ', 'TT Match: ', 'Long Skip1: ', 'Long Skip2: ', '11-Mer WA: ',
           '11-Mer WOA: ',
           '12-Mer WA: ', '12-Mer WOA: ', '13-Mer WA: ', '13-Mer WOA: ', '14-Mer WA: ', '14-Mer WOA: ',
           '15-Mer WA: ', '15-Mer WOA: ']

    # iterate through all file
    for file in os.listdir():
        f = open(file, 'r')
        contents = f.read()
        names.append(file.replace('.txt', ''))
        genes.append(contents)

    for name in names:
        vs = 'Current/' + name
        for currentCat in cat:
            table[0].append(currentCat + vs)

    for gene1 in genes:
        holdNuc = compileDataNuc(gene1)
        print(holdRow)
        for gene2 in genes:
            print(gene1)
            print(gene2)
            holdRow.append(compileDataAlign(gene1, gene2))
            print(holdRow)
        table.append(flatten([names[i], holdNuc, holdRow]))
        print(table)
        holdRow = []
        i += 1

    with open('out.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(table)
