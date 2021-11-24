import dataNavigation as nav
import tools.hashbox as hb
import tools.sequenceFeatures as sf
import tools.smithWaterman as sw
import tools.alignment as al
import numpy as np
import os
import csv

def getRow(seq1, seq2):

    print("Landmark #1: DATA LOADED")

    align1, align2, hang, mismatch, skip = sw.smithWaterman(seq1, seq2, 2, -1, -2, False, True, True, True)
    align1, align2 = sw.fixFront(seq1, seq2, align1, align2)
    align1, align2 = sw.fixBack(seq1, seq2, align1, align2)

    mismatch = np.nan_to_num(mismatch)
    mismatch = mismatch.tolist()

    print("LANDMARK #2: ALIGNMENT COMPLETE")

    # k-mers 11-15
    kmers = []
    for k in range(11, 16):
        kmers.append(hb.getSimilarity(seq1, seq2, k))
        kmers.append(hb.getSimilarityExcludeAlignment(seq1, seq2, k, align1, align2))
    print(f'LANDMARK #3: k-MER SIMILARITY COMPLETE')

    # Nucleotide frequencies
    nfreqs = []
    nfreqs.append(sf.nucleoFreq(seq1))
    nfreqs.append(sf.nucleoFreq(seq2))
    nfreqs.append(sf.dinucleoFreq(seq1))
    nfreqs.append(sf.dinucleoFreq(seq2))
    print("LANDMARK #4: NUCLEOTIDE FREQUENCIES COMPLETE")

    # Complexity
    complex = []
    complex.append(sf.k1Complexity(seq1))
    complex.append(sf.k1Complexity(seq2))
    complex.append(sf.k2Complexity(seq1))
    complex.append(sf.k2Complexity(seq2))
    print("LANDMARK #5: COMPLEXITY COMPLETE")

    ali = al.getPaths(seq1, seq2, 2, -1, -2)
    number, maxAS, avgAS, avgOHt1, avgOHf1,  avgOHb1, avgOHt3, avgOHf2,  avgOHb2, minOHt1, minOHf1, minOHb1, minOHt2, minOHf2, minOHb2, maxOHt1, maxOHf1,  maxOHb1, maxOHt3, maxOHf2,  maxOHb2, amTotalAA, amTotalAC, amTotalAG, amTotalAT, amTotalCA, amTotalCC, amTotalCG, amTotalCT,  amTotalGA, amTotalGC, amTotalGG, amTotalGT,  amTotalTA, amTotalTC, amTotalTG, amTotalTT, avgSkips, maxSkips, avgSD, maxSD  = al.readPaths(ali)

    return flatten([hang, mismatch, skip, kmers, nfreqs, complex, number, maxAS, avgAS, avgOHt1, avgOHf1,  avgOHb1, avgOHt3, avgOHf2,  avgOHb2, minOHt1, minOHf1, minOHb1, minOHt2, minOHf2, minOHb2, maxOHt1, maxOHf1,  maxOHb1, maxOHt3, maxOHf2,  maxOHb2, amTotalAA, amTotalAC, amTotalAG, amTotalAT, amTotalCA, amTotalCC, amTotalCG, amTotalCT,  amTotalGA, amTotalGC, amTotalGG, amTotalGT,  amTotalTA, amTotalTC, amTotalTG, amTotalTT, avgSkips, maxSkips, avgSD, maxSD])

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

    names = [] # Holds the names of sequences
    seqs = [] # Holds the sequence of a corresponding index in names list
    table = [] # Holds all the data that will be written out to the CSV

    # Get all data from the files in the given directory
    for file in os.listdir(directory):

        # Open and read in file contents
        f = open(directory + file, 'r')
        contents = f.read()
        f.close()

        # Add name to names list
        names.append(file.replace('.txt', ''))

        # Add sequence to seqs list
        seqs.append(contents)

    header = ["species_1", "species_2",
              "overhang_1", "overhang_2",
              "fronthang_1", "fronthang_2",
              "backhang_1", "backhang_2",
              "AA_mismatch", "AC_mismatch", "AG_mismatch", "AT_mismatch",
              "CA_mismatch", "CC_mismatch", "CG_mismatch", "CT_mismatch",
              "GA_mismatch", "GC_mismatch", "GG_mismatch", "GT_mismatch",
              "TA_mismatch", "TC_mismatch", "TG_mismatch", "TT_mismatch",
              "long_skip_1", "long_skip_2",
              "11_mer", "11_no_align",
              "12_mer", "12_no_align",
              "13_mer", "13_no_align",
              "14_mer", "14_no_align",
              "15_mer", "15_no_align",
              "A_freq_1", "C_freq_1", "G_freq_1", "T_freq_1",
              "A_freq_2", "C_freq_2", "G_freq_2", "T_freq_2",
              "AA_freq_1", "AC_freq_1", "AG_freq_1", "AT_freq_1",
              "CA_freq_1", "CC_freq_1", "CG_freq_1", "CT_freq_1",
              "GA_freq_1", "GC_freq_1", "GG_freq_1", "GT_freq_1",
              "TA_freq_1", "TC_freq_1", "TG_freq_1", "TT_freq_1",
              "AA_freq_2", "AC_freq_2", "AG_freq_2", "AT_freq_2",
              "CA_freq_2", "CC_freq_2", "CG_freq_2", "CT_freq_2",
              "GA_freq_2", "GC_freq_2", "GG_freq_2", "GT_freq_2",
              "TA_freq_2", "TC_freq_2", "TG_freq_2", "TT_freq_2",
              "k1_1", "k1_2",
              "k2_1", "k2_2",
              "No_Paths",
              "maxAS", "avgAS",
              "avgOHt1", "avgOHf1",  "avgOHb1", "avgOHt3", "avgOHf2",  "avgOHb2",
              "minOHt1", "minOHf1", "minOHb1", "minOHt2", "minOHf2", "minOHb2",
              "maxOHt1", "maxOHf1",  "maxOHb1", "maxOHt3", "maxOHf2",  "maxOHb2",
              "amTotalAA", "amTotalAC", "amTotalAG", "amTotalAT",
              "amTotalCA", "amTotalCC", "amTotalCG", "amTotalCT",
              "amTotalGA", "amTotalGC", "amTotalGG", "amTotalGT",
              "amTotalTA", "amTotalTC", "amTotalTG", "amTotalTT",
              "avgSkips", "maxSkips", "avgSD", "maxSD"]
    table.append(header)

    # Do sequence comparisons for all genes
    for i in range(len(names)): # First gene for comparison
        for j in range(i, len(names)): # Second gene for comparison; formulation skips repeats
            to_add = [names[i], names[j]]
            to_add.append(getRow(seqs[i], seqs[j])) # Do sequence comparison
            to_add.append(nav.getSpeciesDivergence(names[i], names[j])) # Get time of divergence
            table.append(flatten(to_add)) # Add sequence to table

    # Write out data to CSV
    with open('out.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(table)