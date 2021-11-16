import dataNavigation as nav
import tools.hashbox as hb
import tools.sequenceFeatures as sf
import tools.smithWaterman as sw
import numpy as np

def compileData(file1, file2, loc="./data/"):
    seq1 = nav.readSequenceFile(file1, loc)
    seq2 = nav.readSequenceFile(file2, loc)

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

    return flatten([hang, mismatch, skip, kmers, nfreqs, complex])

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