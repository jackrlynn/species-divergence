def hash(seq):
    h = 0
    for i in range(len(seq)):
        if (seq[i] == 'A'):
            h += 0 * (4 ** i)
        elif (seq[i] == 'C'):
            h += 1 * (4 ** i)
        elif (seq[i] == 'G'):
            h += 2 * (4 ** i)
        elif (seq[i] == 'T' or seq[i] == 'U'):
            h += 3 * (4 ** i)
    return h

def dehash(h, kmer_len):
    seq = ""
    i = kmer_len - 1
    while (i > -1):
        if (h // (4**i) == 0):
            seq = 'A' + seq
        elif (h // (4**i) == 1):
            seq = 'C' + seq
        elif (h // (4**i) == 2):
            seq = 'G' + seq
        elif (h // (4**i) == 3):
            seq = 'T' + seq

        h -= (h // (4**i)) * (4 ** i)
        i -= 1

    return seq

def createHashList(seq, k, should_sort):
    hash_lst = []
    for i in range(k-1, len(seq)):
        hash_lst.append(hash(seq[i-k+1:i+1]))
    if should_sort:
        hash_lst.sort()
    return hash_lst

def getSimilarity(seq1, seq2, k):
    lst1 = createHashList(seq1, k, True)
    lst2 = createHashList(seq2, k, True)

    ct = 0
    for kmer1 in lst1:
        for kmer2 in lst2:
            if (kmer1 == kmer2):
                ct += 1

    return ct

def isAligned(align1, align2, pos1, pos2):
    align_pos1 = convertPosToAlignedPos(pos1, align1)
    align_pos2 = convertPosToAlignedPos(pos2, align2)
    if (align_pos1 == align_pos2):
        return True
    else:
        return False

def convertPosToAlignedPos(pos, align):
    new_pos = 0
    i = 0
    while (i != pos):
        if (align[new_pos] == '-'):
            new_pos += 1
        else:
            new_pos += 1
            i += 1
    return new_pos

def getSmilarityExcludeAlignment(seq1, seq2, k, align1, align2):
    lst1 = createHashList(seq1, k, False)
    lst2 = createHashList(seq2, k, False)

    ct = 0
    for i in range(len(lst1)):
        for j in range(len(lst2)):
            kmer1 = lst1[i]
            kmer2 = lst2[j]
            if (kmer1 == kmer2):
                if (not isAligned(align1, align2, i+k-1, j+k-1)):
                    ct += 1

    return ct