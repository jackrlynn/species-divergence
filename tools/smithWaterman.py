import numpy as np
import matplotlib.pyplot as plt

def smithWaterman(seq1, seq2, ms, mp, gp, include_matrix=False, include_hang=False, include_mismatch=False, include_skip=False):
    len1 = len(seq1)
    len2 = len(seq2)

    # Create alignment matrix
    matrix = np.zeros((len1 + 1, len2 + 1))

    # Populate alignment matrix
    for i in range(1, len1 + 1):
        for j in range(1, len2 + 1):
            pos = np.array([matrix[i, j - 1] + gp])
            pos = np.append(pos, [matrix[i - 1, j] + gp])
            if seq1[i - 1] == seq2[j - 1]:
                pos = np.append(
                    pos, [matrix[i - 1][j - 1] + ms])
            else:
                pos = np.append(
                    pos, [matrix[i - 1][j - 1] + mp])
            pos = np.append(pos, [0])
            matrix[i, j] = pos.max()

    start = np.unravel_index(matrix.argmax(), matrix.shape)
    top = ""
    bottom = ""
    position = np.array(start)

    skips1 = 0
    skips2 = 0
    longestSkipstretch1 = 0
    longestSkipstretch2 = 0
    mismatchArray = np.empty(16)
    mismatchArray = np.resize(mismatchArray, [4, 4])

    skipStretch1 = 0
    longestSkipStretch1 = 0
    skipStretch2 = 0
    longestSkipStretch2 = 0

    while matrix[position[0], position[1]] != 0:
        options = np.array([matrix[position[0] - 1, position[1]]])
        options = np.append(options, matrix[position[0], position[1] - 1])
        options = np.append(options, matrix[position[0] - 1, position[1] - 1])

        ##Which squares could it have come from
        legal = [False, False, False]
        if ((options[2] == matrix[position[0], position[1]] - ms) or (
                options[2] == matrix[position[0], position[1]] - mp)):
            legal[2] = True
        if (options[1] == matrix[position[0], position[1]] - gp):
            legal[1] = True
        if (options[0] == matrix[position[0], position[1]] - gp):
            legal[0] = True

        # Which is best
        choices = np.empty(3)
        for x in range(2):
            if legal[x]:
                choices[x] = options[x]
            else:
                choices[x] = 0
        choice = np.argmax(choices)

        ##Record added base
        if (choice == 0):
            top = seq1[position[0] - 1] + top
            bottom = "-" + bottom
            position = [position[0] - 1, position[1]]
            skips2 = skips2 + 1
            skipStretch2 = skipStretch2 + 1
        if (choice == 1):
            bottom = seq2[position[1] - 1] + bottom
            top = "-" + top
            position = [position[0], position[1] - 1]
            skips1 = skips1 + 1
            skipStretch1 = skipStretch1 + 1
        if (choice == 2):
            top = seq1[position[0] - 1] + top
            bottom = seq2[position[1] - 1] + bottom
            position = [position[0] - 1, position[1] - 1]
            skipStretch1 = 0
            skipStretch2 = 0
            if (seq1[position[0] - 1] == "A"): k = 0
            if (seq1[position[0] - 1] == "C"): k = 1
            if (seq1[position[0] - 1] == "G"): k = 2
            if (seq1[position[0] - 1] == "T"): k = 3
            if (seq2[position[1] - 1] == "A"): l = 0
            if (seq2[position[1] - 1] == "C"): l = 1
            if (seq2[position[1] - 1] == "G"): l = 2
            if (seq2[position[1] - 1] == "T"): l = 3
            mismatchArray[k, l] = mismatchArray[k, l] + 1

        if (skipStretch1 > longestSkipStretch1):
            longestSkipStretch1 = skipStretch1
        if (skipStretch2 > longestSkipStretch2):
            longestSkipStretch2 = skipStretch2

    return_arr = [top, bottom]

    if (include_matrix):
        return_arr.append(matrix)

    if (include_hang):
        overhang1 = len(seq1) - (len(top) - skips1)
        overhang2 = len(seq2) - (len(bottom) - skips2)
        fronthang1 = position[0]
        fronthang2 = position[1]
        backhang1 = len(seq1) - start[0]
        backhang2 = len(seq2) - start[1]
        hang = [overhang1, overhang2, fronthang1, fronthang2, backhang1, backhang2]
        return_arr.append(hang)

    if (include_mismatch):
        return_arr.append(mismatchArray)

    if (include_skip):
        skip = [longestSkipStretch1, longestSkipStretch2]
        return_arr.append(skip)

    return return_arr

def addFrontToAlignment(seq, align, index):
    align_save = align
    align_no_dash = align.replace("-", "")
    i = seq.index(align_no_dash)
    if index == True:
        return [seq[:i] + align_save, i]
    else:
        return seq[:i] + align_save

def addBackToAlignment(seq, align, index):
    align_save = align
    align_no_dash = align.replace("-", "")
    i = seq.index(align_no_dash) + len(align_no_dash)
    if index == True:
        return [align_save + seq[i:], i]
    else:
        return align_save + seq[i:]

def fixFront(seq1, seq2, align1, align2):
    [new_seq1, i1] = addFrontToAlignment(seq1, align1, True)
    [new_seq2, i2] = addFrontToAlignment(seq2, align2, True)

    if (i1 < i2):
        new_seq1 = (i2 - i1) * '-' + new_seq1
    elif (i1 > i2):
        new_seq2 = (i1 - i2) * '-' + new_seq2

    return [new_seq1, new_seq2]

def fixBack(seq1, seq2, align1, align2):
    new_seq1 = addBackToAlignment(seq1, align1, False)
    new_seq2 = addBackToAlignment(seq2, align2, False)

    return [new_seq1, new_seq2]