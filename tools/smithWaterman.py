import numpy as np
import matplotlib.pyplot as plt

def smithWaterman(seq1, seq2, ms, mp, gp):
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
        if (choice == 1):
            bottom = seq2[position[1] - 1] + bottom
            top = "-" + top
            position = [position[0], position[1] - 1]
        if (choice == 2):
            top = seq1[position[0] - 1] + top
            bottom = seq2[position[1] - 1] + bottom
            position = [position[0] - 1, position[1] - 1]

    return [top, bottom]

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