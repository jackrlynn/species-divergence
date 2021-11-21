import math


# K1 complexity
def k1Complexity(sequence):
    window = len(sequence)
    CV = nucleoFreq(sequence)
    k1 = 0

    k1 += (1 / window) * (math.log(math.factorial(window) // (
                math.factorial(CV[0]) * math.factorial(CV[1]) * math.factorial(CV[2]) * math.factorial(CV[3])), 4))

    return k1


# K2 complexity
def k2Complexity(sequence):
    window = len(sequence)
    CV = nucleoFreq(sequence)
    K2 = 0

    for i in range(0, 4):
        if CV[i] == 0:
            K2 += 0
        else:
            K2 += (CV[i] / window) * (math.log2(CV[i] / window))

    return -K2


def dinucleoFreq(sequence):
    diNucDict = {"AA": 0, "AC": 0, "AG": 0, "AT": 0, "CA": 0, "CC": 0, "CG": 0, "CT": 0, "TA": 0,
                 "TC": 0, "TG": 0, "TT": 0, "GA": 0, "GC": 0, "GG": 0, "GT": 0}
    diNuc = ""

    for i in range(0, len(sequence), 1):
        if i > 0:
            diNuc = sequence[i - 1] + sequence[i]
            for key in diNucDict:
                if key == diNuc:
                    diNucDict[key] += 1
    print(diNucDict)
    return list(diNucDict.values())


# Nucleotide Frequencies
def nucleoFreq(sequence):
    count = [0, 0, 0, 0]

    for currentNuc in sequence:
        if currentNuc == 'C':
            count[0] += 1
        elif currentNuc == 'A':
            count[1] += 1
        elif currentNuc == 'T':
            count[2] += 1
        elif currentNuc == 'G':
            count[3] += 1
    return count