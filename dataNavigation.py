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