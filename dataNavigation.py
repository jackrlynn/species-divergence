import pandas as pd

# Dictionary that converts scientific names (case sensitive) to abbreviated names
convSciToAbbr = {}

def initializeDataTools():

    # Add abbreviation conversions
    abbr = pd.read_csv("data/species_abbr.csv")
    abbr_lst = abbr["Abbr."]
    sci_lst = abbr["Scientific Name"]
    for i in range(len(abbr_lst)):
        convSciToAbbr[sci_lst[i]] = abbr_lst[i]


if __name__ == '__main__':
    initializeDataTools()