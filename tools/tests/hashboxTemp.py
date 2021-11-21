import tools.hashbox as hb

def main():
    hb.getSimilarityExcludeAlignment("CCCAAACCC", "CCCAAA", 3, "CCCAAACCC", "---AAACCC")

if __name__ == '__main__':
    main()