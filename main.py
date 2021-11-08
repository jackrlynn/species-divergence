import dataNavigation as nav

def main():
    print(nav.getSpeciesDivergence("Cercocebus atys", "Homo sapiens"))
    print(nav.getSpeciesDivergence("Sapajus apella", "Sapajus apella"))
    print(nav.getSpeciesDivergence("Homo sapiens", "Homo sapiens"))

if __name__ == '__main__':
    main()
