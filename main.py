import dataNavigation as nav

def main():
    print(nav.getSpeciesDivergence("Sus scrofa", "Sapajus apella"))
    print(nav.getSpeciesDivergence("Sapajus apella", "Sapajus apella"))
    print(nav.getSpeciesDivergence("Homo sapiens", "Homo sapiens"))

if __name__ == '__main__':
    main()
