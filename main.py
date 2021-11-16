import compileData as cd

def main():
    arr = cd.compileData("Homo sapiens.txt", "Trichosurus vulpecula.txt", "./data/hfe/training/")
    print(arr)

if __name__ == '__main__':
    main()
