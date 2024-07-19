def count(file):
    store = []
    with open(file, "r") as f:
        lines = f.readlines()
        for line in lines:
            if not line.startswith(">"):
                if len(line) != 143:
                    print("not equal 143")
                store.append(len(line))
    # print(store)


count("./split_1.fasta")
