import ete3

species_tree = ete3.Tree()

profiles = "/Users/adriandavin/Desktop/Bioinformatics/SimuLyon/Cedric/ANIKET/TEST3/G/Profiles/Profiles.tsv"
extant_tree_file = "/Users/adriandavin/Desktop/Bioinformatics/SimuLyon/Cedric/ANIKET/TEST3/T/ExtantTree.nwk"

with open(extant_tree_file) as f:
    tree = ete3.Tree(f.readline().strip(), format=1)

leaves = [x.name for x in tree.get_leaves()]


with open(profiles) as f:

    header = f.readline().strip().split("\t")
    entries = [i for i, x in enumerate(header) if x in leaves]

    for line in f:

        handle = line.strip().split()

        fam = handle[0]
        v = sum([1 for x in entries if int(handle[x]) != 0])
        print(fam, str(v))



