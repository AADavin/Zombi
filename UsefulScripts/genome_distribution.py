import os
import ete3
import pandas

mainpath  = "/Users/davin/Desktop/Zombi/BRISBANE/TestPhil"

with open(os.path.join(mainpath, "T","CompleteTree.nwk")) as f:
    stree = ete3.Tree(f.readline().strip(), format=1)

nodes = dict()
leaves = set()

for node in stree.traverse():
    if node.is_leaf():
        nodes[node.name] = {"type" : "leaf"}
        leaves.add(node.name)
    else:
        nodes[node.name] = {"type" : "inner"}
nodes["Initial"] = {"type" : "initial"}

for genome_file in os.listdir(os.path.join(mainpath, "G","Genomes")):
    with open(os.path.join(mainpath, "G","Genomes",genome_file)) as f:
        node = genome_file.split("_")[0]
        g  = f.readlines()
        nodes[node]["size"] = len(g)

with open(os.path.join(mainpath, "G","GenomeParameters.tsv")) as f:

    params = f.readlines()
    duplication_rate = float(params[5].split(":")[-1])
    transfer_rate = float(params[6].split(":")[-1])
    loss_rate = float(params[7].split(":")[-1])
    origination_rate = float(params[10].split(":")[-1])
    d_e = float(params[14].split(":")[-1])
    t_e = float(params[15].split(":")[-1])
    l_e = float(params[16].split(":")[-1])

# Universality

gene_families = dict()

profiles = pandas.read_csv(os.path.join(mainpath, "G","Profiles","Profiles.tsv"), sep = "\t")
leaf_profiles = profiles[list(leaves)]

universality = (leaf_profiles.astype(bool).sum(axis=1))[1:]
universality.to_csv(os.path.join(mainpath,"Universality.tsv"), index =False, header = False)

with open(os.path.join(mainpath, "Summary.tsv"), "w") as f:
    for n in nodes:
        line = "\t".join([n, nodes[n]["type"], str(nodes[n]["size"])]) + "\n"
        f.write(line)

events = {"C":0, "I":0, "D":0, "LT":0, "L":0, "O":0, "F":0, "AT":0}

# Number of events


event_files = os.listdir(os.path.join(mainpath, "G", "Events_per_branch"))

for event_file in event_files:

    with open(os.path.join(mainpath, "G","Events_per_branch", event_file)) as f:

        f.readline()

        for line in f:
            event = line.strip().split("\t")[1]
            events[event] +=1

with open(os.path.join(mainpath, "TotalNumberOfEvents.tsv"), "w") as f:
    f.write("Event\tNumber\n")
    for k,v in events.items():
        line = "\t".join([k, str(v)]) + "\n"
        f.write(line)









