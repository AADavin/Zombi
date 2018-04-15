import os
import ete3

experiment = "FairyTree5"

path_to_whole_tree = "/Users/aadavin/Desktop/simuLyon/Tests/FairyTree5/T/WholeTree.nwk"
path_to_extant_tree = "/Users/aadavin/Desktop/simuLyon/Tests/FairyTree5/T/ExtantTree.nwk"
path_to_events = "/Users/aadavin/Desktop/simuLyon/Tests/FairyTree5/T/Events.tsv"

with open(path_to_whole_tree) as f:
    whole_tree = ete3.Tree(f.readline().strip(), format=1)
    r = whole_tree.get_tree_root()
    r.name = "Root"
    whole_nodes_dict = dict()
    for node in whole_tree.traverse():
        whole_nodes_dict[node.name] = node

with open(path_to_extant_tree) as f:
    extant_tree = ete3.Tree(f.readline().strip(), format=1)
    r = extant_tree.get_tree_root()
    r.name = "Root"
    extant_nodes_dict = dict()
    for node in extant_tree.traverse():
        extant_nodes_dict[node.name] = node

extinct_leaves = list()
time_origin = dict()
time_disappear = dict()



with open(path_to_events) as f:
    f.readline()
    for line in f:
        t, e, nodes = line.strip().split("\t")
        if e == "S":
            p,c1,c2 = nodes.split(";")
            time_origin[c1] = float(t)
            time_origin[c2] = float(t)
            time_disappear[p] = float(t)

        elif e == "E":
            extinct_leaves.append(nodes)

missing_lineages = dict()

for leaf in extinct_leaves:
    extinct_lineages = list()
    mynode = leaf

    while mynode not in extant_nodes_dict:
        extinct_lineages.append(mynode)
        mynode = whole_nodes_dict[mynode].up.name

    for node in extinct_lineages:
        if node not in missing_lineages:
            missing_lineages[node] = mynode


mypath = "/Users/aadavin/Desktop/simuLyon/Tests/FairyTree5/G/Events_per_branch"

events_per_branch = os.listdir(mypath)

transfer_events = list()

for myfile in events_per_branch:
    with open(os.path.join(mypath, myfile)) as f:
        f.readline()
        for line in f:
            t, e, nodes = line.strip().split("\t")
            if e == "AT":
                transfer_events.append((t,e,nodes))

# Now I have read all the transfer events. Now I have to add those lineages to the extant tree
# I need to get the time of leaving the tree and the time to arrive to the tree for each transfer
###

diversity = list()

already_seen = set()

for t, e, n in transfer_events:

    handle = n.split(";")
    donor = handle[0]
    recipient = handle[-2]

    if donor in missing_lineages:

        if missing_lineages[donor] in already_seen:
            continue

        t1 = time_disappear[missing_lineages[donor]]
        t2 = time_origin[recipient]

        already_seen.add(missing_lineages[donor])

    diversity.append((t1, donor, "LT"))
    diversity.append((t2, recipient, "AT"))


###




root = whole_tree.get_tree_root()
extant_leaves = {x.name for x in extant_tree.get_leaves()}
lca = whole_tree.get_common_ancestor(extant_leaves)
lca_t = lca.get_distance(root)

root = extant_tree.get_tree_root()

for node in extant_tree.traverse():
    if node.is_leaf():
        continue
    d = node.get_distance(root) + lca_t
    diversity.append((d, node.name, "S"))

diversity = sorted(diversity, key = lambda x: x[0])

for item in diversity:
    print(item)

with open("/Users/aadavin/Desktop/simuLyon/Tests/PlotsFairyTree/EXPERIMENT_FairyTreeLineages.tsv".replace("EXPERIMENT", experiment), "w") as f:
    counter = 1
    for time, node, event in diversity:
        if event == "LT" or event == "S":
            counter += 1
        if event == "AT":
            counter -= 1
        line = "\t".join([str(time), str(counter)]) + "\n"
        f.write(line)





















