import ete3
lineages_counter = 0
output = list()
experiment = "FairyTree5"


#### Whole Tree

with open("/Users/aadavin/Desktop/simuLyon/Tests/EXPERIMENT/T/Events.tsv".replace("EXPERIMENT", experiment)) as f:
    f.readline()
    for line in f:
        time, e, nodes = line.strip().split("\t")
        if e == "S":
            lineages_counter += 1
        if e == "E":
            lineages_counter -=1
        if e == "F":
            final_time = float(time)
            break
        output.append((time, str(lineages_counter)))

with open("/Users/aadavin/Desktop/simuLyon/Tests/PlotsFairyTree/EXPERIMENT_WholeTreeLineages.tsv".replace("EXPERIMENT", experiment), "w") as f:

    for item in output:
        line = "\t".join(item) + "\n"
        f.write(line)


### Extant Tree

tree = "/Users/aadavin/Desktop/simuLyon/Tests/EXPERIMENT/T/ExtantTree.nwk".replace("EXPERIMENT", experiment)
with open(tree) as f:
    mytree = ete3.Tree(f.readline().strip(), format=1)
extant_leaves = {x.name for x in mytree.get_leaves()}

tree = "/Users/aadavin/Desktop/simuLyon/Tests/EXPERIMENT/T/WholeTree.nwk".replace("EXPERIMENT", experiment)
with open(tree) as f:
    mytree = ete3.Tree(f.readline().strip(), format=1)

root = mytree.get_tree_root()
lca = mytree.get_common_ancestor(extant_leaves)
lca_t = lca.get_distance(root)


print(lca.name)
print(lca_t)


tree = "/Users/aadavin/Desktop/simuLyon/Tests/EXPERIMENT/T/ExtantTree.nwk".replace("EXPERIMENT", experiment)
with open(tree) as f:
    mytree = ete3.Tree(f.readline().strip(), format=1)

root = mytree.get_tree_root()
nodes = list()

for node in mytree.traverse():
    if node.is_leaf():
        continue
    d = node.get_distance(root) + lca_t
    nodes.append((d, node.name))

nodes = sorted(nodes, key = lambda x: x[0])

with open("/Users/aadavin/Desktop/simuLyon/Tests/PlotsFairyTree/EXPERIMENT_ExtantTreeLineages.tsv".replace("EXPERIMENT", experiment), "w") as f:
    counter = 1
    for node in nodes:
        time = node[0]
        counter += 1
        line = "\t".join([str(time), str(counter)]) + "\n"
        f.write(line)

