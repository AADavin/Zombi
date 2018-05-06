import ete3


with open("/Users/adriandavin/Desktop/Bioinformatics/SimuLyon/TEO/proof") as f:
    tree = ete3.Tree(f.readline(),format=1)


for node in tree.traverse():
    if hasattr(node, "l_dist"):
        print(node.name + "\t" + node.l_dist)