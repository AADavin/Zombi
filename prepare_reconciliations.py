import ete3
import random
from io import StringIO

project = ete3.Phyloxml()

# SPECIES TREE

phylo = ete3.phyloxml.PhyloxmlTree(newick="((A,B),C);")
project.add_phylogeny(phylo)
phylo.phyloxml_phylogeny.set_name("Species_Tree")

if len(phylo.children) <= 2:
    phylo.phyloxml_phylogeny.set_rooted("true")
else:
    phylo.phyloxml_phylogeny.set_rooted("false")

project.add_phylogeny(phylo)

#print(project.export())


# GENE TREE

phylo = ete3.phyloxml.PhyloxmlTree(newick="((a,b),c);")

# We add the events:

for node in phylo.traverse():

    node.add_feature("rec", 0)

project.add_phylogeny(phylo)
phylo.phyloxml_phylogeny.set_name("Gene_Tree")

if len(phylo.children) <= 2:
    phylo.phyloxml_phylogeny.set_rooted("true")
else:
    phylo.phyloxml_phylogeny.set_rooted("false")

project.add_phylogeny(phylo)
print(project.export())

