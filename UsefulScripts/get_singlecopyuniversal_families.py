import ete3
import os
import sys
import argparse

def get_singlecopyuniversal_families(extant_tree_path, infolder, cutoff = 1):

    with open(extant_tree_path) as f:
        extant_tree = ete3.Tree(f.readline().strip(), format=1)

    sps = extant_tree.get_leaf_names()
    tnsp = len(sps)


    gene_trees = [x for x in os.listdir(infolder) if "_prunedtree.nwk" in x]

    for gene_tree in gene_trees:

        with open(os.path.join(infolder, gene_tree)) as f:

             t = f.readline().strip()
             if "(" not in t:
                 continue
             else:
                 t = ete3.Tree(t, format=1)

        sps = [x.split("_")[0] for x in t.get_leaf_names()]
        lsps = len(sps)

        if lsps == len(set(sps)) and lsps >= (tnsp * cutoff):
            print(gene_tree)

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("e", type=str, help="ExtantTree")
    parser.add_argument("f", type=str,  help="Folder with gene trees")
    parser.add_argument("--c", type=float, help="Universality cutoff")

    args = parser.parse_args()
    extant_tree, infolder, cutoff = args.e, args.f, args.c


    get_singlecopyuniversal_families(extant_tree, infolder, cutoff)

