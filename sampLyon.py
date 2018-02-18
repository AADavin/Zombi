import ete3
import random
import sys
import os
from globals import *

def prune_species_tree(infolder, outfolder, N):

    treefile = os.path.join(infolder, "WholeTree")

    with open(treefile) as f:
        mytree = ete3.Tree(f.readline().strip(), format=1)

    lineagesfile = os.path.join(infolder, "LineagesInTime.tsv")

    with open(lineagesfile) as f:
        for line in f:
            pass

    leaves_alive = line.strip().split("\t")[1].split(";")  # We pick up the last line of the file
    leaves_sampled = random.sample(leaves_alive, int(len(leaves_alive) * N))

    mytree.prune(leaves_sampled, preserve_branch_length=True)
    with open(os.path.join(outfolder, "SampledTree"),"w") as f:
        f.write(mytree.write(format=1))

def prune_gene_trees(infolder, outfolder):

    fams = os.listdir(os.path.join(infolder, "RawGeneFamilies"))

    with open(os.path.join(outfolder, "SampledTree")) as f:

        mytree = ete3.Tree(f.readline().strip(), format=1)
        leaves_sampled = [x.name for x in mytree.get_leaves()]

    for fam in fams:


        print("Pruning gene family %s" % fam)

        with open(os.path.join(os.path.join(infolder, "RawGeneFamilies"), fam)) as f:
            gf_tree = ete3.Tree(f.readline(), format=1)

        geneleaves_sampled = [x.name for x in gf_tree.get_leaves() if
                              x.name.split("_")[0] in leaves_sampled and x.name.split("_")[1] == "A"]

        if len(geneleaves_sampled) == 0:
            continue

        if len(geneleaves_sampled) < 3: # Ignoring small families
            continue

        else:
            gf_tree.prune(geneleaves_sampled, preserve_branch_length=True)

        with open(os.path.join(outfolder, fam), "w") as f:
            f.write(gf_tree.write(format=1))


if __name__ == "__main__":

    args = sys.argv[1:]


    if args[0] == "T":

        infolder, outfolder, N = args[1:]

        if float(N) <0 or float(N) > 1:
            print("Error, N must be comprised between 0 and 1")
        else:
            os.mkdir(os.path.join(infolder, outfolder))
            prune_species_tree(infolder, os.path.join(infolder, outfolder), float(N))

    elif args[0] == "G":

        infolder, outfolder = args[1:]
        prune_gene_trees(infolder, os.path.join(infolder, outfolder))

    else:
        print("Incorrect usage. Please read the manual. The usual way to run this script is:")
        print("python sampLyon.py fraction_of_lineages_sampled /Output_folder /Sampling_folder")







