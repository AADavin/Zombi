import ete3
import random
import sys
import os
from globals import *

class SamplingManager():

    def __init__(self, N, infolder, outfolder):

        self.infolder = infolder
        self.outfolder = os.path.join(infolder,outfolder)

        treefile = os.path.join(infolder,"WholeTree")

        with open(treefile) as f:
            self.mytree = ete3.Tree(f.readline().strip(), format=1)

        lineagesfile = os.path.join(infolder, "LineagesInTime.tsv")

        with open(lineagesfile) as f:
            for line in f:
                pass

        self.leaves_alive = line.strip().split("\t")[1].split(";") # We pick up the last line of the file

        if int(N) > len(self.leaves_alive):
            print("Error")
            pass

        self.leaves_sampled = random.sample(self.leaves_alive, int(N))

    def prune_species_tree(self):

        self.mytree.prune(self.leaves_sampled, preserve_branch_length=True)
        with open(os.path.join(self.outfolder, "SampledTree"),"w") as f:
            f.write(self.mytree.write(format=1))

    def prune_gene_trees(self):

        gene_families = os.path.join(self.infolder, "Profiles.tsv")

        fams = list()

        with open(gene_families) as f:
            f.readline()
            for line in f:
                fam = line.strip().split("\t")[0]
                fams.append(fam)

        for fam in fams:

            with open(os.path.join(os.path.join(self.infolder, "RawGeneFamilies"), fam)) as f:
                gf_tree = ete3.Tree(f.readline(), format=1)

            geneleaves_sampled = [x.name for x in gf_tree.get_leaves() if
                                  x.name.split("_")[0] in self.leaves_sampled and x.name.split("_")[1] == "A"]

            if len(geneleaves_sampled) == 0:
                continue
            else:
                gf_tree.prune(geneleaves_sampled, preserve_branch_length=True)

            with open(os.path.join(self.outfolder, fam), "w") as f:
                f.write(gf_tree.write(format=1))


if __name__ == "__main__":

    args = sys.argv[1:]
    if len(args) <= 2:
        print("Incorrect usage. Please read the manual. The usual way to run this script is:")
        print("python sampLyon.py N /Output_folder /Sampling_folder")
    else:
        N, infolder, outfolder = args

        SM = SamplingManager(N,infolder,outfolder)

        if os.path.isdir(os.path.join(infolder,outfolder)):
            pass
        else:
            os.mkdir(os.path.join(infolder,outfolder))

        SM.prune_species_tree()
        SM.prune_gene_trees()

