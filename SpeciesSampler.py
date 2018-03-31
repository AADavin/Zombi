import AuxiliarFunctions as af
import argparse
import os
import random
import numpy

class SpeciesSampler():

    def __init__(self, experiment_folder):

        self.experiment_folder = experiment_folder
        self.output_folder = os.path.join(experiment_folder, "SAMPLE")

        if not os.path.isdir(self.output_folder):
            os.mkdir(self.output_folder)

        self.existing_folders = os.listdir(experiment_folder)
        self.species_tree_events = self.read_events(os.path.join(experiment_folder, "T/Events.tsv"))
        self.gene_family_folder = os.path.join(experiment_folder, "G/Gene_families")

    def read_events(self, events_file):

        events = list()

        with open(events_file) as f:
            f.readline()
            for line in f:
                events.append(tuple(line.strip().split("\t")))

        return events


    def prune_tree_events(self, species_to_sample):

        # Now we prune the events
        pruned_events = list()

        for t, e, nodes in self.species_tree_events:

            if e == "F":

                if nodes not in species_to_sample:

                    pruned_events.append((t, "E", nodes))

                else:

                    pruned_events.append((t, e, nodes))

            else:

                pruned_events.append((t, e, nodes))

        return pruned_events

    def prune_gene_events(self, gene_tree_events, species_to_sample):
        # Now we prune the events

        pruned_events = list()
        for t, e, nodes in gene_tree_events:
            if e == "F":
                if nodes.split(";")[0] not in species_to_sample:

                    pruned_events.append((t, "E", nodes))
                else:
                    pruned_events.append((t, e, nodes))
            else:
                pruned_events.append((t, e, nodes))

        return pruned_events


    def cut(self, species_sampled):

        cut_trees = list()
         # Now we just prune the tree

        pruned_events = self.prune_tree_events(species_sampled)
        w, p = af.generate_newick_trees(pruned_events)

        cut_trees.append(("SampleSpeciesTree.nwk", p))

        # p stores the pruned Species Tree

        # Now we cut the gene trees

        events_files = [os.path.join(self.gene_family_folder, x) for x in os.listdir(self.gene_family_folder)]

        for event_file in events_files:
            gf_events = self.read_events(event_file)
            pruned_events = self.prune_gene_events(gf_events, species_sampled)
            w1, sgt = af.generate_gene_tree(pruned_events)
            cut_trees.append((event_file.split("/")[-1], sgt))

        return cut_trees


    def mode_i(self, myfile):

        species_to_sample = set()

        with open(myfile) as f:
            for line in f:
                species_to_sample.add(line.strip())

        trees = self.cut(species_to_sample)
        self.write_sampled_trees(trees)

    def mode_r(self, value):

        # First we get the list of all the extant species

        species_alive = [x[2] for x in self.species_tree_events if x[1] == "F"]
        # Then we sample
        value = float(value)
        if value < 0 or value > 1:
            return 0
        total_species_sampled = int(value * len(species_alive))
        species_to_sample = random.sample(species_alive, total_species_sampled)
        trees = self.cut(species_to_sample)
        self.write_sampled_trees(trees)


    def mode_n(self, value):

        species_alive = [x[2] for x in self.species_tree_events if x[1] == "F"]
        # Then we sample
        species_to_sample = random.sample(species_alive, value)
        trees = self.cut(species_to_sample)
        self.write_sampled_trees(trees)


    def mode_w(self, input):

        myfile, n = input.split(";")

        species = []
        weights = []

        with open(myfile) as f:
            for line in f:
                sp, wt = line.strip().split("\t")
                species.append(sp)
                weights.append(wt)

        species_to_sample = numpy.random.choice(species, n, p=af.normalize(weights))
        trees = self.cut(species_to_sample)
        self.write_sampled_trees(trees)


    def write_sampled_trees(self, trees):

        for name, tree in trees:

            print("Writing %s" % name)
            with open(os.path.join(self.output_folder,name), "w") as f:
                f.write(tree)




if __name__ == "__main__":


    parser = argparse.ArgumentParser()

    parser.add_argument("mode", type=str, choices=["i","r","n","w"], help="Mode")
    parser.add_argument("input", type=str, help="Input for SpeciesSampler. If the mode is i, file with species to sample. If the mode is r, number between 0 and 1 with the probability of sampling a lineage. If the mode is n, the number of lineages to preserve. If the mode is w, .tsv file with all the species in the extant trees next to the probabilities of sampling each lineage")
    parser.add_argument("output", type=str, help="Name of the experiment folder")

    args = parser.parse_args()

    mode, input, experiment_folder =  args.mode, args.input, args.output

    ss = SpeciesSampler(experiment_folder)

    if mode == "i":

        ss.mode_i(input)

    elif mode == "r":

        ss.mode_r(input)

    elif mode == "n":

        ss.mode_n(input)

    elif mode == "w":


        ss.mode_w(input)



