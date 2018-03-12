from SpeciesTreeSimulator import SpeciesTreeGenerator
from GenomeSimulator import GenomeSimulator
import AuxiliarFunctions as af
import argparse
import os
import sys

class simuLyon():

    def __init__(self):

        self.tree_parameters = dict()
        self.genome_parameters = dict()
        self.sequence_parameters = dict()


    def T(self, parameters_file, experiment_folder, advanced_mode):

        tree_folder = os.path.join(experiment_folder, "T")

        parameters = af.prepare_species_tree_parameters(af.read_parameters(parameters_file))
        stg = SpeciesTreeGenerator(parameters)

        run_counter = 0
        success = False

        while success == False and run_counter <= 50:
            run_counter+=1
            print("Computing Species Tree. Trial number %s" % str(run_counter))
            if advanced_mode == "0":
                success = stg.run()
            if advanced_mode == "a":
                success = stg.run()
            if advanced_mode == "b":
                success = stg.run()
            if advanced_mode == "i":
                success = stg.run()
            if advanced_mode == "p":
                success = stg.run_p()

        if run_counter >= 50:
            print("Aborting computation of the Species Tree. Please use other speciation and extinction rates!")
            return 0

        whole_tree_file = os.path.join(tree_folder, "WholeTree.nwk")
        extant_tree_file = os.path.join(tree_folder, "ExtantTree.nwk")
        events_file = os.path.join(tree_folder, "Events.tsv")
        stg.generate_whole_tree()
        stg.generate_extant_tree()
        stg.write_whole_tree(whole_tree_file)
        stg.write_extant_tree(extant_tree_file)
        stg.write_events_file(events_file)


    def G(self,parameters_file, experiment_folder):

        parameters = af.prepare_genome_parameters(af.read_parameters(parameters_file))
        events_file = os.path.join(experiment_folder, "T/Events.tsv")
        genome_folder = os.path.join(experiment_folder, "G")

        genomes_folder = os.path.join(genome_folder, "Genomes")
        gene_families_folder = os.path.join(genome_folder, "Gene_families")
        gene_trees_folder = os.path.join(genome_folder, "Gene_trees")
        events_per_branch_folder = os.path.join(genome_folder, "Events_per_branch")
        profiles_folder = os.path.join(genome_folder, "Profiles")

        gss = GenomeSimulator(parameters, events_file)

        gss.run()

        print("Writing Genomes")
        gss.write_genomes(genomes_folder)
        print("Writing Profiles")
        gss.write_profiles(profiles_folder)
        print("Writing Gene Families")
        gss.write_gene_family_events(gene_families_folder)
        print("Writing Events Per Branch")
        gss.write_events_per_branch(events_per_branch_folder)
        print("Writing Gene Trees")
        gss.write_gene_trees(gene_trees_folder)


    def S(self, parameters_file, experiment_folder):

        import pyvolve
        from ete3 import Tree

        genome_folder = os.path.join(experiment_folder, "Genomes")
        gene_trees_folder = os.path.join(genome_folder, "GeneTrees")
        sequences_folder = os.path.join(experiment_folder, "Sequences")

        if not os.path.isdir(sequences_folder):
            os.mkdir(sequences_folder)

        self.read_parameters(parameters_file, self.sequence_parameters)
        size = int(self.sequence_parameters["SEQUENCE_SIZE"])
        sequence = self.sequence_parameters["SEQUENCE"]
        sequence_type = "The type of sequence be either 'nucleotide', 'amino-acid' or 'codon'"
        assert sequence in ['nucleotide', 'amino-acid', 'codon'], sequence_type

        if sequence == 'nucleotide':
            nucleotides = ['A', 'C', 'G', 'T']
            state_freqs = []
            custom_mu = {}

            for source in nucleotides:
                state_freqs.append(float(self.sequence_parameters[source]))
                for target in nucleotides:
                    if source != target:
                        pair = source + target
                        custom_mu[pair] = float(self.sequence_parameters[pair])
            assert abs(sum(state_freqs) - 1) < 1e-6, "Equilibrium frequencies of nucleotides must sum to 1.0"

            model = pyvolve.Model("nucleotide", {"mu": custom_mu, "state_freqs": state_freqs})
        elif sequence == 'amino-acid':

            model = pyvolve.Model(self.sequence_parameters['AA_MODEL'])
        else:
            codon_params = {}
            for param in ["ALPHA", "BETA", "KAPPA"]:
                codon_params[param.lower()] = float(self.sequence_parameters[param])

            model = pyvolve.Model(self.sequence_parameters['CODON_MODEL'], codon_params, neutral_scaling=True)

        for tree_file in os.listdir(gene_trees_folder):
            tree_path = os.path.join(gene_trees_folder, tree_file)
            ete_tree = Tree(tree_path)
            if len(ete_tree) != 1:
                tree = pyvolve.read_tree(tree=ete_tree.write(format=5))
                partition = pyvolve.Partition(models=model, size=size)
                evolver = pyvolve.Evolver(tree=tree, partitions=partition)
                fasta_file = tree_file.replace(".txt", "_") + sequence + ".fasta"
                evolver(seqfile=os.path.join(sequences_folder, fasta_file), ratefile=None, infofile=None)


if __name__ == "__main__":


    parser = argparse.ArgumentParser()
    parser.add_argument("mode", type=str, choices=["T","Ti","Ta","Tb","Tp","S","G"], help="Mode")
    parser.add_argument("params",  type=str, help="Parameters file")
    parser.add_argument("output", type=str, help="Name of the experiment folder")

    args = parser.parse_args()

    mode, parameters_file, experiment_folder =  args.mode, args.params, args.output

    if len(mode) == 1:
        main_mode = mode[0]
        advanced_mode = "0"

    elif len(mode) == 2:
        main_mode = mode[0]
        advanced_mode = mode[1]

    else:
        print ("Incorrect value for mode")

    SL = simuLyon()

    if main_mode == "T":

        if not os.path.isdir(experiment_folder):
            os.mkdir(experiment_folder)

        if not os.path.isdir(os.path.join(experiment_folder,"T")):
            os.mkdir(os.path.join(experiment_folder,"T"))

        SL.T(parameters_file, experiment_folder, advanced_mode)

    elif main_mode == "G":

        genome_folder = os.path.join(experiment_folder, "G")

        if not os.path.isdir(genome_folder):
            os.mkdir(genome_folder)

        SL.G(parameters_file, experiment_folder, advanced_mode)

    elif main_mode == "S":

        SL.obtain_sequences(parameters_file, experiment_folder, advanced_mode)


    else:
        print("Incorrect usage. Please select a mode: T, G or S")

