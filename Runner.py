from SpeciesTreeSimulator import SpeciesTreeGenerator
from GenomeSimulator import GenomeSimulator
import AuxiliarFunctions as af
import os
import sys

class simuLyon():

    def __init__(self):

        self.tree_parameters = dict()
        self.genome_parameters = dict()
        self.sequence_parameters = dict()


    def T(self, parameters_file, experiment_folder):

        parameters = af.prepare_species_tree_parameters(af.read_parameters(parameters_file))
        stg = SpeciesTreeGenerator(parameters)

        stg.run()

        whole_tree_file = os.path.join(experiment_folder,"WholeTree.newk")
        extant_tree_file = os.path.join(experiment_folder, "ExtantTree.newk")
        events_file = os.path.join(experiment_folder, "Events.tsv")

        stg.write_whole_tree(whole_tree_file)
        stg.write_extant_tree(extant_tree_file)
        stg.write_events_file(events_file)


    def F(self,parameters_file, experiment_folder):
        pass

    def G(self,parameters_file, experiment_folder):

        parameters = af.prepare_genome_parameters(af.read_parameters(parameters_file))
        events_file = os.path.join(experiment_folder, "Events.tsv")

        genome_folder = os.path.join(experiment_folder, "Genomes")
        gene_families_folder = os.path.join(experiment_folder, "Gene_families")
        events_per_branch_folder = os.path.join(experiment_folder, "Events_per_branch")

        gss = GenomeSimulator(parameters, events_file)

        gss.run()

        print("Writing Genomes")
        gss.write_genomes(genome_folder)
        print("Writing Gene Families")
        gss.write_gene_family_events(gene_families_folder)
        #gss.write_gene_trees("/Users/adriandavin/Desktop/Bioinformatics/SimuLyon/Cedric/Gillespie/All/Gene_trees/")
        print("Writing Events Per Branch")
        gss.write_events_per_branch(events_per_branch_folder)



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

    args = sys.argv[1:]

    if len(args) != 3:
        print("Incorrect usage. Please read the manual. The usual way to run this script is:")
        print("python simuLyon.py T Parameters_file.tsv /Output_folder")
    else:
        mode, parameters_file, experiment_folder = args

        SL = simuLyon()

        if mode == "T":

            if not os.path.isdir(experiment_folder):
                os.mkdir(experiment_folder)

            SL.T(parameters_file, experiment_folder)

        elif mode == "G":

            genome_folder = os.path.join(experiment_folder, "Genomes")

            if not os.path.isdir(genome_folder):
                os.mkdir(genome_folder)

            SL.G(parameters_file, experiment_folder)

        elif mode == "F":

            SL.F()

        elif mode == "S":

            SL.obtain_sequences(parameters_file, experiment_folder)

        else:
            print("Incorrect usage. Please select a mode: T, G, F or S")

