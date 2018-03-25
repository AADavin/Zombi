import pyvolve
import os
import ete3
import AuxiliarFunctions as af

class SequenceSimulator():

    def __init__(self, parameters):

        self.parameters = parameters
        self.size = self.parameters["SEQUENCE_SIZE"]
        self.sequence = self.parameters["SEQUENCE"]

        sequence_type = "The type of sequence be either 'nucleotide', 'amino-acid' or 'codon'"
        assert self.sequence in ['nucleotide', 'amino-acid', 'codon'], sequence_type

        if self.sequence == 'nucleotide':
            self.model = self.get_nucleotide_model()
        elif self.sequence == 'amino-acid':
            self.model = self.get_aminoacid_model()
        elif self.sequence == "codon":
            self.model = self.get_codon_model()

    def run(self, tree_file, sequences_folder):

        with open(tree_file) as f:
            my_tree = ete3.Tree(f.readline().strip(), format=1)

        tree = pyvolve.read_tree(tree=my_tree.write(format=5))
        partition = pyvolve.Partition(models=self.model, size=self.size)
        evolver = pyvolve.Evolver(tree=tree, partitions=partition)
        fasta_file = tree_file.split("/")[-1].replace("_wholetree.nwk", "_") + self.sequence + ".fasta"
        evolver(seqfile=os.path.join(sequences_folder, fasta_file), ratefile=None, infofile=None)

    def run_b(self, tree_file, sequences_folder):

        with open(tree_file) as f:
            my_tree = ete3.Tree(f.readline().strip(), format=1)

        root = my_tree.get_tree_root()
        root.name = "Root"

        for node in my_tree.traverse():
            node.dist = node.dist * self.branch_rates[node.name.split("_")[0]]

        tree = pyvolve.read_tree(tree=my_tree.write(format=5))
        partition = pyvolve.Partition(models=self.model, size=self.size)
        evolver = pyvolve.Evolver(tree=tree, partitions=partition)
        fasta_file = tree_file.split("/")[-1].replace("_wholetree.nwk", "_") + self.sequence + ".fasta"
        evolver(seqfile=os.path.join(sequences_folder, fasta_file), ratefile=None, infofile=None)

    def obtain_rates_multiplier(self, tree_file):

        with open(tree_file) as f:
            mytree = ete3.Tree(f.readline().strip(), format=1)


        root = mytree.get_tree_root()
        root.name = "Root"

        self.branch_rates = dict()

        for node in mytree.traverse():
            self.branch_rates[node.name] = af.obtain_value(self.parameters["RATE_MULTIPLIERS"])

    def write_rates_tree(self, inputtree_file, output_tree):

        with open(inputtree_file) as f:
            mytree = ete3.Tree(f.readline().strip(), format=1)


        root = mytree.get_tree_root()
        root.name = "Root"

        for node in mytree.traverse():
            node.dist = self.branch_rates[node.name] * node.dist

        with open(output_tree, "w") as f:
            f.write(mytree.write(format=1))


    def write_rates(self, rates_file):

        with open(rates_file, "w") as f:

            line = "\t".join(["lineage", "multiplier"]) + "\n"
            f.write(line)

            for lineage, value in self.branch_rates.items():

                line = "\t".join(map(str,[lineage, value])) + "\n"
                f.write(line)


    def get_nucleotide_model(self):

        nucleotides = ['A', 'C', 'G', 'T']
        state_freqs = []
        custom_mu = {}

        for source in nucleotides:
            state_freqs.append(float(self.parameters[source]))
            for target in nucleotides:
                if source != target:
                    pair = source + target
                    custom_mu[pair] = float(self.parameters[pair])

        assert abs(sum(state_freqs) - 1) < 1e-6, "Equilibrium frequencies of nucleotides must sum to 1.0"
        return pyvolve.Model("nucleotide", {"mu": custom_mu, "state_freqs": state_freqs})

    def get_aminoacid_model(self):

        return pyvolve.Model(self.parameters['AA_MODEL'])

    def get_codon_model(self):
        codon_params = {}
        for param in ["ALPHA", "BETA", "KAPPA"]:
            codon_params[param.lower()] = float(self.parameters[param])
        return pyvolve.Model(self.parameters['CODON_MODEL'], codon_params, neutral_scaling=True)





