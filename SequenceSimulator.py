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

            line = f.readline().strip()
            if "(" not in line or line == ";":
                return None
            else:
                my_tree = ete3.Tree(line, format=1)

        tree = pyvolve.read_tree(tree=my_tree.write(format=5))
        partition = pyvolve.Partition(models=self.model, size=self.size)
        evolver = pyvolve.Evolver(tree=tree, partitions=partition)
        fasta_file = tree_file.split("/")[-1].replace("_wholetree.nwk", "_whole") + ".fasta"
        evolver(seqfile=os.path.join(sequences_folder, fasta_file), ratefile=None, infofile=None)

    def run_u(self, tree_file, sequences_folder):

        with open(tree_file) as f:
            line = f.readline().strip()
            if "(" not in line or line == ";":
                return None
            else:
                my_tree = ete3.Tree(line, format=1)

        root = my_tree.get_tree_root()
        root.name = "Root"

        # in this case we need to read the multipliers
        # First we apply the multipliers per family
        # Second, the multipliers per species tree branch

        gf_multiplier = self.gf_multipliers[tree_file.split("_")[1].split("/")[-1]]

        for node in my_tree.traverse():
            node.dist = node.dist * gf_multiplier * self.st_multipliers[node.name.split("_")[0]]


        tree = pyvolve.read_tree(tree=my_tree.write(format=5))
        partition = pyvolve.Partition(models=self.model, size=self.size)
        evolver = pyvolve.Evolver(tree=tree, partitions=partition)
        fasta_file = tree_file.split("/")[-1].replace("_wholetree.nwk", "_") +  "whole.fasta"
        evolver(seqfile=os.path.join(sequences_folder, fasta_file), ratefile=None, infofile=None)

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


    def obtain_rates_multipliers(self, gt_file, st_file):

        self.gf_multipliers = dict()
        self.st_multipliers = dict()

        with open(gt_file) as f:
            f.readline()
            for line in f:
                fm, m = line.strip().split("\t")
                self.gf_multipliers[fm] = float(m)

        with open(st_file) as f:
            f.readline()
            for line in f:
                clade, m = line.strip().split("\t")
                self.st_multipliers[clade] = float(m)

    def write_rates_sttree(self, whole_tree, rates_tree):

        with open(whole_tree) as f:
            whole_tree = ete3.Tree(f.readline().strip(), format=1)

        r = whole_tree.get_tree_root()
        r.name = "Root"

        for n in whole_tree.traverse():
            n.dist *= self.st_multipliers[n.name]

        with open(rates_tree, "w") as f:
            f.write(whole_tree.write(format=1))













