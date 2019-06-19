import pyvolve
import os
import ete3
import numpy
import AuxiliarFunctions as af



class SequenceSimulator():

    def __init__(self, parameters):

        self.parameters = parameters

        if self.parameters["SEED"] != 0:


            random.seed(parameters["SEED"])
            numpy.random.seed(parameters["SEED"])

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

        tree = pyvolve.read_tree(tree=my_tree.write(format=5), scale_tree = self.parameters["SCALING"])
        name_mapping = self.get_mapping_internal_names(tree, my_tree)
        partition = pyvolve.Partition(models=self.model, size=self.size)
        evolver = pyvolve.Evolver(tree=tree, partitions=partition)
        fasta_file = tree_file.split("/")[-1].replace("_completetree.nwk", "_complete") + ".fasta"
        evolver(seqfile=os.path.join(sequences_folder, fasta_file), ratefile=None, infofile=None, write_anc=True)

        # Correct the names
        self.correct_names(os.path.join(sequences_folder, fasta_file), name_mapping)

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

        tree = pyvolve.read_tree(tree=my_tree.write(format=5), scale_tree = self.parameters["SCALING"])
        name_mapping = self.get_mapping_internal_names(tree, my_tree)
        partition = pyvolve.Partition(models=self.model, size=self.size)
        evolver = pyvolve.Evolver(tree=tree, partitions=partition)
        fasta_file = tree_file.split("/")[-1].replace("_completetree.nwk", "_") +  "complete.fasta"
        evolver(seqfile=os.path.join(sequences_folder, fasta_file), ratefile=None, infofile=None, write_anc=True)
        # Correct the names
        self.correct_names(os.path.join(sequences_folder, fasta_file), name_mapping)

    def run_f(self, tree_file, gene_length, sequences_folder):

        if self.parameters["SEQUENCE"] != "codon":
            self.model = self.get_codon_model()

        with open(tree_file) as f:

            line = f.readline().strip()
            if "(" not in line or line == ";":
                self.simulate_single_sequence(line.replace(";",""),gene_length, tree_file, sequences_folder)
                return None
            else:
                my_tree = ete3.Tree(line, format=1)
                tree = pyvolve.read_tree(tree=my_tree.write(format=5), scale_tree = self.parameters["SCALING"])
                name_mapping = self.get_mapping_internal_names(tree, my_tree)
                partition = pyvolve.Partition(models=self.model, size=gene_length)
                evolver = pyvolve.Evolver(tree=tree, partitions=partition)
                fasta_file = tree_file.split("/")[-1].replace("_completetree.nwk", "_complete") + ".fasta"
                evolver(seqfile=os.path.join(sequences_folder, fasta_file), ratefile=None, infofile=None, write_anc=True)
                self.correct_names(os.path.join(sequences_folder, fasta_file), name_mapping)


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

    def write_rates_sttree(self, complete_tree, rates_tree):

        with open(complete_tree) as f:
            complete_tree = ete3.Tree(f.readline().strip(), format=1)
        r = complete_tree.get_tree_root()
        r.name = "Root"
        for n in complete_tree.traverse():
            n.dist *= self.st_multipliers[n.name]
        with open(rates_tree, "w") as f:
            f.write(complete_tree.write(format=1))

    def retrieve_sequences(self, name, gf, sequences_folder):

        for n,s in af.fasta_reader(os.path.join(sequences_folder, gf + "_complete.fasta")):
            if n[1:] == name:
                return s
        return None

    def retrieve_orientation(self, species, gene_name, lengths_folder):

        with open(os.path.join(lengths_folder, species + "_GENOME.tsv")) as f:
            f.readline()
            for line in f:
                h = line.strip().split("\t")
                orientation = h[2]
                gf = h[1]
                id = h[3]
                if gene_name == gf + "_" + id:
                    return orientation
        return None

    def simulate_single_sequence(self, name, gene_length, tree_file, sequences_folder):

        my_tree = "(A:1,B:1);".replace("A",name)
        tree = pyvolve.read_tree(tree=my_tree)
        partition = pyvolve.Partition(models=self.model, size=gene_length)
        evolver = pyvolve.Evolver(tree=tree, partitions=partition)

        fasta_file = tree_file.split("/")[-1].replace("_completetree.nwk", "_complete") + ".fasta"
        evolver(seqfile=os.path.join(sequences_folder, fasta_file), ratefile=None, infofile=None, write_anc=True)

        # Select single sequence

        entries = list()

        for n, v in af.fasta_reader(os.path.join(sequences_folder, fasta_file)):
            if n[1:] != name:
                continue
            else:
                entries.append((n,v))
        af.fasta_writer(os.path.join(sequences_folder, fasta_file), entries)


    def generate_intergenic_sequences(self, l):

        return("".join(numpy.random.choice(["A", "T", "C", "G"], l)))

    def get_mapping_internal_names(self, pytree, ettree):

        pyvolvemap = dict()

        def traverse(root):
            if root:
                if len(root.children) == 0:
                    return None
                else:
                    traverse(root.children[0])
                    traverse(root.children[1])
                    pyvolvemap[root.children[0].name + "+" + root.children[1].name] = root.name
        traverse(pytree)
        good_mapping = dict()

        etroot = ettree.get_tree_root().name
        good_mapping["myroot"] = etroot

        for n in ettree.traverse(strategy="postorder"):
            if not n.is_leaf():
                c1, c2 = n.get_children()
                n1 = c1.name + "+" + c2.name
                n2 = c2.name + "+" + c1.name
                if n1 in pyvolvemap:
                    good_mapping[pyvolvemap[n1]] = n.name
                    n.name = pyvolvemap[n1]
                if n2 in pyvolvemap:
                    good_mapping[pyvolvemap[n2]] = n.name
                    n.name = pyvolvemap[n2]


        return good_mapping

    def correct_names(self,fasta_file, good_mapping):

        entries = list()

        for n,v in af.fasta_reader(fasta_file):
            if "root" in n:
                entries.append((">"+good_mapping["myroot"], v))
            elif "internal" in n:
                entries.append((">"+good_mapping[n[1:]], v))
            else:
                entries.append((n, v))

        af.fasta_writer(fasta_file, entries)
















