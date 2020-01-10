import pyvolve
import os
import ete3
import numpy
import random
import itertools
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

        gf_multiplier = self.gf_multipliers[tree_file.split("_")[-2].split("/")[-1]]

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
                
                
    def run_s(self, tree_file, gene_length, sequences_folder):

        # 
        
        pass


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
        
    def _read_events_file(self, events_file):

        events = list()
        with open(events_file) as f:
            f.readline()
            for line in f:
                handle = line.strip().split("\t")
                if handle[1] == "SE" or handle[1] == "SS":
                    continue
                events.append(handle)
        return events
    
    def get_time_to_next_event(self, n, events):

        total = 0.0
        for __ in range(n):
            total += sum(events)

        if total == 0:
            return 1000000000000000 # We sent an arbitrarily big number. Probably not the most elegant thing to do
        else:
            time = numpy.random.exponential(1/total)
            return time
    
    def simulate_shifts(self, events_file):
        
        sr = float(self.parameters["SHIFT_SUBSTITUTION_RATE"])
        cats = int(self.parameters["SHIFT_CATEGORIES"])
        hbase_rate = self.parameters["BASE_RATE"]
       
        # We get the categories for speciations
        
        distribution, value = hbase_rate.split(":")
        value = float(value)
        
        if distribution == "g":
            substitution_rates = af.discretize(value,cats,"gamma")        
        elif distribution == "l":
            substitution_rates = af.discretize(value,cats,"lognorm")        
        else:
            print("Unrecognized distribution. Please, use g or l")
            return False
        
        self.substitution_rates = substitution_rates
        
        self.tree_events = self._read_events_file(events_file)
        self.shift_events = list()
        
        current_species_tree_event = 0
        current_time = 0.0
        all_species_tree_events = len(self.tree_events)
        
        self.category_position = dict()
        self.category_position["Root"] = int(len(substitution_rates)/2)     
        
        self.branchwise_rates = dict()
        self.branchwise_rates["Root"] = list()
        self.branchwise_rates["Root"].append((0, "S", substitution_rates[self.category_position["Root"]]))
        
        # Second, we compute the time to the next event:
        
        elapsed_time = 0.0
        self.active_genomes = set()
        self.active_genomes.add("Root")

        while current_species_tree_event < all_species_tree_events:

            time_of_next_species_tree_event, event, nodes = self.tree_events[current_species_tree_event]
            time_of_next_species_tree_event = float(time_of_next_species_tree_event)

            if self.parameters["VERBOSE"] == 1:
                print("Simulating shifts. Time %s" % str(current_time))
            time_to_next_genome_event = self.get_time_to_next_event(len(self.active_genomes), [sr])
            elapsed_time = float(current_time) - elapsed_time
            if time_to_next_genome_event + current_time >= float(time_of_next_species_tree_event):
                current_species_tree_event +=1
                current_time = time_of_next_species_tree_event
                if event == "S":
                    sp,c1,c2 = nodes.split(";")
                    # First we keep track of the active and inactive genomes
                    self.active_genomes.discard(sp)
                    self.active_genomes.add(c1)
                    self.active_genomes.add(c2)
                    
                    self.category_position[c1] = self.category_position[sp]
                    self.category_position[c2] = self.category_position[sp]
                    
                    
                    self.branchwise_rates[sp].append((current_time, "ES", substitution_rates[self.category_position[sp]]))
                    
                    self.branchwise_rates[c1] = list()
                    self.branchwise_rates[c1].append((current_time, "S", substitution_rates[self.category_position[c1]]))
                    
                    self.branchwise_rates[c2] = list()
                    self.branchwise_rates[c2].append((current_time, "S", substitution_rates[self.category_position[c2]]))
                                        
                elif event == "E":                    
                    self.active_genomes.discard(nodes)
                    self.branchwise_rates[nodes].append((current_time, "E", substitution_rates[self.category_position[nodes]]))
                                        
                elif event == "F":                    
                    self.branchwise_rates[nodes].append((current_time, "F", substitution_rates[self.category_position[nodes]]))
                    
            else:
                
                current_time += time_to_next_genome_event
                # A shift event occurs in a randomly selected lineage
                # First we select that lineage
                
                lineage = random.choice(list(self.active_genomes))
                
                # The substitution rate changes                
                    
                cat = self.category_position[lineage]                  
                oldcat = cat
                if cat == cats-1:         # Meaning that we are in the border
                    direction = numpy.random.choice([-1,0])
                elif cat == 0:
                    direction = numpy.random.choice([0,1])            
                else:
                    direction = numpy.random.choice([-1,1], p = [0.5,0.5])

                p = self.category_position[lineage]                    
                self.category_position[lineage] = p + direction

                new_rate = substitution_rates[self.category_position[lineage]]                                        
                
                self.shift_events.append((current_time, "SR",  lineage + ";" + str(substitution_rates[oldcat]) + "->" + str(new_rate)))
                self.branchwise_rates[lineage].append((current_time, "SS", substitution_rates[self.category_position[lineage]]))
                
              
    def write_events(self, events_file):

        header = ["TIME","EVENT","NODES"]
        header = "\t".join(map(str, header)) + "\n"

        with open(events_file, "w") as f:
            f.write(header)
            for item in self.shift_events:
                line = "\t".join(map(str,item)) + "\n"
                f.write(line)
                
    def write_substitution_scaled_stree(self, complete_tree, extant_tree, substitution_scaled_complete_tree_file, 
                                        substitution_scaled_extant_tree_file, branchwise_file):

        with open(complete_tree) as f:
            complete_tree = ete3.Tree(f.readline().strip(), format=1)
        
        r = complete_tree.get_tree_root()
        r.name = "Root"
        
        
        with open(extant_tree) as f:            
            extant_tree= ete3.Tree(f.readline().strip(), format=1)
            er = extant_tree.get_tree_root()
        
        extant_sps = {x.name for x in extant_tree.get_leaves()}
        
        # We transform the branchwise rates into intervals:
        
        self.eff_multiplier = dict()
        
        for node, vls in self.branchwise_rates.items():
        
            self.eff_multiplier[node] = 0
            
            # There are as many intervals as n - 1
            
            # We get the total time the branch exists (first event and last event)
            
            t1, *_ = vls[0]
            t2, *_ = vls[-1]
            
            tt = float(t2) - float(t1)            
       
            for vl1, vl2 in zip(vls, vls[1:]):             
                t1, e1, sr1 = vl1
                t2, e2, sr2 = vl2
                t = float(t2 - t1) / tt                
                self.eff_multiplier[node] += t * float(sr1)
                  
        for n in complete_tree.traverse():
            n.dist *= self.eff_multiplier[n.name]
            
        with open(substitution_scaled_complete_tree_file, "w") as f:
            f.write(complete_tree.write(format=1))  
        
        with open(substitution_scaled_extant_tree_file, "w") as f:
        
            f.write(self.quick_pruner(complete_tree, extant_sps, er.name).write(format=1, format_root_node = True))
        
            
        with open(branchwise_file, "w") as f:
            for node, vls in self.branchwise_rates.items(): 
                line = node + "\t" + "\t".join([";".join([str(x) for x in vl]) for vl in vls]) + "\n"
                f.write(line)
                
    def quick_pruner(self, complete_tree, extant_sps, initial_node):       
        
        initial_node = complete_tree&initial_node
                
        # Now, we need to assign to every branch a number (0:dead, 1:preserved, 2:extant)
        
        for n in initial_node.iter_descendants("postorder"):
            if n.is_leaf(): 
                if n.name in extant_sps:
                    n.add_feature("state",2)
                else:                    
                    n.add_feature("state",0)
                
            else:                
                
                c1,c2 = n.get_children()
                
                if c1.state == 2 and c2.state == 2:
                    n.add_feature("state", 2)
                if c1.state == 1 and c2.state == 1:
                    n.add_feature("state", 2)
                if c1.state == 0 and c2.state == 0:
                    n.add_feature("state", 0)
                if (c1.state == 1 and c2.state == 0) or (c2.state == 1 and c1.state == 0):
                    n.add_feature("state", 1)
                if (c1.state == 2 and c2.state == 0) or (c2.state == 2 and c1.state == 0):
                    n.add_feature("state", 1)
                if (c1.state == 2 and c2.state == 1) or (c2.state == 2 and c1.state == 1):
                    n.add_feature("state", 2)
        
        initial_node.state = 2 
        
        # Now we create modify the extant tree
        
        carrying = 0
        n2dist = dict()
        
        for n in initial_node.traverse("preorder"):            
            if n.state == 1:
                carrying += n.dist                
                n.delete(prevent_nondicotomic=True, preserve_branch_length=False)
            elif n.state == 2:
                n2dist[n.name] = n.dist + carrying                  
                carrying = 0                  
            elif n.state == 0:
                n.delete(prevent_nondicotomic=True, preserve_branch_length=True)
        
        for n in initial_node.traverse("preorder"):
            n.dist = n2dist[n.name]
        
        return initial_node                
                
    def write_effective_gtree(self, complete_gtree, events_gtree):
        
        # Parser of the effective lengths
        
        def parse_eff(eff):            
            it = iter(eff)
            for x in it:
                yield (x, next(it), next(it))
        
        # We read the tree
        
        with open(complete_gtree) as f:
            line = f.readline().strip()
            if "(" not in line or line == ";":
                # the tree is too small
                return None
            else:
                gtree = ete3.Tree(line, format=1)
        
        # We get all the nodes
        
        #all_nodes = {n.name:0 for n in gtree.traverse()}
        
        # We read the events and the beginning point of a node and the ending
        all_nodes = dict()
        
        with open(events_gtree) as f:              
            
            vls = f.readlines()[1:]            
            vls = [x.strip().split("\t") for x in vls if x.split("\t")[1] in ["S","O","D","T","L","E","F"]]
            
            for t, event, nodes in vls:                  
                
                if event == "O":
                    all_nodes[nodes.strip() + "_1"] = [t]                    
                else:
                    nodes = nodes.split(";")                
                    it = iter(nodes)                
                    for x in it:
                        mnode = x + "_" + next(it)
                        if mnode not in all_nodes:
                            all_nodes[mnode] = [t]
                        else:
                            all_nodes[mnode].append(t)
        node2eff = dict()                
        for node, time in all_nodes.items():
            
            eff = []
            o_t, e_t = time # Origin time, ending time
            o_t, e_t = float(o_t), float(e_t)
            tt = e_t - o_t
            sp_n = node.split("_")[0]    
                        
            # First we get the first substitution rate
            
            vls = self.branchwise_rates[node.split("_")[0]]       
            
            origin_found = False
            ending_found = False
                        
            for vl in vls:                                             
                t, e, sr = vl                
                if origin_found == False and float(t) <= o_t:
                    eff.append(o_t)
                    eff.append(sr)
                    origin_found = True   
                    
                if origin_found == True and e_t <= float(t):                    
                    eff.append(e_t)      
                    break
                    
                elif origin_found == True and e == "SS":
                    # Switch in the substitution rate
                    eff.append(t)
                    eff.append(t)
                    eff.append(sr)
                    
                    
            #print(sp_n, o_t, e_t, self.branchwise_rates[sp_n], eff)
            t_eff = 0
            for x in parse_eff(eff):
                
                t1, cat, t2 = x
                t_eff += ((t2 - t1) / tt) * cat
        
            node2eff[node] = t_eff
        # We multiply the tree
        
        for n in gtree.traverse():
            if n.name == "Root":
                name = "Root_1"
            else:
                name = n.name
            n.dist *= node2eff[name]
                
        return gtree.write(format=1, format_root_node=True)
        
        
    def write_categories(self, categories_file):        
        
        with open(categories_file, "w") as f:

            f.write("SUBSTITUTION_RATE_CATEGORIES\n")
            f.write("\t".join(map(str,self.substitution_rates))+"\n")

