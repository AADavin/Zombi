import AuxiliarFunctions as af
import ete3
import numpy
import copy
import random
import os
import networkx as nx
import ReconciledTree as RT
from itertools import cycle
from functools import reduce

# from GenomeClasses import GeneFamily, Gene, Intergene, CircularChromosome, LinearChromosome, Genome

class GenomeSimulator():

    def __init__(self, parameters, events_file):

        self.parameters = parameters

        mseed = self.parameters["SEED"]
        if mseed != 0:
            random.seed(parameters["SEED"])
            numpy.random.seed(parameters["SEED"])

        self.tree_events = self._read_events_file(events_file)
        self.distances_to_start = self._read_distances_to_start(events_file) # Only useful when computing assortative transfers
        self.complete_tree = self._read_tree(events_file.replace("Events.tsv", "CompleteTree.nwk"))

        self.all_genomes = dict()
        self.all_gene_families = dict()

        self.gene_families_counter = 0
        self.active_genomes = set()

        if self.parameters["RATE_FILE"] != "False":
            if self.parameters["SCALE_RATES"] == "True":
                self.crown_length = self._read_crown_length(events_file.replace("Events", "Lengths"))
                self.empirical_rates = af.read_empirical_rates(rates_file=self.parameters["RATE_FILE"], scale_rates= self.crown_length)
            else:
                self.empirical_rates = af.read_empirical_rates(rates_file=self.parameters["RATE_FILE"])


    def write_genomes(self, genome_folder, intergenic_sequences = False):

        if not os.path.isdir(genome_folder):
            os.mkdir(genome_folder)

        for genome_name,genome in self.all_genomes.items():

            with open(os.path.join(genome_folder, genome_name + "_GENOME.tsv"), "w") as f:

                header = ["POSITION", "GENE_FAMILY", "ORIENTATION", "GENE_ID"]
                header = "\t".join(map(str, header)) + "\n"
                f.write(header)

                for chromosome in genome:
                    for index, gene in enumerate(chromosome):

                        line = [index, gene.gene_family, gene.orientation, gene.gene_id]
                        line = "\t".join(map(str,line)) +"\n"
                        f.write(line)

            if intergenic_sequences == True:

                with open(os.path.join(genome_folder, genome_name + "_LENGTHS.tsv"), "w") as f:

                    header = ["POSITION", "IDENTITY", "LENGTH"]
                    header = "\t".join(map(str, header)) + "\n"
                    f.write(header)

                    for chromosome in genome:
                        i = 0
                        for j, gene in enumerate(chromosome.genes):
                            line = [i, "G(" + str(gene.gene_family) + "_" + str(gene.gene_id) + ")", str(gene.length)]
                            line = "\t".join(map(str, line)) + "\n"
                            f.write(line)
                            i += 1
                            line = [i, "I", str(chromosome.intergenes[j].length)]
                            line = "\t".join(map(str, line)) + "\n"
                            f.write(line)
                            i += 1



    def write_gene_family_lengths(self, genome_folder):

        with open(os.path.join(genome_folder, "GeneFamily_lengths.tsv"), "w") as f:
            header = ["GENE_FAMILY", "LENGTH"]
            header = "\t".join(map(str, header)) + "\n"
            f.write(header)

            for gene_family_name, gene_family in self.all_gene_families.items():
                line = "\t".join([gene_family_name, str(gene_family.length)]) + "\n"
                f.write(line)


    def write_gene_family_events(self, gene_family_events_folder):

        if not os.path.isdir(gene_family_events_folder):
            os.mkdir(gene_family_events_folder)

        for gene_family_name, gene_family in self.all_gene_families.items():

            with open(os.path.join(gene_family_events_folder, gene_family_name + "_events.tsv"),"w") as f:

                header = ["TIME","EVENT","NODES"]
                header = "\t".join(map(str, header)) + "\n"

                f.write(header)

                for time, event, nodes in gene_family.events:
                    line = [time, event, nodes]
                    line = "\t".join(map(str, line)) + "\n"
                    f.write(line)

    def write_gene_trees(self, gene_tree_folder, gene_trees = True, reconciliations = False):

        if not os.path.isdir(gene_tree_folder):
            os.mkdir(gene_tree_folder)

        for gene_family_name, gene_family in self.all_gene_families.items():

            complete_tree, pruned_tree, rec_tree = gene_family.generate_tree()

            if gene_trees == True:

                with open(os.path.join(gene_tree_folder, gene_family_name + "_completetree.nwk"), "w") as f:
                    f.write(complete_tree)
                if pruned_tree != None:
                    with open(os.path.join(gene_tree_folder, gene_family_name + "_prunedtree.nwk"), "w") as f:
                        f.write(pruned_tree)
            if reconciliations == True:
                with open(os.path.join(gene_tree_folder, gene_family_name + "_rec.xml"), "w") as f:
                    f.write(rec_tree)

    def write_events_per_branch(self, events_per_branch_folder, scale, scaled_file, events_file):

        if not os.path.isdir(events_per_branch_folder):
            os.mkdir(events_per_branch_folder)

        events_per_branch = dict()

        for gene_family_name, gene_family in self.all_gene_families.items():

            for time, event, nodes in gene_family.events:

                name = nodes.split(";")[0]

                if name not in events_per_branch:
                    events_per_branch[name] = list()

                if event == "S" or event == "E" or event == "F":
                    continue

                elif event == "T":
                    donor = name
                    recipient = nodes.split(";")[4]

                    handle = nodes.split(";")
                    gene_names = list(map(lambda x: gene_family_name + "_" + x, [handle[1], handle[3], handle[5]]))
                    new_nodes = ";".join([donor, gene_names[0], donor, gene_names[1], recipient, gene_names[2]])


                    if donor not in events_per_branch:
                        events_per_branch[donor] = list()

                    events_per_branch[donor].append((time, "LT", new_nodes))

                    if recipient not in events_per_branch:
                        events_per_branch[recipient] = list()

                    events_per_branch[recipient].append((time, "AT", new_nodes))

                elif event == "D":

                    handle = nodes.split(";")
                    new_nodes = ";".join(map(lambda x: gene_family_name + "_" + x,[handle[1],handle[3],handle[5]]))
                    events_per_branch[name].append((time, event, new_nodes))

                else:

                    gene_id = nodes.split(";")[-1]
                    events_per_branch[name].append((time, event, gene_family_name + "_" + gene_id))


        for name, events in events_per_branch.items():

            with open(os.path.join(events_per_branch_folder, name + "_branchevents.tsv"), "w") as f:

                header = ["TIME", "EVENT", "NODES"]
                header = "\t".join(map(str, header)) + "\n"

                f.write(header)

                for time, event, nodes in sorted(events, key = lambda x: float(x[0])):

                    line = [str(time), event, nodes]
                    line = "\t".join(line) + "\n"
                    f.write(line)
        
            if scale != 0: # Only working if Species Tree has been scaled too the same distance!

                # First I read where the root is:

                with open(scaled_file) as f:
                    f.readline()
                    for l in f:
                        t, event, nodes = l.strip().split("\t")
                        if float(t) == 0:
                            eroot = nodes.split(";")[0]


                # Second, I read the total length of the tree

                with open(events_file) as f:                
                    for l in f:                    
                        t, event, nodes = l.strip().split("\t")
                        node1 = nodes.split(";")[0]
                        if node1 == eroot:
                            beginning_time = float(t)
                        continue
                    t, event, nodes = l.strip().split("\t")

                    totaltime = float(t) - beginning_time

                    mfactor = scale / totaltime


                with open(os.path.join(events_per_branch_folder, name + "_brancheventsscaled.tsv"), "w") as f:

                    header = ["TIME", "EVENT", "NODES"]
                    header = "\t".join(map(str, header)) + "\n"

                    f.write(header)

                    for time, event, nodes in sorted(events, key = lambda x: float(x[0])):
                        time = (float(time) - beginning_time) * mfactor
                        line = [str(time), event, nodes]
                        line = "\t".join(line) + "\n"
                        f.write(line)

            
        
        
        
        
        

    def write_profiles(self, profiles_folder):

        if not os.path.isdir(profiles_folder):
            os.mkdir(profiles_folder)


        genome_names = [x for x in self.all_genomes.keys()]
        gene_family_names = [str(x) for x in self.all_gene_families.keys()]

        # For clarity, I start with Initial Genome

        genome_names[0], genome_names[1] = genome_names[1], genome_names[0]

        data = list()
        data.append(["GENOME"] + gene_family_names)


        for genome_name in genome_names:
            line = dict()

            genome = self.all_genomes[genome_name]

            for gene_family_name in gene_family_names:
                if gene_family_name not in line:
                    line[gene_family_name] = 0

            for chromosome in genome:
                for index, gene in enumerate(chromosome):
                    line[gene.gene_family] += 1

            data.append([genome_name] + [str(line[fm]) for fm in gene_family_names])

        # We will transpose the data

        mlenx = len(data)
        mleny = len(data[0])

        with open(os.path.join(profiles_folder, "Profiles.tsv"), "w") as f:
            for i in range(mleny):
                line = list()
                for j in range(mlenx):
                    line.append(str(data[j][i]))
                line = "\t".join(line) + "\n"
                f.write(line)


    def write_interactomes(self, genome_folder):

        if not os.path.isdir(genome_folder):
            os.mkdir(genome_folder)

        for genome_name, genome in self.all_genomes.items():

            if not hasattr(genome, "interactome"):
                continue

            with open(os.path.join(genome_folder, genome_name + "_INTERACTOME.tsv"), "w") as f:

                header = ["GENE_1", "GENE_2"]
                header = "\t".join(map(str, header)) + "\n"
                f.write(header)

                for g1, g2 in genome.interactome.edges:
                    f.write("\t".join([g1, g2]) + "\n")

    def write_family_rates(self, genome_folder):

        with open(os.path.join(genome_folder, "Family_rates.tsv"), "w") as f:
            header = ["GENE_FAMILY", "D", "T", "L"]
            header = "\t".join(map(str, header)) + "\n"
            f.write(header)

            for gene_family_name, gene_family in self.all_gene_families.items():

                d = gene_family.rates["DUPLICATION"]
                t = gene_family.rates["TRANSFER"]
                l = gene_family.rates["LOSS"]

                f.write("\t".join(map(str,[gene_family_name, d,t,l])) + "\n")


    def _read_events_file(self, events_file):

        events = list()
        with open(events_file) as f:
            f.readline()
            for line in f:
                handle = line.strip().split("\t")
                events.append(handle)
        return events

    def _read_distances_to_start(self, events_file):

        # This function could be fusion with the function above

        distances_to_start = dict()

        with open(events_file) as f:
            f.readline()
            for line in f:
                time, _, nodes = line.strip().split("\t")
                n = nodes.split(";")[0]
                distances_to_start[n] = float(time)
        return distances_to_start


    def _read_crown_length(self, length_file):

        with open(length_file) as f:

            cl = float(f.readlines()[-1].strip().split("\t")[-1])

        return cl

    def _read_tree(self, tree_file):

        with open(tree_file) as f:
            t = ete3.Tree(f.readline().strip(), format=1)
        return t

    def return_new_identifiers_for_segment(self, segment):

        new_identifiers = list()

        for gene in segment:
            gf = gene.gene_family
            new_id = self.all_gene_families[gf].obtain_new_gene_id()
            new_identifiers.append(new_id)

        return new_identifiers

    def fill_genome(self, intergenic_sequences = False, family_rates = False, interactome = False):

        genome = Genome()
        genome.species = "Root"
        time = 0

        initial_genome_size = self.parameters["INITIAL_GENOME_SIZE"].split(";")
        shape = "C"

        for n_genes in initial_genome_size:

            if shape == "L":
                chromosome = LinearChromosome()
                chromosome.shape = "L"
            elif shape == "C":
                chromosome = CircularChromosome()
                chromosome.shape = "C"

            if intergenic_sequences == True:
                chromosome.has_intergenes = True
                mean_length = int(self.parameters["INTERGENE_LENGTH"])
                intergene_lengths = [int(x * mean_length * int(n_genes)) for x in
                                     af.sample_from_dirichlet(int(n_genes))]

                for i in range(int(n_genes)):
                    intergenic_sequence = Intergene()
                    intergenic_sequence.length = intergene_lengths[i]
                    chromosome.intergenes.append(intergenic_sequence)

            for i in range(int(n_genes)):

                # We fill the chromosomes and we create also the gene families
                if family_rates == True and self.parameters["RATE_FILE"] == "False":
                        gene, gene_family = self.make_origination(genome.species, time, family_mode=True)

                elif family_rates == True and self.parameters["RATE_FILE"] != "False":

                        gene, gene_family = self.make_origination(genome.species, time, family_mode=True,
                                                                  empirical_rates=True)
                else:
                    gene, gene_family = self.make_origination(genome.species, time)
                initial_gene = copy.deepcopy(gene)
                initial_gene.species = "Initial"
                gene_family.genes.append(initial_gene)
                chromosome.genes.append(gene)
                self.all_gene_families[str(self.gene_families_counter)] = gene_family
                if intergenic_sequences == True:
                    gene.length = int(af.obtain_value(self.parameters["GENE_LENGTH"]))
            if intergenic_sequences == True:
                chromosome.obtain_flankings()
                chromosome.obtain_locations()
            genome.chromosomes.append(chromosome)

            if interactome == True:
                genome.create_interactome()


        return genome


    def run(self):

        d = af.obtain_value(self.parameters["DUPLICATION"])
        t = af.obtain_value(self.parameters["TRANSFER"])
        l = af.obtain_value(self.parameters["LOSS"])
        i = af.obtain_value(self.parameters["INVERSION"])
        c = af.obtain_value(self.parameters["TRANSPOSITION"])

        o = af.obtain_value(self.parameters["ORIGINATION"])

        # First we prepare the first genome

        genome = self.fill_genome()

        self.active_genomes.add(genome.species)
        self.all_genomes["Root"] = genome

        # We add the original genome too

        self.all_genomes["Initial"] = copy.deepcopy(genome)

        current_species_tree_event = 0
        current_time = 0.0
        all_species_tree_events = len(self.tree_events)
        # Second, we compute the time to the next event:

        elapsed_time = 0.0

        while current_species_tree_event < all_species_tree_events:

            time_of_next_species_tree_event, event, nodes = self.tree_events[current_species_tree_event]
            time_of_next_species_tree_event = float(time_of_next_species_tree_event)

            if self.parameters["VERBOSE"] == 1:
                print("Simulating genomes. Time %s" % str(current_time))

            time_to_next_genome_event = self.get_time_to_next_event(len(self.active_genomes), [d, t, l, i, c, o])

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

                    # Second, we speciate the genomes

                    genome_c1, genome_c2 = self.make_speciation(sp, c1, c2, current_time)

                    self.all_genomes[c1] = genome_c1
                    self.all_genomes[c2] = genome_c2

                elif event == "E":
                    self.make_extinction(nodes, current_time)
                    self.active_genomes.discard(nodes)

                elif event == "F":
                    self.make_end(current_time)
                    break

            else:

                current_time += time_to_next_genome_event
                self.evolve_genomes(d, t, l, i, c, o, current_time)

    def run_i(self):

        # Interactome mode

        d = af.obtain_value(self.parameters["DUPLICATION"])
        t = af.obtain_value(self.parameters["TRANSFER"])
        l = af.obtain_value(self.parameters["LOSS"])
        i = af.obtain_value(self.parameters["INVERSION"])
        c = af.obtain_value(self.parameters["TRANSPOSITION"])
        o = af.obtain_value(self.parameters["ORIGINATION"])
        rm = af.obtain_value(self.parameters["REMOVE"])
        rw = af.obtain_value(self.parameters["REWIRE"])

        # First we prepare the first genome

        genome = self.fill_genome(interactome = True)

        # We prepare to important dicts in this mode

        self.active_genomes.add(genome.species)
        self.all_genomes["Root"] = genome

        # We add the initial genome too

        self.all_genomes["Initial"] = copy.deepcopy(genome)

        current_species_tree_event = 0
        current_time = 0.0
        all_species_tree_events = len(self.tree_events)

        # Second, we compute the time to the next event:

        elapsed_time = 0.0

        while current_species_tree_event < all_species_tree_events:

            time_of_next_species_tree_event, event, nodes = self.tree_events[current_species_tree_event]
            time_of_next_species_tree_event = float(time_of_next_species_tree_event)

            if self.parameters["VERBOSE"] == 1:
                print("Simulating genomes. Time %s" % str(current_time))

            time_to_next_genome_event = self.get_time_to_next_event(len(self.active_genomes), [d, t, l, i, c, o, rm, rw])

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

                    # Second, we speciate the genomes

                    genome_c1, genome_c2 = self.make_speciation(sp, c1, c2, current_time)

                    self.all_genomes[c1] = genome_c1
                    self.all_genomes[c2] = genome_c2


                elif event == "E":
                    self.make_extinction(nodes, current_time)
                    self.active_genomes.discard(nodes)

                elif event == "F":
                    self.make_end(current_time)
                    break

            else:

                current_time += time_to_next_genome_event
                self.evolve_genomes_i(d, t, l, i, c, o, rm, rw, current_time)

    def run_m(self):

        # First we prepare the first genome

        genome = self.fill_genome(family_rates=True)

        # We prepare to important dicts in this mode

        self.active_genomes.add(genome.species)
        self.all_genomes["Root"] = genome

        # We add the initial genome too

        self.all_genomes["Initial"] = copy.deepcopy(genome)

        current_species_tree_event = 0
        current_time = 0.0
        all_species_tree_events = len(self.tree_events)

        # Second, we compute the time to the next event:

        elapsed_time = 0.0


        while current_species_tree_event < all_species_tree_events:

            time_of_next_species_tree_event, event, nodes = self.tree_events[current_species_tree_event]
            time_of_next_species_tree_event = float(time_of_next_species_tree_event)

            if self.parameters["VERBOSE"] == 1:
                print("Simulating genomes. Time %s" % str(current_time))

            time_to_next_genome_event = self.get_time_to_next_event_family_mode()

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

                    # Second, we speciate the genomes

                    genome_c1, genome_c2 = self.make_speciation(sp, c1, c2, current_time)

                    self.all_genomes[c1] = genome_c1
                    self.all_genomes[c2] = genome_c2


                elif event == "E":
                    self.make_extinction(nodes, current_time)
                    self.active_genomes.discard(nodes)

                elif event == "F":
                    self.make_end(current_time)
                    break

            else:

                current_time += time_to_next_genome_event
                self.evolve_genomes_m( current_time)


    def run_u(self):

        genome = self.fill_genome()
        self.active_genomes.add(genome.species)
        self.all_genomes["Root"] = genome

        # We add the original genome too

        self.all_genomes["Initial"] = copy.deepcopy(genome)

        current_species_tree_event = 0
        current_time = 0.0
        all_species_tree_events = len(self.tree_events)

        elapsed_time = 0.0

        while current_species_tree_event < all_species_tree_events:

            time_of_next_species_tree_event, event, nodes = self.tree_events[current_species_tree_event]
            time_of_next_species_tree_event = float(time_of_next_species_tree_event)

            if self.parameters["VERBOSE"] == 1:
                print("Simulating genomes. Time %s" % str(current_time))

            time_to_next_genome_event = self.get_time_to_next_event_advanced_modes()

            elapsed_time = float(current_time) - elapsed_time
            if time_to_next_genome_event + current_time >= float(time_of_next_species_tree_event):
                current_species_tree_event += 1
                current_time = time_of_next_species_tree_event
                if event == "S":

                    sp, c1, c2 = nodes.split(";")

                    # First we keep track of the active and inactive genomes
                    self.active_genomes.discard(sp)
                    self.active_genomes.add(c1)
                    self.active_genomes.add(c2)

                    # Second, we speciate the genomes

                    genome_c1, genome_c2 = self.make_speciation(sp, c1, c2, current_time)
                    self.all_genomes[c1] = genome_c1
                    self.all_genomes[c2] = genome_c2

                elif event == "E":
                    self.make_extinction(nodes, current_time)
                    self.active_genomes.discard(nodes)

                elif event == "F":
                    self.make_end(current_time)
                    break

            else:

                current_time += time_to_next_genome_event
                self.advanced_evolve_genomes(current_time)

    def run_f(self):

        d = af.obtain_value(self.parameters["DUPLICATION"])
        t = af.obtain_value(self.parameters["TRANSFER"])
        l = af.obtain_value(self.parameters["LOSS"])
        i = af.obtain_value(self.parameters["INVERSION"])
        c = af.obtain_value(self.parameters["TRANSPOSITION"])
        o = af.obtain_value(self.parameters["ORIGINATION"])

        # First we prepare the first genome

        genome = self.fill_genome(intergenic_sequences=True)

        ## These two lines are important for this mode

        for chromosome in genome:
            chromosome.obtain_flankings()
            chromosome.obtain_locations()

        self.active_genomes.add(genome.species)
        self.all_genomes["Root"] = genome

        # We add the original genome too

        self.all_genomes["Initial"] = copy.deepcopy(genome)

        current_species_tree_event = 0
        current_time = 0.0
        all_species_tree_events = len(self.tree_events)
        # Second, we compute the time to the next event:

        elapsed_time = 0.0

        while current_species_tree_event < all_species_tree_events:

            time_of_next_species_tree_event, event, nodes = self.tree_events[current_species_tree_event]
            time_of_next_species_tree_event = float(time_of_next_species_tree_event)

            if self.parameters["VERBOSE"] == 1:
                print("Simulating genomes. Time %s" % str(current_time))

            time_to_next_genome_event = self.get_time_to_next_event(len(self.active_genomes), [d, t, l, i, c, o])

            elapsed_time = float(current_time) - elapsed_time

            if time_to_next_genome_event + current_time >= float(time_of_next_species_tree_event):

                current_species_tree_event += 1
                current_time = time_of_next_species_tree_event

                if event == "S":

                    sp, c1, c2 = nodes.split(";")

                    # First we keep track of the active and inactive genomes

                    self.active_genomes.discard(sp)
                    self.active_genomes.add(c1)
                    self.active_genomes.add(c2)

                    # Second, we speciate the genomes

                    genome_c1, genome_c2 = self.make_speciation(sp, c1, c2, current_time)

                    self.all_genomes[c1] = genome_c1
                    self.all_genomes[c2] = genome_c2

                elif event == "E":
                    self.make_extinction(nodes, current_time)
                    self.active_genomes.discard(nodes)

                elif event == "F":
                    self.make_end(current_time)
                    break

            else:

                current_time += time_to_next_genome_event
                self.advanced_evolve_genomes_f(d, t, l, i, c, o, current_time)

    def generate_new_rates(self):

        d = af.obtain_value(self.parameters["DUPLICATION"])
        t = af.obtain_value(self.parameters["TRANSFER"])
        l = af.obtain_value(self.parameters["LOSS"])
        i = af.obtain_value(self.parameters["INVERSION"])
        p = af.obtain_value(self.parameters["TRANSPOSITION"])
        o = af.obtain_value(self.parameters["ORIGINATION"])

        return d,t,l,i,p,o

    def generate_empirical_rates(self):

        mlen = len(self.empirical_rates)
        d,t,l = self.empirical_rates[numpy.random.randint(mlen)]

        return d,t,l

    def read_rates(self, rates_folder):

        self.branch_event_rates = dict()
        self.branch_extension_rates = dict()
        self.transfer_rates = dict()

        with open(os.path.join(rates_folder, "Event_rates.tsv")) as f:
            f.readline()
            for line in f:
                sp, d, t, l, i, c, o,  = line.split("\t")
                self.branch_event_rates[sp] = tuple([float(x) for x in (d, t, l, i, c, o)])

        with open(os.path.join(rates_folder, "Extension_rates.tsv")) as f:
            f.readline()
            for line in f:
                sp, d, t, l, i, c,  = line.split("\t")
                self.branch_extension_rates[sp] = tuple([x for x in (d, t, l, i, c)])

        with open(os.path.join(rates_folder, "Transfer_rates.tsv")) as f:
            f.readline()
            for line in f:
                dn, rc, wt  = line.split("\t")
                if dn not in self.transfer_rates:
                    self.transfer_rates[dn] = dict()
                if rc not in self.transfer_rates[dn]:
                    self.transfer_rates[dn][rc] = 0.0

                self.transfer_rates[dn][rc] = float(wt)


    def choose_event(self, duplication, transfer, loss, inversion, transposition, origination):

        draw = numpy.random.choice(["D", "T", "L", "I", "P", "O"], 1,
                                   p=af.normalize([duplication, transfer, loss, inversion, transposition, origination]))
        return draw

    def choose_event_i(self, duplication, transfer, loss, inversion, transposition, origination, remove, rewire):

        draw = numpy.random.choice(["D", "T", "L", "I", "P", "O", "RM", "RW"], 1,
                                   p=af.normalize([duplication, transfer, loss, inversion, transposition, origination, remove, rewire]))
        return draw

    def choose_recipient(self, lineages_alive, donor):
        possible_recipients = [x for x in lineages_alive if x != donor]
        if len(possible_recipients) > 1:
            recipient = random.choice(possible_recipients)
            return recipient
        else:
            return None

    def evolve_genomes(self, duplication, transfer, loss, inversion, transposition, origination, time):

        lineage = random.choice(list(self.active_genomes))
        event = self.choose_event(duplication, transfer, loss, inversion, transposition, origination)

        if event == "D":
            d_e = self.parameters["DUPLICATION_EXTENSION"]
            self.make_duplication(d_e, lineage, time)
            return "D", lineage

        elif event == "T":

            t_e = self.parameters["TRANSFER_EXTENSION"]

            possible_recipients = [x for x in self.active_genomes if x != lineage]

            if len(possible_recipients) > 0:

                donor = lineage

                # We choose a recipient

                if self.parameters["ASSORTATIVE_TRANSFER"] == "True":
                    recipient = self.choose_assortative_recipient(time, possible_recipients, donor)
                    if recipient == None:
                        return None
                else:
                    recipient = random.choice(possible_recipients)

                self.make_transfer(t_e, donor, recipient, time)
                return "T", donor+"->"+recipient

            else:
                return None

        elif event == "L":

            l_e = self.parameters["LOSS_EXTENSION"]

            self.make_loss(l_e, lineage, time)
            return "L", lineage

        elif event == "I":
            i_e = self.parameters["INVERSION_EXTENSION"]
            self.make_inversion(i_e, lineage, time)
            return "I", lineage

        elif event == "P":
            c_e = self.parameters["TRANSPOSITION_EXTENSION"]
            self.make_transposition(c_e, lineage, time)
            return "P",lineage

        elif event == "O":

            gene, gene_family = self.make_origination(lineage, time)
            chromosome = self.all_genomes[lineage].select_random_chromosome()
            position = chromosome.select_random_position()
            segment = [gene]
            chromosome.insert_segment(position, segment)
            return "O", lineage

    def evolve_genomes_i(self, duplication, transfer, loss, inversion, transposition, origination, remove, rewire, time):

        d_e = self.parameters["DUPLICATION_EXTENSION"]
        t_e = self.parameters["TRANSFER_EXTENSION"]
        l_e = self.parameters["LOSS_EXTENSION"]
        i_e = self.parameters["INVERSION_EXTENSION"]
        c_e = self.parameters["TRANSPOSITION_EXTENSION"]

        lineage = random.choice(list(self.active_genomes))
        event = self.choose_event_i(duplication, transfer, loss, inversion, transposition, origination, remove, rewire)

        if event == "D":

            self.make_duplication_interactome(d_e, lineage, time)
            return "D", lineage

        elif event == "T":

            # We choose a recipient

            possible_recipients = [x for x in self.active_genomes if x != lineage]

            if len(possible_recipients) > 0:

                if self.parameters["ASSORTATIVE_TRANSFER"] == "True":
                    recipient = self.choose_assortative_recipient(time, possible_recipients, donor)
                    if recipient == None:
                        return None
                else:
                    recipient = random.choice(possible_recipients)

                donor = lineage
                self.make_transfer_interactome(t_e, donor, recipient, time)
                return "T", donor + "->" + recipient

            else:
                return None

        elif event == "L":

            self.make_loss_interactome(l_e, lineage, time)
            return "L", lineage

        elif event == "I":
            self.make_inversion(i_e, lineage, time)
            return "I", lineage

        elif event == "P":
            self.make_transposition(c_e, lineage, time)
            return "P", lineage

        elif event == "O":

            gene, gene_family = self.make_origination(lineage, time)
            chromosome = self.all_genomes[lineage].select_random_chromosome()
            position = chromosome.select_random_position()
            segment = [gene]
            chromosome.insert_segment(position, segment)

            # We need to insert the gene in the interactome too, with preferential attachment

            interactome = self.all_genomes[lineage].interactome

            node_degrees = [d + 1 for n, d in interactome.degree()]
            choice = numpy.random.choice(interactome.nodes, 1, p=af.normalize(node_degrees))[0]

            interactome.add_node(str(gene))
            interactome.add_edge(str(gene), choice)

            return "O", lineage

        elif event == "RM":

            self.make_remove_edge(lineage, time)
            return "RM", lineage


        elif event == "RW":

            self.make_rewiring_edge(lineage, time)

            return "RW", lineage


    def evolve_genomes_m(self, time):

        d_e = self.parameters["DUPLICATION_EXTENSION"]
        t_e = self.parameters["TRANSFER_EXTENSION"]
        l_e = self.parameters["LOSS_EXTENSION"]
        i_e = self.parameters["INVERSION_EXTENSION"]
        c_e = self.parameters["TRANSPOSITION_EXTENSION"]

        ####


        ####

        mactive_genomes = list(self.active_genomes)
        mweights = list()
        for genome in mactive_genomes:
            lineage_weight = 0
            for chromosome in self.all_genomes[genome]:
                for gene in chromosome:
                    for r,vl in self.all_gene_families[gene.gene_family].rates.items():
                        lineage_weight += vl
            mweights.append(lineage_weight)

        lineage = numpy.random.choice(mactive_genomes, 1, p=af.normalize(mweights))[0]

        d, t, l, i, p, o = 0, 0, 0, 0, 0, 0

        for chromosome in self.all_genomes[lineage]:
            for gene in chromosome:
                d += self.all_gene_families[gene.gene_family].rates["DUPLICATION"]
                t += self.all_gene_families[gene.gene_family].rates["TRANSFER"]
                l += self.all_gene_families[gene.gene_family].rates["LOSS"]

        i += af.obtain_value((self.parameters["INVERSION"]))
        p += af.obtain_value((self.parameters["TRANSPOSITION"]))
        o += af.obtain_value((self.parameters["ORIGINATION"]))

        #print(d,t,l,i,p,o)

        event = self.choose_event(d, t, l, i, p, o)

        ####

        if event == "D":

            self.make_duplication(d_e, lineage, time, family_mode=True)
            return "D", lineage

        elif event == "T":

            # We choose a recipient

            possible_recipients = [x for x in self.active_genomes if x != lineage]

            if len(possible_recipients) > 0:

                if self.parameters["ASSORTATIVE_TRANSFER"] == "True":
                    recipient = self.choose_assortative_recipient(time, possible_recipients, donor)
                    if recipient == None:
                        return None
                else:
                    recipient = random.choice(possible_recipients)

                donor = lineage
                self.make_transfer(t_e, donor, recipient, time, family_mode = True)
                return "T", donor + "->" + recipient

            else:
                return None

        elif event == "L":
            self.make_loss(l_e, lineage, time, family_mode = True)
            return "L", lineage

        elif event == "I":
            self.make_inversion(i_e, lineage, time)
            return "I", lineage

        elif event == "P":
            self.make_transposition(c_e, lineage, time)
            return "P", lineage

        elif event == "O":


            if  self.parameters["RATE_FILE"] == "False":
                gene, gene_family = self.make_origination(lineage, time, family_mode=True)

            elif  self.parameters["RATE_FILE"] != "False":
                gene, gene_family = self.make_origination(lineage, time, family_mode=True,
                                                          empirical_rates=True)
            chromosome = self.all_genomes[lineage].select_random_chromosome()
            position = chromosome.select_random_position()
            segment = [gene]
            chromosome.insert_segment(position, segment)

            return "O", lineage


    def advanced_evolve_genomes(self, time):

        active_genomes = list(self.active_genomes)
        lineage = numpy.random.choice(active_genomes, 1, p=af.normalize(
            [sum(self.branch_event_rates[x]) for x in active_genomes]))[0]

        d,t,l,i,c,o = self.branch_event_rates[lineage]

        event = self.choose_event(d,t,l,i,c,o)

        d_e, t_e, l_e, i_e, c_e = self.branch_extension_rates[lineage]

        if event == "D":
            self.make_duplication(d_e, lineage, time)
            return "D", lineage

        elif event == "T":

            # We choose a recipient

            possible_recipients = [x for x in self.active_genomes if x != lineage]

            if len(possible_recipients) > 0:

                recipient = self.choose_advanced_recipient(possible_recipients, lineage)
                if recipient != None:
                    donor = lineage
                    self.make_transfer(t_e, donor, recipient, time)
                    return "T", donor+"->"+recipient

            else:
                return None

        elif event == "L":

            self.make_loss(l_e, lineage, time)
            return "L", lineage

        elif event == "I":
            self.make_inversion(i_e, lineage, time)
            return "I", lineage

        elif event == "P":
            self.make_transposition(c_e, lineage, time)
            return "P",lineage

        elif event == "O":

            gene, gene_family = self.make_origination(lineage, time)

            chromosome = self.all_genomes[lineage].select_random_chromosome()
            position = chromosome.select_random_position()
            segment = [gene]
            chromosome.insert_segment(position, segment)

            return "O", lineage

    def advanced_evolve_genomes_f(self, duplication, transfer, loss, inversion, transposition, origination, time):


        d_e = int(af.obtain_value(self.parameters["DUPLICATION_EXTENSION"]))
        t_e = int(af.obtain_value(self.parameters["TRANSFER_EXTENSION"]))
        l_e = int(af.obtain_value(self.parameters["LOSS_EXTENSION"]))
        i_e = int(af.obtain_value(self.parameters["INVERSION_EXTENSION"]))
        c_e = int(af.obtain_value(self.parameters["TRANSPOSITION_EXTENSION"]))

        #print(d_e, t_e, l_e, i_e, c_e)

        mean_gene_length = int()

        distribution  = self.parameters["GENE_LENGTH"].split(":")[0]

        if distribution == "f" or distribution == "n":
            mean_gene_length = int(self.parameters["GENE_LENGTH"].split(":")[1].split(";")[0])
        elif distribution == "u":
            u1,u0 = self.parameters["GENE_LENGTH"].split(":")[1].split(";")
            mean_gene_length = (int(u1) - int(u0))/2
        else:
            print("Error, please switch the distribution type por the gene length")
            return 0

        mean_intergene_length = int(self.parameters["INTERGENE_LENGTH"])
        multiplier = 1.0 / (mean_gene_length + mean_intergene_length)


        lineage = random.choice(list(self.active_genomes))
        event = self.choose_event(duplication, transfer, loss, inversion, transposition, origination)

        for chromosome in self.all_genomes[lineage]:
            chromosome.obtain_flankings()
            chromosome.obtain_locations()

        if event == "D":

            r = self.select_advanced_length(lineage, d_e * multiplier)
            if r == None:
                return None
            else:
                c1, c2, d = r
                self.make_duplication_within_intergene(c1, c2, d, lineage, time)

            return "D", lineage

        elif event == "T":


            # We choose a recipient

            possible_recipients = [x for x in self.active_genomes if x != lineage]

            if len(possible_recipients) > 0:

                if self.parameters["ASSORTATIVE_TRANSFER"] == "True":
                    recipient = self.choose_assortative_recipient(time, possible_recipients, donor)
                    if recipient == None:
                        return None
                else:
                    recipient = random.choice(possible_recipients)

                donor = lineage

                r = self.select_advanced_length(lineage, t_e * multiplier)

                if r == None:
                    return None
                else:
                    c1, c2, d = r
                    self.make_transfer_intergenic(c1, c2, d, donor, recipient, time)

                return "T", donor + "->" + recipient

            else:
                return None

        elif event == "L":

            r = self.select_advanced_length(lineage, l_e * multiplier)

            if r == None:
                return None
            else:
                c1, c2, d = r
                pseudo = False
                if numpy.random.uniform(0,1) <= float(self.parameters["PSEUDOGENIZATION"]):
                    pseudo = True
                self.make_loss_intergenic(c1, c2, d, lineage, time, pseudo)

            return "L", lineage

        elif event == "I":

            r = self.select_advanced_length(lineage, i_e * multiplier)

            if r == None:
                return None
            else:
                c1, c2, d = r
                self.make_inversion_intergenic(c1, c2, d, lineage, time)
            return "I", lineage

        elif event == "P":

            r = self.select_advanced_length(lineage, c_e * multiplier)
            if r == None:
                return None
            else:
                c1, c2, d = r
                self.make_transposition_intergenic(c1, c2, d, lineage, time)

            return "P", lineage

        elif event == "O":

            gene, gene_family = self.make_origination(lineage, time)
            chromosome = self.all_genomes[lineage].select_random_chromosome()
            intergene_coordinate = chromosome.select_random_coordinate_in_intergenic_regions()
            location = chromosome.return_location_by_coordinate(intergene_coordinate, within_intergene=True)
            chromosome.insert_gene_within_intergene(intergene_coordinate, location, gene)

            return "O", lineage


    def get_time_to_next_event(self, n, events):

        total = 0.0
        for __ in range(n):
            total += sum(events)

        if total == 0:
            return 1000000000000000 # We sent an arbitrarily big number. Probably not the most elegant thing to do
        else:
            time = numpy.random.exponential(1/total)
            return time

    def get_time_to_next_event_advanced_modes(self):
        # To obtain the time to next event in case that we have different rates per branch
        total = 0.0

        for lineage in self.active_genomes:
            total += sum(self.branch_event_rates[lineage])

        if total == 0:
            return 1000000000000000 # We sent an arbitrarily big number. Probably not the most elegant thing to do
        time = numpy.random.exponential(1 / total)
        return time

    def get_time_to_next_event_family_mode(self):

        total = 0.0

        for lineage in self.active_genomes:
            for chromosome in self.all_genomes[lineage]:
                for gene in chromosome:
                    for r,vl in self.all_gene_families[gene.gene_family].rates.items():
                        total += vl

        total_active =  len(self.active_genomes)

        total += af.obtain_value((self.parameters["INVERSION"])) * total_active
        total += af.obtain_value((self.parameters["TRANSPOSITION"])) * total_active
        total += af.obtain_value((self.parameters["ORIGINATION"])) * total_active

        if total == 0:
            return 1000000000000000 # We sent an arbitrarily big number. Probably not the most elegant thing to do
        time = numpy.random.exponential(1 / total)

        return time

    def increase_distances(self, time_to_next_event, active_lineages):

        for node in active_lineages:
            node.dist += time_to_next_event

    def make_origination(self, species_tree_node, time, family_mode = False, empirical_rates = False):

        self.gene_families_counter += 1
        gene_family_id = str(self.gene_families_counter)

        gene = Gene()
        gene.determine_orientation()

        gene.gene_family = str(self.gene_families_counter)
        gene.species = species_tree_node

        gene_family = GeneFamily(gene_family_id, time)
        gene_family.length = int(af.obtain_value(self.parameters["GENE_LENGTH"]))
        gene.length = gene_family.length

        gene_family.genes.append(gene)
        gene.gene_id = gene_family.obtain_new_gene_id()

        self.all_gene_families[gene_family_id] = gene_family
        self.all_gene_families[gene.gene_family].register_event(str(time), "O", species_tree_node)

        if family_mode == True and empirical_rates == False:

            d, t,l, _, _, _ = self.generate_new_rates()
            gene_family.rates["DUPLICATION"] = d
            gene_family.rates["TRANSFER"] = t
            gene_family.rates["LOSS"] = l

        elif family_mode == True and empirical_rates == True:
            d, t, l = self.generate_empirical_rates()
            gene_family.rates["DUPLICATION"] = d
            gene_family.rates["TRANSFER"] = t
            gene_family.rates["LOSS"] = l

        return gene, gene_family

    def make_speciation(self, sp, c1, c2, time):

        # This function receives a genome and the names of the two branching lineages of the species node

        genome_sp = self.all_genomes[sp]

        genome1 = Genome()
        genome2 = Genome()

        if hasattr(genome_sp, 'interactome'):
            genome1.interactome = copy.deepcopy(genome_sp.interactome)
            genome2.interactome = copy.deepcopy(genome_sp.interactome)

            new_names_1 = dict()
            new_names_2 = dict()

        for chromosome in genome_sp:

            shape = chromosome.shape

            if shape == "C":
                ch1 = CircularChromosome()
                ch2 = CircularChromosome()
                ch1.shape = "C"
                ch2.shape = "C"

            elif shape == "L":
                ch1 = LinearChromosome()
                ch2 = LinearChromosome()
                ch1.shape = "L"
                ch2.shape = "L"

            genome1.chromosomes.append(ch1)
            genome2.chromosomes.append(ch2)

            if chromosome.has_intergenes:
                ch1.has_intergenes = True
                ch2.has_intergenes = True

            for gene in chromosome:

                new_id1 = self.return_new_identifiers_for_segment([gene])
                new_id2 = self.return_new_identifiers_for_segment([gene])

                new_gene1 = af.copy_segment([Gene()], new_id1)[0]
                new_gene2 = af.copy_segment([Gene()], new_id2)[0]

                new_gene1.species = c1
                new_gene2.species = c2

                new_gene1.orientation = gene.orientation
                new_gene2.orientation = gene.orientation

                new_gene1.length = gene.length
                new_gene2.length = gene.length

                gene_family = self.all_gene_families[gene.gene_family]
                gene_family.genes.append(new_gene1)
                gene_family.genes.append(new_gene2)

                new_gene1.gene_family = gene.gene_family
                new_gene2.gene_family = gene.gene_family

                ch1.genes.append(new_gene1)
                ch2.genes.append(new_gene2)

                if hasattr(genome_sp, 'interactome'):

                    new_names_1[str(gene)] = str(new_gene1)
                    new_names_2[str(gene)] = str(new_gene2)

                gene.active = False

                # The code for the node is:
                # 1. Branch of the species tree that splits
                # 2. Id of the gene that is split
                # 3. Branch of the species tree first child
                # 4. Id of the gene that goes to first child
                # 5. Branch of the species tree second child
                # 6. Id of the gene that goes to second child

                nodes = [sp,
                         gene.gene_id,
                         c1,
                         new_gene1.gene_id,
                         c2,
                         new_gene2.gene_id
                         ]

                self.all_gene_families[gene.gene_family].register_event(str(time), "S", ";".join(map(str,nodes)))

            for intergene in chromosome.intergenes:

                new_intergene1 = Intergene()
                new_intergene2 = Intergene()
                new_intergene1.length = intergene.length
                new_intergene2.length = intergene.length
                ch1.intergenes.append(new_intergene1)
                ch2.intergenes.append(new_intergene2)


        genome1.update_genome_species(c1)
        genome2.update_genome_species(c2)

        if hasattr(genome_sp, 'interactome'):
            genome1.interactome = nx.relabel_nodes(genome1.interactome, new_names_1)
            genome2.interactome = nx.relabel_nodes(genome2.interactome, new_names_2)
            #nx.relabel_nodes(genome1.interactome, new_names_1)
            #nx.relabel_nodes(genome2.interactome, new_names_2)

        return genome1, genome2

    def make_extinction(self, sp, time):

        # We have to inactivate all the genes

        genome = self.all_genomes[sp]

        for chromosome in genome:
            for gene in chromosome:
                gene.active = False
                self.all_gene_families[gene.gene_family].register_event(str(time), "E", ";".join(map(str,[sp, gene.gene_id])))

    def make_end(self, time):

        for genome_name in self.active_genomes:
            genome = self.all_genomes[genome_name]
            for chromosome in genome:
                for gene in chromosome:
                    gene.active = False
                    self.all_gene_families[gene.gene_family].register_event(str(time), "F", ";".join(
                        map(str, [genome.species, gene.gene_id])))


    def make_duplication(self, p, lineage, time, family_mode = False):

        chromosome = self.all_genomes[lineage].select_random_chromosome()

        if family_mode == True:
            affected_genes = chromosome.obtain_affected_genes_accounting_for_family_rates(p, self.all_gene_families, "DUPLICATION")
        else:
            affected_genes = chromosome.obtain_affected_genes(p)

        segment = chromosome.obtain_segment(affected_genes)

        new_identifiers1 = self.return_new_identifiers_for_segment(segment)
        new_identifiers2 = self.return_new_identifiers_for_segment(segment)

        # Now we create two segments

        copied_segment1 = af.copy_segment(segment, new_identifiers1)
        copied_segment2 = af.copy_segment(segment, new_identifiers2)

        # We insert the two new segments after the last position of the old segment

        chromosome.insert_segment(affected_genes[-1], copied_segment1 + copied_segment2)

        # And we remove the old segment

        chromosome.remove_segment(segment)

        # We have to register in the affected gene families that there has been a duplication

        for i, gene in enumerate(segment):

            nodes = [gene.species,
                     gene.gene_id,
                     copied_segment1[i].species,
                     copied_segment1[i].gene_id,
                     copied_segment2[i].species,
                     copied_segment2[i].gene_id]

            gene.active = False

            # We add the genes to the list of genes in the gene family
            gene_family = gene.gene_family

            self.all_gene_families[gene_family].genes.append(copied_segment1[i])
            self.all_gene_families[gene_family].genes.append(copied_segment2[i])
            self.all_gene_families[gene_family].register_event(time, "D", ";".join(map(str, nodes)))

    def make_duplication_interactome(self, p, lineage, time):

        chromosome = self.all_genomes[lineage].select_random_chromosome()
        affected_genes = chromosome.obtain_affected_genes(p)
        segment = chromosome.obtain_segment(affected_genes)

        new_identifiers1 = self.return_new_identifiers_for_segment(segment)
        new_identifiers2 = self.return_new_identifiers_for_segment(segment)

        # Now we create two segments

        copied_segment1 = af.copy_segment(segment, new_identifiers1)
        copied_segment2 = af.copy_segment(segment, new_identifiers2)

        # We insert the two new segments after the last position of the old segment

        chromosome.insert_segment(affected_genes[-1], copied_segment1 + copied_segment2)

        # And we remove the old segment

        chromosome.remove_segment(segment)

        # We have to register in the affected gene families that there has been a duplication

        # We need to create two dicts to change the names of the interactome

        new_genes_1 = dict()

        for i, gene in enumerate(segment):

            nodes = [gene.species,
                     gene.gene_id,
                     copied_segment1[i].species,
                     copied_segment1[i].gene_id,
                     copied_segment2[i].species,
                     copied_segment2[i].gene_id]

            gene.active = False

            gene_family = gene.gene_family

            self.all_gene_families[gene_family].genes.append(copied_segment1[i])
            self.all_gene_families[gene_family].genes.append(copied_segment2[i])

            self.all_gene_families[gene.gene_family].register_event(time, "D", ";".join(map(str, nodes)))

            # We update the new gene name (which by default is going to be the copied_segment1)

            gene1 = copied_segment1[i]
            gene2 = copied_segment2[i]

            new_genes_1[str(gene)] = str(gene1)

            self.all_genomes[lineage].interactome = nx.relabel_nodes(self.all_genomes[lineage].interactome,
                                                                     new_genes_1)

            # WE ADD THE NEW NODE

            self.all_genomes[lineage].interactome.add_node(str(gene2))

            # We distribute node depending on the parameter PROPORTION

            # If p is 1, all the links go to the first node
            # If p is 0, all the links go to the second node
            # If p is 0.5, equal repartition.

            PROPORTION = 0.5

            n_edges_to_old_node = int(PROPORTION * len(self.all_genomes[lineage].interactome.edges(str(gene1))))
            myedges = list(self.all_genomes[lineage].interactome.edges(str(gene1)))
            random.shuffle(myedges)
            edges_to_new_node = myedges[n_edges_to_old_node:]
            edges_to_add_to_new_node = [(str(gene2), x[1]) for x in edges_to_new_node]
            edges_to_remove_to_old_node = edges_to_new_node

            # Now, I have to remove the edges

            self.all_genomes[lineage].interactome.remove_edges_from(edges_to_remove_to_old_node)
            self.all_genomes[lineage].interactome.add_edges_from(edges_to_add_to_new_node)


    def make_duplication_within_intergene(self, c1, c2, d, lineage, time):

        chromosome = self.all_genomes[lineage].select_random_chromosome()
        r = chromosome.return_affected_region(c1, c2, d)

        if r == None:
            return None

        else:
            r1, r2, r3, r4 = r
            segment = chromosome.obtain_segment(r1)
            intergene_segment = chromosome.obtain_intergenic_segment(r2[1:])

        new_identifiers1 = self.return_new_identifiers_for_segment(segment)
        new_identifiers2 = self.return_new_identifiers_for_segment(segment)

        # We duplicate the genes

        new_segment_1 = af.copy_segment(segment, new_identifiers1)
        new_segment_2 = af.copy_segment(segment, new_identifiers2)

        # And the intergenes

        new_intergene_segment_1 = [copy.deepcopy(chromosome.intergenes[x]) for x in r2[1:]]
        new_intergene_segment_2 = [copy.deepcopy(chromosome.intergenes[x]) for x in r2[1:]]

        scar1 = new_intergene_segment_1[-1]
        scar2 = new_intergene_segment_2[-1]

        new_segment = new_segment_1 + new_segment_2
        new_intergene_segment = new_intergene_segment_1 + new_intergene_segment_2

        ###
        ###

        position = r1[-1] + 1

        for i, gene in enumerate(new_segment):
            chromosome.genes.insert(position + i, gene)
        for i, intergene in enumerate(new_intergene_segment):
            chromosome.intergenes.insert(position + i, intergene)

        # We remove the old copies:

        chromosome.remove_segment(segment)
        chromosome.remove_intersegment(intergene_segment)

        # We adjust the new intergenes lengths

        if d == "left":
            r3, r4 = r4, r3

        scar1.length = r3[1] + r4[0]
        scar2.length = r4[1]

        for i, gene in enumerate(segment):
            nodes = [gene.species,
                     gene.gene_id,
                     new_segment_1[i].species,
                     new_segment_1[i].gene_id,
                     new_segment_2[i].species,
                     new_segment_2[i].gene_id]

            gene.active = False

            gene_family = gene.gene_family

            self.all_gene_families[gene_family].genes.append(new_segment_1[i])
            self.all_gene_families[gene_family].genes.append(new_segment_2[i])

            self.all_gene_families[gene.gene_family].register_event(time, "D", ";".join(map(str, nodes)))


    def choose_assortative_recipient(self, time, possible_recipients, donor):

        alpha = self.parameters["ALPHA"]
        weights = list()

        mdonor = self.complete_tree&donor

        for recipient in possible_recipients:
            mrecipient = self.complete_tree&recipient
            ca = self.complete_tree.get_common_ancestor(mrecipient, mdonor).name
            x1 = self.distances_to_start[ca]
            td =  time - x1
            weights.append(td)

        beta = min(alpha * af.normalize(weights))
        val = (alpha * af.normalize(weights)) - beta
        pvector = af.normalize(numpy.exp(-val))

        draw = numpy.random.choice(possible_recipients, 1, p=pvector)[0]

        return draw

    def choose_advanced_recipient(self, possible_recipients, donor):

        weights = list()

        for recipient in possible_recipients:
            weights.append(self.transfer_rates[donor][recipient])

        if sum(weights) == 0:
            return None


        draw = numpy.random.choice(possible_recipients, 1, p=af.normalize(weights))[0]

        return draw

    def make_transfer(self, p, donor, recipient, time, family_mode = False):

        chromosome1 = self.all_genomes[donor].select_random_chromosome()

        if family_mode == True:
            affected_genes = chromosome1.obtain_affected_genes_accounting_for_family_rates(p, self.all_gene_families, "TRANSFER")
        else:
            affected_genes = chromosome1.obtain_affected_genes(p)

        segment = chromosome1.obtain_segment(affected_genes)
        new_identifiers1 = self.return_new_identifiers_for_segment(segment)
        new_identifiers2 = self.return_new_identifiers_for_segment(segment)

        inverted = False

        # Now we create two segments

        copied_segment1 = af.copy_segment(segment, new_identifiers1)
        copied_segment2 = af.copy_segment(segment, new_identifiers2)

        # We insert the first segment (leaving transfer) in the same position than the previous segment
        # We do this just to change the identifiers of the numbers

        chromosome1.insert_segment(affected_genes[0], copied_segment1)

        # And we remove the old segment

        chromosome1.remove_segment(segment)

        # Now we insert the transfer segment in the recipient genome in one of the homologous position.

        if numpy.random.uniform(0,1) <= self.parameters["REPLACEMENT_TRANSFER"]:

            possible_positions = list()

            for chromosome in self.all_genomes[recipient]:
                for direction, positions in chromosome.get_homologous_position(segment):
                    possible_positions.append((direction, positions, chromosome))

            if len(possible_positions) != 0:

                direction, positions, chromosome2 = random.choice(possible_positions)


                if direction == "F":

                    # I replace gene by gene the segment

                    i = 0

                    for position in positions:

                        # First I inactivate the gene

                        gene = chromosome2.genes[position]
                        gene.active = False
                        self.all_gene_families[gene.gene_family].register_event(time, "L", ";".join(
                            map(str, [recipient, gene.gene_id])))

                        # And then I replace

                        chromosome2.genes[position] = copied_segment2[i]

                        i += 1

                elif direction == "B":
                    # I invert the segment and I replace gene by gene the segment
                    inverted = True

                    copied_segment2 = copied_segment2[::-1]

                    for gene in copied_segment2:
                        gene.change_sense()

                    i = 0
                    for position in positions:
                        # First I inactivate the gene
                        gene = chromosome2.genes[position]
                        gene.active = False
                        self.all_gene_families[gene.gene_family].register_event(time, "L", ";".join(
                            map(str, [recipient, gene.gene_id])))

                        # And then I replace

                        chromosome2.genes[position] = copied_segment2[i]


                        i += 1
            else:

                # Normal transfers
                chromosome2 = self.all_genomes[recipient].select_random_chromosome()
                position = chromosome2.select_random_position()
                chromosome2.insert_segment(position, copied_segment2)
        else:
            # Normal transfer
            chromosome2 = self.all_genomes[recipient].select_random_chromosome()
            position = chromosome2.select_random_position()
            chromosome2.insert_segment(position, copied_segment2)

        # We have to register in the affected gene families that there has been a transfer event

        if inverted == True:
            # We invert again to store the event
            copied_segment2 = copied_segment2[::-1]

        for i, gene in enumerate(segment):

            gene.active = False

            # The code for the node is:
            # 1. Branch of the species tree for the donor genome
            # 2. Id of the gene that is transferred
            # 3. Id of the gene that remains in the donor genome
            # 4. Branch of the species tree for the recipient genome
            # 5. Id of the new gene arriving

            copied_segment1[i].species = donor
            copied_segment2[i].species = recipient

            nodes = [gene.species,
                     gene.gene_id,
                     copied_segment1[i].species,
                     copied_segment1[i].gene_id,
                     copied_segment2[i].species,
                     copied_segment2[i].gene_id]

            self.all_gene_families[gene.gene_family].register_event(time, "T", ";".join(map(str,nodes)))

    def make_transfer_interactome(self, p, donor, recipient, time):

        chromosome1 = self.all_genomes[donor].select_random_chromosome()
        affected_genes = chromosome1.obtain_affected_genes(p)
        segment = chromosome1.obtain_segment(affected_genes)

        new_identifiers1 = self.return_new_identifiers_for_segment(segment)
        new_identifiers2 = self.return_new_identifiers_for_segment(segment)

        inverted = False

        is_replacement_transfer = False
        replaced_genes = list()

        # Now we create two segments

        copied_segment1 = af.copy_segment(segment, new_identifiers1)
        copied_segment2 = af.copy_segment(segment, new_identifiers2)

        # We insert the first segment (leaving transfer) in the same position than the previous segment
        # We do this just to change the identifiers of the numbers

        chromosome1.insert_segment(affected_genes[0], copied_segment1)

        # And we remove the old segment

        chromosome1.remove_segment(segment)

        # Now we insert the transfer segment in the recipient genome in one of the homologous position.

        if numpy.random.uniform(0,1) <= self.parameters["REPLACEMENT_TRANSFER"]:

            possible_positions = list()

            for chromosome in self.all_genomes[recipient]:
                for direction, positions in chromosome.get_homologous_position(segment):
                    possible_positions.append((direction, positions, chromosome))

            if len(possible_positions) != 0:

                direction, positions, chromosome2 = random.choice(possible_positions)


                if direction == "F":

                    # I replace gene by gene the segment

                    i = 0

                    for position in positions:

                        # First I inactivate the gene

                        gene = chromosome2.genes[position]
                        gene.active = False


                        self.all_gene_families[gene.gene_family].register_event(time, "L", ";".join(
                            map(str, [recipient, gene.gene_id])))
                        # And then I replace
                        chromosome2.genes[position] = copied_segment2[i]
                        i += 1
                        is_replacement_transfer = True

                        replaced_genes.append((gene, copied_segment2[i]))

                elif direction == "B":
                    # I invert the segment and I replace gene by gene the segment
                    inverted = True
                    copied_segment2 = copied_segment2[::-1]
                    for gene in copied_segment2:
                        gene.change_sense()

                    i = 0
                    for position in positions:
                        # First I inactivate the gene
                        gene = chromosome2.genes[position]
                        gene.active = False

                        replaced_genes.append(gene)

                        self.all_gene_families[gene.gene_family].register_event(time, "L", ";".join(
                            map(str, [recipient, gene.gene_id])))
                        # And then I replace
                        chromosome2.genes[position] = copied_segment2[i]
                        i += 1
                        is_replacement_transfer = True

                        replaced_genes.append((gene, copied_segment2[i]))

            else:

                # Normal transfers
                chromosome2 = self.all_genomes[recipient].select_random_chromosome()
                position = chromosome2.select_random_position()
                chromosome2.insert_segment(position, copied_segment2)
        else:
            # Normal transfer
            chromosome2 = self.all_genomes[recipient].select_random_chromosome()
            position = chromosome2.select_random_position()
            chromosome2.insert_segment(position, copied_segment2)

        # We have to register in the affected gene families that there has been a transfer event

        if inverted == True:
            # We invert again to store the event
            copied_segment2 = copied_segment2[::-1]


        ## We register the event

        new_names_1 = dict()
        new_names_2 = dict()

        for i, gene in enumerate(segment):

            gene.active = False

            # The code for the node is:
            # 1. Branch of the species tree for the donor genome
            # 2. Id of the gene that is transferred
            # 3. Id of the gene that remains in the donor genome
            # 4. Branch of the species tree for the recipient genome
            # 5. Id of the new gene arriving

            copied_segment1[i].species = donor
            copied_segment2[i].species = recipient

            nodes = [gene.species,
                     gene.gene_id,
                     copied_segment1[i].species,
                     copied_segment1[i].gene_id,
                     copied_segment2[i].species,
                     copied_segment2[i].gene_id]

            self.all_gene_families[gene.gene_family].register_event(time, "T", ";".join(map(str,nodes)))

            new_names_1[str(gene)] = str(copied_segment1[i])

            ## We update the interactome

            # First we update the interactome in the donor lineage

            self.all_genomes[donor].interactome = nx.relabel_nodes(self.all_genomes[donor].interactome, new_names_1)

            # Second we update the interactome in the recipient lineage

            if is_replacement_transfer == True:

                ## All links pass now to the new gene
                ## The old gene gets no links

                new_names_2 = {str(n1):str(n2) for n1,n2 in replaced_genes}
                #new_names_2[str(gene)] = str(copied_segment2[i])
                self.all_genomes[recipient].interactome = nx.relabel_nodes(self.all_genomes[recipient].interactome, new_names_2)


            else:

                # It is not a replacement transfer. Preferential attachment

                node_degrees = [d + 1 for n, d in self.all_genomes[recipient].interactome.degree()]
                choice = numpy.random.choice(self.all_genomes[recipient].interactome.nodes, 1, p=af.normalize(node_degrees))[0]
                self.all_genomes[recipient].interactome.add_node(str(copied_segment2[i]))
                self.all_genomes[recipient].interactome.add_edge(str(copied_segment2[i]), choice)



    def make_transfer_intergenic(self, c1, c2, d, donor, recipient, time):

        chromosome1 = self.all_genomes[donor].select_random_chromosome()
        r = chromosome1.return_affected_region(c1, c2, d)

        if r == None:
            return None

        else:
            r1, r2, r3, r4 = r
            segment = chromosome1.obtain_segment(r1)
            
        new_identifiers1 = self.return_new_identifiers_for_segment(segment)
        new_identifiers2 = self.return_new_identifiers_for_segment(segment)

        # Now we create two segments

        copied_segment1 = af.copy_segment(segment, new_identifiers1)
        copied_segment2 = af.copy_segment(segment, new_identifiers2)

        new_intergene_segment = [copy.deepcopy(chromosome1.intergenes[x]) for x in r2[1:]]

        # We insert the first segment (leaving transfer) in the same position than the previous segment
        # We do this just to change the identifiers of the numbers

        # We insert in the same place

        position = r1[-1] + 1

        for i, gene in enumerate(copied_segment1):
            chromosome1.genes.insert(position + i, gene)

        # We remove the old copies:

        chromosome1.remove_segment(segment)

        # Normal transfer

        chromosome2 = self.all_genomes[recipient].select_random_chromosome()
        chromosome2.obtain_flankings()
        chromosome2.obtain_locations()
        intergene_coordinate = chromosome2.select_random_coordinate_in_intergenic_regions()
        l = chromosome2.return_location_by_coordinate(intergene_coordinate, within_intergene=True)



        position = int(l[4]) + 1

        
        for i, gene in enumerate(copied_segment2):
            chromosome2.genes.insert(position + i, gene)
        for i, intergene in enumerate(new_intergene_segment):
            chromosome2.intergenes.insert(position + i, intergene)

        cut_position = (intergene_coordinate - l[2], l[3] - intergene_coordinate)

        scar1 = chromosome2.intergenes[int(l[4])]
        scar2 = chromosome2.intergenes[position + i]

        if d == "left":
            r3,r4 = r4,r3

        scar1.length = r3[1] + cut_position[0]
        scar2.length = r4[0] + cut_position[1]



        # We have to register in the affected gene families that there has been a transfer event


        for i, gene in enumerate(segment):

            gene.active = False

            # The code for the node is:
            # 1. Branch of the species tree for the donor genome
            # 2. Id of the gene that is transferred
            # 3. Id of the gene that remains in the donor genome
            # 4. Branch of the species tree for the recipient genome
            # 5. Id of the new gene arriving

            copied_segment1[i].species = donor
            copied_segment2[i].species = recipient

            nodes = [gene.species,
                     gene.gene_id,
                     copied_segment1[i].species,
                     copied_segment1[i].gene_id,
                     copied_segment2[i].species,
                     copied_segment2[i].gene_id]

            self.all_gene_families[gene.gene_family].register_event(time, "T", ";".join(map(str, nodes)))


    def make_loss(self, p, lineage, time, family_mode = False):

        chromosome = self.all_genomes[lineage].select_random_chromosome()
        if family_mode == True:
            affected_genes = chromosome.obtain_affected_genes_accounting_for_family_rates(p, self.all_gene_families, "LOSS")
        else:
            affected_genes = chromosome.obtain_affected_genes(p)
        segment = chromosome.obtain_segment(affected_genes)

        # Now we check we are not under the minimum size

        if len(chromosome) - len(affected_genes) <= self.parameters["MIN_GENOME_SIZE"]:
            return 0

        chromosome.remove_segment(segment)

        # We have to register in the affected gene families that there has been as loss
        # All genes affected must be returned

        for gene in segment:
            gene.active = False
            self.all_gene_families[gene.gene_family].register_event(time, "L", ";".join(map(str,[lineage, gene.gene_id])))

    def make_loss_intergenic(self, c1, c2, d, lineage, time, pseudo =False):

        chromosome = self.all_genomes[lineage].select_random_chromosome()
        r = chromosome.return_affected_region(c1, c2, d)

        if r == None:
            return None

        else:
            r1, r2, r3, r4 = r

            segment = chromosome.obtain_segment(r1)
            intergene_segment = chromosome.obtain_intergenic_segment(r2[1:])

            scar1 = chromosome.intergenes[r2[0]]

            # Now we remove the genes

            for gene in segment:
                chromosome.genes.remove(gene)

            # Now we remove the intergenes

            for intergene in intergene_segment:
                chromosome.intergenes.remove(intergene)

            # We modify the length of the scar:

            if d == "left":
                r3, r4 = r4, r3

            if pseudo == True:

                # We need to add the lenght of the genes removed
                scar1.length = sum(r3) + sum(r4) \
                               + sum([x.length for x in segment]) \
                               + sum([x.length for x in intergene_segment[:-1]])
            else:

                scar1.length = r3[0] + r4[1]

        # We have to register in the affected gene families that there has been as loss
        # All genes affected must be returned

        for gene in segment:
            gene.active = False
            self.all_gene_families[gene.gene_family].register_event(time, "L", ";".join(map(str,[lineage, gene.gene_id])))

    def make_loss_family_mode(self, p, lineage, time):

        chromosome = self.all_genomes[lineage].select_random_chromosome()

        affected_genes = chromosome.obtain_affected_genes_accounting_for_family_rates(p, interactome)
        segment = chromosome.obtain_segment(affected_genes)

        # Now we check we are not under the minimum size

        if len(chromosome) - len(affected_genes) <= self.parameters["MIN_GENOME_SIZE"]:
            return 0

        chromosome.remove_segment(segment)

        # We have to register in the affected gene families that there has been as loss
        # All genes affected must be returned

        for gene in segment:
            gene.active = False
            self.all_gene_families[gene.gene_family].register_event(time, "L", ";".join(map(str,[lineage, gene.gene_id])))

    def make_loss_interactome(self, p, lineage, time):

        interactome = self.all_genomes[lineage].interactome
        chromosome = self.all_genomes[lineage].select_random_chromosome()

        affected_genes = chromosome.obtain_affected_genes_accounting_for_connectedness(p, interactome)
        segment = chromosome.obtain_segment(affected_genes)

        # Now we check we are not under the minimum size

        if len(chromosome) - len(affected_genes) <= 0:
            return 0

        chromosome.remove_segment(segment)

        # We have to register in the affected gene families that there has been as loss
        # All genes affected must be returned

        for gene in segment:
            gene.active = False
            self.all_gene_families[gene.gene_family].register_event(time, "L", ";".join(map(str,[lineage, gene.gene_id])))
            # We remove from the connectome

            interactome.remove_node(str(gene))



    def make_inversion(self, p, lineage, time):

        chromosome = self.all_genomes[lineage].select_random_chromosome()
        affected_genes = chromosome.obtain_affected_genes(p)
        segment = chromosome.obtain_segment(affected_genes)
        chromosome.invert_segment(affected_genes)

        for i, gene in enumerate(segment):
            self.all_gene_families[gene.gene_family].register_event(str(time), "I", ";".join(map(str,[lineage, gene.gene_id])))

    def make_inversion_intergenic(self, c1, c2, d, lineage, time):
        
        chromosome = self.all_genomes[lineage].select_random_chromosome()
        r = chromosome.return_affected_region(c1, c2, d)

        if r== None:

            return None
        else:

            r1, r2, r3, r4 = r

            segment = chromosome.obtain_segment(r1)
            chromosome.invert_segment(r1)

            if d == "left":
                r3, r4 = r4, r3

            scar1 = chromosome.intergenes[r2[0]]
            scar2 = chromosome.intergenes[r2[-1]]

            scar1.length = r3[0] + r4[0]
            scar2.length = r4[1] + r3[1]

            for i, gene in enumerate(segment):
                self.all_gene_families[gene.gene_family].register_event(str(time), "I", ";".join(map(str,[lineage, gene.gene_id])))

    def make_transposition(self, p, lineage, time):

        chromosome = self.all_genomes[lineage].select_random_chromosome()
        affected_genes = chromosome.obtain_affected_genes(p)
        segment = chromosome.obtain_segment(affected_genes)
        chromosome.cut_and_paste(affected_genes)

        for i, gene in enumerate(segment):
            self.all_gene_families[gene.gene_family].register_event(str(time), "P", ";".join(map(str,[lineage, gene.gene_id])))

    def make_transposition_intergenic(self, c1, c2, d, lineage, time):

        chromosome = self.all_genomes[lineage].select_random_chromosome()
        r = chromosome.return_affected_region(c1, c2, d)

        if r == None:
            return None

        else:
            r1, r2, r3, r4 = r

        success = False
        counter = 0
        while success == False and counter < 100:
            counter += 1
            c3 = chromosome.select_random_coordinate_in_intergenic_regions()
            l3 = chromosome.return_location_by_coordinate(c3, within_intergene=True)

            tc3_1, tc3_2, sc3_1, sc3_2, p, t = l3
            tc3_1, tc3_2, sc3_1, sc3_2, p = map(int, (tc3_1, tc3_2, sc3_1, sc3_2, p))

            if p not in r2:
                success = True

        if success == False:
            return None

        segment = chromosome.obtain_segment(r1)
        intergene_segment = chromosome.obtain_intergenic_segment(r2[1:])

        scar1 = chromosome.intergenes[r2[0]]
        scar2 = chromosome.intergenes[p]
        scar3 = chromosome.intergenes[r2[-1]]

        new_segment = list()
        new_intergene_segment = list()

        # If we insert in the intergene i, the gene must occupy the position i - 1
        # We store it for reference

        left_gene = chromosome.genes[p]

        # Now we pop the genes

        for gene in segment:
            new_segment.append(chromosome.genes.pop(chromosome.genes.index(gene)))

        # And now we insert the genes at the right of the gene we saved before

        position = chromosome.genes.index(left_gene) + 1

        for i, gene in enumerate(new_segment):
            chromosome.genes.insert(position + i, gene)

        # We move the intergene on the right also

        # We save the position for insertion

        left_intergene = chromosome.intergenes[p]

        for intergene in intergene_segment:
            new_intergene_segment.append(chromosome.intergenes.pop(chromosome.intergenes.index(intergene)))

        # And now we insert the genes at the right of the gene we saved before

        position = chromosome.intergenes.index(left_intergene) + 1

        for i, intergene in enumerate(new_intergene_segment):
            chromosome.intergenes.insert(position + i, intergene)

        # Finally, we modify the segments so that they have the right length

        r5 = (c3 - sc3_1, sc3_2 - c3)

        if d == "left":
            r3, r4 = r4, r3

        scar1.length = r3[0] + r4[1]
        scar2.length = r3[1] + r5[0]
        scar3.length = r4[0] + r5[1]

        for i, gene in enumerate(segment):
            self.all_gene_families[gene.gene_family].register_event(str(time), "P", ";".join(map(str,[lineage, gene.gene_id])))


    def make_rewiring_edge(self, lineage, time):

        chromosome = self.all_genomes[lineage].select_random_chromosome()
        position = chromosome.select_random_position()
        interactome = self.all_genomes[lineage].interactome
        normalized_weights = af.normalize([d + 1 for n, d in interactome.degree()])
        n1 = chromosome.genes[position]
        n2 = numpy.random.choice(interactome.nodes, 1, p=normalized_weights)[0]

        while (str(n1) == str(n2)):
            n2 = numpy.random.choice(interactome.nodes, 1, p=normalized_weights)[0]

        self.all_genomes[lineage].interactome.add_edge(str(n1), str(n2))
        self.all_gene_families[n1.gene_family].register_event(str(time), "RW", ";".join(map(str, [lineage, n1, n2])))

    def make_remove_edge(self, lineage, time):

        myedges = list(self.all_genomes[lineage].interactome.edges())

        if len(myedges) == 0:
            return None

        myedge = myedges[random.randint(0, len(myedges) - 1)]
        self.all_genomes[lineage].interactome.remove_edge(*myedge)

        self.all_gene_families[myedge[0].split("_")[0]].register_event(str(time), "RM", ";".join([lineage, myedge[0], myedge[1]]))


    def get_gene_family_tree(self):

        if len(self.gene_family["Gene_tree"].get_leaves()) < 3:
            return "None"
        else:
            return self.gene_family["Gene_tree"].write(format=1)


    def select_advanced_length(self, lineage, p):

        chromosome = self.all_genomes[lineage].select_random_chromosome()
        total_genome_length = chromosome.map_of_locations[-1][1]
        success = False

        counter = 0

        while counter <= 100 and success == False:

            counter += 1

            sc1 = chromosome.select_random_coordinate_in_intergenic_regions()
            tc1 = chromosome.return_total_coordinate_from_specific_coordinate(sc1, "I")
            d = numpy.random.choice(("left", "right"), p=[0.5, 0.5])

            ## CHANGE

            extension = numpy.random.geometric(p)
            #extension = numpy.random.randint(1000000+1000)


            if d == "right":

                if tc1 + extension >= total_genome_length:
                    tc2 = total_genome_length - (extension - tc1)
                    if tc2 < tc1:
                        success = True
                    else:
                        # The event covers the whole genome
                        pass
                else:
                    tc2 = tc1 + extension
                    success = False

            elif d == "left":

                if tc1 - extension <= 0:
                    tc2 = total_genome_length - extension - (0 - tc1)
                    if tc1 < tc2:
                        success = True
                    else:
                        # The event covers the whole genome
                        success = False
                else:
                    tc2 = tc1 - extension
                    success = True

            if success == True and tc2 >= 0 and tc2 <= total_genome_length:

                sc2 = chromosome.return_specific_coordinate_from_total_coordinate(tc2)
                if sc2 == None:
                    success = False
                else:
                    return sc1, sc2, d
        return None

##################

###### GENOME CLASSES

##################

##################

##################

##################

class GeneFamily():

    def __init__(self, identifier, time):

        self.identifier = identifier
        self.origin = time

        self.genes = list()
        self.events = list()
        self.event_counter = 0  # Each time that the family is modified in any form, we have to update the event counter
        self.gene_ids_counter = 0

        self.length = 0

        self.rates = dict() # Only in Gm mode


    def register_event(self, time, event, genes):

        self.events.append((time, event, genes))

    def generate_tree(self):


        def find_descendant(surviving_nodes, node):

            found = 0
            mynode = surviving_nodes[node]["descendant"]

            while found == 0:

                if surviving_nodes[mynode]["state"] == 1:
                    found = 1
                else:
                    mynode = surviving_nodes[mynode]["descendant"]

            return mynode

        # Eric's algorithm

        # First we will iterate the events from the end

        events = self.events

        surviving_nodes = dict()
        times = dict()

        family_size = 0

        for current_time, event, nodes in events[::-1]:

            if event == "F":

                nodename = nodes.replace(";","_")
                times[nodename] = float(current_time)
                surviving_nodes[nodename] = {"state": 1, "descendant": "None"}

                family_size += 1

            elif event == "E" or event == "L":

                nodename = nodes.replace(";", "_")

                times[nodename] = float(current_time)
                surviving_nodes[nodename] = {"state": 0, "descendant": "None"}

            elif event == "S" or event == "D" or event == "T":

                p, g0, c1, g1, c2, g2 = nodes.split(";")

                pnodename = p + "_" + g0
                c1nodename = c1 + "_" + g1
                c2nodename = c2 + "_" + g2

                times[pnodename] = float(current_time)

                if surviving_nodes[c1nodename]["state"] == 1 and surviving_nodes[c2nodename]["state"] == 1:

                    surviving_nodes[pnodename] = {"state": 1, "descendant": c1nodename + ";" + c2nodename}

                elif surviving_nodes[c1nodename]["state"] == 0 and surviving_nodes[c2nodename]["state"] == 0:

                    surviving_nodes[pnodename] = {"state": 0, "descendant": "None"}

                elif surviving_nodes[c1nodename]["state"] == -1 and surviving_nodes[c2nodename]["state"] == -1:

                    mynode1 = find_descendant(surviving_nodes, c1nodename)
                    mynode2 = find_descendant(surviving_nodes, c2nodename)

                    surviving_nodes[pnodename] = {"state": 1, "descendant": mynode1 + ";" + mynode2}

                elif surviving_nodes[c1nodename]["state"] == 1 and surviving_nodes[c2nodename]["state"] == 0:

                    surviving_nodes[pnodename] = {"state": -1, "descendant": c1nodename}

                elif surviving_nodes[c1nodename]["state"] == 0 and surviving_nodes[c2nodename]["state"] == 1:

                    surviving_nodes[pnodename] = {"state": -1, "descendant": c2nodename}

                elif surviving_nodes[c1nodename]["state"] == 1 and surviving_nodes[c2nodename]["state"] == -1:

                    mynode = find_descendant(surviving_nodes, c2nodename)
                    surviving_nodes[pnodename] = {"state": 1, "descendant": c1nodename + ";" + mynode}

                elif surviving_nodes[c1nodename]["state"] == -1 and surviving_nodes[c2nodename]["state"] == 1:

                    mynode = find_descendant(surviving_nodes, c1nodename)
                    surviving_nodes[pnodename] = {"state": 1, "descendant": mynode + ";" + c2nodename}

                elif surviving_nodes[c1nodename]["state"] == -1 and surviving_nodes[c2nodename]["state"] == 0:

                    mynode = find_descendant(surviving_nodes, c1nodename)
                    surviving_nodes[pnodename] = {"state": -1, "descendant": mynode}

                elif surviving_nodes[c1nodename]["state"] == 0 and surviving_nodes[c2nodename]["state"] == -1:

                    mynode = find_descendant(surviving_nodes, c2nodename)
                    surviving_nodes[pnodename] = {"state": -1, "descendant": mynode}

        extanttree = RT.ReconciledTree()
        completetree = RT.ReconciledTree()

        eroot = extanttree.get_tree_root()
        eroot.name = ""

        wquick_nodes = dict()
        equick_nodes = dict()

        for i, values in enumerate(events):

            current_time, event, nodes = values

            if event == "O":

                wroot = completetree.get_tree_root()
                wroot.name = nodes + "_1"
                wquick_nodes[wroot.name] = wroot

            if event == "L" or event == "E":

                p, g0 = nodes.split(";")
                pnodename = p + "_" + g0
                mynode = wquick_nodes[pnodename]
                e = RT.RecEvent("L", p, int(float(current_time)))
                mynode.addEvent(e, append=True)

            if event == "F":

                p, g0 = nodes.split(";")
                pnodename = p + "_" + g0
                mynode = wquick_nodes[pnodename]
                e = RT.RecEvent("P", p, int(float(current_time)))
                mynode.addEvent(e, append=True)

            if event == "S" or event == "D" or event == "T":

                p, g0, c1, g1, c2, g2 = nodes.split(";")
                pnodename = p + "_" + g0
                c1nodename = c1 + "_" + g1
                c2nodename = c2 + "_" + g2

                mynode = wquick_nodes[pnodename]
                myc1 = mynode.add_child()
                myc2 = mynode.add_child()
                myc1.name = c1nodename
                myc2.name = c2nodename
                myc1.dist = times[c1nodename] - times[pnodename]
                myc2.dist = times[c2nodename] - times[pnodename]

                wquick_nodes[c1nodename] = myc1
                wquick_nodes[c2nodename] = myc2

                state = surviving_nodes[pnodename]["state"]

                ### Now we add the reconciled events

                e = RT.RecEvent(event, p, int(float(current_time)))
                mynode.addEvent(e, append=True)

                if state == 1:  # Now the extant tree

                    c1name, c2name = surviving_nodes[pnodename]["descendant"].split(";")

                    if eroot.name == "":
                        eroot.name = pnodename
                        equick_nodes[pnodename] = eroot

                    mynode = equick_nodes[pnodename]

                    myc1 = mynode.add_child()
                    myc2 = mynode.add_child()

                    myc1.name = c1name
                    myc2.name = c2name

                    myc1.dist = times[c1name] - times[pnodename]
                    myc2.dist = times[c2name] - times[pnodename]

                    equick_nodes[c1name] = myc1
                    equick_nodes[c2name] = myc2

        if family_size == 0:

            extanttree = ";"

        elif family_size == 1:


            extanttree = [k for k, v in surviving_nodes.items() if v["state"] == 1 and v["descendant"] == "None"][
                             0] + ";"

        else:

            extanttree = extanttree.write(format=1, format_root_node=True)


        rec = completetree.getTreeRecPhyloXML()

        if len(completetree) == 0:
            completetree = ";"
        elif len(completetree) == 1:
            completetree = completetree.get_leaves()[0].name + ";"
        else:
            completetree = completetree.write(format=1, format_root_node=True)


        return completetree, extanttree, rec


    def generate_oldtree(self):

        tree = ete3.Tree()

        current_time, event, nodes = self.events[0]

        sp = tree.get_tree_root()
        sp.name = nodes + "_1"
        sp.add_feature("is_active", True)

        elapsed_time = float(current_time)

        for current_time, event, nodes in self.events[1:]:

            elapsed_time = float(current_time) - elapsed_time
            active_nodes = [x for x in tree.get_leaves() if x.is_active == True]
            for node in active_nodes:
                node.dist += elapsed_time
            elapsed_time = float(current_time)

            if event == "S":

                sp, gp, c1, g1, c2, g2 = nodes.split(";")

                myname = sp + "_" + gp
                mynode = tree & myname
                mynode.is_active = False

                gc1 = mynode.add_child(dist=0)
                gc1.name = c1 + "_" + g1
                gc1.add_feature("is_active", True)

                gc2 = mynode.add_child(dist=0)
                gc2.name = c2 + "_" + g2
                gc2.add_feature("is_active", True)

            elif event == "E":
                sp, gp = nodes.split(";")
                myname = sp + "_" + gp
                mynode = tree & myname
                mynode.is_active = False

            elif event == "L":
                sp, gp = nodes.split(";")
                myname = sp + "_" + gp
                mynode = tree & myname
                mynode.is_active = False

            elif event == "D":

                sp, gp, c1, g1, c2, g2 = nodes.split(";")
                myname = sp + "_" + gp
                mynode = tree & myname

                mynode.is_active = False

                gc1 = mynode.add_child(dist=0)
                gc1.name = c1 + "_" + g1
                gc1.add_feature("is_active", True)

                gc2 = mynode.add_child(dist=0)
                gc2.name = c2 + "_" + g2
                gc2.add_feature("is_active", True)

            elif event == "T":
                sp, gp, c1, g1, c2, g2 = nodes.split(";")

                myname = sp + "_" + gp

                mynode = tree & myname
                mynode.is_active = False

                gc1 = mynode.add_child(dist=0)
                gc1.name = c1 + "_" + g1
                gc1.add_feature("is_active", True)

                gc2 = mynode.add_child(dist=0)
                gc2.name = c2 + "_" + g2
                gc2.add_feature("is_active", True)

            elif event == "F":
                break

        complete_tree = tree.write(format=1, format_root_node=True)
        active_nodes = [x for x in tree.get_leaves() if x.is_active == True]

        if len(active_nodes) < 3:
            pruned_tree = None

        else:
            tree.prune(active_nodes, preserve_branch_length=True)
            pruned_tree = tree.write(format=1, format_root_node=True)

        return complete_tree, pruned_tree

    def obtain_new_gene_id(self):
        self.gene_ids_counter += 1
        return self.gene_ids_counter

    def __str__(self):

        return "GeneFamily_" + str(self.identifier) + ";" + ";".join([str(x) for x in self.genes])

        #return ";".join(map(str, self.genes))

    def __len__(self):

        return len([x for x in self.genes if x.active == True])

    def __iter__(self):
        for gene in self.genes:
            yield gene




class Gene():

    def __init__(self):

        self.active = True
        self.orientation = ""
        self.gene_family = ""
        self.gene_id = ""
        self.sequence = ""
        self.species = ""
        self.importance = 0
        self.length = 0
        self.total_flanking = 0
        self.specific_flanking = 0

    def determine_orientation(self):

        if numpy.random.binomial(1,0.5):
            self.orientation = "+"

        else:
            self.orientation = "-"

    def change_sense(self):

        if self.orientation == "+":
            self.orientation = "-"
        elif self.orientation == "-":
            self.orientation = "+"

    def __str__(self):

        myname = "_".join(map(str, (self.species, self.gene_family, self.gene_id)))
        #myname = "_".join(map(str, (self.gene_family, self.orientation)))
        #myname = str(self.gene_family) + "_" + str(self.gene_id)
        #myname = "_".join(map(str, (self.gene_family, self.length)))
        return myname


class Intergene():

    def __init__(self):

        self.length = 0
        self.total_flanking = 0
        self.specific_flanking = 0
        self.id = 0 # Only for debugging purposes

    def __str__(self):

        return "(" + str(self.length) + ")"
        #return "I_" + str(self.length)
        #return "I_" + str(self.id) + "_" + str(self.length)


class Chromosome():

    def __init__(self):

        self.has_intergenes = False
        self.intergenes = list()
        self.genes = list()
        self.shape = ""
        self.length = 0

        self.total_locations = list()
        self.map_of_locations = list()

        self.total_rates = 0

    def obtain_total_itergenic_length(self):

        total_length = 0
        for intergene in self.intergenes:
            total_length += intergene.length
        return total_length

    def select_random_position(self):

        return numpy.random.randint(len(self.genes))

    def obtain_flankings(self):

        if self.has_intergenes:

            self.genes[0].total_flanking = (0, self.genes[0].length)
            self.genes[0].specific_flanking = (0, self.genes[0].length)
            self.intergenes[0].total_flanking = (self.genes[0].total_flanking[1], self.genes[0].total_flanking[1] + self.intergenes[0].length)
            self.intergenes[0].specific_flanking = (0, self.intergenes[0].length)

            for i in range(len(self.genes)):

                if i == 0:
                    continue

                lb = self.intergenes[i-1].total_flanking[1]
                ub = lb + self.genes[i].length

                lbg = self.genes[i-1].specific_flanking[1] + 1
                ubg = lbg + self.genes[i].length

                self.genes[i].total_flanking = (lb, ub)
                self.genes[i].specific_flanking = (lbg, ubg)

                lbi = self.intergenes[i - 1].specific_flanking[1] + 1
                ubi = lbi + self.intergenes[i].length

                self.intergenes[i].total_flanking = (ub, ub + self.intergenes[i].length)
                self.intergenes[i].specific_flanking = (lbi, ubi)

            #self.intergenes[i].total_flanking = (ub, 0)

    def obtain_locations(self):

        self.map_of_locations = list()

        # The structure of total location is:
        # tc1, tc2, Specific coordinate 1, specific coordinate 2, position, Gene/intergene

        for i in range(len(self.genes)):

            tc1 = self.genes[i].total_flanking[0]
            tc2 = self.genes[i].total_flanking[1]
            sc1 = self.genes[i].specific_flanking[0]
            sc2 = self.genes[i].specific_flanking[1]
            self.map_of_locations.append((tc1, tc2, sc1, sc2, str(i), "G"))

            tc1 = self.intergenes[i].total_flanking[0]
            tc2 = self.intergenes[i].total_flanking[1]
            sc1 = self.intergenes[i].specific_flanking[0]
            sc2 = self.intergenes[i].specific_flanking[1]
            self.map_of_locations.append((tc1, tc2, sc1, sc2, str(i), "I"))


    def select_random_coordinate_in_intergenic_regions(self):

        # We weight the position by the length of the region

        t = sum([x.length for x in self.intergenes]) + len(self.intergenes) - 1
        return random.randint(0, int(t))

    def return_total_coordinate_from_specific_coordinate(self, c, type = "I", debug = False):

        tc = None
        for r in self.map_of_locations:
            tc1, tc2, spc1, spc2, sp, t = r
            if debug == True:
                print(r)
            if t != type:
                continue
            if c >= spc1 and c <= spc2:
                distance_to_lower_bound = c - spc1
                tc = tc1 + distance_to_lower_bound
        return tc

    def return_specific_coordinate_from_total_coordinate(self, c, debug = False):


        sc = None
        for r in self.map_of_locations:
            tc1, tc2, spc1, spc2, sp, t = r
            if debug == True:
                print(r)
            if t == "I" and c >= tc1 and c <= tc2:

                distance_to_lower_bound = c - tc1
                sc = spc1 + distance_to_lower_bound

                if debug == True:
                    print("PRINTING R")
                    print(r)

        return sc

    def return_location_by_coordinate(self, c, within_intergene = False):

        ### Returns
        ### 1 and 2 Limits for total coordinates
        ### 3 and 4 Limits for specific coordinates
        ### 5  Position in the list of genes or itergenes
        ### 6  Intergene (I) or Gene(G)

        if within_intergene == False:

            for l in self.map_of_locations:
                tc1, tc2, spc1, spc2, sp, t = l
                if c >= tc1 and c <= tc2:
                    return l
        else:

            for l in self.map_of_locations:
                tc1, tc2, spc1, spc2, sp, t = l
                if t != "I":
                    continue
                if c >= spc1 and c <= spc2:
                    return l

    def return_affected_region(self, c1, c2, direction):

        # It returns a tuple
        # 1. List of the position of the genes affected. ALWAYS FROM LEFT TO RIGHT
        # 2. List of the position of the intergenes affected. ALWAYS FROM LEFT TO RIGHT
        # 3. Tuple with left and right cuts of first intergene.
        # #  Watch out, the fact of calling it left or right can be confusing.
        # 4. Tuple with left and right cuts of last intergene. Same note than above

        l1 = self.return_location_by_coordinate(c1, within_intergene=True)
        l2 = self.return_location_by_coordinate(c2, within_intergene=True)

        tc1_1, tc1_2, sc1_1, sc1_2, p1, t1, = l1
        tc2_1, tc2_2, sc2_1, sc2_2, p2, t2 = l2

        p1 = int(p1)
        p2 = int(p2)

        affected_genes = list()
        affected_intergenes = list()

        t_length = len(self.intergenes)

        left_limits = None
        right_limits = None

        if c1 == c2:
            return None

        elif p1 == p2:
            return None

        elif c1 < c2 and direction == "right":

            affected_genes = [i + 1 for i in range(p1, p2)]
            affected_intergenes = [i for i in range(p1,p2 + 1)]
            left_limits = (c1 - sc1_1, sc1_2 - c1)
            right_limits = (c2 - sc2_1, sc2_2 - c2)

        elif c1 > c2 and direction == "right":

            affected_genes = [i + 1 for i in range(p1, t_length - 1)]
            affected_genes += [i for i in range(0, p2 + 1)]

            affected_intergenes = [i for i in range(p1, t_length)]
            affected_intergenes += [i for i in range(0, p2 + 1)]

            left_limits = (c1 - sc1_1, sc1_2 - c1)
            right_limits = (c2 - sc2_1, sc2_2 - c2)

        elif c1 > c2 and direction == "left":

            affected_genes = [i for i in range(p1, p2, - 1)]
            affected_intergenes = [i for i in range(p1, p2 - 1, -1)]
            left_limits = (c1 - sc1_1, sc1_2 - c1)
            right_limits = (c2 - sc2_1, sc2_2 - c2)

            affected_genes.reverse()
            affected_intergenes.reverse()

        elif c1 < c2 and direction == "left":

            affected_genes = [i for i in range(p1, -1, - 1)]
            affected_genes += [i for i in range(t_length - 1, p2, -1)]

            affected_intergenes = [i for i in range(p1, -1, -1)]
            affected_intergenes += [i for i in range(t_length - 1, p2 - 1, -1)]

            affected_genes.reverse()
            affected_intergenes.reverse()

            left_limits = (c1 - sc1_1, sc1_2 - c1)
            right_limits = (c2 - sc2_1, sc2_2 - c2)

        return (affected_genes, affected_intergenes, left_limits, right_limits)


    def select_random_length(self, p):

        return int(af.obtain_value(p))

    def return_rates(self):

        # NOT WORKING FOR NOW

        self.total_rates = 0.0

        for gene in self.genes:
            d = gene.gene_family.rates["DUPLICATION"] ## GENE.GENE_FAMILY SHOULD POINT TO A GENE FAMILY OBJECT, NOT A STR
            t = gene.gene_family.rates["TRANSFER"]
            l = gene.gene_family.rates["LOSS"]

            self.total_rates += d + t +l

        return self.total_rates



    def __len__(self):

        # Watch out!! This is probably no the safest thing to do
        if not self.has_intergenes:
            return len(self.genes)
        else:
            return self.length

    def __str__(self):

        if self.has_intergenes == True:

            #return ";".join(["CHROMOSOME"] + [str(self.genes[i])+";"+str(self.intergenes[i]) for i in range(len(self.genes))])
            return "".join(["CHROMOSOME: "] + [str(self.genes[i])+str(self.intergenes[i]) for i in range(len(self.genes))])

        else:

            return ";".join(["CHROMOSOME"] + [str(gene) for gene in self.genes])

    def __iter__(self):

        for x in self.genes:
            yield x


class CircularChromosome(Chromosome):

    def __init__(self):
        super().__init__()

    def obtain_segment(self, affected_genes):

        segment = [self.genes[x] for x in affected_genes]

        return segment

    def obtain_intergenic_segment(self, affected_intergenes):

        segment = [self.intergenes[x] for x in affected_intergenes]

        return segment

    def remove_segment(self, segment):

        for gene in segment:
            self.genes.remove(gene)

    def remove_intersegment(self, intersegment):

        for intergene in intersegment:
            self.intergenes.remove(intergene)

    def insert_segment(self, position, segment):

        for i, x in enumerate(segment):
            self.genes.insert(position + i, x)

    def invert_segment(self, affected_genes):

        segment = [self.genes[x] for x in affected_genes]

        reversed_segment = segment[::-1]

        for gene in reversed_segment:
            gene.change_sense()

        for i,x in enumerate(affected_genes):
            self.genes[x] = reversed_segment[i]


    def cut_and_paste(self, affected_genes):

        segment = [self.genes[x] for x in affected_genes]
        new_segment = list()

        if len(segment) == len(self.genes):
            return 0

        for gene in segment:
            new_segment.append(self.genes.pop(self.genes.index(gene)))

        position = self.select_random_position()
        for i, gene in enumerate(new_segment):
            self.genes.insert(position + i, gene)


    def obtain_affected_genes(self, p_extension):

        # Returns the index list of the affected genes

        position = self.select_random_position()
        length = self.select_random_length(p_extension)
        total_length = len(self.genes)
        affected_genes = list()

        if length >= total_length:
            affected_genes = [x for x in range(total_length)]
            return affected_genes

        for i in range(position, position + length):
            if i >= total_length:
                affected_genes.append(i - total_length)
            else:
                affected_genes.append(i)
        return affected_genes

    def obtain_affected_genes_accounting_for_family_rates(self, p_extension, gene_families, mrate):

        # In this first version, length is 1. For a more advanced version, I should extent the interactome model
        # Returns N genes accounting for the family rates


        gene2rate = {gene: gene_families[gene.gene_family].rates[mrate] for gene in self.genes}

        if p_extension == 1:

            norm = af.normalize([vl for x, vl in gene2rate.items()])
            mgenes = [i for i in range(len(self.genes))]
            affected_genes = numpy.random.choice(mgenes,size = 1,p=norm)
            return affected_genes

        else:

            # If there is an extension

            length = self.select_random_length(p_extension)
            total_length = len(self.genes)

            if length >= total_length:
                affected_genes = [x for x in range(total_length)]
                return affected_genes
            else:

                all_weights = list()

                # If the extension is shorter than the whole genome length

                for each_start, gene in enumerate(self.genes):

                    position = each_start
                    affected_genes = list()

                    for i in range(position, position + length):
                        if i >= total_length:
                            affected_genes.append(i - total_length)
                        else:
                            affected_genes.append(i)

                    all_weights.append(reduce(lambda x, y: x * y, [gene2rate[self.genes[x]] for x in affected_genes]))
                #print(all_weights)
                position = numpy.random.choice([i for i, g in enumerate(self.genes)], 1, p=af.normalize(all_weights))[0]
                affected_genes = list()

                # Returns the index list of the affected genes

                for i in range(position, position + length):
                    if i >= total_length:
                        affected_genes.append(i - total_length)
                    else:
                        affected_genes.append(i)


                return affected_genes


    def obtain_affected_genes_accounting_for_connectedness(self, p_extension, interactome):

        # Returns N genes accounting for the inverse of the connectedness

        node_degrees = {n:d for n,d in interactome.degree()}

        corrected_node_degrees = list()
        all_weights = list()

        for gene in self.genes:

            corrected_node_degrees.append(1/(node_degrees[str(gene)] + 1))

        length = self.select_random_length(p_extension)
        total_length = len(self.genes)

        for each_start, gene in enumerate(self.genes):
            position = each_start
            affected_genes = list()

            if length >= total_length:
                # We select the whole genome
                affected_genes = [x for x in range(total_length)]
                return affected_genes

            else:

                for i in range(position, position + length):
                    if i >= total_length:
                        affected_genes.append(i - total_length)
                    else:
                        affected_genes.append(i)

                # We obtain the total weight of this option
                # This means, multiplying all the weights if we start in a given position

                all_weights.append(reduce(lambda x, y: x * y, [corrected_node_degrees[x] for x in affected_genes]))

        position = numpy.random.choice([i for i,g in enumerate(self.genes)], 1, p=af.normalize(all_weights))[0]
        affected_genes = list()

        # Returns the index list of the affected genes

        for i in range(position, position + length):
            if i >= total_length:
                affected_genes.append(i - total_length)
            else:
                affected_genes.append(i)
        return affected_genes


    def get_homologous_position(self, segment):

        homologous = list()

        segment_length = len(segment)
        genes_length = len(self.genes)

        genes = [x.gene_family + "_" + x.orientation for x in self.genes]
        mysegment = [x.gene_family + "_" + x.orientation for x in segment]


        # First we traverse the genome forwards

        for i, gene in enumerate(genes):

            length_counter = 0
            positions = list()

            name_gene_in_genome = gene
            name_gene_in_segment = mysegment[0]

            if name_gene_in_genome == name_gene_in_segment:

                positions.append(i)
                length_counter += 1

                for j, x in enumerate(mysegment):
                    if length_counter == segment_length:
                        homologous.append(("F", tuple(positions)))
                        break
                    if 1 + i + j >= genes_length:
                        if genes[(i + j + 1) - genes_length] == mysegment[j + 1]:
                            positions.append(i+j+1 - genes_length)
                            length_counter += 1
                        else:
                            break
                    else:
                        if genes[i + j + 1] == mysegment[j + 1]:
                            positions.append(i + j + 1)
                            length_counter += 1
                        else:
                            break

        # Second we traverse the genome backwards

        inverted_segment = list()

        for gene in mysegment[::-1]:

            inverted_segment.append(gene.replace("+","A").replace("-", "+").replace("A", "-"))

        for i, gene in enumerate(genes):

            positions = list()
            length_counter = 0
            name_gene_in_segment = inverted_segment[0]

            if genes[i] == name_gene_in_segment:

                positions.append(i)
                length_counter += 1

                for j, x in enumerate(inverted_segment):
                    if length_counter == segment_length:
                        homologous.append(("B", tuple(positions)))
                        break
                    if 1 + i + j >= genes_length:
                        if genes[(i + j + 1) - genes_length] == inverted_segment[j + 1]:
                            positions.append(i + j + 1 - genes_length)
                            length_counter += 1
                        else:
                            break
                    else:
                        if genes[i + j + 1] == inverted_segment[j + 1]:
                            positions.append(i + j + 1)
                            length_counter += 1
                        else:
                            break

        return homologous

    ### From here, functions related to intergenic regions

    def remove_segment_with_intergenic(self, segment):

        for gene in segment:
            self.genes.remove(gene)

    def insert_gene_within_intergene(self, coordinate, location, gene):

        tc1, tc2, sc1, sc2, position, t = location

        # Convert to ints:

        tc1, tc2, sc1, sc2, position = map(int,(tc1,tc2,sc1,sc2,position))

        # The first part is easier - We simply add the gene to the list of genes

        self.genes.insert(position + 1, gene)

        # The second part is cutting the intergene, obtaining the distances

        left_limits = coordinate - sc1
        right_limits = sc2 - coordinate

        # Now we insert the new intergene in the position i + 1

        intergene = Intergene()
        intergene.length = right_limits
        self.intergenes[position].length = left_limits
        self.intergenes.insert(position + 1, intergene)


class LinearChromosome(Chromosome):

    pass

class Genome():

    def __init__(self):

        self.species = ""
        self.chromosomes = list()

    def start_genome(self, input):

        for size, shape in input:

            if shape == "L":
                self.chromosomes.append(LinearChromosome(size))
            elif shape == "C":
                self.chromosomes.append(CircularChromosome(size))

    def select_random_chromosome(self):

        # I have to weight by the length of each chromosome

        #chromosome = numpy.random.choice(self.chromosomes, 1, p=af.normalize([len(x) for x in self.chromosomes]))[0]

        # So far, only one chromosome per genome, I can safely return the first chromosome

        return self.chromosomes[0]

    def update_genome_species(self, species):

        self.species = species

        for ch in self.chromosomes:

            for gene in ch:
                gene.species = species

    def create_interactome(self, network_model = "BA"):

        import networkx as nx
        import random

        if network_model == "BA":
            self.interactome  = nx.barabasi_albert_graph(len(self.chromosomes[0]), 1)
        else:
            self.interactome = nx.barabasi_albert_graph(len(self.chromosomes[0]), 1)

        ## Need to shuffle the nodes!!

        randomly_ordered_genes = list(self.chromosomes[0].genes)
        random.shuffle(randomly_ordered_genes)

        self.interactome = nx.relabel_nodes(self.interactome, {i:str(n) for i,n in enumerate(randomly_ordered_genes)})

    def __str__(self):

        return ";".join(["GENOME:"] + [str(x) for x in self.chromosomes])

    def __iter__(self):
        for chromosome in self.chromosomes:
            yield chromosome
