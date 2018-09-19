import AuxiliarFunctions as af
import ete3
import numpy
import copy
import random
import os

from GenomeClasses import GeneFamily, Gene, Intergene, CircularChromosome, LinearChromosome, Genome

class GenomeSimulator():

    def __init__(self, parameters, events_file):

        self.parameters = parameters
        self.tree_events = self._read_events_file(events_file)

        self.all_genomes = dict()
        self.all_gene_families = dict()

        self.gene_families_counter = 0
        self.active_genomes = set()


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

    def write_events_per_branch(self, events_per_branch_folder):

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


    def write_profiles(self, profiles_folder):

        if not os.path.isdir(profiles_folder):
            os.mkdir(profiles_folder)


        genome_names = [x for x in self.all_genomes.keys()]

        with open(os.path.join(profiles_folder, "Profiles.tsv"), "w") as f:

            header = ["FAMILY"] + genome_names
            header = "\t".join(map(str, header)) + "\n"

            f.write(header)

            for gene_family_name, gene_family in self.all_gene_families.items():

                if self.parameters["VERBOSE"] == 1:
                    print("Writing profile for family %s" % str(gene_family_name))

                line = ["Fam" + gene_family_name]
                for genome in genome_names:
                    n = 0
                    for gene in gene_family:
                        if gene.species == genome:
                            n+=1
                    line.append(str(n))

                line = "\t".join(line)+"\n"

                f.write(line)

    def _read_events_file(self, events_file):

        events = list()
        with open(events_file) as f:
            f.readline()
            for line in f:
                handle = line.strip().split("\t")
                events.append(handle)
        return events

    def return_new_identifiers_for_segment(self, segment):

        new_identifiers = list()

        for gene in segment:
            gf = gene.gene_family
            new_id = self.all_gene_families[gf].obtain_new_gene_id()
            new_identifiers.append(new_id)

        return new_identifiers

    def fill_genome(self, intergenic_sequences = False):

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

            for i in range(int(n_genes)):

                # We fill the chromosomes and we create also the gene families

                gene, gene_family = self.make_origination(genome.species, time)

                chromosome.genes.append(gene)
                self.all_gene_families[str(self.gene_families_counter)] = gene_family

                if intergenic_sequences == True:

                    gene.length = int(af.obtain_value(self.parameters["GENE_SIZE"]))
                    intergenic_sequence = Intergene()
                    intergenic_sequence.length = int(af.obtain_value(self.parameters["INTERGENIC_LENGTH"]))
                    chromosome.intergenes.append(intergenic_sequence)

            if intergenic_sequences == True:
                chromosome.obtain_flankings()
                chromosome.obtain_locations()

            genome.chromosomes.append(chromosome)
        return genome


    def run(self):

        d = af.obtain_value(self.parameters["DUPLICATION"])
        t = af.obtain_value(self.parameters["TRANSFER"])
        l = af.obtain_value(self.parameters["LOSS"])
        i = af.obtain_value(self.parameters["INVERSION"])
        c = af.obtain_value(self.parameters["TRANSLOCATION"])

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

            time_to_next_genome_event = self.get_time_to_next_event(len(self.active_genomes), d, t, l, i, c, o)

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

    def run_s(self):

        d = af.obtain_value(self.parameters["DUPLICATION"])
        t = af.obtain_value(self.parameters["TRANSFER"])
        l = af.obtain_value(self.parameters["LOSS"])
        i = af.obtain_value(self.parameters["INVERSION"])
        c = af.obtain_value(self.parameters["TRANSLOCATION"])
        o = af.obtain_value(self.parameters["ORIGINATION"])

        # First we prepare the first genome

        genome = self.fill_genome()

        # We add the importance of the genes

        for gene in genome:
            gene.importance = af.obtain_value(self.parameters["GENE_IMPORTANCE"])

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

            time_to_next_genome_event = self.get_time_to_next_event(len(self.active_genomes), d, t, l, i, c, o)

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
                self.evolve_genomes_selection(d, t, l, i, c, o, current_time)

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

    def run_i(self):

        d = af.obtain_value(self.parameters["DUPLICATION"])
        t = af.obtain_value(self.parameters["TRANSFER"])
        l = af.obtain_value(self.parameters["LOSS"])
        i = af.obtain_value(self.parameters["INVERSION"])
        c = af.obtain_value(self.parameters["TRANSLOCATION"])
        o = af.obtain_value(self.parameters["ORIGINATION"])

        # First we prepare the first genome

        genome = self.fill_genome(intergenic_sequences=True)

        ## These two lines are important for this mode

        for chromosome in genome:
            chromosome.obtain_flankings()

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

            time_to_next_genome_event = self.get_time_to_next_event(len(self.active_genomes), d, t, l, i, c, o)

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
                self.advanced_evolve_genomes_i(d, t, l, i, c, o, current_time)

    def generate_new_rates(self):

        d = af.obtain_value(self.parameters["DUPLICATION"])
        t = af.obtain_value(self.parameters["TRANSFER"])
        l = af.obtain_value(self.parameters["LOSS"])
        i = af.obtain_value(self.parameters["INVERSION"])
        c = af.obtain_value(self.parameters["TRANSLOCATION"])
        o = af.obtain_value(self.parameters["ORIGINATION"])

        return d,t,l,i,c,o

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
                self.branch_extension_rates[sp] = tuple([float(x) for x in (d, t, l, i, c)])

        with open(os.path.join(rates_folder, "Transfer_rates.tsv")) as f:
            f.readline()
            for line in f:
                dn, rc, wt  = line.split("\t")
                if dn not in self.transfer_rates:
                    self.transfer_rates[dn] = dict()
                if rc not in self.transfer_rates[dn]:
                    self.transfer_rates[dn][rc] = 0.0

                self.transfer_rates[dn][rc] = float(wt)




    def choose_event(self, duplication, transfer, loss, inversion, translocation, origination):

        draw = numpy.random.choice(["D", "T", "L", "I", "C", "O"], 1,
                                   p=af.normalize([duplication, transfer, loss, inversion, translocation, origination]))
        return draw

    def choose_recipient(self, lineages_alive, donor):
        possible_recipients = [x for x in lineages_alive if x != donor]
        if len(possible_recipients) > 1:
            recipient = random.choice(possible_recipients)
            return recipient
        else:
            return None

    def evolve_genomes(self, duplication, transfer, loss, inversion, translocation, origination, time):

        d_e = af.obtain_value(self.parameters["DUPLICATION_EXTENSION"])
        t_e = af.obtain_value(self.parameters["TRANSFER_EXTENSION"])
        l_e = af.obtain_value(self.parameters["LOSS_EXTENSION"])
        i_e = af.obtain_value(self.parameters["INVERSION_EXTENSION"])
        c_e = af.obtain_value(self.parameters["TRANSLOCATION_EXTENSION"])

        lineage = random.choice(list(self.active_genomes))
        event = self.choose_event(duplication, transfer, loss, inversion, translocation, origination)

        if event == "D":
            self.make_duplication(d_e, lineage, time)
            return "D", lineage

        elif event == "T":

            # We choose a recipient

            possible_recipients = [x for x in self.active_genomes if x != lineage]

            if len(possible_recipients) > 0:


                recipient = random.choice(possible_recipients)
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

        elif event == "C":
            self.make_translocation(c_e, lineage, time)
            return "C",lineage

        elif event == "O":

            gene, gene_family = self.make_origination(lineage, time)

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

        elif event == "C":
            self.make_translocation(c_e, lineage, time)
            return "C",lineage

        elif event == "O":

            gene, gene_family = self.make_origination(lineage, time)

            chromosome = self.all_genomes[lineage].select_random_chromosome()
            position = chromosome.select_random_position()
            segment = [gene]
            chromosome.insert_segment(position, segment)

            return "O", lineage

    def advanced_evolve_genomes_i(self, duplication, transfer, loss, inversion, translocation, origination, time):

        extension_multiplier = self.parameters["EXTENSION_MULTIPLIER"]

        d_e = af.obtain_value(self.parameters["DUPLICATION_EXTENSION"])
        t_e = af.obtain_value(self.parameters["TRANSFER_EXTENSION"])
        l_e = af.obtain_value(self.parameters["LOSS_EXTENSION"])
        i_e = af.obtain_value(self.parameters["INVERSION_EXTENSION"])
        c_e = af.obtain_value(self.parameters["TRANSLOCATION_EXTENSION"])

        lineage = random.choice(list(self.active_genomes))
        event = self.choose_event(duplication, transfer, loss, inversion, translocation, origination)

        if event == "D":
            self.make_duplication(d_e, lineage, time)
            return "D", lineage

        elif event == "T":

            # We choose a recipient

            possible_recipients = [x for x in self.active_genomes if x != lineage]

            if len(possible_recipients) > 0:

                recipient = random.choice(possible_recipients)
                donor = lineage
                self.make_transfer(t_e, donor, recipient, time)
                return "T", donor + "->" + recipient

            else:
                return None

        elif event == "L":

            self.make_loss(l_e, lineage, time)
            return "L", lineage

        elif event == "I":

            self.make_inversion(i_e, lineage, time)
            return "I", lineage

        elif event == "C":
            self.make_translocation(c_e, lineage, time)
            return "C", lineage

        elif event == "O":

            gene, gene_family = self.make_origination(lineage, time)
            chromosome = self.all_genomes[lineage].select_random_chromosome()
            chromosome.obtain_flankings()
            chromosome.obtain_locations()
            intergene_coordinate = chromosome.select_random_coordinate_in_intergenic_regions()
            location = chromosome.return_location_by_coordinate(intergene_coordinate, within_intergene=True)
            chromosome.insert_gene_within_intergene(intergene_coordinate, location, gene)

            return "O", lineage



    def get_time_to_next_event(self, n, d, t, l ,i , c, o):

        total = 0.0
        for j in range(n):
            total += sum((d, t, l, i, c, o))

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

    def increase_distances(self, time_to_next_event, active_lineages):

        for node in active_lineages:
            node.dist += time_to_next_event


    def make_origination(self, species_tree_node, time):

        self.gene_families_counter += 1
        gene_family_id = str(self.gene_families_counter)

        gene = Gene()
        gene.determine_orientation()

        gene.gene_family = str(self.gene_families_counter)
        gene.species = species_tree_node

        gene_family = GeneFamily(gene_family_id, time)
        gene_family.length = int(af.obtain_value(self.parameters["GENE_SIZE"]))
        gene.length = gene_family.length

        gene_family.genes.append(gene)
        gene.gene_id = gene_family.obtain_new_gene_id()

        self.all_gene_families[gene_family_id] = gene_family
        self.all_gene_families[gene.gene_family].register_event(str(time), "O", species_tree_node)

        return gene, gene_family

    def make_speciation(self, sp, c1, c2, time):

        # This function receives a genome and the names of the two branching lineages of the species node

        genome_sp = self.all_genomes[sp]

        genome1 = Genome()
        genome2 = Genome()

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

                #new_gene1.selection_coefficient = gene.selection_coefficient
                #new_gene2.selection_coefficient = gene.selection_coefficient

                gene_family = self.all_gene_families[gene.gene_family]
                gene_family.genes.append(new_gene1)
                gene_family.genes.append(new_gene2)

                new_gene1.gene_family = gene.gene_family
                new_gene2.gene_family = gene.gene_family

                ch1.genes.append(new_gene1)
                ch2.genes.append(new_gene2)

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


    def make_duplication(self, p, lineage, time):

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

        for i, gene in enumerate(segment):

            nodes = [gene.species,
                     gene.gene_id,
                     copied_segment1[i].species,
                     copied_segment1[i].gene_id,
                     copied_segment2[i].species,
                     copied_segment2[i].gene_id]

            gene.active = False

            self.all_gene_families[gene.gene_family].register_event(time, "D", ";".join(map(str, nodes)))

    def choose_precise_distance_recipient(self, time, possible_recipients, donor):

        ### Deprecated

        weights = list()
        mydonor = self.distances_to_root[donor][0]

        for recipient in possible_recipients:

            myrecipient = self.distances_to_root[recipient][0]
            phylo_d = mydonor.get_distance(myrecipient)
            td = phylo_d + (2 * time) - self.distances_to_root[donor][1] - self.distances_to_root[recipient][1]
            weights.append(td)

        draw = numpy.random.choice(possible_recipients, 1, p=af.normalize(af.inverse(weights)))[0]

        return draw

    def choose_advanced_recipient(self, possible_recipients, donor):

        weights = list()

        for recipient in possible_recipients:
            weights.append(self.transfer_rates[donor][recipient])

        if sum(weights) == 0:
            return None


        draw = numpy.random.choice(possible_recipients, 1, p=af.normalize(weights))[0]

        return draw

    def make_transfer(self, p, donor, recipient, time):

        chromosome1 = self.all_genomes[donor].select_random_chromosome()
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


    def make_loss(self, p, lineage, time):

        chromosome = self.all_genomes[lineage].select_random_chromosome()
        affected_genes = chromosome.obtain_affected_genes(p)
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

    def make_loss_intergenic(self, p, lineage, time):

        chromosome = self.all_genomes[lineage].select_random_chromosome()
        affected_genes = chromosome.obtain_affected_genes(p)
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


    def make_loss_selection_mode(self, p, lineage, time):

        chromosome = self.all_genomes[lineage].select_random_chromosome()
        affected_genes = chromosome.obtain_affected_genes(p)
        segment = chromosome.obtain_segment(affected_genes)

        # Now we check we are not under the minimum size

        if len(chromosome) - len(affected_genes) <= 0:
            return 0

        importance = sum([gene.importance for gene in segment])
        if importance >= 0:
            if numpy.random.exponential():
                pass

            return 0

        chromosome.remove_segment(segment)

        # We have to register in the affected gene families that there has been as loss
        # All genes affected must be returned

        for gene in segment:
            gene.active = False
            self.all_gene_families[gene.gene_family].register_event(time, "L", ";".join(map(str,[lineage, gene.gene_id])))

    def make_inversion(self, p, lineage, time):

        chromosome = self.all_genomes[lineage].select_random_chromosome()
        affected_genes = chromosome.obtain_affected_genes(p)
        segment = chromosome.obtain_segment(affected_genes)
        chromosome.invert_segment(affected_genes)

        for i, gene in enumerate(segment):
            self.all_gene_families[gene.gene_family].register_event(str(time), "I", ";".join(map(str,[lineage, gene.gene_id])))

    def make_inversion_intergenic(self, p, lineage, time):

        chromosome = self.all_genomes[lineage].select_random_chromosome()
        affected_genes = chromosome.obtain_affected_genes(p)
        segment = chromosome.obtain_segment(affected_genes)
        chromosome.invert_segment(affected_genes)

        for i, gene in enumerate(segment):
            self.all_gene_families[gene.gene_family].register_event(str(time), "I", ";".join(map(str,[lineage, gene.gene_id])))

    def make_translocation(self, p, lineage, time):

        chromosome = self.all_genomes[lineage].select_random_chromosome()
        affected_genes = chromosome.obtain_affected_genes(p)
        segment = chromosome.obtain_segment(affected_genes)
        chromosome.cut_and_paste(affected_genes)

        for i, gene in enumerate(segment):
            self.all_gene_families[gene.gene_family].register_event(str(time), "C", ";".join(map(str,[lineage, gene.gene_id])))

    def get_gene_family_tree(self):

        if len(self.gene_family["Gene_tree"].get_leaves()) < 3:
            return "None"
        else:
            return self.gene_family["Gene_tree"].write(format=1)