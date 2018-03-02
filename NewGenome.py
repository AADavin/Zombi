import AuxiliarFunctions as af
import numpy
import copy
import ete3
import random
import os

numpy.random.seed(245)
random.seed(10)


def read_parameters(parameters_file):

    parameters = dict()

    with open(parameters_file) as f:
        for line in f:
            if line[0] == "#" or line == "\n":
                continue
            if "\t" in line:
                parameter, value = line.strip().split("\t")
                parameters[parameter] = value
            elif " " in line:
                parameter, value = line.strip().split(" ")
                parameters[parameter] = value

    return parameters


def obtain_value(value):

    handle = value.split(":")[0]

    if handle[0] == "f":
        # Fixed value
        value =  float(handle[1])

    if handle[0] == "n":
        # normal distribution
        value =  float(handle[1])

    if handle[0] == "l":
        # lognormal distribution
        value =  float(handle[1])

    return value


def prepare_species_tree_parameters(parameters):

    for parameter, value in parameters:

        if parameter == "SPECIATION":
            parameters[parameter] = obtain_value(value)
        if parameter == "EXTINCTION":
            parameters[parameter] = obtain_value(value)

def prepare_parameters(parameters):

    for parameter, value in parameters:

        if parameter == "DUPLICATION" or parameter == "TRANSFER" or parameter == "LOSSES" or \
                parameter == "INVERSION" or parameter == "TRANSLOCATION" or parameter == "ORIGINATION":
            parameters[parameter] = obtain_value(value)

        if parameter == "DUPLICATION_EXTENSION" or parameter == "TRANSFER_EXTENSION" \
                or parameter == "LOSSES_EXTENSION" or parameter == "INVERSION_EXTENSION" or \
                parameter == "TRANSLOCATION_EXTENSION" or parameter == "ORIGINATION_EXTENSION":

            parameters[parameter] = obtain_value(value)

def generate_events(tree_file): # I have to round distances

    events = []
    lineage_counter = 0

    with open(tree_file) as f:

        tree = ete3.Tree(f.readline().strip(), format=1)

    root = tree.get_tree_root()
    root.name = "Root"
    nodes = [(x, int(round(x.get_distance(root)))) for x in tree.traverse()]

    # We get the time for present time

    total_time = int(round(root.get_farthest_leaf()[1]))

    # We order the events # This can be more efficient. But who cares
    sorted(nodes, key = lambda x: x[1])
    for node, time in nodes:
        lineage_counter += 1
        if not node.is_leaf():
            node.name = "n" + str(lineage_counter)

    for node, time in nodes:

        if node.is_leaf():

            if time == total_time:
                event = "F"
            elif time != total_time:
                event = "E"
            else:
                event = "ERROR"
            events.append((str(time), event, node.name))
        else:
            c1, c2 = node.get_children()
            event = "S"
            events.append((str(time), event, ";".join((node.name,c1.name,c2.name))))

    for event in events:
        print(event)


def read_parameters(parameters_file):

    parameters = dict()

    with open(parameters_file) as f:
        for line in f:
            line.strip().split("\t")

    return parameters


def copy_segment(segment, new_identifiers):

    new_segment = list()

    for i,gene in enumerate(segment):
        new_gene = copy.deepcopy(gene)
        new_gene.gene_id = new_identifiers[i]
        new_segment.append(new_gene)

    return new_segment


class GeneFamily():

    def __init__(self, identifier, time):

        self.identifier = identifier
        self.origin = time

        self.genes = list()
        self.events = list()
        self.event_counter = 0  # Each time that the family is modified in any form, we have to update the event counter
        self.gene_ids_counter = 0

    def register_event(self, time, event, genes):

        self.events.append((time, event, genes))

    def generate_tree(self):

        tree = ete3.Tree()

        sp = tree.get_tree_root()
        sp.name = "Root_1"
        sp.add_feature("is_active", True)

        elapsed_time = 0

        for current_time, event, nodes in self.events:

            elapsed_time = float(current_time) - elapsed_time
            active_nodes = [x for x in tree.get_leaves() if x.is_active == True]
            for node in active_nodes:
                node.dist += elapsed_time
            elapsed_time = float(current_time)

            if event == "S":

                sp, gp, c1, g1, c2, g2 = nodes.split(";")

                myname = sp + "_" + gp
                mynode = tree&myname
                mynode.is_active = False

                gc1 = mynode.add_child(dist=0)
                gc1.name = c1 + "_" + g1
                gc1.add_feature("is_active",True)

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

        whole_tree = tree.write(format=1)
        active_nodes = [x for x in tree.get_leaves() if x.is_active == True]

        if len(active_nodes) < 3:
            pruned_tree = None

        else:
            tree.prune(active_nodes, preserve_branch_length=True)
            pruned_tree = tree.write(format=1)

        return whole_tree, pruned_tree

    def obtain_new_gene_id(self):
        self.gene_ids_counter += 1
        return self.gene_ids_counter

    def __str__(self):

        return "GeneFamily_" + str(self.identifier) + ";" + ";".join([str(x) for x in self.genes])

        #return ";".join(map(str, self.genes))

    def __len__(self):

        return len([x for x in self.genes if x.active == True])

class Gene():

    def __init__(self):

        self.active = True
        self.orientation = ""
        self.gene_family = ""
        self.gene_id = ""
        self.sequence = ""
        self.species = ""
        self.selection_coefficient = ""

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
        #myname = "_".join(map(str,(self.genome, self.orientation, self.gene_family, self.gene_id)))
        #myname = "_".join(map(str, (self.genome, self.gene_family, self.gene_id)))
        myname = "_".join(map(str, (self.orientation, self.gene_family)))
        return myname

class Chromosome():

    def __init__(self):

        self.genes = list()
        self.intergene_distances = list()
        self.shape = ""

    def select_random_position(self):

        return numpy.random.randint(len(self.genes))

    def select_random_length(self, p):

        return numpy.random.geometric(p)

    def __len__(self):

        return len(self.genes)

    def __str__(self):

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

    def remove_segment(self, segment):

        for gene in segment:
            self.genes.remove(gene)

    def insert_segment(self, position, segment):

        for i, x in enumerate(segment):
            self.genes.insert(position + i, x)

    def invert_segment(self, affected_genes):

        new_genes = list()
        inverted_genes = list()

        for i in range(len(self.genes)):
            if i not in affected_genes:
                new_genes.append(self.genes[i])
        for i in affected_genes:
            inverted_genes.append(self.genes[i])
        for gene in inverted_genes:
            gene.change_sense()
        self.genes = new_genes + inverted_genes[::-1]

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



    def get_homologous_position(self, segment):

        # This function returns if there is an homologous position to the segment transferred

        for gene in self.genes:

            pass




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

        chromosome = numpy.random.choice(self.chromosomes, 1, p=af.normalize([len(x) for x in self.chromosomes]))[0]
        return chromosome

    def update_genome_species(self, species):

        self.species = species

        for ch in self.chromosomes:

            for gene in ch:
                gene.species = species


    def __str__(self):

        return ";".join(["GENOME:"] + [str(x) for x in self.chromosomes])

    def __iter__(self):
        for chromosome in self.chromosomes:
            yield chromosome


class GenomeSimulator():

    def __init__(self, events_file):

        self.parameters = dict()

        self.parameters["DUPLICATION_E"] = 0.2
        self.parameters["LOSS_E"] = 0.8
        self.parameters["TRANSFER_E"] = 0.8
        self.parameters["INVERSION_E"] = 0.1
        self.parameters["TRANSLOCATION_E"] = 0.1
        self.parameters["GENOME_FILE"] = "/Users/adriandavin/Desktop/Bioinformatics/SimuLyon/Cedric/Gillespie/Genome_file.tsv"

        self.tree_events = self._read_events_file(events_file)

        self.all_genomes = dict()
        self.all_gene_families = dict()

        self.gene_families_counter = 0

        self.active_genomes = set()


    def write_genomes(self, genome_folder):

        if not os.path.isdir(genome_folder):
            os.mkdir(genome_folder)

        for genome_name,genome in self.all_genomes.items():
            # Open file
            with open(os.path.join(genome_folder,genome_name + "_GENOME.tsv"),"w") as f:

                header = ["POSITION","GENE_FAMILY","ORIENTATION","GENE_ID"]
                header = "\t".join(map(str, header)) + "\n"

                f.write(header)

                for chromosome in genome:
                    for index, gene in enumerate(chromosome):
                        #print(chromosome) I have to add the chromosome id
                        line = [index, gene.gene_family, gene.orientation, gene.gene_id]
                        line = "\t".join(map(str,line)) +"\n"
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

    def write_gene_trees(self, gene_tree_folder):

        if not os.path.isdir(gene_tree_folder):
            os.mkdir(gene_tree_folder)

        for gene_family_name, gene_family in self.all_gene_families.items():

            whole_tree, pruned_tree = gene_family.generate_tree()

            with open(os.path.join(gene_tree_folder, gene_family_name + "_wholetree.nwk"),"w") as f:

                f.write(whole_tree)

            if pruned_tree != None:

                with open(os.path.join(gene_tree_folder, gene_family_name + "_prunedtree.nwk"),"w") as f:

                    f.write(pruned_tree)

    def write_events_per_branch(self, events_per_branch_folder):

        if not os.path.isdir(events_per_branch_folder):
            os.mkdir(events_per_branch_folder)

        events_per_branch = dict()

        for gene_family_name, gene_family in self.all_gene_families.items():

            for time, event, nodes in gene_family.events:

                name = nodes.split(";")[0]

                if event == "S" or event == "E" or event == "F":
                    continue

                if name not in events_per_branch:
                    events_per_branch[name] = list()

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

                    if name not in events_per_branch:
                        events_per_branch[name] = list()

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


    def _read_events_file(self, events_file):

        events = list()
        with open(events_file) as f:
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

    def _FillGenome(self, genome_file, all_gene_families):

        genome = Genome()
        genome.species = "Root"
        time = 0

        with open(genome_file) as f:

            for line in f:

                n_genes, shape = line.strip().split("\t")

                # We create a chromosome for each entry

                if shape == "L":
                    chromosome = LinearChromosome()
                    chromosome.shape = "L"
                elif shape == "C":
                    chromosome = CircularChromosome()
                    chromosome.shape = "C"

                for i in range(int(n_genes)):
                    # We fill the chromosomes and we create also the gene families

                    gene, gene_family = self.make_origination(genome.species, time)
                    chromosome.genes.append(gene)
                    all_gene_families[str(self.gene_families_counter)] = gene_family

                genome.chromosomes.append(chromosome)

        return genome

    def run(self):

        d, t, l, i, c, o = (4, 4, 0, 0, 0, 0)

        # First we prepare the first genome

        genome = self._FillGenome(self.parameters["GENOME_FILE"], self.all_gene_families)
        self.active_genomes.add(genome.species)
        self.all_genomes["Root"] = genome

        current_species_tree_event = 0
        current_time = 0.0
        all_species_tree_events = len(self.tree_events)
        # Second, we compute the time to the next event:

        elapsed_time = 0.0

        while current_species_tree_event < all_species_tree_events:

            time_of_next_species_tree_event, event, nodes = self.tree_events[current_species_tree_event]
            time_of_next_species_tree_event = float(time_of_next_species_tree_event)

            # We check that we are not exactly in the same span of time WRITE THIS!!

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

            else:

                current_time += time_to_next_genome_event

                self.evolve_genomes(d, t, l, i, c, o, current_time)


    def verbose_run(self):

        # Only for debugging purposes

        detailed_genomes = list()

        d, t, l, i, c, o = (4, 4, 0, 0, 0, 0)

        # First we prepare the first genome

        genome = self._FillGenome(self.parameters["GENOME_FILE"], self.all_gene_families)
        self.active_genomes.add(genome.species)
        self.all_genomes["Root"] = genome

        detailed_genomes.append((0, genome.species, copy.deepcopy(genome)))

        current_species_tree_event = 0
        current_time = 0.0
        all_species_tree_events = len(self.tree_events)
        # Second, we compute the time to the next event:

        elapsed_time = 0.0

        while current_species_tree_event < all_species_tree_events:

            time_of_next_species_tree_event, event, nodes = self.tree_events[current_species_tree_event]
            time_of_next_species_tree_event = float(time_of_next_species_tree_event)

            # We check that we are not exactly in the same span of time WRITE THIS!!

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

                    detailed_genomes.append((current_time, c1, copy.deepcopy(genome_c1)))
                    detailed_genomes.append((current_time, c2,  copy.deepcopy(genome_c2)))


                elif event == "E":
                    self.make_extinction(nodes, current_time)
                    self.active_genomes.discard(nodes)

                    detailed_genomes.append((current_time, nodes, copy.deepcopy(self.all_genomes[nodes])))

                elif event == "F":
                    self.make_end(current_time)

            else:

                current_time += time_to_next_genome_event
                genome_event, affected_linage = self.evolve_genomes(d, t, l, i, c, o, current_time)
                if genome_event == "T":

                    donor, recipient = affected_linage.split("->")

                    detailed_genomes.append((current_time, donor + ";" + "LT",
                                             copy.deepcopy(self.all_genomes[donor])))

                    detailed_genomes.append((current_time, recipient + ";" + "AT",
                                             copy.deepcopy(self.all_genomes[recipient])))

                else:
                    detailed_genomes.append((current_time, affected_linage + ";" + genome_event,
                                             copy.deepcopy(self.all_genomes[affected_linage])))


        for time, lineage, genome in detailed_genomes:
            print(time, lineage, genome)


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

        total_probability_of_event = duplication + transfer + loss + inversion + translocation + origination

        d_e = float(self.parameters["DUPLICATION_E"])
        t_e = float(self.parameters["TRANSFER_E"])
        l_e = float(self.parameters["LOSS_E"])
        i_e = float(self.parameters["INVERSION_E"])
        c_e = float(self.parameters["TRANSLOCATION_E"])

        if numpy.random.uniform(0, 1) <= total_probability_of_event:  # An event takes place

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
                return "O", lineage
                pass
                #gene, gene_family = make_origination(1,lineage)

    def get_time_to_next_event(self, n, d, t, l ,i , c, o):

        total = 0.0
        for j in range(n):
            total += sum((d,t,l,i,c,o))
        time = numpy.random.exponential(1/total)
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
        gene.genome = species_tree_node

        gene_family = GeneFamily(gene_family_id, time)
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

            for gene in chromosome:

                new_id1 = self.return_new_identifiers_for_segment([gene])
                new_id2 = self.return_new_identifiers_for_segment([gene])

                new_gene1 = copy_segment([Gene()], new_id1)[0]
                new_gene2 = copy_segment([Gene()], new_id2)[0]

                new_gene1.genome = c1
                new_gene2.genome = c2

                new_gene1.orientation = gene.orientation
                new_gene2.orientation = gene.orientation

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

        copied_segment1 = copy_segment(segment, new_identifiers1)
        copied_segment2 = copy_segment(segment, new_identifiers2)

        # We insert the two new segments after the last position of the old segment

        chromosome.insert_segment(affected_genes[-1], copied_segment1 + copied_segment2)

        # And we remove the old segment

        chromosome.remove_segment(segment)

        # We have to register in the affected gene families that there has been a duplication

        for i, gene in enumerate(segment):

            nodes = [gene.genome,
                     gene.gene_id,
                     copied_segment1[i].genome,
                     copied_segment1[i].gene_id,
                     copied_segment2[i].genome,
                     copied_segment2[i].gene_id]

            gene.active = False

            self.all_gene_families[gene.gene_family].register_event(time, "D", ";".join(map(str, nodes)))

    def make_transfer(self, p, donor, recipient, time):

        chromosome1 = self.all_genomes[donor].select_random_chromosome()
        chromosome2 = self.all_genomes[recipient].select_random_chromosome()

        affected_genes = chromosome1.obtain_affected_genes(p)
        segment = chromosome1.obtain_segment(affected_genes)

        new_identifiers1 = self.return_new_identifiers_for_segment(segment)
        new_identifiers2 = self.return_new_identifiers_for_segment(segment)

        # Now we create two segments

        copied_segment1 = copy_segment(segment, new_identifiers1)
        copied_segment2 = copy_segment(segment, new_identifiers2)

        # We insert the first segment (leaving transfer) in the same position than the previous segment
        # We do this just to change the identifiers of the numbers

        chromosome1.insert_segment(affected_genes[0], copied_segment1)

        # And we remove the old segment

        chromosome1.remove_segment(segment)

        # Now we insert the transfer segment in the recipient genome

        position = chromosome2.select_random_position()

        chromosome2.insert_segment(position, copied_segment2)

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

            self.all_gene_families[gene.gene_family].register_event(time, "T", ";".join(map(str,nodes)))

    def make_loss(self, p, lineage, time):

        chromosome = self.all_genomes[lineage].select_random_chromosome()
        affected_genes = chromosome.obtain_affected_genes(p)
        segment = chromosome.obtain_segment(affected_genes)
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


class SpeciesTreeGenerator():

    def __init__(self):

        self.whole_species_tree = ete3.Tree()

        self.lineages_counter = 0
        self.events = list()
        self.events.append(("0", "S", "Root;n1;n2"))

        self._start_tree(self.whole_species_tree)
        self.parameters = dict()

        self.parameters["SPECIATION"] = 1.2
        self.parameters["EXTINCTION"] = 0.2
        self.parameters["TOTAL_TIME"] = 6
        self.parameters["STOPPING_RULE"] = 1
        self.parameters["TOTAL_LINEAGES"] = 5

    def _start_tree(self, tree):

        tree.name = "Root"
        tree.add_feature("is_alive", True)

        c1 = tree.add_child(dist=0)
        self.lineages_counter += 1
        c1.name = "n" + str(self.lineages_counter)
        c1.add_feature("is_alive", True)  # One indicates that the species is alive at present time

        c2 = tree.add_child(dist=0)
        self.lineages_counter += 1
        c2.name = "n" + str(self.lineages_counter)
        c2.add_feature("is_alive", True)  # One indicates that the species is alive at present time

        tree.is_alive = False

    def generate_precomputed_tree(self, events):

        new_tree = ete3.Tree()
        self._start_tree(new_tree)

    def generate_new_tree(self):

        speciation = float(self.parameters["SPECIATION"])
        extinction = float(self.parameters["EXTINCTION"])
        stopping_rule = int(self.parameters["STOPPING_RULE"])
        total_time = float(self.parameters["TOTAL_TIME"])
        total_lineages = float(self.parameters["TOTAL_LINEAGES"])

        time = 0
        success = False

        while True:

            print(time)

            lineages_alive = [x for x in self.whole_species_tree.get_leaves() if x.is_alive == True]

            time_to_next_event = self.get_time_to_next_event(len(lineages_alive), speciation, extinction)

            if stopping_rule == 0 and time + time_to_next_event >= total_time:

                self.increase_distances(total_time - time, lineages_alive)
                self.events.append((total_time, "F", "None"))

                break
            elif stopping_rule == 1 and len(lineages_alive) >= total_lineages:

                self.increase_distances(time_to_next_event, lineages_alive)
                self.events.append((time+time_to_next_event, "F", "None"))

                break

            elif len(lineages_alive) == 0:

                self.increase_distances(time_to_next_event, lineages_alive)
                print("All dead")
                break

            else:
                # In this case we do the normal the computation

                time += time_to_next_event

                self.increase_distances(time_to_next_event, lineages_alive)

                event = self.choose_event(speciation, extinction)
                lineage = random.choice(lineages_alive)

                if event == "S":
                    self._get_speciated(lineage, time)

                elif event == "E":
                    self._get_extinct(lineage, time)

    def increase_distances(self, time, lineages_alive):

        for lineage in lineages_alive:
            lineage.dist += time

    def get_time_to_next_event(self, n, s, e):

        total = 0.0
        for i in range(n):
            total += s
            total += e
        time = numpy.random.exponential(1/total)
        return time

    def _get_speciated(self, lineage, time):

        c1 = lineage.add_child(dist=0)
        self.lineages_counter += 1
        c1.name = "n" + str(self.lineages_counter)
        c1.add_feature("is_alive", True)

        c2 = lineage.add_child(dist=0)
        self.lineages_counter += 1
        c2.name = "n" + str(self.lineages_counter)
        c2.add_feature("is_alive", True)

        self.events.append((time,"S", ";".join((lineage.name,c1.name,c2.name))))  # Store the event
        lineage.is_alive = False

    def _get_extinct(self, lineage, time):

        lineage.is_alive = False

        self.events.append((time, "E", lineage.name))  # Store the event

    def get_extant_tree(self):

        leaves_alive = [x.name for x in self.whole_species_tree.get_leaves() if x.is_alive == True]

        if len(leaves_alive) < 3:
            return 0

        self.extant_species_tree = ete3.Tree(self.whole_species_tree.write(format=1), format=1)
        self.extant_species_tree.prune(leaves_alive, preserve_branch_length=True)

    def choose_event(self, speciation, extinction):

        if numpy.random.uniform(0, 1) <= (speciation / (speciation + extinction)):
            return "S"
        else:
            return "E"

    def write_events_file(self, events_file):

        with open(events_file, "w") as f:
            for item in self.events:
                line = "\t".join(map(str,item)) + "\n"
                f.write(line)

    def write_tree(self, tree_file):

        with open(tree_file, "w") as f:
            f.write(self.whole_species_tree.write(format=1))


#GenomeEvolver()

#stg = SpeciesTreeGenerator()
#stg.generate_new_tree()
#print(stg.whole_species_tree.write(format=1))

myevents = "/Users/adriandavin/Desktop/Bioinformatics/SimuLyon/Cedric/Gillespie/Events.tsv"
#mytree = "/Users/adriandavin/Desktop/Bioinformatics/SimuLyon/Cedric/Gillespie/WholeTree"
tree_file = "/Users/adriandavin/Desktop/Bioinformatics/SimuLyon/Cedric/Gillespie/CyanoTree"
#stg.write_events_file(myevents)
#stg.write_tree(mytree)

#species_tree_events = list()
#species_counter = 0

#gss = GenomeSimulator(myevents)
#gss.verbose_run()


p = read_parameters("/Users/adriandavin/PycharmProjects/SimuLYON/Parameters/SpeciesTreeParameters.tsv")
p = prepare_parameters(p)



#gss.write_genomes("/Users/adriandavin/Desktop/Bioinformatics/SimuLyon/Cedric/Gillespie/All/Genomes/")
#gss.write_gene_family_events("/Users/adriandavin/Desktop/Bioinformatics/SimuLyon/Cedric/Gillespie/All/Gene_families/")
#gss.write_gene_trees("/Users/adriandavin/Desktop/Bioinformatics/SimuLyon/Cedric/Gillespie/All/Gene_trees/")
#gss.write_events_per_branch("/Users/adriandavin/Desktop/Bioinformatics/SimuLyon/Cedric/Gillespie/All/Events_per_branch")

#for item in gss.all_gene_families["3"].events:
#    print(item)

#print(gss.all_gene_families["3"].generate_tree())
#generate_events(tree_file)

# Genome evolution in continuous time

## When I come back:

# 1. Fix function make_duplication so that it takes lineage instead of chromosome 10 min # Done
# 2. Fix inversions and translocations 1 h # I think it is done
# 3. Function to prune the trees 10 m - Done
# 4. Write nice input and output 2 h - Almost done
# 5. Write replacement transfers 1 h



# I shouldn't continue with the next points till having al the other things solved

# 6. Transfers proportional to distance 2h
# 7. Variable rates for the species tree 8h
# 8. Variable rates for genomes 16h
# 9. Documenting software 8h
# 10. Adding sequences evolution 16h
# 11. Adding linear chromosomes 8h
# 12. Adding intergenic distances 8h


# Write output
# Replacement transfers
# Linear chromosomes
# Fissions and fusions of chromosomes
# Transfers proportional to distance
# Variable rates of speciation and extinction
# Variable rates of evolution in gene families
# Adding simulation of sequences
# Testing different pieces of software
# Simulation of independent gene families (?)

