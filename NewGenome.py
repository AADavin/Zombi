import AuxiliarFunctions as af
import numpy
import copy
import sys
import ete3
import random

input = [(5, "C"), (10, "C")]
all_gene_families = dict()
species_tree_events = list()
species_counter = 0

def copy_segment(segment, new_identifiers):

    new_segment = list()

    for i,gene in enumerate(segment):
        new_gene = copy.deepcopy(gene)
        new_gene.gene_id = new_identifiers[i]
        new_segment.append(new_gene)

    return new_segment

def invert_segment(segment):

    for gene in segment:
        gene.change_sense()


def return_new_identifiers_for_segment(segment, all_gene_families):

    new_identifiers = list()

    for gene in segment:
        gf = gene.gene_family
        new_id = all_gene_families[gf].obtain_new_gene_id()
        new_identifiers.append(new_id)

    return new_identifiers


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
        pass

    def obtain_new_gene_id(self):
        self.gene_ids_counter += 1
        return self.gene_ids_counter

    def __str__(self):
        return "GeneFamily_" + str(self.identifier)

        #return ";".join(map(str, self.genes))

    def __len__(self):

        return len([x for x in self.genes if gene.active == True])

class Gene():

    def __init__(self):

        self.active = True
        self.orientation = ""
        self.gene_family = ""
        self.gene_id = ""
        self.sequence = ""
        self.genome = ""

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
        myname = "_".join(map(str,(self.genome, self.orientation, self.gene_family, self.gene_id)))
        return myname

class Chromosome():

    def __init__(self, shape):

        self.genes = list()
        self.intergene_distances = list()
        self.shape = shape
        self.species = ""

    def select_random_position(self):

        return numpy.random.randint(len(self.genes))

    def select_random_length(self, p):

        return numpy.random.geometric(p)

    def __len__(self):

        return len(self.genes)

    def __str__(self):

        return ";".join([str(gene) for gene in self.genes])

    def __iter__(self):

        for x in self.genes:
            yield x


class CircularChromosome(Chromosome):

    def obtain_segment(self, affected_genes):

        segment = [self.genes[x] for x in affected_genes]

        return segment

    def remove_segment(self, segment):

        for gene in segment:
            self.genes.remove(gene)

    def insert_segment(self, position, segment):

        for i, x in enumerate(segment):

            self.genes.insert(position + i, x)

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

    def get_homologous_position(self):
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

        chromosome = numpy.random.choice(self.chromosomes, 1, p=af.normalize([len(x) for x in self.chromosomes]))
        return chromosome

    def __str__(self):

        return ";".join([str(x) for x in self.chromosomes])

    def __iter__(self):
        for chromosome in self.chromosomes:
            yield chromosome

def make_speciation(sp, c1, c2):

    genome1 = Genome()
    genome2 = Genome()

    genome1.species = c1.name
    genome2.species = c2.name

    for chromosome in sp.genome:

        shape = chromosome.shape

        ch1 = Chromosome(shape)
        ch2 = Chromosome(shape)

        genome1.chromosomes.append(ch1)
        genome2.chromosomes.append(ch2)

        for gene in chromosome:

            new_id1 = return_new_identifiers_for_segment([gene], all_gene_families)
            new_id2 = return_new_identifiers_for_segment([gene], all_gene_families)

            new_gene1 = copy_segment([Gene()], new_id1)[0]
            new_gene2 = copy_segment([Gene()], new_id2)[0]

            new_gene1.genome = c1.name
            new_gene2.genome = c2.name

            new_gene1.orientation = gene.orientation
            new_gene2.orientation = gene.orientation

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

            nodes = [sp.name,
                     gene.gene_id,
                     c1.name,
                     new_gene1.gene_id,
                     c2.name,
                     new_gene2.gene_id
                     ]

            all_gene_families[gene.gene_family].register_event(0, "S", ";".join(nodes))

    c1.add_feature("genome", genome1)
    c2.add_feature("genome", genome2)

def make_extinction(sp):

    # We have to inactivate all the genes

    genome = sp.genome

    for chromosome in genome:
        for gene in chromosome:
            gene.active = False

def make_loss(p, chromosome):

    affected_genes = chromosome.obtain_affected_genes(p)
    segment = chromosome.obtain_segment(affected_genes)
    chromosome.remove_segment(segment)

    # We have to register in the affected gene families that there has been as loss
    # All genes affected must be returned

    for gene in segment:
        gene.active = False
        all_gene_families[gene.gene_family].register_event(0, "L", gene.gene_id)


def make_duplication(p, chromosome):

    affected_genes = chromosome.obtain_affected_genes(p)
    segment = chromosome.obtain_segment(affected_genes)

    new_identifiers1 = return_new_identifiers_for_segment(segment, all_gene_families)
    new_identifiers2 = return_new_identifiers_for_segment(segment, all_gene_families)

    # Now we create two segments

    copied_segment1 = copy_segment(segment, new_identifiers1)
    copied_segment2 = copy_segment(segment, new_identifiers2)

    # We insert the two new segments after the last position of the old segment

    chromosome.insert_segment(affected_genes[-1], copied_segment1 + copied_segment2)

    # And we remove the old segment

    chromosome.remove_segment(segment)

    # We have to register in the affected gene families that there has been a duplication

    for i,gene in enumerate(segment):
        gene.active = False
        all_gene_families[gene.gene_family].register_event(0, "D", ";".join(
            (gene.gene_id, new_identifiers1[i], new_identifiers2[i])))

def make_inversion(p, chromosome):

    affected_genes = chromosome.obtain_affected_genes(p)
    segment = chromosome.obtain_segment(affected_genes)

    new_identifiers = return_new_identifiers_for_segment(segment, all_gene_families)
    new_segment = copy_segment(segment, new_identifiers)
    invert_segment(new_segment)

    chromosome.insert_segment(affected_genes[0], new_segment[::-1])
    chromosome.remove_segment(segment)

    for i,gene in enumerate(segment):
        gene.active = False
        all_gene_families[gene.gene_family].register_event(0, "I", str(gene.gene_id) + ";" + str(new_identifiers[i]))


def make_translocation(p, chromosome):

    affected_genes = chromosome.obtain_affected_genes(p)
    segment = chromosome.obtain_segment(affected_genes)

    new_identifiers = return_new_identifiers_for_segment(segment, all_gene_families)
    new_segment = copy_segment(segment, new_identifiers)

    # For now all translocations keep the orientation of the DNA

    position = chromosome.select_random_position()

    # Watch out! Right now translocations can occur within the sequence (resulting in no changes at all)

    chromosome.insert_segment(position, new_segment)
    chromosome.remove_segment(segment)

    for i, gene in enumerate(segment):
        gene.active = False
        all_gene_families[gene.gene_family].register_event(0, "C",
                                                           str(gene.gene_id) + ";" + str(new_identifiers[i]))

def make_origination(gene_family_id, species_tree_node):

    gene = Gene()
    gene.determine_orientation()
    gene.gene_family = gene_family_id
    gene.genome = species_tree_node

    gene_family_id += 1
    gene_family = GeneFamily(gene_family_id, 0)
    gene_family.genes.append(gene)

    gene.gene_id = gene_family.obtain_new_gene_id()

    all_gene_families[gene_family_id] = gene_family

    return gene, gene_family

def make_transfer(p, chromosome1, chromosome2):

    affected_genes = chromosome1.obtain_affected_genes(p)
    segment = chromosome1.obtain_segment(affected_genes)

    new_identifiers1 = return_new_identifiers_for_segment(segment, all_gene_families)
    new_identifiers2 = return_new_identifiers_for_segment(segment, all_gene_families)

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

        nodes = [chromosome1.species,
                 segment[i].gene_id,
                 new_identifiers1[i],
                 chromosome2.species,
                 new_identifiers2[i],
                 ]

        all_gene_families[gene.gene_family].register_event(0, "T", ";".join(nodes))


# First we create genes and their gene families
# A gene belongs both to a gene family and to a chromosome
# So first you create chromosomes and then they create gene families

def FillGenome(gene_family, genome, input):
    pass



def GenomeEvolver():

    tree = ete3.Tree()
    genome = Genome()
    tree.name = "Root"

    sp = tree.get_tree_root()

    gene_family_id = 0

    for n_genes, shape in input:

        # We create a chromosome for each entry

        if shape == "L":
            chromosome = LinearChromosome(shape)
        elif shape == "C":
            chromosome = CircularChromosome(shape)

        chromosome.species = tree.name

        for i in range(n_genes):

            # We fill the chromosomes and we create also the gene families

            gene, gene_family = make_origination(gene_family_id, chromosome.species)
            chromosome.genes.append(gene)
            all_gene_families[gene_family_id] = gene_family

            gene_family_id += 1

        genome.chromosomes.append(chromosome)

class GenomeSimulator():

    def __init__(self, events):

        self.parameters = dict()

        self.parameters["DUPLICATION_E"] = 0.1
        self.parameters["LOSS_E"] = 0.1
        self.parameters["TRANSFER_E"] = 0.1
        self.parameters["INVERSION_E"] = 0.1
        self.parameters["TRANSLOCATION_E"] = 0.1

        self.species_tree = ete3.Tree()
        self.tree_events = events
        self._start_tree()

    def _start_tree(self):

        gnm = Genome()
        root = self.species_tree.get_tree_root()
        root.name = "Root"
        root.add_feature("genome",gnm)
        root.add_feature("is_alive", True)

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

    def evolve_genomes(self, lineages_active, duplication, transfer, loss, inversion, translocation, origination):

        total_probability_of_event = duplication + transfer + loss + inversion + translocation + origination

        d_e = float(self.parameters["DUPLICATION_E"])
        t_e = float(self.parameters["TRANSFER_E"])
        l_e = float(self.parameters["LOSS_E"])
        i_e = float(self.parameters["INVERSION_E"])
        c_e = float(self.parameters["TRANSLOCATION_E"])

        lineage = random.choice(lineages_active)
        chromosome = lineage.genome.select_random_chromosome()

        if numpy.random.uniform(0, 1) <= total_probability_of_event:  # An event takes place

            event = self.choose_event(duplication, transfer, loss, inversion, translocation, origination)

            if event == "D":
                make_duplication(d_e, chromosome)

            elif event == "T":

                # We choose a recipient

                possible_recipients = [x for x in lineages_active if x != lineage]

                if len(possible_recipients) > 0:

                    recipient = random.choice(possible_recipients)
                    chromosome_r = recipient.genome.select_random_chromosome()
                    make_transfer(t_e, chromosome, chromosome_r)

            elif event == "L":
                make_loss(l_e, chromosome)

            elif event == "I":
                make_inversion(i_e, chromosome)

            elif event == "C":
                make_translocation(c_e, chromosome)

            elif event == "O":
                pass
                #gene, gene_family = make_origination(1,lineage)

    def run(self):

        d, t, l, i, c, o = (0.1, 0.1, 0.1, 0.1, 0.1, 0.1)

        event_counter = 0
        time = 0

        time_to_next_st_event, next_st_event, nodes_in_event = self.tree_events[event_counter]

        while True:

            lineages_alive = [x for x in self.species_tree.get_leaves() if x.is_alive == True]
            time_to_next_event = self.get_time_to_next_event(len(lineages_alive), d, t, l, i, c, o)

            if time + float(time_to_next_event) >= float(time_to_next_st_event):

                time_to_next_st_event, next_st_event, nodes_in_event = self.tree_events[event_counter]
                event_counter += 1

                if next_st_event == "S":
                    self._get_speciated(nodes_in_event)

                elif next_st_event == "E":
                    self._get_extinct(nodes_in_event)

                elif next_st_event == "F":
                    pass

            else:

                time += float(time_to_next_event)
                event = self.choose_event(d,t,l,i,c,o)
                self.evolve_genomes(lineages_alive, d,t,l,i,c,o)
                self.increase_distances(float(time_to_next_event), lineages_alive)

    def get_time_to_next_event(self, n, d, t, l ,i , c, o):

        total = 0.0
        for j in range(n):
            total += sum((d,t,l,i,c,o))
        time = numpy.random.exponential(1/total)
        return time

    def increase_distances(self, time_to_next_event, active_lineages):

        for node in active_lineages:
            node.dist += time_to_next_event

    def _get_speciated(self, nodes_in_event):

        # First we take care of the species tree

        sp, c1, c2 = nodes_in_event.split(";")

        st_sp = self.species_tree&sp

        sc1 = st_sp.add_child(dist=0)
        sc1.name = c1
        sc1.add_feature("is_alive", True)

        sc2 = st_sp.add_child(dist=0)
        sc2.name = c2
        sc2.add_feature("is_alive", True)

        st_sp.is_alive = False

        # Then, we take care of the genomes affected by the speciation event

        make_speciation(st_sp, sc1, sc2)

    def _get_extinct(self, nodes_in_event):

        sp_st = self.species_tree&nodes_in_event
        sp_st.is_alive = False
        make_extinction(sp_st)

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
        self.parameters["TOTAL_LINEAGES"] = 100

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

    def write_log(self):

        for item in self.events:
            print(item)


#GenomeEvolver()

stg = SpeciesTreeGenerator()
stg.generate_new_tree()
print(stg.whole_species_tree.write(format=1))
events = stg.events

gss = GenomeSimulator(events)
gss.run()




#chromosome.remove_segment(segment)

#for gene in segment:
#    gene.active = False
#    all_gene_families[gene.gene_family].register_loss(0, gene)

#print(chromosome)


if __name__ == "__main__":

    args = sys.argv[1:]
    if len(args) != 3:
        print("Incorrect usage. Please read the manual. The usual way to run this script is:")
        print("python simuLyon.py T Parameters_file.tsv /Output_folder")
