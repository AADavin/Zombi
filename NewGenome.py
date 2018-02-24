import AuxiliarFunctions as af
import numpy
import copy

import ete3

from ete3 import Tree

input = [(10, "C"), (10, "C")]
all_gene_families = dict()
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

        if numpy.random.randint(1):
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

        chromosome = numpy.random.choice([self.chromosomes], 1, p=af.normalize([len(x) for x in self.chromosomes]))

        return chromosome


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

        genome1.chromosomes.append(Chromosome(shape))
        genome2.chromosomes.append(Chromosome(shape))

        for gene in chromosome:

            new_id1 = return_new_identifiers_for_segment(gene, all_gene_families)
            new_id2 = return_new_identifiers_for_segment(gene, all_gene_families)

            new_gene1 = copy_segment(Gene(), new_id1)
            new_gene2 = copy_segment(Gene(), new_id2)

            genome1.chromosomes.genes.append(new_gene1)
            genome2.chromosomes.genes.append(new_gene2)

            gene.active = False

    return genome1, genome2




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
    gene_family_id += 1
    gene.gene_family = gene_family_id

    gene_family = GeneFamily(gene_family_id, 0)

    gene.determine_orientation()
    gene.gene_family = gene_family_id
    gene.gene_id = gene_family.obtain_new_gene_id()
    gene.genome = species_tree_node
    gene_family.genes.append(gene)
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

def GenomeEvolver():

    tree = ete3.Tree()
    genome = Genome()
    tree.name = "Root"
    tree.add_feature("Genome", genome)
    sp = tree.get_tree_root()

    gene_family_id = 0

    for n_genes, shape in input:

        # We create a chromosome for each entry

        if shape == "L":
            chromosome = LinearChromosome(shape)
        elif shape == "C":
            chromosome = CircularChromosome(shape)

        chromosome.species = "Root"

        for i in range(n_genes):

            # We fill the chromosomes and we create also the gene families

            gene_family_id += 1
            gene, gene_family = make_origination(gene_family_id, chromosome.species)
            chromosome.genes.append(gene)
            all_gene_families[gene_family_id] = gene_family


    c1 = tree.add_child(dist=0)

    c1.name = "n1"
    c2 = tree.add_child(dist=0)
    c2.name = "n2"

    g1,g2 = make_speciation(sp,c1,c2)


## Testing






GenomeEvolver()



#chromosome.remove_segment(segment)

#for gene in segment:
#    gene.active = False
#    all_gene_families[gene.gene_family].register_loss(0, gene)

#print(chromosome)


