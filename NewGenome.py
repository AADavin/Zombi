import numpy
import copy

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

    def __init__(self, identifier):

        self.identifier = identifier
        self.genes = list()
        self.events = list()
        self.event_counter = 0  # Each time that the family is modified in any form, we have to update the event counter
        self.gene_ids_counter = 0

    def register_loss(self, time, genes):

        self.events.append((time, "L", genes))

    def register_duplication(self, time, genes):

        self.events.append((time, "D", genes))

    def register_inversion(self, time, genes):

        self.events.append((time, "I", genes))

    def register_translocation(self, time, genes):

        self.events.append((time, "C", genes))

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

    def select_random_position(self):

        return numpy.random.randint(len(self.genes))

    def select_random_length(self, p):

        return numpy.random.geometric(p)

    def __len__(self):

        return len(self.genes)

    def __str__(self):

        return ";".join([str(gene) for gene in self.genes])

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

        # I have to ponderate by the length of each chromosome

        return numpy.random.randint(len(self.genes))






input = [(10, "C"), (10, "C")]

# First we create genes and their gene families
# A gene belongs both to a gene family and to a chromosome
# So first you create chromosomes and then they create gene families

gene_families_counter = 0
all_gene_families = dict()
genome = Genome()
gf_name_id = 0

for n_genes, shape in input:

    # We create a chromosome for each entry

    if shape == "L":
        chromosome = LinearChromosome(shape)
    elif shape == "C":
        chromosome = CircularChromosome(shape)

    for i in range(n_genes):
        # We fill the chromosomes and we create also the gene families

        gf_name_id += 1
        gene_family = GeneFamily(gf_name_id)

        gene = Gene()
        gene.determine_orientation()
        gene.gene_family = gf_name_id

        gene.gene_id = gene_family.obtain_new_gene_id()
        gene.genome = "Root"

        gene_family.genes.append(gene)
        chromosome.genes.append(gene)
        all_gene_families[gf_name_id] = gene_family



def make_loss(p):

    affected_genes = chromosome.obtain_affected_genes(p)
    segment = chromosome.obtain_segment(affected_genes)
    chromosome.remove_segment(segment)

    # We have to register in the affected gene families that there has been as loss
    # All genes affected must be returned

    for gene in segment:
        gene.active = False
        all_gene_families[gene.gene_family].register_loss(0, gene.gene_id)


def make_duplication(p):

    affected_genes = chromosome.obtain_affected_genes(p)
    segment = chromosome.obtain_segment(affected_genes)

    new_identifiers1 = return_new_identifiers_for_segment(segment, all_gene_families)
    new_identifiers2 = return_new_identifiers_for_segment(segment, all_gene_families)

    # Now we create two segments

    copied_segment1 = copy_segment(segment, new_identifiers1)
    copied_segment2 = copy_segment(segment, new_identifiers2)

    # We remove the old segment

    # We insert the two new segments after the last position of the old segment

    chromosome.insert_segment(affected_genes[-1], copied_segment1 + copied_segment2)

    # And we remove the old segment

    chromosome.remove_segment(segment)

    # We have to register in the affected gene families that there has been a duplication

    for i,gene in enumerate(segment):
        gene.active = False
        all_gene_families[gene.gene_family].register_duplication(0,(gene.gene_id, new_identifiers1[i], new_identifiers2[i]))

def make_inversion(p):

    affected_genes = chromosome.obtain_affected_genes(p)
    segment = chromosome.obtain_segment(affected_genes)

    new_identifiers = return_new_identifiers_for_segment(segment, all_gene_families)
    new_segment = copy_segment(segment, new_identifiers)
    invert_segment(new_segment)

    chromosome.insert_segment(affected_genes[0], new_segment[::-1])
    chromosome.remove_segment(segment)

    for i,gene in enumerate(segment):
        gene.active = False
        all_gene_families[gene.gene_family].register_inversion(0, str(gene.gene_id) + ";" + str(new_identifiers[i]))


def make_translocation():
    
    affected_genes = chromosome.obtain_affected_genes(p)
    segment = chromosome.obtain_segment(affected_genes)

    new_identifiers = return_new_identifiers_for_segment(segment, all_gene_families)
    new_segment = copy_segment(segment, new_identifiers)
    invert_segment(new_segment)

    chromosome.insert_segment(affected_genes[0], new_segment[::-1])
    chromosome.remove_segment(segment)

    for i, gene in enumerate(segment):
        gene.active = False
        all_gene_families[gene.gene_family].register_inversion(0, str(gene.gene_id) + ";" + str(new_identifiers[i]))


def make_origination():
    pass


def make_transfer():
    pass



print(chromosome)
make_inversion(0.4)
print(chromosome)










#chromosome.remove_segment(segment)

#for gene in segment:
#    gene.active = False
#    all_gene_families[gene.gene_family].register_loss(0, gene)

#print(chromosome)



'''    
# THIS IS WORKING

for x in all_gene_families:
    print(str(all_gene_families[x]))

print("NOW WE CHANGE")

chromosome.genes[3].gene_id = "4"

for x in all_gene_families:
    print(str(all_gene_families[x]))


print(chromosome)
















#genome.start_genome(input)

'''
