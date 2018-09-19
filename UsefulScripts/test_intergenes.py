import numpy
import GenomeClasses as GC
import AuxiliarFunctions as af


def fill_genome(intergenic_sequences = False):

        genome = GC.Genome()
        genome.species = "Root"
        time = 0

        initial_genome_size = [3]
        shape = "C"

        for n_genes in initial_genome_size:

            if shape == "L":
                chromosome = GC.LinearChromosome()
                chromosome.shape = "L"
            elif shape == "C":
                chromosome = GC.CircularChromosome()
                chromosome.shape = "C"

            if intergenic_sequences == True:
                chromosome.has_intergenes = True

            for i in range(int(n_genes)):

                # We fill the chromosomes and we create also the gene families

                gene = GC.Gene()
                gene.length = 3
                gene.gene_family = i
                gene.gene_id = 1
                chromosome.genes.append(gene)

                if intergenic_sequences == True:
                    intergenic_sequence = GC.Intergene()
                    intergenic_sequence.length = 3
                    chromosome.intergenes.append(intergenic_sequence)

            genome.chromosomes.append(chromosome)

        return chromosome

g = fill_genome(intergenic_sequences=True)
print(g)
g.obtain_flankings()
g.obtain_locations()

for i in g.map_of_locations:
    print(i)

def test_inversion(g, c1, c2, d):
    r1, r2, r3, r4 = g.return_affected_region(c1, c2, d)
    g.invert_segment(r1)
    left_intergene = g.intergenes[r2[0]]
    right_intergene = g.intergenes[r2[-1]]
    left_intergene.length -= r3[1]
    left_intergene.length += r4[0]
    right_intergene.length -= r4[0]
    right_intergene.length += r3[1]


def select_lengthy(p, extension_multiplier):

    counter = 0
    total_genome_length = g.map_of_locations[-1][1]
    success = False

    while counter <= 100 and success == False:

        counter += 1
        sc1 = g.select_random_coordinate_in_intergenic_regions()
        tc1 = g.return_total_coordinate_from_specific_coordinate(sc1)
        d = numpy.random.choice(("left", "right"), p=[0.5, 0.5])
        extension = numpy.random.geometric(p) * extension_multiplier

        if d == "right":
            if tc1 + extension >= total_genome_length:
                tc2 = extension - (total_genome_length - tc1)
                if tc2 < tc1:
                    success = True
                else:
                    # The event covers the whole genome
                    pass
            else:
                tc2 = tc1 + extension
                success = True

        elif d == "left":

            if tc1 - extension <= 0:
                tc2 = total_genome_length - extension - (0 - tc1)
                if tc1 < tc2:
                    success = True
                else:
                    # The event covers the whole genome
                    pass
            else:
                tc2 = tc1 - extension
                success = True

        if success == True and tc2 >= 0 and tc2 <= total_genome_length:

            sc2 = g.return_specific_coordinate_from_total_coordinate(tc2)

            if sc2 == None:
                success = False
            else:
                return sc1, sc2, d

    return None







p = 0.5
extension_multiplier = 6

sc1, sc2, d = select_lengthy(p, extension_multiplier)
print(sc1, sc2, d)
print(g)
test_inversion(g, sc1, sc2, d)
print(g)



#test_inversion(g, 3, 11, "left")

#print(g)








def test2():
    for j in range(1000):
        c1 = g.select_random_coordinate_in_intergenic_regions()
        l1 = g.return_location_by_coordinate(c1, within_intergene=True)
        gene = GC.Gene()
        gene.length = 3
        gene.gene_family = "A"+str(j)
        g.insert_gene_within_intergene(c1, l1, gene)
        g.obtain_flankings()
        g.obtain_locations()
    a=0
    for i in g.intergenes:
        a += i.length

    print(a)



def test1():

    c1 = g.select_random_coordinate_in_intergenic_regions()
    c2 = g.select_random_coordinate_in_intergenic_regions()
    #c1 = 2
    #c2 = 6


    d = numpy.random.choice(("left","right"),p=[0.5,0.5])
    print(d)

    r = g.return_affected_segment(c1, c2, d)
    if r != None:
        r1, r2, r3, r4 = r

        print(r1)
        print(r2)
        print(r3)
        print(r4)
    else:
        print("Ignoring event")

