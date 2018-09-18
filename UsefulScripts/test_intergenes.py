import numpy
import GenomeClasses as GC
import AuxiliarFunctions as af


def fill_genome(intergenic_sequences = False):

        genome = GC.Genome()
        genome.species = "Root"
        time = 0

        initial_genome_size = [2]
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
                gene.length = 4
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

g.obtain_flankings()
g.obtain_locations()
for i in g.map_of_locations:
    print(i)





c1 = g.select_random_coordinate_in_intergenic_regions()
c2 = g.select_random_coordinate_in_intergenic_regions()
l1 = g.return_location_by_coordinate(c1, within_intergene=True)
l2 = g.return_location_by_coordinate(c2, within_intergene=True)
print(c1)
print(l1)
print("....")
print(c2)
print(l2)
print("....")

d = numpy.random.choice(("left","right"),p=[0.5,0.5])
print(d)
print(g.affected_segment(l1, l2, d))





