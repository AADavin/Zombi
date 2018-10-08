import copy
import GenomeClasses as GC
import AuxiliarFunctions as af



def return_new_identifiers_for_segment(segment):
    iden = 0
    new_identifiers = list()

    for gene in segment:
        iden += 1
        gf = gene.gene_family
        new_id = iden
        new_identifiers.append(new_id)

    return new_identifiers

def fill_genome(intergenic_sequences=False):

    genome = GC.Genome()
    genome.species = "Root"

    initial_genome_size = [4]
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
                intergenic_sequence.id = i
                chromosome.intergenes.append(intergenic_sequence)
        genome.chromosomes.append(chromosome)

    return chromosome

def make_transfer_intergenic(c1, c2, d, chr_donor, chr_recipient):

        r = chr_donor.return_affected_region(c1, c2, d)

        if r == None:
            return None

        else:
            r1, r2, r3, r4 = r
            segment = chr_donor.obtain_segment(r1)

        new_identifiers1 = return_new_identifiers_for_segment(segment)
        new_identifiers2 = return_new_identifiers_for_segment(segment)

        # Now we create two segments

        copied_segment1 = af.copy_segment(segment, new_identifiers1)
        copied_segment2 = af.copy_segment(segment, new_identifiers2)

        new_intergene_segment = [copy.deepcopy(chr_donor.intergenes[x]) for x in r2[1:]]


        # We insert the first segment (leaving transfer) in the same position than the previous segment
        # We do this just to change the identifiers of the numbers

        # We insert in the same place

        position = r1[-1] + 1

        for i, gene in enumerate(copied_segment1):
            chr_donor.genes.insert(position + i, gene)

        # We remove the old copies:

        chr_donor.remove_segment(segment)

        # Normal transfer

        intergene_coordinate = chr_recipient.select_random_coordinate_in_intergenic_regions()
        intergene_coordinate = 9
        l = chr_recipient.return_location_by_coordinate(intergene_coordinate, within_intergene=True)
        position = int(l[4]) + 1

        for i, gene in enumerate(copied_segment2):
            chr_recipient.genes.insert(position + i, gene)
        for i, intergene in enumerate(new_intergene_segment):
            chr_recipient.intergenes.insert(position + i, intergene)

        cut_position = (intergene_coordinate - l[2], l[3] - intergene_coordinate)

        scar1 = chr_recipient.intergenes[int(l[4])]
        scar2 = chr_recipient.intergenes[position + i]

        if d == "left":
            r3,r4 = r4,r3

        scar1.length = r3[1] + cut_position[0]
        scar2.length = r4[0] + cut_position[1]


        return str(chr_recipient)





        # We have to register in the affected gene families that there has been a transfer event

def test_intergene(a,b,d):
    chr_donor = fill_genome(intergenic_sequences=True)
    chr_donor.obtain_flankings()
    chr_donor.obtain_locations()
    chr_recipient = fill_genome(intergenic_sequences=True)
    chr_recipient.obtain_flankings()
    chr_recipient.obtain_locations()
    r = make_transfer_intergenic(a, b, d, chr_donor, chr_recipient)
    return r


a = 14
b = 5

r1 = test_intergene(a,b,"right")
r2 = test_intergene(b,a,"left")

if r1 == r2:
    print("Success")


#make_transfer_intergenic(2, 5, "right", chr_donor, chr_recipient)
#make_transfer_intergenic(2, 5, "right", chr_donor, chr_recipient)

