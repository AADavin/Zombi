

def test_TreeGenerator():

    parameters_file = "/Users/adriandavin/PycharmProjects/SimuLYON/Parameters_file.tsv"
    tg = TreeGenerator(parameters_file)
    tg.generate_tree_mode_0()
    outpath = "/Users/adriandavin/Desktop/Bioinformatics/SimuLyon/NEWTEST"
    tg.store_log(outpath)

def test_GeneFamilySimulator():

    parameters_file = "/Users/adriandavin/PycharmProjects/SimuLYON/Parameters_file.tsv"
    events_file = "/Users/adriandavin/Desktop/Bioinformatics/SimuLyon/NEWTEST/SpeciesTreeEvents.tsv"
    lineages_file = "/Users/adriandavin/Desktop/Bioinformatics/SimuLyon/NEWTEST/LineagesInTime.tsv"

    gf = GeneFamilySimulator(parameters_file,events_file, lineages_file)
    gf.origination("Root",0,0)
    gf.run(0)
    gf.origination("6", 800, 1)
    gf.run(1)

    print(gf.gene_families[0]["Gene_tree"].write(format=1))
    print(gf.gene_families[1]["Gene_tree"].write(format=1))

def test_FamilyOriginator():

     fo = FamilyOriginator("/Users/adriandavin/Desktop/Bioinformatics/SimuLyon/NEWTEST/WholeTree")
     for i in range (500):
         print(fo.create_families(2))

def test_GenomeSimulator():

    parameters_file = "/Users/adriandavin/PycharmProjects/SimuLYON/Parameters_file.tsv"
    events_file = "/Users/adriandavin/Desktop/Bioinformatics/SimuLyon/NEWTEST/SpeciesTreeEvents.tsv"
    lineages_file = "/Users/adriandavin/Desktop/Bioinformatics/SimuLyon/NEWTEST/LineagesInTime.tsv"

    gf = GeneFamilySimulator(parameters_file, events_file, lineages_file)
    fo = FamilyOriginator("/Users/adriandavin/Desktop/Bioinformatics/SimuLyon/NEWTEST/WholeTree")

    for i in range(20000):
        node, time = fo.create_families(2)
        gf.origination(node, time, i)
        gf.run(i)

        print(node, time, i)
        #print(gf.gene_families[i]["Gene_tree"].write(format=1))

def complete_test():

    #test_TreeGenerator()
    #test_GeneFamilySimulator()
    #test_FamilyOriginator()
    test_GenomeSimulator()

#complete_test()
#import numpy
#print(numpy.random.exponential(100))

def get_homologous_position(segment, genes):

    segment_length = len(segment)
    positions = list()
    genes_length = len(genes)

    # First we traverse the genome forwards

    for i, gene in enumerate(genes):

        #name_gene_in_genome = gene.orientation + "_" + gene.gene_family
        #name_gene_in_segment = segment[0].orientation + "_" + segment[0].gene_family

        length_counter = 0

        name_gene_in_genome = gene
        name_gene_in_segment = segment[0]

        if name_gene_in_genome == name_gene_in_segment:
            length_counter += 1
            for j, x in enumerate(segment):
                if length_counter == segment_length:
                    positions.append(("F",i))
                    break
                if 1 + i + j >= genes_length:
                    if genes[(i + j + 1) - genes_length] == segment[j + 1]:
                        length_counter += 1
                    else:
                        break
                else:
                     if genes[i + j + 1] == segment[j + 1]:
                          length_counter += 1
                     else:
                         break

    # Second we traverse the genome backwards

    inverted_genome = list()

    for gene in genes[::-1]:
        if "+" in gene:
            name_gene_in_genome = gene.replace("+","-")
        elif "-" in gene:
            name_gene_in_genome = gene.replace("-", "+")
        inverted_genome.append(name_gene_in_genome)

    for i, gene in enumerate(inverted_genome):

        # name_gene_in_genome = gene.orientation + "_" + gene.gene_family
        # name_gene_in_segment = segment[0].orientation + "_" + segment[0].gene_family

        length_counter = 0

        name_gene_in_segment = segment[0]

        if inverted_genome[i] == name_gene_in_segment:
            length_counter += 1
            for j, x in enumerate(segment):
                if length_counter == segment_length:
                    positions.append(("B", i))
                    break
                if 1 + i + j >= genes_length:
                    if inverted_genome[(i + j + 1) - genes_length] == segment[j + 1]:
                        length_counter += 1
                    else:
                        break
                else:
                    if inverted_genome[i + j + 1] == segment[j + 1]:
                        length_counter += 1
                    else:
                        break

    return positions

# numpy.random.seed(245)
# random.seed(10)
# GenomeEvolver()
# stg = SpeciesTreeGenerator()
# stg.generate_new_tree()
# print(stg.whole_species_tree.write(format=1))

myevents = "/Users/adriandavin/Desktop/Bioinformatics/SimuLyon/Cedric/Gillespie/Events.tsv"
# mytree = "/Users/adriandavin/Desktop/Bioinformatics/SimuLyon/Cedric/Gillespie/WholeTree"
tree_file = "/Users/adriandavin/Desktop/Bioinformatics/SimuLyon/Cedric/Gillespie/CyanoTree"
# stg.write_events_file(myevents)
# stg.write_tree(mytree)

# species_tree_events = list()
# species_counter = 0

gss = GenomeSimulator(myevents)
gss.verbose_run()

# p = read_parameters("/Users/adriandavin/PycharmProjects/SimuLYON/Parameters/SpeciesTreeParameters.tsv")
# p = prepare_parameters(p)

# gss.write_genomes("/Users/adriandavin/Desktop/Bioinformatics/SimuLyon/Cedric/Gillespie/All/Genomes/")
# gss.write_gene_family_events("/Users/adriandavin/Desktop/Bioinformatics/SimuLyon/Cedric/Gillespie/All/Gene_families/")
gss.write_gene_trees("/Users/adriandavin/Desktop/Bioinformatics/SimuLyon/Cedric/Gillespie/All/Gene_trees/")
# gss.write_events_per_branch("/Users/adriandavin/Desktop/Bioinformatics/SimuLyon/Cedric/Gillespie/All/Events_per_branch")

# for item in gss.all_gene_families["3"].events:
#    print(item)

# print(gss.all_gene_families["3"].generate_tree())
# generate_events(tree_file)
