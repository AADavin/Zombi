from TreeGenerator import TreeGenerator
from GenomeGenerator import GeneFamilySimulator
from GenomeGenerator import FamilyOriginator

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

import numpy

print(numpy.random.randint(1))





