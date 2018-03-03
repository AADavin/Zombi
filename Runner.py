import AuxiliarFunctions as af
import numpy
import copy
import ete3
import random
import os
from SimulatorClasses import GenomeSimulator
#numpy.random.seed(245)
#random.seed(10)
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

gss = GenomeSimulator(myevents)
gss.verbose_run()


#p = read_parameters("/Users/adriandavin/PycharmProjects/SimuLYON/Parameters/SpeciesTreeParameters.tsv")
#p = prepare_parameters(p)



#gss.write_genomes("/Users/adriandavin/Desktop/Bioinformatics/SimuLyon/Cedric/Gillespie/All/Genomes/")
#gss.write_gene_family_events("/Users/adriandavin/Desktop/Bioinformatics/SimuLyon/Cedric/Gillespie/All/Gene_families/")
gss.write_gene_trees("/Users/adriandavin/Desktop/Bioinformatics/SimuLyon/Cedric/Gillespie/All/Gene_trees/")
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

