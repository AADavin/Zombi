import ete3
import numpy
import os
import random
import math
import sys
import copy
# Some auxiliar functions



def obtain_distances(whole_tree, time_counter, candidates):

    with open(whole_tree) as f:
        mytree = ete3.Tree(f.readline().strip())

    myroot = mytree.get_tree_root()

    for n1 in candidates:
        mynode1 = mytree & n1
        for n2 in candidates:
            mynode2 = mytree & n2
            if n1==n2:
                continue
            else:
                d = mynode1.get_distance(mynode2)
                e1 = mynode1.get_distance(myroot) - time_counter
                e2 = mynode2.get_distance(myroot) - time_counter
                td = d - e1 - e2
    return(td)


def normalize(array):
    total = numpy.sum(array)
    return (array/total)

def inverse(array):
    transformed_array = 1./numpy.array(array)
    return (transformed_array)

def logtransform(array):
    transformed_array = numpy.log(array)
    return (transformed_array)

def calculate_mean_genome_size(myfile):
    genomes = list()
    with open(myfile) as f:
        header = f.readline().strip().split("\t")[1:]
        for node in header:
            genomes.append(0)
        for line in f:
            handle = line.strip().split("\t")[1:]
            for i,gf in enumerate(handle):
                genomes[i] += int(gf)
    return(numpy.array(genomes).mean())

def divide_by_time_increase(array):
    transformed_array = numpy.log(array)
    return (transformed_array)

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

        if parameters == "ROOT_GENOME":
            parameters[parameter] = value.split(";")

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



'''

npd = normalize(pd)

# Divide all values by TIME_INCREASE
# Apply logarithms
# Inverse constants
# Normalize
# Sample

print(normalize((inverse(pd))))

draw = numpy.random.choice(list_of_candidates, 1, p = normalize((inverse(pd))))
print(draw)

'''