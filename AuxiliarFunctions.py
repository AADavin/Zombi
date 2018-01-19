import ete3
import numpy
import os
import random
import math
import sys

# Some auxiliar functions

TIME_INCREASE = 0.001
TOTAL_TIME = 1

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

def generate_tree(events_file):

    active_lineages = dict()

    tree_events = dict()

    with open(events_file) as f:
        f.readline()
        for line in f:
            time, dn, ln, cld = line.strip().split("\t")
            if int(time) not in tree_events:
                tree_events[int(time)] = list()
            tree_events[int(time)].append((dn, ln, cld))

    mytree = ete3.Tree()
    mytree.name = "Root"
    mytree.add_feature("is_alive", True)

    for time_counter in range(int(TOTAL_TIME / TIME_INCREASE)):
        active_lineages[time_counter] = list()
        lineages_alive = [x for x in mytree.get_leaves() if x.is_alive == True]
        active_lineages[time_counter] = lineages_alive

        for lineage in lineages_alive:
            lineage.dist += TIME_INCREASE

        if time_counter in tree_events: # In this case the gene follows the species tree

            for event, snode, children in tree_events[time_counter]:
                if event == "EX":
                    mynode = mytree&snode
                    mynode.is_alive = False
                elif event == "SP":
                    mynode = mytree&snode
                    sc1, sc2 = children.split("+")
                    c1 = mynode.add_child(dist=0)
                    c1.name = sc1
                    c1.add_feature("is_alive", True)
                    c2 = mynode.add_child(dist=0)
                    c2.name = sc2
                    c2.add_feature("is_alive", True)
                    mynode.is_alive = False

    print(mytree.write(format=1))

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