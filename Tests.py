import numpy
import random


def get_homologous_position(segment, genes):

    segment_length = len(segment)
    positions = list()
    genes_length = len(genes)

    # First we traverse the genome forwards

    for i, gene in enumerate(genes):

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

import numpy

def select_random_position(genes):

    return numpy.random.randint(len(genes))

def select_random_length(p):

    return numpy.random.geometric(p)

def obtain_affected_genes(genes, p_extension):

    # Returns the index list of the affected genes

    position = select_random_position(genes)
    length = select_random_length(p_extension)
    print(length, position)
    total_length = len(genes)
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

def invert_segment(genes, affected_genes):

    segment = [genes[x] for x in affected_genes]

    reversed_segment = segment[::-1]

    for i,x in enumerate(affected_genes):
        genes[x] = reversed_segment[i]

    return genes

def cut_and_paste(genes, affected_genes):

    segment = [genes[x] for x in affected_genes]
    new_segment = list()

    if len(segment) == len(genes):
        return 0

    for gene in segment:
        new_segment.append(genes.pop(genes.index(gene)))

    position = select_random_position()
    for i, gene in enumerate(new_segment):
        genes.insert(position + i, gene)

import ete3

def efficient_generator(events):

    # Eric's algorithm

    # First we will iterate the events from the end

    surviving_nodes = dict()
    times = dict()

    for current_time, event, nodes in events[::-1]:

        if event == "F":

            times[nodes] = float(current_time)
            surviving_nodes[nodes] = {"state": 1, "descendant": "None"}

        elif event == "E":

            times[nodes] = float(current_time)
            surviving_nodes[nodes] = {"state": 0,  "descendant": "None"}

        elif event == "S":

            p, c1, c2 = nodes.split(";")

            times[p] = float(current_time)

            if surviving_nodes[c1]["state"] == 1 and surviving_nodes [c2]["state"] == 1:

                surviving_nodes[p] = {"state": 1, "descendant": c1 + ";" + c2}

            elif surviving_nodes[c1]["state"] == 0 and surviving_nodes [c2]["state"] == 0:

                surviving_nodes[p] = {"state": 0,  "descendant": "None"}

            elif surviving_nodes[c1]["state"] == -1 and surviving_nodes[c2]["state"] == -1:

                mynode1 = find_descendant(surviving_nodes, c1)
                mynode2 = find_descendant(surviving_nodes, c2)

                surviving_nodes[p] = {"state": 1,  "descendant": mynode1 + ";" + mynode2}


            elif surviving_nodes[c1]["state"] == 1 and surviving_nodes[c2]["state"] == 0:

                surviving_nodes[p] = {"state": -1, "descendant": c1}

            elif surviving_nodes[c1]["state"] == 0 and surviving_nodes[c2]["state"] == 1:

                surviving_nodes[p] = {"state": -1, "descendant": c2}


            elif surviving_nodes[c1]["state"] == 1 and surviving_nodes[c2]["state"] == -1:

                mynode = find_descendant(surviving_nodes, c2)
                surviving_nodes[p] = {"state": 1, "descendant": c1 + ";" + mynode}

            elif surviving_nodes[c1]["state"] == -1 and surviving_nodes[c2]["state"] == 1:

                mynode = find_descendant(surviving_nodes, c1)
                surviving_nodes[p] = {"state": 1,  "descendant": mynode + ";" + c2}


            elif surviving_nodes[c1]["state"] == -1 and surviving_nodes[c2]["state"] == 0:

                mynode = find_descendant(surviving_nodes, c1)
                surviving_nodes[p] = {"state": -1, "descendant": mynode}

            elif surviving_nodes[c1]["state"] == 0 and surviving_nodes[c2]["state"] == -1:

                mynode = find_descendant(surviving_nodes, c2)
                surviving_nodes[p] = {"state": -1, "descendant": mynode}

    extanttree = ete3.Tree()
    wholetree = ete3.Tree()
    eroot = extanttree.get_tree_root()
    eroot.name = ""
    wroot = wholetree.get_tree_root()
    wroot.name = "Root"

    t = (len(events))

    wquick_nodes = dict()
    equick_nodes = dict()

    wquick_nodes["Root"] = wroot


    for i, values in enumerate(events):

        current_time, event, nodes = values

        if event == "S":

            p, c1, c2 = nodes.split(";")

            mynode = wquick_nodes[p]
            myc1 = mynode.add_child()
            myc2 = mynode.add_child()
            myc1.name = c1
            myc2.name = c2
            myc1.dist = times[c1] - times[p]
            myc2.dist = times[c2] - times[p]

            wquick_nodes[c1] = myc1
            wquick_nodes[c2] = myc2

            state = surviving_nodes[p]["state"]

            if state == 1: # Now the extant tree

                c1name, c2name = surviving_nodes[p]["descendant"].split(";")

                if eroot.name == "":
                    eroot.name = p
                    equick_nodes[p] = eroot

                mynode = equick_nodes[p]

                myc1 = mynode.add_child()
                myc2 = mynode.add_child()

                myc1.name = c1name
                myc2.name = c2name

                myc1.dist = times[c1name] - times[p]
                myc2.dist = times[c2name] - times[p]

                equick_nodes[c1name] = myc1
                equick_nodes[c2name] = myc2


    return wholetree.write(format=1), extanttree.write(format=1)


def find_descendant(surviving_nodes, node):

    found = 0
    mynode = surviving_nodes[node]["descendant"]

    while found == 0:

        if surviving_nodes[mynode]["state"] == 1:
            found = 1
        else:
            mynode = surviving_nodes[mynode]["descendant"]

    return mynode

def inheritance_of_alleles():

    population_1 = 10000
    population_2 = 600

    aleles = [2000,4000,1000,500,2500]
    new_aleles = []

    # You should shuffle the list

    for i, alele in enumerate(aleles[1:]):

        p = alele/population_1
        q = ((population_1 - alele) / population_1)
        variance = numpy.sqrt(p * q * population_1)
        total_alleles_pop2 = sum(new_aleles)

        print(alele, p, q, variance)

        if total_alleles_pop2 >= population_2:
            break
        else:
            new_value = int(numpy.random.normal(p * population_2,variance))
            if new_value <= 0:
                new_value = 0
            new_aleles.append(new_value)

    total_alleles_pop2 = sum(new_aleles)
    dif = abs(total_alleles_pop2 - population_2)
    new_aleles.insert(0, dif)
    print(new_aleles)
    print(sum(new_aleles))


import ete3

def return_vector_of_distances(time):

    distances_to_root = dict()
    distances = dict()

    with open("/Users/adriandavin/Desktop/Simulator/Tests/Test8/T/WholeTree.nwk") as f:

        mytree = ete3.Tree(f.readline().strip(), format=1)
        root = mytree.get_tree_root()
        root.name = "Root"
        for node in mytree.traverse():
            if node.is_root():
                continue
            distances[node.name] = (node.up.get_distance(root), node.get_distance(root))
            distances_to_root[node.name] = node.get_distance(root)

    alive_lineages = list()
    for k,v in distances.items():
        if time <= v[1] and time >= v[0]:
            alive_lineages.append(k)

    corrected_distances = dict()

    for node1 in alive_lineages:
        mynode1 = mytree & node1
        for node2 in alive_lineages:
            if node1 == node2:
                continue
            mynode2 = mytree & node2
            phylo_d = mynode1.get_distance(mynode2)
            td = phylo_d + (2 * time) - distances_to_root[node1] - distances_to_root[node2]
            print(node1,node2, td)



return_vector_of_distances(2.2)





