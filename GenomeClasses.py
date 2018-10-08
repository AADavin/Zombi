import AuxiliarFunctions as af
import numpy
import ete3
import ReconciledTree as RT
import random
from itertools import cycle

class GeneFamily():

    def __init__(self, identifier, time):

        self.identifier = identifier
        self.origin = time

        self.genes = list()
        self.events = list()
        self.event_counter = 0  # Each time that the family is modified in any form, we have to update the event counter
        self.gene_ids_counter = 0

        self.length = 0

    def register_event(self, time, event, genes):

        self.events.append((time, event, genes))

    def generate_tree(self):


        def find_descendant(surviving_nodes, node):

            found = 0
            mynode = surviving_nodes[node]["descendant"]

            while found == 0:

                if surviving_nodes[mynode]["state"] == 1:
                    found = 1
                else:
                    mynode = surviving_nodes[mynode]["descendant"]

            return mynode

        # Eric's algorithm

        # First we will iterate the events from the end

        events = self.events

        surviving_nodes = dict()
        times = dict()

        family_size = 0

        for current_time, event, nodes in events[::-1]:

            if event == "F":

                nodename = nodes.replace(";","_")
                times[nodename] = float(current_time)
                surviving_nodes[nodename] = {"state": 1, "descendant": "None"}

                family_size += 1

            elif event == "E" or event == "L":

                nodename = nodes.replace(";", "_")

                times[nodename] = float(current_time)
                surviving_nodes[nodename] = {"state": 0, "descendant": "None"}

            elif event == "S" or event == "D" or event == "T":

                p, g0, c1, g1, c2, g2 = nodes.split(";")

                pnodename = p + "_" + g0
                c1nodename = c1 + "_" + g1
                c2nodename = c2 + "_" + g2

                times[pnodename] = float(current_time)

                if surviving_nodes[c1nodename]["state"] == 1 and surviving_nodes[c2nodename]["state"] == 1:

                    surviving_nodes[pnodename] = {"state": 1, "descendant": c1nodename + ";" + c2nodename}

                elif surviving_nodes[c1nodename]["state"] == 0 and surviving_nodes[c2nodename]["state"] == 0:

                    surviving_nodes[pnodename] = {"state": 0, "descendant": "None"}

                elif surviving_nodes[c1nodename]["state"] == -1 and surviving_nodes[c2nodename]["state"] == -1:

                    mynode1 = find_descendant(surviving_nodes, c1nodename)
                    mynode2 = find_descendant(surviving_nodes, c2nodename)

                    surviving_nodes[pnodename] = {"state": 1, "descendant": mynode1 + ";" + mynode2}

                elif surviving_nodes[c1nodename]["state"] == 1 and surviving_nodes[c2nodename]["state"] == 0:

                    surviving_nodes[pnodename] = {"state": -1, "descendant": c1nodename}

                elif surviving_nodes[c1nodename]["state"] == 0 and surviving_nodes[c2nodename]["state"] == 1:

                    surviving_nodes[pnodename] = {"state": -1, "descendant": c2nodename}

                elif surviving_nodes[c1nodename]["state"] == 1 and surviving_nodes[c2nodename]["state"] == -1:

                    mynode = find_descendant(surviving_nodes, c2nodename)
                    surviving_nodes[pnodename] = {"state": 1, "descendant": c1nodename + ";" + mynode}

                elif surviving_nodes[c1nodename]["state"] == -1 and surviving_nodes[c2nodename]["state"] == 1:

                    mynode = find_descendant(surviving_nodes, c1nodename)
                    surviving_nodes[pnodename] = {"state": 1, "descendant": mynode + ";" + c2nodename}

                elif surviving_nodes[c1nodename]["state"] == -1 and surviving_nodes[c2nodename]["state"] == 0:

                    mynode = find_descendant(surviving_nodes, c1nodename)
                    surviving_nodes[pnodename] = {"state": -1, "descendant": mynode}

                elif surviving_nodes[c1nodename]["state"] == 0 and surviving_nodes[c2nodename]["state"] == -1:

                    mynode = find_descendant(surviving_nodes, c2nodename)
                    surviving_nodes[pnodename] = {"state": -1, "descendant": mynode}

        extanttree = RT.ReconciledTree()
        completetree = RT.ReconciledTree()

        eroot = extanttree.get_tree_root()
        eroot.name = ""

        wquick_nodes = dict()
        equick_nodes = dict()

        for i, values in enumerate(events):

            current_time, event, nodes = values

            if event == "O":

                wroot = completetree.get_tree_root()
                wroot.name = nodes + "_1"
                wquick_nodes[wroot.name] = wroot

            if event == "L" or event == "E":

                p, g0 = nodes.split(";")
                pnodename = p + "_" + g0
                mynode = wquick_nodes[pnodename]
                e = RT.RecEvent("L", p, int(float(current_time)))
                mynode.addEvent(e, append=True)

            if event == "F":

                p, g0 = nodes.split(";")
                pnodename = p + "_" + g0
                mynode = wquick_nodes[pnodename]
                e = RT.RecEvent("C", p, int(float(current_time)))
                mynode.addEvent(e, append=True)

            if event == "S" or event == "D" or event == "T":

                p, g0, c1, g1, c2, g2 = nodes.split(";")
                pnodename = p + "_" + g0
                c1nodename = c1 + "_" + g1
                c2nodename = c2 + "_" + g2

                mynode = wquick_nodes[pnodename]
                myc1 = mynode.add_child()
                myc2 = mynode.add_child()
                myc1.name = c1nodename
                myc2.name = c2nodename
                myc1.dist = times[c1nodename] - times[pnodename]
                myc2.dist = times[c2nodename] - times[pnodename]

                wquick_nodes[c1nodename] = myc1
                wquick_nodes[c2nodename] = myc2

                state = surviving_nodes[pnodename]["state"]

                ### Now we add the reconciled events

                e = RT.RecEvent(event, p, int(float(current_time)))
                mynode.addEvent(e, append=True)

                if state == 1:  # Now the extant tree

                    c1name, c2name = surviving_nodes[pnodename]["descendant"].split(";")

                    if eroot.name == "":
                        eroot.name = pnodename
                        equick_nodes[pnodename] = eroot

                    mynode = equick_nodes[pnodename]

                    myc1 = mynode.add_child()
                    myc2 = mynode.add_child()

                    myc1.name = c1name
                    myc2.name = c2name

                    myc1.dist = times[c1name] - times[pnodename]
                    myc2.dist = times[c2name] - times[pnodename]

                    equick_nodes[c1name] = myc1
                    equick_nodes[c2name] = myc2

        if family_size == 0:

            extanttree = ";"

        elif family_size == 1:


            extanttree = [k for k, v in surviving_nodes.items() if v["state"] == 1 and v["descendant"] == "None"][
                             0] + ";"

        else:

            extanttree = extanttree.write(format=1, format_root_node=True)


        rec = completetree.getTreeRecPhyloXML()

        if len(completetree) == 0:
            completetree = ";"
        elif len(completetree) == 1:
            completetree = completetree.get_leaves()[0].name + ";"
        else:
            completetree = completetree.write(format=1, format_root_node=True)


        return completetree, extanttree, rec


    def generate_oldtree(self):

        tree = ete3.Tree()

        current_time, event, nodes = self.events[0]

        sp = tree.get_tree_root()
        sp.name = nodes + "_1"
        sp.add_feature("is_active", True)

        elapsed_time = float(current_time)

        for current_time, event, nodes in self.events[1:]:

            elapsed_time = float(current_time) - elapsed_time
            active_nodes = [x for x in tree.get_leaves() if x.is_active == True]
            for node in active_nodes:
                node.dist += elapsed_time
            elapsed_time = float(current_time)

            if event == "S":

                sp, gp, c1, g1, c2, g2 = nodes.split(";")

                myname = sp + "_" + gp
                mynode = tree & myname
                mynode.is_active = False

                gc1 = mynode.add_child(dist=0)
                gc1.name = c1 + "_" + g1
                gc1.add_feature("is_active", True)

                gc2 = mynode.add_child(dist=0)
                gc2.name = c2 + "_" + g2
                gc2.add_feature("is_active", True)

            elif event == "E":
                sp, gp = nodes.split(";")
                myname = sp + "_" + gp
                mynode = tree & myname
                mynode.is_active = False

            elif event == "L":
                sp, gp = nodes.split(";")
                myname = sp + "_" + gp
                mynode = tree & myname
                mynode.is_active = False

            elif event == "D":

                sp, gp, c1, g1, c2, g2 = nodes.split(";")
                myname = sp + "_" + gp
                mynode = tree & myname

                mynode.is_active = False

                gc1 = mynode.add_child(dist=0)
                gc1.name = c1 + "_" + g1
                gc1.add_feature("is_active", True)

                gc2 = mynode.add_child(dist=0)
                gc2.name = c2 + "_" + g2
                gc2.add_feature("is_active", True)

            elif event == "T":
                sp, gp, c1, g1, c2, g2 = nodes.split(";")

                myname = sp + "_" + gp

                mynode = tree & myname
                mynode.is_active = False

                gc1 = mynode.add_child(dist=0)
                gc1.name = c1 + "_" + g1
                gc1.add_feature("is_active", True)

                gc2 = mynode.add_child(dist=0)
                gc2.name = c2 + "_" + g2
                gc2.add_feature("is_active", True)

            elif event == "F":
                break

        complete_tree = tree.write(format=1, format_root_node=True)
        active_nodes = [x for x in tree.get_leaves() if x.is_active == True]

        if len(active_nodes) < 3:
            pruned_tree = None

        else:
            tree.prune(active_nodes, preserve_branch_length=True)
            pruned_tree = tree.write(format=1, format_root_node=True)

        return complete_tree, pruned_tree

    def obtain_new_gene_id(self):
        self.gene_ids_counter += 1
        return self.gene_ids_counter

    def __str__(self):

        return "GeneFamily_" + str(self.identifier) + ";" + ";".join([str(x) for x in self.genes])

        #return ";".join(map(str, self.genes))

    def __len__(self):

        return len([x for x in self.genes if x.active == True])

    def __iter__(self):
        for gene in self.genes:
            yield gene


class Gene():

    def __init__(self):

        self.active = True
        self.orientation = ""
        self.gene_family = ""
        self.gene_id = ""
        self.sequence = ""
        self.species = ""
        self.importance = 0
        self.length = 0
        self.total_flanking = 0
        self.specific_flanking = 0

    def determine_orientation(self):

        if numpy.random.binomial(1,0.5):
            self.orientation = "+"

        else:
            self.orientation = "-"

    def change_sense(self):

        if self.orientation == "+":
            self.orientation = "-"
        elif self.orientation == "-":
            self.orientation = "+"

    def __str__(self):
        #myname = "_".join(map(str,(self.genome, self.orientation, self.gene_family, self.gene_id)))
        #myname = "_".join(map(str, (self.species, self.gene_family, self.gene_id)))
        #myname = "_".join(map(str, (self.gene_family, self.orientation)))
        myname = str(self.gene_family)
        #myname = "_".join(map(str, (self.gene_family, self.length)))
        return myname


class Intergene():

    def __init__(self):

        self.length = 0
        self.total_flanking = 0
        self.specific_flanking = 0
        self.id = 0 # Only for debugging purposes

    def __str__(self):

        return "(" + str(self.length) + ")"
        #return "I_" + str(self.length)
        #return "I_" + str(self.id) + "_" + str(self.length)


class Chromosome():

    def __init__(self):

        self.has_intergenes = False
        self.intergenes = list()
        self.genes = list()
        self.shape = ""
        self.length = 0

        self.total_locations = list()
        self.map_of_locations = list()

    def obtain_total_itergenic_length(self):

        total_length = 0
        for intergene in self.intergenes:
            total_length += intergene.length
        return total_length


    def select_random_position(self):

        return numpy.random.randint(len(self.genes))

    def obtain_flankings(self):

        if self.has_intergenes:

            self.genes[0].total_flanking = (0, self.genes[0].length)
            self.genes[0].specific_flanking = (0, self.genes[0].length)
            self.intergenes[0].total_flanking = (self.genes[0].total_flanking[1], self.genes[0].total_flanking[1] + self.intergenes[0].length)
            self.intergenes[0].specific_flanking = (0, self.intergenes[0].length)

            for i in range(len(self.genes)):

                if i == 0:
                    continue

                lb = self.intergenes[i-1].total_flanking[1]
                ub = lb + self.genes[i].length

                lbg = self.genes[i-1].specific_flanking[1] + 1
                ubg = lbg + self.genes[i].length

                self.genes[i].total_flanking = (lb, ub)
                self.genes[i].specific_flanking = (lbg, ubg)

                lbi = self.intergenes[i - 1].specific_flanking[1] + 1
                ubi = lbi + self.intergenes[i].length

                self.intergenes[i].total_flanking = (ub, ub + self.intergenes[i].length)
                self.intergenes[i].specific_flanking = (lbi, ubi)

            #self.intergenes[i].total_flanking = (ub, 0)

    def obtain_locations(self):

        self.map_of_locations = list()

        # The structure of total location is:
        # tc1, tc2, Specific coordinate 1, specific coordinate 2, position, Gene/intergene

        for i in range(len(self.genes)):

            tc1 = self.genes[i].total_flanking[0]
            tc2 = self.genes[i].total_flanking[1]
            sc1 = self.genes[i].specific_flanking[0]
            sc2 = self.genes[i].specific_flanking[1]
            self.map_of_locations.append((tc1, tc2, sc1, sc2, str(i), "G"))

            tc1 = self.intergenes[i].total_flanking[0]
            tc2 = self.intergenes[i].total_flanking[1]
            sc1 = self.intergenes[i].specific_flanking[0]
            sc2 = self.intergenes[i].specific_flanking[1]
            self.map_of_locations.append((tc1, tc2, sc1, sc2, str(i), "I"))


    def select_random_coordinate_in_intergenic_regions(self):

        # We weight the position by the length of the region

        t = sum([x.length for x in self.intergenes]) + len(self.intergenes) - 1

        return random.randint(0, t)

    def return_total_coordinate_from_specific_coordinate(self, c, type = "I", debug = False):

        tc = None
        for r in self.map_of_locations:
            tc1, tc2, spc1, spc2, sp, t = r
            if debug == True:
                print(r)
            if t != type:
                continue
            if c >= spc1 and c <= spc2:
                distance_to_lower_bound = c - spc1
                tc = tc1 + distance_to_lower_bound
        return tc

    def return_specific_coordinate_from_total_coordinate(self, c, debug = False):


        sc = None
        for r in self.map_of_locations:
            tc1, tc2, spc1, spc2, sp, t = r
            if debug == True:
                print(r)
            if t == "I" and c >= tc1 and c <= tc2:

                distance_to_lower_bound = c - tc1
                sc = spc1 + distance_to_lower_bound

                if debug == True:
                    print("PRINTING R")
                    print(r)

        return sc

    def return_location_by_coordinate(self, c, within_intergene = False):

        if within_intergene == False:

            for l in self.map_of_locations:
                tc1, tc2, spc1, spc2, sp, t = l
                if c >= tc1 and c <= tc2:
                    return l
        else:

            for l in self.map_of_locations:
                tc1, tc2, spc1, spc2, sp, t = l
                if t != "I":
                    continue
                if c >= spc1 and c <= spc2:
                    return l

    def return_affected_region(self, c1, c2, direction):

        # It returns a tuple
        # 1. List of the position of the genes affected. ALWAYS FROM LEFT TO RIGHT
        # 2. List of the position of the intergenes affected. ALWAYS FROM LEFT TO RIGHT
        # 3. Tuple with left and right cuts of first intergene.
        # #  Watch out, the fact of calling it left or right can be confusing.
        # 4. Tuple with left and right cuts of last intergene. Same note than above

        l1 = self.return_location_by_coordinate(c1, within_intergene=True)
        l2 = self.return_location_by_coordinate(c2, within_intergene=True)

        tc1_1, tc1_2, sc1_1, sc1_2, p1, t1, = l1
        tc2_1, tc2_2, sc2_1, sc2_2, p2, t2 = l2

        p1 = int(p1)
        p2 = int(p2)

        affected_genes = list()
        affected_intergenes = list()

        t_length = len(self.intergenes)

        left_limits = None
        right_limits = None

        if c1 == c2:
            return None

        elif p1 == p2:
            return None

        elif c1 < c2 and direction == "right":

            affected_genes = [i + 1 for i in range(p1, p2)]
            affected_intergenes = [i for i in range(p1,p2 + 1)]
            left_limits = (c1 - sc1_1, sc1_2 - c1)
            right_limits = (c2 - sc2_1, sc2_2 - c2)

        elif c1 > c2 and direction == "right":

            affected_genes = [i + 1 for i in range(p1, t_length - 1)]
            affected_genes += [i for i in range(0, p2 + 1)]

            affected_intergenes = [i for i in range(p1, t_length)]
            affected_intergenes += [i for i in range(0, p2 + 1)]

            left_limits = (c1 - sc1_1, sc1_2 - c1)
            right_limits = (c2 - sc2_1, sc2_2 - c2)

        elif c1 > c2 and direction == "left":

            affected_genes = [i for i in range(p1, p2, - 1)]
            affected_intergenes = [i for i in range(p1, p2 - 1, -1)]
            left_limits = (c1 - sc1_1, sc1_2 - c1)
            right_limits = (c2 - sc2_1, sc2_2 - c2)

            affected_genes.reverse()
            affected_intergenes.reverse()

        elif c1 < c2 and direction == "left":

            affected_genes = [i for i in range(p1, -1, - 1)]
            affected_genes += [i for i in range(t_length - 1, p2, -1)]

            affected_intergenes = [i for i in range(p1, -1, -1)]
            affected_intergenes += [i for i in range(t_length - 1, p2 - 1, -1)]

            affected_genes.reverse()
            affected_intergenes.reverse()

            left_limits = (c1 - sc1_1, sc1_2 - c1)
            right_limits = (c2 - sc2_1, sc2_2 - c2)

        return (affected_genes, affected_intergenes, left_limits, right_limits)


    def select_random_length(self, p):

        # Watch out here

        if not self.has_intergenes:
            return numpy.random.geometric(p)
        else:
            return numpy.random.geometric(p)

    def __len__(self):

        # Watch out!! This is probably no the safest thing to do
        if not self.has_intergenes:
            return len(self.genes)
        else:
            return self.length

    def __str__(self):

        if self.has_intergenes == True:

            #return ";".join(["CHROMOSOME"] + [str(self.genes[i])+";"+str(self.intergenes[i]) for i in range(len(self.genes))])
            return "".join(["CHROMOSOME: "] + [str(self.genes[i])+str(self.intergenes[i]) for i in range(len(self.genes))])

        else:

            return ";".join(["CHROMOSOME"] + [str(gene) for gene in self.genes])

    def __iter__(self):

        for x in self.genes:
            yield x


class CircularChromosome(Chromosome):

    def __init__(self):
        super().__init__()

    def obtain_segment(self, affected_genes):

        segment = [self.genes[x] for x in affected_genes]

        return segment

    def obtain_intergenic_segment(self, affected_intergenes):

        segment = [self.intergenes[x] for x in affected_intergenes]

        return segment

    def remove_segment(self, segment):

        for gene in segment:
            self.genes.remove(gene)

    def remove_intersegment(self, intersegment):

        for intergene in intersegment:
            self.intergenes.remove(intergene)

    def insert_segment(self, position, segment):

        for i, x in enumerate(segment):
            self.genes.insert(position + i, x)

    def invert_segment(self, affected_genes):

        segment = [self.genes[x] for x in affected_genes]

        reversed_segment = segment[::-1]

        for gene in reversed_segment:
            gene.change_sense()

        for i,x in enumerate(affected_genes):
            self.genes[x] = reversed_segment[i]


    def cut_and_paste(self, affected_genes):

        segment = [self.genes[x] for x in affected_genes]
        new_segment = list()

        if len(segment) == len(self.genes):
            return 0

        for gene in segment:
            new_segment.append(self.genes.pop(self.genes.index(gene)))

        position = self.select_random_position()
        for i, gene in enumerate(new_segment):
            self.genes.insert(position + i, gene)


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



    def get_homologous_position(self, segment):

        homologous = list()

        segment_length = len(segment)
        genes_length = len(self.genes)

        genes = [x.gene_family + "_" + x.orientation for x in self.genes]
        mysegment = [x.gene_family + "_" + x.orientation for x in segment]


        # First we traverse the genome forwards

        for i, gene in enumerate(genes):

            length_counter = 0
            positions = list()

            name_gene_in_genome = gene
            name_gene_in_segment = mysegment[0]

            if name_gene_in_genome == name_gene_in_segment:

                positions.append(i)
                length_counter += 1

                for j, x in enumerate(mysegment):
                    if length_counter == segment_length:
                        homologous.append(("F", tuple(positions)))
                        break
                    if 1 + i + j >= genes_length:
                        if genes[(i + j + 1) - genes_length] == mysegment[j + 1]:
                            positions.append(i+j+1 - genes_length)
                            length_counter += 1
                        else:
                            break
                    else:
                        if genes[i + j + 1] == mysegment[j + 1]:
                            positions.append(i + j + 1)
                            length_counter += 1
                        else:
                            break

        # Second we traverse the genome backwards

        inverted_segment = list()

        for gene in mysegment[::-1]:

            inverted_segment.append(gene.replace("+","A").replace("-", "+").replace("A", "-"))

        for i, gene in enumerate(genes):

            positions = list()
            length_counter = 0
            name_gene_in_segment = inverted_segment[0]

            if genes[i] == name_gene_in_segment:

                positions.append(i)
                length_counter += 1

                for j, x in enumerate(inverted_segment):
                    if length_counter == segment_length:
                        homologous.append(("B", tuple(positions)))
                        break
                    if 1 + i + j >= genes_length:
                        if genes[(i + j + 1) - genes_length] == inverted_segment[j + 1]:
                            positions.append(i + j + 1 - genes_length)
                            length_counter += 1
                        else:
                            break
                    else:
                        if genes[i + j + 1] == inverted_segment[j + 1]:
                            positions.append(i + j + 1)
                            length_counter += 1
                        else:
                            break

        return homologous

    ### From here, functions related to intergenic regions

    def remove_segment_with_intergenic(self, segment):

        for gene in segment:
            self.genes.remove(gene)

    def insert_gene_within_intergene(self, coordinate, location, gene):

        tc1, tc2, sc1, sc2, position, t = location

        # Convert to ints:

        tc1, tc2, sc1, sc2, position = map(int,(tc1,tc2,sc1,sc2,position))

        # The first part is easier - We simply add the gene to the list of genes

        self.genes.insert(position + 1, gene)

        # The second part is cutting the intergene, obtaining the distances

        left_limits = coordinate - sc1
        right_limits = sc2 - coordinate

        # Now we insert the new intergene in the position i + 1

        intergene = Intergene()
        intergene.length = right_limits
        self.intergenes[position].length = left_limits
        self.intergenes.insert(position + 1, intergene)


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

        # I have to weight by the length of each chromosome

        chromosome = numpy.random.choice(self.chromosomes, 1, p=af.normalize([len(x) for x in self.chromosomes]))[0]
        return chromosome

    def update_genome_species(self, species):

        self.species = species

        for ch in self.chromosomes:

            for gene in ch:
                gene.species = species

    def __str__(self):

        return ";".join(["GENOME:"] + [str(x) for x in self.chromosomes])

    def __iter__(self):
        for chromosome in self.chromosomes:
            yield chromosome
