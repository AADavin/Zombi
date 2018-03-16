import AuxiliarFunctions as af
import numpy
import ete3

class GeneFamily():

    def __init__(self, identifier, time):

        self.identifier = identifier
        self.origin = time

        self.genes = list()
        self.events = list()
        self.event_counter = 0  # Each time that the family is modified in any form, we have to update the event counter
        self.gene_ids_counter = 0

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

        for current_time, event, nodes in events[::-1]:

            if event == "F":

                nodename = nodes.replace(";","_")

                times[nodename] = float(current_time)
                surviving_nodes[nodename] = {"state": 1, "descendant": "None"}

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

        extanttree = ete3.Tree()
        wholetree = ete3.Tree()
        eroot = extanttree.get_tree_root()
        eroot.name = ""

        wquick_nodes = dict()
        equick_nodes = dict()

        for i, values in enumerate(events):

            current_time, event, nodes = values

            if event == "O":

                wroot = wholetree.get_tree_root()
                wroot.name = nodes + "_1"
                wquick_nodes[wroot.name] = wroot

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

        return wholetree.write(format=1), extanttree.write(format=1)


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

        whole_tree = tree.write(format=1)
        active_nodes = [x for x in tree.get_leaves() if x.is_active == True]

        if len(active_nodes) < 3:
            pruned_tree = None

        else:
            tree.prune(active_nodes, preserve_branch_length=True)
            pruned_tree = tree.write(format=1)

        return whole_tree, pruned_tree

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
        myname = "_".join(map(str, (self.gene_family, self.orientation)))
        return myname

class Chromosome():

    def __init__(self):

        self.genes = list()
        self.intergene_distances = list()
        self.shape = ""

    def select_random_position(self):

        return numpy.random.randint(len(self.genes))

    def select_random_length(self, p):

        return numpy.random.geometric(p)

    def __len__(self):

        return len(self.genes)

    def __str__(self):

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

    def remove_segment(self, segment):

        for gene in segment:
            self.genes.remove(gene)

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
