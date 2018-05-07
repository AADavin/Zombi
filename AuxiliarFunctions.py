import ete3
import numpy
import copy

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

    handle = value.split(":")

    if handle[0] == "f":
        # Fixed value
        value =  float(handle[1])

    elif handle[0] == "n":
        # normal distribution
        params = handle[1].split(";")
        value = abs(numpy.random.normal(float(params[0]), float(params[1])))

    elif handle[0] == "l":
        # lognormal distribution
        params = handle[1].split(";")
        value = abs(numpy.random.lognormal(float(params[0]), float(params[1])))

    elif handle[0] == "u":
        # uniform distribution
        params = handle[1].split(";")
        value = abs(numpy.random.uniform(float(params[0]), float(params[1])))

    return value


def prepare_sequence_parameters(parameters):

    for parameter, value in parameters.items():

        if parameter == "SEQUENCE_SIZE" or parameter == "VERBOSE":
            parameters[parameter] = int(value)

        if parameter == "SCALING":
            parameters[parameter] = float(value)

    return parameters

def prepare_species_tree_parameters(parameters):

    for parameter, value in parameters.items():

        if parameter == "TURNOVER":
            parameters[parameter] = obtain_value(value)

        if parameter == "TOTAL_TIME":
            parameters[parameter] = float(value)

        if parameter == "POPULATION_SIZES":
            parameters[parameter] = [tuple([int(j) for j in x.split("-")]) for x in value.split(";")]

        if parameter == "SPECIES_EVOLUTION_MODE" or parameter == "N_LINEAGES" or parameter == "MIN_LINEAGES" \
                or parameter == "TOTAL_LINEAGES" or parameter == "STOPPING_RULE" or parameter == "MAX_LINEAGES"\
                or parameter == "VERBOSE":
            parameters[parameter] = int(value)

    return parameters


def fasta_reader(fasta_file):

    with open(fasta_file) as f:

        seq = ""
        for line in f:
            if ">" == line[0]:
                if seq != "":
                    yield header, seq
                    header = line.strip()
                    seq = ""
                else:
                    header = line.strip()
                    seq = ""
            else:
                seq += line.strip()

        yield header, seq

def fasta_writer(outfile, entries):

    x = 80
    with open(outfile, "w") as f:
        for h, seq in entries:
            f.write(h + "\n")
            lines = [seq[i: i + x] for i in range(0, len(seq), x)]
            for line in lines:
                f.write(line +"\n")


def prepare_genome_parameters(parameters):

    for parameter, value in parameters.items():

        #if parameter == "DUPLICATION_EXTENSION" or parameter == "TRANSFER_EXTENSION" \
        #        or parameter == "LOSS_EXTENSION" or parameter == "INVERSION_EXTENSION" or \
        #        parameter == "TRANSLOCATION_EXTENSION" or parameter == "ORIGINATION_EXTENSION":

        #    parameters[parameter] = obtain_value(value)

        if parameter == "ROOT_GENOME":
            parameters[parameter] = value.split(";")

        if parameter == "REPLACEMENT_TRANSFER":
            parameters[parameter] = float(value)

        if parameter == "PROFILES" or parameter == "EVENTS_PER_BRANCH" or parameter == "GENE_TREES" \
                or parameter == "PRUNE_TREES" or parameter == "TRANSFER_PREFERENCE" or parameter == "RECONCILED_TREES" or parameter =="VERBOSE":

            parameters[parameter] = int(value)

    return parameters

def generate_events(tree_file):

    events = []

    with open(tree_file) as f:
        tree = ete3.Tree(f.readline().strip(), format=1)

    root = tree.get_tree_root()
    root.name = "Root"
    total_time = root.get_farthest_leaf()[1]

    # There might be slighlty variations in the branch length that we have to account for. So all nodes
    # that are at 0.1% distance of the farthes leaf will be considered to be alive

    error_margin = total_time * 0.001

    nodes = list()

    for node in tree.traverse():

        node_dist = node.get_distance(root)
        if node.is_leaf():
            if  total_time <= node_dist + error_margin  and total_time >= node_dist - error_margin:
                nodes.append((node, "A", node_dist))
            else:
                nodes.append((node, "E", node_dist))
        else:
            nodes.append((node, "S", node_dist))

    # We order the events # This can be more efficient. But who cares

    nodes = sorted(nodes, key = lambda x: x[2])

    for node, estate, time in nodes:

        if estate == "A":
            events.append((str(time), "F", node.name))

        elif estate == "E":
            events.append((str(time), "E", node.name))
        elif estate == "S":
            c1, c2 = node.get_children()
            events.append((str(time), "S", ";".join((node.name, c1.name, c2.name))))

    return events


def copy_segment(segment, new_identifiers):

    new_segment = list()

    for i,gene in enumerate(segment):
        new_gene = copy.deepcopy(gene)
        new_gene.gene_id = new_identifiers[i]
        new_segment.append(new_gene)

    return new_segment


def return_vector_of_distances(self, tree_file):

    self.distances_to_root = dict()

    with open(tree_file) as f:

        self.mytree = ete3.Tree(f.readline().strip(), format=1)
        root = self.mytree.get_tree_root()
        root.name = "Root"
        for node in self.mytree.traverse():
            if node.is_root():
                continue

            self.distances_to_root[node.name] = (node, node.get_distance(root))

def choose_advanced_recipient(self, time, alive_lineages, donor):

    # Chooses and advanced recipient according to the logarithm of the phylogenetic distance

    possible_recipients = list()
    weights = list()

    mydonor = self.distances_to_root[donor][0]

    for recipient in alive_lineages:

        if donor == recipient:
            continue

        myrecipient = self.distances_to_root[recipient][0]
        phylo_d = mydonor.get_distance(myrecipient)

        td = phylo_d + (2 * time) - self.distances_to_root[donor][1] - self.distances_to_root[recipient][1]

        possible_recipients.append(recipient)
        weights.append(td)

    draw = numpy.random.choice(possible_recipients, 1, p= normalize(weights))


def generate_newick_trees(events):

    ### THIS IS FOR GENERATING SPECIES TREES, AND MOST LIKELY REDUNDANT

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

    surviving_nodes = dict()
    times = dict()

    for current_time, event, nodes in events[::-1]:


        if event == "F":

            times[nodes] = float(current_time)
            surviving_nodes[nodes] = {"state": 1, "descendant": "None"}

        elif event == "E":

            times[nodes] = float(current_time)
            surviving_nodes[nodes] = {"state": 0, "descendant": "None"}

        elif event == "S":

            p, c1, c2 = nodes.split(";")

            times[p] = float(current_time)

            if surviving_nodes[c1]["state"] == 1 and surviving_nodes[c2]["state"] == 1:

                surviving_nodes[p] = {"state": 1, "descendant": c1 + ";" + c2}

            elif surviving_nodes[c1]["state"] == 0 and surviving_nodes[c2]["state"] == 0:

                surviving_nodes[p] = {"state": 0, "descendant": "None"}

            elif surviving_nodes[c1]["state"] == -1 and surviving_nodes[c2]["state"] == -1:

                mynode1 = find_descendant(surviving_nodes, c1)
                mynode2 = find_descendant(surviving_nodes, c2)

                surviving_nodes[p] = {"state": 1, "descendant": mynode1 + ";" + mynode2}


            elif surviving_nodes[c1]["state"] == 1 and surviving_nodes[c2]["state"] == 0:

                surviving_nodes[p] = {"state": -1, "descendant": c1}

            elif surviving_nodes[c1]["state"] == 0 and surviving_nodes[c2]["state"] == 1:

                surviving_nodes[p] = {"state": -1, "descendant": c2}


            elif surviving_nodes[c1]["state"] == 1 and surviving_nodes[c2]["state"] == -1:

                mynode = find_descendant(surviving_nodes, c2)
                surviving_nodes[p] = {"state": 1, "descendant": c1 + ";" + mynode}

            elif surviving_nodes[c1]["state"] == -1 and surviving_nodes[c2]["state"] == 1:

                mynode = find_descendant(surviving_nodes, c1)
                surviving_nodes[p] = {"state": 1, "descendant": mynode + ";" + c2}


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

            if state == 1:  # Now the extant tree

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



def generate_gene_tree(events):

 # THIS IS FOR GENERATING GENE TREES


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


    if len(wholetree) == 0:
        wholetree = ";"

    elif len(wholetree) == 1:
        wholetree =  wholetree.get_leaves()[0].name + ";"

    else:
        wholetree = wholetree.write(format=1)

    if len(extanttree) == 0:
        extanttree = ";"

    elif len(extanttree) == 1:
        extanttree = extanttree.get_leaves()[0].name + ";"

    else:
        extanttree = extanttree.write(format=1)

    return wholetree, extanttree


def write_pruned_sequences(tree_file, fasta_folder):
    with open(tree_file) as f:
        line = f.readline().strip()
        if "(" not in line or line == ";":
            return None
        else:
            my_tree = ete3.Tree(line, format=1)

    surviving_nodes = {x.name for x in my_tree.get_leaves()}
    file_name = tree_file.split("/")[-1].split("_")[0]
    entries = fasta_reader(fasta_folder + "/" + file_name + "_whole.fasta")

    clean_entries = list()
    for h, seq in entries:
        if h[1:] in surviving_nodes:
            clean_entries.append((h, seq))

    fasta_writer(fasta_folder + "/" + file_name + "_pruned.fasta", clean_entries)


def write_sampled_sequences(tree_file, infasta_folder, outfasta_folder):

    with open(tree_file) as f:
        line = f.readline().strip()
        if "(" not in line or line == ";":
            return None
        else:
            my_tree = ete3.Tree(line, format=1)
    surviving_nodes = {x.name for x in my_tree.get_leaves()}
    file_name = tree_file.split("/")[-1].split("_")[0]
    entries = fasta_reader(infasta_folder + "/" + file_name + "_whole.fasta")

    clean_entries = list()
    for h, seq in entries:
        if h[1:] in surviving_nodes:
            clean_entries.append((h, seq))

    fasta_writer(outfasta_folder + "/" + file_name + "_sampled.fasta", clean_entries)
