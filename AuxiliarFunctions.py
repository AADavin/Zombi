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

        if parameter == "SEQUENCE_SIZE":
            parameters[parameter] = int(value)

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
                or parameter == "TOTAL_LINEAGES" or parameter == "STOPPING_RULE" or parameter == "MAX_LINEAGES":
            parameters[parameter] = int(value)

    return parameters


def prepare_genome_parameters(parameters):

    for parameter, value in parameters.items():

        if parameter == "P_ESSENTIAL_GENE":
            parameters[parameter] = obtain_value(value)

        if parameter == "DUPLICATION_EXTENSION" or parameter == "TRANSFER_EXTENSION" \
                or parameter == "LOSS_EXTENSION" or parameter == "INVERSION_EXTENSION" or \
                parameter == "TRANSLOCATION_EXTENSION" or parameter == "ORIGINATION_EXTENSION":

            parameters[parameter] = obtain_value(value)

        if parameter == "ROOT_GENOME":
            parameters[parameter] = value.split(";")

        if parameter == "REPLACEMENT_TRANSFER":
            parameters[parameter] = float(value)

        if parameter == "PROFILES" or parameter == "EVENTS_PER_BRANCH" or parameter == "GENE_TREES" \
                or parameter == "PRUNE_TREES" or parameter == "TRANSFER_PREFERENCE":

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

    draw = numpy.random.choice(possible_recipients, 1, p=af.normalize(weights))



