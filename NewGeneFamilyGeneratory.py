import ete3
import os
import numpy
import random
import AuxiliarFunctions as af
from RatesManager import GeneEvolutionRates
import copy


class UnorderedGenome():

    def __init__(self):

        self.genes = list()

    def start_genome(self, gf):

        name = "_".join(["Root", "+", str(gf), "1"])
        self.genes.append(name)

    def select_random_gene(self):

        return numpy.randint(len(self.genes))

    def duplicate_gene(self, species_node, homologous, time, mychoice):

        cb, sense, gf, cp = self.genes[mychoice].split("_")

        homologous[gf]["Copies"] += 1
        name1 = cb + "_" + sense + "_" + gf + "_" + str(homologous[gf]["Copies"])
        homologous[gf]["Copies"] += 1
        name2 = cb + "_" + sense + "_" + gf + "_" + str(homologous[gf]["Copies"])

        self.genes.append(name1)
        self.genes.append(name2)

        homologous[gf]["Events"].append(
            ("D", time, species_node + ";" + cp + ";" + name1.split("_")[3] + ";" + name2.split("_")[3]))


    def loss_gene(self, species_node, homologous, time, mychoice):

        cb, sense, gf, cp = self.genes[mychoice].split("_")
        homologous[gf]["Events"].append(("L", time, species_node + ";" + cp))
        self.genes.remove(self.genes[mychoice])

    def insert_gene(self, species_node, homologous, time, mygene):

        self.genes.append(mygene)

    def update_homologous(self, species_node, homologous):

        for i,gene in enumerate(self.genes):
            cb,sense,gf,id = gene.split("_")
            homologous[gf]["Copies"] += 1
            new_id = homologous[gf]["Copies"]
            self.genes[i] = "_".join((species_node,sense,gf,str(new_id)))

        return self.genes

    def write_genome(self, genome_file):

        with open(genome_file, "w") as f:

            f.write("\t".join(("Position", "Gene_family","Orientation","Id"))+"\n")

            for i,gene in enumerate(self.genes):

                cb, sense, gf, id = gene.split("_")
                line = "\t".join((str(i), gf,sense,id)) + "\n"
                f.write(line)

class GeneFamilySimulator():

    def __init__(self, parameters_file, events_file, lineages_in_time_file):

        self.homologous = dict()
        self.parameters = dict()
        self.species_tree = ete3.Tree()

        self.read_parameters(parameters_file)

        if self.parameters["SEED"] != "0":
            SEED = int(self.parameters["SEED"])
            random.seed(SEED)
            numpy.random.seed(SEED)

        self.gf_number = int(self.parameters["STEM_FAMILIES"])

        self.rm = GeneEvolutionRates(self.parameters)
        self.tree_events = dict()
        self._start_tree()

        with open(events_file) as f:
            f.readline()
            for line in f:
                time, dn, ln, cld = line.strip().split("\t")
                if dn == "END":
                    self.total_time = int(time)
                    continue
                elif int(time) not in self.tree_events:
                    self.tree_events[int(time)] = list()
                self.tree_events[int(time)].append((dn,ln,cld))

        self.lineages_in_time = dict()

        with open(lineages_in_time_file) as f:
            f.readline()
            for line in f:
                k,v = line.strip().split("\t")
                self.lineages_in_time[int(k)] = v.split(";")

    def read_parameters(self, parameters_file):

        with open(parameters_file) as f:
            for line in f:
                if line[0] == "#":
                    continue
                parameter, value = line.strip().split("\t")
                self.parameters[parameter] = value

    def _start_tree(self):

        gnm = UnorderedGenome()
        gnm.start_genome(1)
        root = self.species_tree.get_tree_root()
        root.name = "Root"
        root.add_feature("Genome",gnm)
        root.add_feature("Is_alive", True)

        for i in range(1, int(self.parameters["STEM_FAMILIES"])+1):
            self.homologous[str(i)] = dict()
            self.homologous[str(i)]["Copies"] = 1
            self.homologous[str(i)]["Events"] = list()
            self.homologous[str(i)]["Events"].append(("O", str(0), "Root" + ";" + "1"))

    def choose_event(self, duplication, transfer, loss):

        draw = numpy.random.choice(["D", "T", "L"], 1, p=af.normalize([duplication, transfer, loss]))
        return draw

    def choose_recipient(self, time_counter, donor):

        possible_recipients = [x for x in self.lineages_in_time[time_counter] if x != donor]
        if len(possible_recipients) > 1:
            recipient = random.choice(possible_recipients)
            return recipient
        else:
            return None

    def evolve_genomes(self, duplication, transfer, loss, time_counter):

        total_probability_of_event = duplication + transfer + loss

        active_genomes = [x for x in self.species_tree.get_leaves() if x.is_alive == True]

        random.shuffle(active_genomes)

        for node in active_genomes:

            genome = node.Genome
            snode = node.name

            mygenes = list(genome.genes)
            random.shuffle(mygenes)

            for i,gene in mygenes:

                if numpy.random.uniform(0, 1) <= total_probability_of_event:  # An event takes place

                    event = self.choose_event(duplication, transfer, loss)

                    if event == "D":

                        genome.duplicate_gene(snode, homologous, time_counter, gene)

                    elif event == "T":
                        recipient = self.choose_recipient(time_counter, node.name, 0)
                        if recipient == None:
                            continue

                        a = genome.select_random

                        old_segment = list()
                        new_segment = list()

                        for i in a:

                            cb, sense, gf, id = genome.genes[i].split("_")

                            self.homologous[gf]["Copies"] += 1
                            name1  = snode + "_" + sense + "_" + gf + "_" + str(self.homologous[gf]["Copies"])
                            old_segment.append(name1)

                            self.homologous[gf]["Copies"] += 1
                            name2 = recipient + "_" + sense + "_" + gf + "_" + str(self.homologous[gf]["Copies"])
                            new_segment.append(name2)

                            self.homologous[gf]["Events"].append(
                                ("T", str(time_counter), snode + ";" + id + ";" + name1.split("_")[3] + ";" + recipient + ";" + name2.split("_")[3]))

                        # Now I have prepared the two segments. First I am going to update de donor segment

                        elements_to_remove = [genome.genes[x] for x in a]

                        for element in elements_to_remove:
                            genome.genes.remove(element)

                        position = a[0]

                        for i, x in enumerate(old_segment):
                            genome.genes.insert(position + i + 1, x)

                        # Then I update the receptor segment

                            my_recipient = self.species_tree&recipient
                            recipient_genome = my_recipient.Genome
                            p = recipient_genome.select_random_position()
                            recipient_genome.insert_segment(p, new_segment)

                    elif event == "L":

                        continue

                        # We have to check that the minimal size has not been attained

    def run(self, genome_folder):

        duplication, transfer, loss = self.rm.mode_0()

        for time_counter in range(int(self.total_time)):

            if time_counter % 100 == 0:
                print("Simulating gene family evolution %s".format() % str('{:.1%}'.format(time_counter/int(self.total_time))))

            if time_counter in self.tree_events:

                for event, snode, children in self.tree_events[time_counter]:

                    if event == "EX":
                        self.get_extinct(time_counter, snode)

                    elif event == "SP":

                        # First we write the ancestral genome

                        mynode = self.species_tree&snode

                        sc1, sc2 = children.split("+")
                        self.get_speciated(time_counter, snode, sc1, sc2)

            self.evolve_genomes(duplication, transfer, loss, time_counter)
            self.increase_distances()


    def increase_distances(self):

        active_lineages = [x for x in self.species_tree.get_leaves() if x.is_alive == True]

        for node in active_lineages:
            node.dist += 1

    def get_extinct(self, time, sp):

        sp = self.species_tree&sp
        sp.is_alive = False
        parent_genome = sp.Genome

        for i, gene in enumerate(parent_genome.genes):
            cb, sense, gf, id = gene.split("_")
            self.homologous[gf]["Events"].append(("E", time, id))

    def get_speciated(self, time, sp, c1, c2):

        sp = self.species_tree&sp
        parent_genome = sp.Genome

        sc1 = sp.add_child(dist=0)
        sc1.name = c1
        sc1.add_feature("is_alive", True)
        sc1.add_feature("Genome", copy.deepcopy(parent_genome))
        genes_affected_1 = sc1.Genome.update_homologous(sc1.name, self.homologous)

        sc2 = sp.add_child(dist=0)
        sc2.name = c2
        sc2.add_feature("is_alive", True)
        sc2.add_feature("Genome", copy.deepcopy(parent_genome))
        genes_affected_2 = sc2.Genome.update_homologous(sc2.name, self.homologous)

        for i, gene in enumerate(parent_genome.genes):

            cb, sense, gf, id = gene.split("_")
            speciation_event = ";".join(
                (sp.name, id, sc1.name, genes_affected_1[i].split("_")[3], sc2.name, genes_affected_2[i].split("_")[3]))
            self.homologous[gf]["Events"].append(("S", time, speciation_event))

        sp.is_alive = False

    def get_gene_family_tree(self):

        if len(self.gene_family["Gene_tree"].get_leaves()) < 3:
            return "None"
        else:
            return self.gene_family["Gene_tree"].write(format=1)

    def output_profile(self, profile, which):

        with open(profile, "r+") as f:
            header = f.readline().strip().split("\t")[1:]
            myline = list()
            myline.append(str(self.gene_family["Name"]))

            for name in header:

                myline.append(str(self.gene_family[which][name]))

            myline = "\t".join(myline) +"\n"
            f.write(myline)

    def write_transfers(self, transfers_file):

        with open(transfers_file, "a") as f:
            for dn, rc in self.gene_family["Transfers"]:
                line = self.gene_family["Name"] + "\t" + dn + "\t" + rc + "\n"
                f.write(line)

    def write_log(self, events_per_family_folder, gene_trees_folder):

        # we write the events per gene family

        for gf in self.homologous:

            events_file = os.path.join(events_per_family_folder, self.parameters["PREFIX"] + str(gf) + "_events.tsv")

            with open(events_file, "w") as f:

                f.write("Event\tTime\tNode\n")

                for time,event,node in self.homologous[gf]["Events"]:

                    line = "\t".join(map(str,[time, event, node])) + "\n"
                    f.write(line)

        if self.parameters["OUTPUT_GENETREES"] == "1":

            events_files = os.listdir(events_per_family_folder)
            for event_file in events_files:
                genetree_file = os.path.join(gene_trees_folder, event_file.replace("_events.tsv","_genetree.txt"))
                self.write_gene_tree(os.path.join(events_per_family_folder, event_file), genetree_file)


    def write_gene_tree(self, events_file, genetree_file):

        print("Computing gene tree %s" % genetree_file.split("/")[-1])

        events = dict()

        with open(events_file) as f:
            f.readline()

            for line in f:
                event, time, nodes = line.strip().split("\t")
                time = int(time)
                if time not in events:
                    events[time] = list()
                events[time].append((event, nodes))

                if event == "O":
                    origin_time = time
                    origin_node = nodes.split(";")[0]

        gene_tree = ete3.Tree()
        root = gene_tree.get_tree_root()
        root.name = "1"
        root.add_feature("is_alive", True)
        root.add_feature("current_branch", origin_node)

        for time in range(origin_time, self.total_time):

            if time in events:

                for event, nodes in events[time]:

                    if event == "S":

                        sp,parent,spc1,c1,spc2,c2 = nodes.split(";")

                        n = gene_tree&parent
                        n.is_alive = False
                        n.current_branch = sp

                        gc1 = n.add_child()
                        gc1.name = c1
                        gc1.add_feature("is_alive", True)
                        gc1.add_feature("current_branch", spc1)
                        gc1.dist = 0

                        gc2 = n.add_child()
                        gc2.name = c2
                        gc2.add_feature("is_alive", True)
                        gc2.add_feature("current_branch", spc2)
                        gc2.dist = 0

                    elif event == "E":

                        n = gene_tree&nodes
                        n.is_alive = False

                    elif event == "D":

                        sp, parent, c1, c2 = nodes.split(";")

                        n = gene_tree & parent
                        n.is_alive = False

                        gc1 = n.add_child()
                        gc1.name = c1
                        gc1.add_feature("is_alive", True)
                        gc1.add_feature("current_branch", sp)

                        gc2 = n.add_child()
                        gc2.name = c2
                        gc2.add_feature("is_alive", True)
                        gc2.add_feature("current_branch", sp)

                    elif event == "T":

                        sp, parent, c1, recipient, c2 = nodes.split(";")

                        n = gene_tree & parent
                        n.is_alive = False

                        gc1 = n.add_child()
                        gc1.name = c1
                        gc1.add_feature("is_alive", True)
                        gc1.add_feature("current_branch", sp)

                        gc2 = n.add_child()
                        gc2.name = c2
                        gc2.add_feature("is_alive", True)
                        gc2.add_feature("current_branch", recipient)

                    elif event == "L":

                        mynode = nodes.split(";")[1]
                        n = gene_tree&mynode
                        n.is_alive = False

            active_lineages = [x for x in gene_tree.get_leaves() if x.is_alive == True]

            for node in active_lineages:
                node.dist+=1

        # Now we add the leaves name

        for n in gene_tree.get_leaves():
            n.name = n.current_branch + "_" + n.name

        with open(genetree_file.replace("_genetree","_wholegenetree"), "w") as f:
            f.write(gene_tree.write(format=1))

        if self.parameters["PRUNE_GENETREES"] == "1":
            print("Pruning %s" % genetree_file.split("/")[-1])
            active_lineages = [x for x in gene_tree.get_leaves() if x.is_alive == True]

            if len(active_lineages) >= 3:

                gene_tree.prune(active_lineages, preserve_branch_length=True)
                with open(genetree_file.replace("_genetree", "_extantgenetree"), "w") as f:
                    f.write(gene_tree.write(format=1))


def testing():

    params = "/Users/adriandavin/Desktop/Bioinformatics/SimuLyon/Cedric/GeneFamilyParameters.tsv"
    events_file = "/Users/adriandavin/Desktop/Bioinformatics/SimuLyon/Cedric/Example1/SpeciesTreeEvents.tsv"
    lineages_intime = "/Users/adriandavin/Desktop/Bioinformatics/SimuLyon/Cedric/Example1/LineagesInTime.tsv"

    gfs = GeneFamilySimulator(params, events_file, lineages_intime)
    gfs.run("/Users/adriandavin/Desktop/Bioinformatics/SimuLyon/Cedric/TEST")
    gfs.write_log("/Users/adriandavin/Desktop/Bioinformatics/SimuLyon/Cedric/TEST", "/Users/adriandavin/Desktop/Bioinformatics/SimuLyon/Cedric/TEST")


testing()