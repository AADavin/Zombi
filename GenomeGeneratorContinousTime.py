import ete3
import os
import numpy
import random
import AuxiliarFunctions as af
from RatesManager import GenomeEvolutionRates
from Genome import Genome
import copy

class GenomeSimulator():

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

        self.rm = GenomeEvolutionRates(self.parameters)
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

        gnm = Genome()
        gnm.start_genome(int(self.parameters["STEM_FAMILIES"]))
        root = self.species_tree.get_tree_root()
        root.name = "Root"
        root.add_feature("Genome",gnm)
        root.add_feature("Is_alive", True)

        for i in range(1, int(self.parameters["STEM_FAMILIES"])+1):
            self.homologous[str(i)] = dict()
            self.homologous[str(i)]["Copies"] = 1
            self.homologous[str(i)]["Events"] = list()
            self.homologous[str(i)]["Events"].append(("O", str(0), "Root" + ";" + "1"))

    def choose_event(self, duplication, transfer, loss, inversion, translocation, origination):

        draw = numpy.random.choice(["D", "T", "L", "I", "C", "O"], 1, p=af.normalize([duplication, transfer, loss, inversion, translocation, origination]))

        return draw

    def choose_recipient(self, time_counter, donor, strategy):

        possible_recipients = [x for x in self.lineages_in_time[time_counter] if x != donor]
        if len(possible_recipients) > 1:
            recipient = random.choice(possible_recipients)
            return recipient
        else:
            return None

    def get_time_to_next_event(self, n, d, t, l, i, c, o):

        total = 0.0
        for j in range(n):
            total += d + t + l + i + c + o
        time = numpy.random.exponential(1/total)
        return time

    def evolve_genomes(self, duplication, transfer, loss, inversion, translocation, origination, time_counter):

        total_probability_of_event = duplication + transfer + loss + inversion + translocation + origination

        d_e = float(self.parameters["DUPLICATION_E"])
        t_e = float(self.parameters["TRANSFER_E"])
        l_e = float(self.parameters["LOSS_E"])
        i_e = float(self.parameters["INVERSION_E"])
        c_e = float(self.parameters["TRANSLOCATION_E"])

        active_genomes = [x for x in self.species_tree.get_leaves() if x.is_alive == True]
        random.shuffle(active_genomes)

        for node in active_genomes:

            genome = node.Genome
            snode = node.name

            if numpy.random.uniform(0, 1) <= total_probability_of_event:  # An event takes place

                event = self.choose_event(duplication, transfer, loss, inversion, translocation, origination)

                if event == "D":

                    a = genome.obtain_affected_genes(d_e)
                    genome.duplicate_segment(snode, self.homologous, time_counter, a)

                elif event == "T":
                    recipient = self.choose_recipient(time_counter, node.name, 0)
                    if recipient == None:
                        continue

                    a = genome.obtain_affected_genes(t_e)

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

                    a = genome.obtain_affected_genes(l_e)

                    # We have to check that the minimal size has not been attained

                    if (len(genome.genes) - len(a)) <= int(self.parameters["MIN_GENOME_SIZE"]):
                        pass
                    else:
                        genome.loss_segment(snode, self.homologous, time_counter, a)

                elif event == "I":


                    a = genome.obtain_affected_genes(i_e)
                    genome.invert_segment(snode, self.homologous,time_counter, a)

                elif event == "C":
                    a = genome.obtain_affected_genes(c_e)
                    genome.translocate_segment(snode, self.homologous, time_counter, a)

                elif event == "O":

                    self.gf_number += 1

                    if numpy.random.randint(2) == 0:
                        sense = "+"
                    else:
                        sense = "-"

                    mygf = str(self.gf_number)
                    self.homologous[mygf] = dict()
                    self.homologous[mygf]["Copies"] = 1
                    self.homologous[mygf]["Events"] = list()
                    self.homologous[mygf]["Events"].append(("O", str(time_counter),  snode + ";" + "1"))

                    segment = ["_".join((snode,sense,mygf,"1"))]

                    p = genome.select_random_position()
                    genome.insert_segment(p, segment)

    def run(self, genome_folder):

        next_event_in_species_tree = 0

        duplication, transfer, loss, inversion, translocation, origination = self.rm.mode_0()

        finished = False

        while finished == False

            time = self.get_time_to_next_event()

            if time >= next_event_in_species_tree:
                print("Do the event in the species tree")

                if event == "EX":
                    self.get_extinct(time_counter, snode)

                elif event == "SP":

                    # First we write the ancestral genome

                    mynode = self.species_tree & snode
                    mynode.Genome.write_genome(os.path.join(genome_folder, snode + "_GENOME.tsv"))

                    sc1, sc2 = children.split("+")
                    self.get_speciated(time_counter, snode, sc1, sc2)
            else:

                self.evolve_genomes(duplication, transfer, loss, inversion, translocation, origination, time_counter)
                self.increase_distances(time)



    def increase_distances(self, time):
        active_lineages = [x for x in self.species_tree.get_leaves() if x.is_alive == True]
        for node in active_lineages:
            node.dist += time

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

    def complete_gene_family_information(self):

        # Add root length

        myroot = self.gene_family["Gene_tree"].get_tree_root()
        one_leaf, distance_to_present = myroot.get_farthest_leaf()
        myroot.dist = (self.total_time - distance_to_present) - self.gene_family["Origin_time"]

        # Add names

        genetree = self.gene_family["Gene_tree"]
        names = dict()

        for leaf in genetree.get_leaves():

            leaf.name = leaf.current_branch
            myname = leaf.name.split("_")[0]

            if myname in names:
                names[myname] += 1
            else:
                names[myname] = 1

            if leaf.is_alive:
                leaf.name = myname + "_" + "A" + "_" + str(names[leaf.name])
            else:
                leaf.name = myname + "_" + "E" + "_" + str(names[leaf.name])

        for leaf in genetree.get_leaves():
            if not leaf.is_alive:
                continue
            if leaf.current_branch not in self.gene_family["Profile"]:
                self.gene_family["Profile"][leaf.current_branch] = 0
            self.gene_family["Profile"][leaf.current_branch] += 1

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

    def write_log(self, gene_trees_folder, events_per_family_folder, complete_genomes_folder, events_per_branch_folder):

        # First we write the genomes

        for n in self.species_tree.get_leaves():
            print("Writing %s genome" % n.name)
            n.Genome.write_genome(os.path.join(complete_genomes_folder, n.name + "_GENOME.tsv"))

        # Second, we write the events per gene family

        for gf in self.homologous:

            events_file = os.path.join(events_per_family_folder, self.parameters["PREFIX"] + str(gf) + "_events.tsv")

            with open(events_file, "w") as f:

                f.write("Event\tTime\tNode\n")

                for time,event,node in self.homologous[gf]["Events"]:

                    line = "\t".join(map(str,[time, event, node])) + "\n"
                    f.write(line)

        # Thirds, we write the gene trees

        if self.parameters["OUTPUT_GENETREES"] == "1":

            events_files = os.listdir(events_per_family_folder)

            for event_file in events_files:

                genetree_file = os.path.join(gene_trees_folder, event_file.replace("_events.tsv","_genetree.txt"))

                self.write_gene_tree(os.path.join(events_per_family_folder, event_file), genetree_file)

        # Fourth, we write the events per branch

        self.write_events_per_branch(events_per_family_folder, events_per_branch_folder)

    def write_events_per_branch(self, events_per_family_folder, events_per_branch_folder):

        fam_events = [os.path.join(events_per_family_folder, x) for x in os.listdir(events_per_family_folder) if
                      "_events.tsv" in x]

        all_events = dict()

        for fam_event in fam_events:

            fam_name = fam_event.split("/")[-1].split("_")[0].replace(self.parameters["PREFIX"],"")

            with open(fam_event) as f:

                f.readline()

                for line in f:

                    event, time, nodes = line.strip().split("\t")

                    if event == "S" or event == "E":
                        continue

                    elif event == "T":

                        leaving_branch, node_gt, remaining_node_gt, arriving_branch, arriving_gt = nodes.split(";")

                        if leaving_branch not in all_events:
                            all_events[leaving_branch] = dict()
                        if int(time) not in all_events[leaving_branch]:
                            all_events[leaving_branch][int(time)] = list()

                        all_events[leaving_branch][int(time)].append(("LT", nodes, fam_name))

                        if arriving_branch not in all_events:
                            all_events[arriving_branch] = dict()
                        if int(time) not in all_events[arriving_branch]:
                            all_events[arriving_branch][int(time)] = list()

                        all_events[arriving_branch][int(time)].append(("AT", nodes, fam_name))

                    else:
                        branch = nodes.split(";")[0]
                        if branch not in all_events:
                            all_events[branch] = dict()
                        if int(time) not in all_events[branch]:
                            all_events[branch][int(time)] = list()
                        all_events[branch][int(time)].append((event, nodes, fam_name))

        for branch in all_events:

            print("Writing events per branch %s" % branch)

            with open(os.path.join(events_per_branch_folder, branch + "_branchevents.tsv"), "w") as f:

                f.write("\t".join(("Time","Event","Nodes", "Gene_family")) + "\n")

                for i in range(self.total_time + 1):

                    if i in all_events[branch]:

                        for event, nodes, fam_name in all_events[branch][i]:

                            if event == "LT" or event == "AT":
                                mynodes = nodes

                            else:
                                mynodes = ";".join(nodes.split(";")[1:])
                            f.write("\t".join((str(i), event, mynodes, fam_name))+"\n")


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

