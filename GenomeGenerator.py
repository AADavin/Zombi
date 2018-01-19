import ete3
import os
import AuxiliarFunctions as af
from RatesManager import GeneEvolutionRates
import numpy
import random
from globals import *

class FamilyOriginator():

    def __init__(self, tree_file, events_file):

        # We need the branch length AND the time of origination of the species

        self.branch_origin = dict()

        with open(events_file) as f:
            f.readline()
            for line in f:
                time, event, node, clades = line.strip().split("\t")
                if event == "SP":
                    c1, c2 = clades.split("+")
                    self.branch_origin[c1] = int(time)
                    self.branch_origin[c2] = int(time)

        with open(tree_file) as f:
            self.whole_tree = ete3.Tree(f.readline().strip(),format=1)

        self.branch_length = dict()

        for node in self.whole_tree.iter_descendants():

            self.branch_length[node.name] = int(node.dist / TIME_INCREASE)

        self.vector_names = [x for x in self.branch_length.keys()]
        self.vector_lengths = [self.branch_length[x] for x in self.vector_names]



    def create_families(self, stem_length =  0, families_in_stem = 0):

        # We select first a branch

        branch = numpy.random.choice(self.vector_names, 1, p=af.normalize(self.vector_lengths))[0]
        while self.branch_length[branch] == 1:
            branch = numpy.random.choice(self.vector_names, 1, p=af.normalize(self.vector_lengths))[0]

        # We select from an uniform distribution a position for the branch

        time_in_branch = numpy.random.randint(1, self.branch_length[branch]) # Families cannot

        # We give the absolute time

        return branch, self.branch_length[branch] + self.branch_origin[branch]


class GeneFamilySimulator():

    def __init__(self, parameters_file, events_file, lineages_in_time_file):

        self.parameters = dict()

        with open(parameters_file) as f:
            for line in f:
                parameter, value = line.strip().split("\t")
                self.parameters[parameter] = value

        self.RM = GeneEvolutionRates(self.parameters)
        self.tree_events = dict()
        self.rates = self.RM.mode_0()


        with open(events_file) as f:
            f.readline()
            for line in f:
                time, dn, ln, cld = line.strip().split("\t")
                if int(time) not in self.tree_events:
                    self.tree_events[int(time)] = list()
                self.tree_events[int(time)].append((dn,ln,cld))

        self.lineages_in_time = dict()

        with open(lineages_in_time_file) as f:
            f.readline()
            for line in f:
                k,v = line.strip().split("\t")
                self.lineages_in_time[int(k)] = v.split(";")

    def choose_event(self, duplication, transfer, loss):
        draw = numpy.random.choice(["D","T","L"], 1, p=af.normalize([duplication, transfer, loss]))
        return draw

    def choose_recipient(self, time_counter, donor, strategy):

        possible_recipients = [x for x in self.lineages_in_time[time_counter] if x != donor]
        if len(possible_recipients) > 1:  # The transfer can go to any other leave, but it cannot be transfer
            recipient = random.choice(possible_recipients)
            return recipient
        else:
            return None

    def evolve_gene_family(self, duplication, transfer, loss, time_counter):

        total_probability_of_event = (duplication + transfer + loss) * TIME_INCREASE

        active_lineages_gt = [x for x in self.gene_family["Gene_tree"].get_leaves() if
                              x.is_alive == True]

        for g_node in active_lineages_gt:

            if g_node.dist == 0:
                continue

            elif numpy.random.uniform(0, 1) <= total_probability_of_event:  # An event takes place

                event = self.choose_event(duplication, transfer, loss)

                if event == "D":
                    self.get_duplicated(g_node, time_counter)

                elif event == "T":

                    recipient = self.choose_recipient(time_counter, g_node.current_branch, 0)
                    if recipient != None:
                        self.get_transferred(g_node, recipient, time_counter)

                elif event == "L":
                    self.get_lost(g_node, time_counter)

    def run_mode_0(self):

        # If rates are global we don't have to compute new each time rates
        # We resort to the global rates already computed

        duplication, transfer, loss = self.rates

        origin_time = self.gene_family["Origin_time"]

        for time_counter in range(origin_time, int(TOTAL_TIME / TIME_INCREASE)):

            self.increase_distances()

            active_branches = {x.current_branch for x in self.gene_family["Gene_tree"].get_leaves()}

            if time_counter in self.tree_events:

                for event, snode, children in self.tree_events[time_counter]:

                    if snode in active_branches:

                        if event == "EX":
                            self.get_extinct(snode)

                        elif event == "SP":
                            sc1, sc2 = children.split("+")
                            self.gene_tree_speciation(snode, sc1, sc2, time_counter)

            self.evolve_gene_family(duplication, transfer, loss, time_counter)

    def run_mode_1(self):

        # Rates are family wise. This is control by the rate manager
        # I must keep track of the rates computed with a new type of file

        duplication, transfer, loss = self.RM.mode_0()

        origin_time = self.gene_family["Origin_time"]

        for time_counter in range(origin_time, int(TOTAL_TIME / TIME_INCREASE)):

            self.increase_distances()

            active_branches = {x.current_branch for x in self.gene_family["Gene_tree"].get_leaves()}

            if time_counter in self.tree_events:

                for event, snode, children in self.tree_events[time_counter]:

                    if snode in active_branches:

                        if event == "EX":
                            self.get_extinct(snode)

                        elif event == "SP":
                            sc1, sc2 = children.split("+")
                            self.gene_tree_speciation(snode, sc1, sc2, time_counter)

            self.evolve_gene_family(duplication, transfer, loss, time_counter)

    def run_mode_3(self):

        # Rates are lineage wise. This is control by the rate manager
        pass



    def increase_distances(self):

        active_lineages_gt = [x for x in self.gene_family["Gene_tree"].get_leaves() if
                                  x.is_alive == True]

        for node in active_lineages_gt:
            node.dist += TIME_INCREASE

    def origination(self, branch, time_counter, name):

        self.gene_family = dict()
        self.gene_family["Name"] = name
        self.gene_family["Origin_time"] = time_counter
        self.gene_family["Events"] = list()
        self.gene_family["Gene_tree"] = ete3.Tree()
        self.gene_family["Gene_tree"].add_feature("is_alive", True)
        self.gene_family["Gene_tree"].add_feature("current_branch", branch)
        self.gene_family["Gene_tree_extant"] = ete3.Tree()
        self.gene_family["Events"].append(("Origination",time_counter * TIME_INCREASE, branch))
        self.gene_family["is_alive"] = True # To avoid too large families
        self.gene_family["Profile"] = dict()

        if time_counter == 0:
            self.gene_family["Profile"]["Root"] = 1

        # To store the events

        self.gene_family["Duplications"] = dict()
        self.gene_family["LeavingTransfers"] = dict()
        self.gene_family["ArrivingTransfers"] = dict()
        self.gene_family["Losses"] = dict()

    def get_duplicated(self, leaf, time_counter):

        self.gene_family["Events"].append(
            ("Duplication", time_counter * TIME_INCREASE, leaf.current_branch))

        if leaf.current_branch not in self.gene_family["Duplications"]:
            self.gene_family["Duplications"][leaf.current_branch] = 0

        self.gene_family["Duplications"][leaf.current_branch] +=1

        gc1 = leaf.add_child(dist=0)
        gc1.add_feature("is_alive", True)
        gc1.add_feature("current_branch", leaf.current_branch)

        gc2 = leaf.add_child(dist=0)
        gc2.add_feature("is_alive", True)
        gc2.add_feature("current_branch", leaf.current_branch)

        leaf.name = "D"
        leaf.is_alive = False

    def get_lost(self, leaf, time_counter):

        self.gene_family["Events"].append(
            ("Loss", time_counter * TIME_INCREASE, leaf.current_branch))

        if leaf.current_branch not in self.gene_family["Losses"]:
            self.gene_family["Losses"][leaf.current_branch] = 0

        self.gene_family["Losses"][leaf.current_branch] +=1

        leaf.is_alive = False
        leaf.name = "L"

    def get_transferred(self, leaf, recipient, time_counter):

        rt = float(self.parameters["REPLACEMENT_T"])

        self.gene_family["Events"].append(
            ("LeavingTransfer", time_counter * TIME_INCREASE, leaf.current_branch))
        self.gene_family["Events"].append(
            ("ArrivingTransfer", time_counter * TIME_INCREASE, recipient))

        if leaf.current_branch not in self.gene_family["LeavingTransfers"]:
            self.gene_family["LeavingTransfers"][leaf.current_branch] = 0

        if leaf.current_branch not in self.gene_family["ArrivingTransfers"]:
            self.gene_family["ArrivingTransfers"][leaf.current_branch] = 0

        self.gene_family["LeavingTransfers"][leaf.current_branch] += 1
        self.gene_family["ArrivingTransfers"][leaf.current_branch] += 1

        leaf.name = "TRANSFER"

        gc1 = leaf.add_child(dist=0)
        gc1.add_feature("is_alive", True)
        gc1.add_feature("current_branch", recipient)
        gc1.name = recipient

        gc2 = leaf.add_child(dist=0)
        gc2.add_feature("is_alive", True)
        gc2.add_feature("current_branch", leaf.current_branch)
        gc2.name = leaf.current_branch

        leaf.is_alive = False

        other_copies = [x for x in self.gene_family["Gene_tree"].get_leaves() if x.current_branch == recipient and x != 0 ]

        if len(other_copies) != 0 and numpy.random.uniform(0,1) < rt:

            # So far, this is like a normal transfer, but now I have to kill one of the extant genes in that branch
            replaced_lineage = random.choice(other_copies)
            replaced_lineage.is_alive = False


    def get_extinct(self, sp_branch):

        g_leaves = self.gene_family["Gene_tree"].get_leaves()

        self.gene_family["Profile"][sp_branch] = 0

        for g_leaf in g_leaves:

            if g_leaf.current_branch == sp_branch and g_leaf.is_alive == True:

                self.gene_family["Profile"][sp_branch] += 1

                g_leaf.is_alive = False

    def gene_tree_speciation(self, sp, c1, c2):

        g_leaves = self.gene_family["Gene_tree"].get_leaves()

        self.gene_family["Profile"][sp] = 0

        for g_leaf in g_leaves:

            if g_leaf.current_branch == sp and g_leaf.is_alive == True:

                self.gene_family["Profile"][sp] += 1

                gc1 = g_leaf.add_child(dist=0)
                gc1.add_feature("is_alive", True)
                gc1.add_feature("current_branch", c1)

                gc2 = g_leaf.add_child(dist=0)
                gc2.add_feature("is_alive", True)
                gc2.add_feature("current_branch", c2)

                g_leaf.name = "SP" + sp
                g_leaf.is_alive = False

    def complete_gene_family_information(self):

        # Add root length

        myroot = self.gene_family["Gene_tree"].get_tree_root()
        one_leaf, distance_to_present = myroot.get_farthest_leaf()
        myroot.dist = (TOTAL_TIME - distance_to_present) - self.gene_family["Origin_time"]

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

        return self.gene_family["Gene_tree"].write(format=1)

    def output_profile(self, profile, which):

        with open(profile, "r+") as f:
            header = f.readline().strip().split("\t")[1:]
            myline = list()
            myline.append(str(self.gene_family["Name"]))

            for name in header:

                try: myline.append(str(self.gene_family[which][name]))
                except: myline.append(str(0))

            myline = "\t".join(myline) +"\n"
            f.write(myline)


