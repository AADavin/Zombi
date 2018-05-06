from old_functions.RatesManager import SpeciesEvolutionRates
import ete3
import numpy
import os
import random
import math

class TreeGenerator():

    def __init__(self, parameters_file):

        # First we determine the time of the massive extinction

        self.parameters = dict()
        self._read_parameters(parameters_file)

        self.RM = SpeciesEvolutionRates(self.parameters)

        self.tree_sizes = dict()
        self.lineages_counter = 0
        self.events = dict()
        self.lineages_in_time = dict()

        self.whole_species_tree = ete3.Tree() # All the species, including extinct lineages
        self._start_tree(self.whole_species_tree)

        self.extant_species_tree = ete3.Tree() # Only alive species
        self._start_tree(self.extant_species_tree)

        self.max_lineages = int(self.parameters["MAX_LINEAGES"])


    def _read_parameters(self, parameters_file):

        with open(parameters_file) as f:
            for line in f:
                if line[0] == "#":
                    continue
                if line[0] == "\n":
                    continue
                parameter, value = line.strip().split("\t")
                self.parameters[parameter] = value

        if self.parameters["SEED"] != "0":
            SEED = int(self.parameters["SEED"])
            random.seed(SEED)
            numpy.random.seed(SEED)

    def _start_tree(self, tree):

        tree.name = "Root"
        tree.add_feature("is_alive",True)
        tree.add_feature("active_gene_families", list())

        c1 = tree.add_child(dist=0)
        self.lineages_counter += 1
        c1.name = "n" + str(self.lineages_counter)
        c1.add_feature("is_alive", True)  # One indicates that the species is alive at present time
        c1.add_feature("active_gene_families", list())

        c2 = tree.add_child(dist=0)
        self.lineages_counter += 1
        c2.name = "n" + str(self.lineages_counter)
        c2.add_feature("is_alive", True)  # One indicates that the species is alive at present time
        c2.add_feature("active_gene_families", list())

        tree.is_alive = False
        self.events[0] = []
        self.events[0].append(("SP", "Root", "n1+n2"))

    def new_tree_generator(self):

        speciation, extinction = self.RM.mode_0()

        total_probability_of_event = speciation + extinction

        if total_probability_of_event >1:
            print("Have a look at your rates, they are too high")
            print("The probability of an event in each step is higher than one")

        stopping_rule = int(self.parameters["STOPPING_RULE"])
        n_lineages = int(self.parameters["N_LINEAGES"])
        total_time = int(self.parameters["TOTAL_TIME"])

        #min_lineages = int(self.parameters["MIN_LINEAGES"])
        #max_lineages = int(self.parameters["MAX_LINEAGES"])

        time_counter = 0

        success = False

        while True:

            time_counter += 1

            lineages_alive = [x for x in self.whole_species_tree.get_leaves() if x.is_alive == True]
            all_lineages = len(self.whole_species_tree.get_leaves())
            dead_lineages = all_lineages - len(lineages_alive)

            if stopping_rule == 0:
                if time_counter == total_time:
                    if time_counter not in self.events:
                        self.events[time_counter] = []
                    self.events[time_counter].append(("END", None, None))
                    success = True
                    break

            elif stopping_rule == 1:
                if all_lineages >= n_lineages:
                    if time_counter not in self.events:
                        self.events[time_counter] = []
                    self.events[time_counter].append(("END", None, None))
                    success = True
                    break

            elif stopping_rule == 2:
                if dead_lineages >= n_lineages:
                    if time_counter not in self.events:
                        self.events[time_counter] = []
                    self.events[time_counter].append(("END", None, None))
                    success = True
                    break

            elif stopping_rule == 3:
                if len(lineages_alive) >= n_lineages:
                    if time_counter not in self.events:
                        self.events[time_counter] = []
                    self.events[time_counter].append(("END", None, None))
                    success = True
                    break

            if len(lineages_alive) == 0:
                print("Simulation interrupted at time %s" % time_counter)
                print("All lineages are dead")
                if time_counter not in self.events:
                    self.events[time_counter] = []
                self.events[time_counter].append(("END", None, None))

                success = False
                return success
                break

            self.lineages_in_time[time_counter] = [x.name for x in lineages_alive]

            for lineage in lineages_alive:
                lineage.dist += 1

                if numpy.random.uniform(0, 1) <= total_probability_of_event:

                    event = self.choose_event(speciation, extinction)
                    if event == "SP":
                        self._get_speciated(lineage, time_counter)
                    elif event == "EX":
                        self._get_extinct(lineage, time_counter)

        # We add one more unit of time to avoide branches with length == 0

        self.lineages_in_time[time_counter] = [x.name for x in lineages_alive]

        lineages_alive = [x for x in self.whole_species_tree.get_leaves() if x.is_alive == True]
        for lineage in lineages_alive:
            lineage.dist += 1


        return success




    def generate_tree_mode_1(self):

        # THIS MUST BE FIXED AND CORRECTED

        # Autocorrelated speciation and extinction rates

        my_speciation, my_extinction = self.RM.mode_1()

        my_root = self.whole_species_tree.get_tree_root()
        my_root.add_feature("speciation_multiplier", 1)
        my_root.add_feature("extinction_multiplier", 1)

        c1, c2 = my_root.get_children()

        c1.add_feature("speciation_multiplier", 1.0)
        c1.add_feature("extinction_multiplier", 1.0)
        c2.add_feature("speciation_multiplier", 1.0)
        c2.add_feature("extinction_multiplier", 1.0)

        c1.speciation_multiplier, c1.extinction_multiplier = self.RM.mode_1(
            my_root.speciation_multiplier, my_root.extinction_multiplier)
        c2.speciation_multiplier, c2.extinction_multiplier = self.RM.mode_1(
            my_root.speciation_multiplier, my_root.extinction_multiplier)

        for time_counter in range(int(TOTAL_TIME / TIME_INCREASE)):

            lineages_alive = [x for x in self.whole_species_tree.get_leaves() if x.is_alive == True]

            print(time_counter, len(lineages_alive))

            self.lineages_in_time[time_counter] = [x.name for x in lineages_alive]

            for lineage in lineages_alive:
                lineage.dist += TIME_INCREASE

                # Now, the probability of extinction is modified according to the multiplier of my father

                speciation = my_speciation * lineage.up.speciation_multiplier
                extinction = my_extinction * lineage.up.extinction_multiplier

                total_probability_of_event = (speciation + extinction) * TIME_INCREASE

                if numpy.random.uniform(0, 1) <= total_probability_of_event:

                    event = self.choose_event(speciation, extinction)
                    if event == "SP":

                        self._get_speciated(lineage, time_counter, multipliers=True)
                        c1,c2 = lineage.get_children()
                        c1.speciation_multiplier, c1.extinction_multiplier = self.RM.mode_1(
                            lineage.speciation_multiplier, lineage.extinction_multiplier)
                        c2.speciation_multiplier, c2.extinction_multiplier = self.RM.mode_1(
                            lineage.speciation_multiplier, lineage.extinction_multiplier)

                    elif event == "EX":

                        self._get_extinct(lineage, time_counter)


        # To see the extinction and speciation multipliers:

        #for n in self.whole_species_tree.traverse():
        #    n.name = n.name + "_" + str(round(n.speciation_multiplier,3)) + "_" + str(round(n.extinction_multiplier,3))

    def generate_tree_mode_2(self):

        # Branch correlated speciation and extinction rates

        my_root = self.whole_species_tree.get_tree_root()

        c1, c2 = my_root.get_children()

        c1.add_feature("speciation_rate", 1.0)
        c1.add_feature("extinction_rate", 1.0)
        c2.add_feature("speciation_rate", 1.0)
        c2.add_feature("extinction_rate", 1.0)

        c1.speciation_rate, c1.extinction_rate = self.RM.mode_2()
        c2.speciation_rate, c2.extinction_rate = self.RM.mode_2()

        for time_counter in range(int(TOTAL_TIME / TIME_INCREASE)):

            lineages_alive = [x for x in self.whole_species_tree.get_leaves() if x.is_alive == True]

            print(time_counter, len(lineages_alive))

            self.lineages_in_time[time_counter] = [x.name for x in lineages_alive]

            for lineage in lineages_alive:
                lineage.dist += TIME_INCREASE

                # Now, the probability of extinction is modified according to the multiplier of my father

                total_probability_of_event = (lineage.speciation_rate + lineage.extinction_rate) * TIME_INCREASE

                if numpy.random.uniform(0, 1) <= total_probability_of_event:

                    event = self.choose_event(lineage.speciation_rate, lineage.extinction_rate)
                    if event == "SP":

                        self._get_speciated(lineage, time_counter, multipliers=True)
                        c1,c2 = lineage.get_children()
                        c1.speciation_rate, c1.extinction_rate = self.RM.mode_2()
                        c2.speciation_rate, c2.extinction_rate = self.RM.mode_2()

                    elif event == "EX":

                        self._get_extinct(lineage, time_counter)


        # To see the extinction and speciation multipliers:

        #for n in self.whole_species_tree.iter_descendants():
        #    n.name = n.name + "_" + str(round(n.speciation_rate,3)) + "_" + str(round(n.extinction_rate,3))


    def generate_tree_mode_3(self):

        # Uncorrelated speciation and extinction rates

        for time_counter in range(int(TOTAL_TIME / TIME_INCREASE)):

            lineages_alive = [x for x in self.whole_species_tree.get_leaves() if x.is_alive == True]

            print(time_counter, len(lineages_alive))

            self.lineages_in_time[time_counter] = [x.name for x in lineages_alive]

            for lineage in lineages_alive:
                lineage.dist += TIME_INCREASE

                speciation_rate, extinction_rate = self.RM.mode_3()

                total_probability_of_event = (speciation_rate + extinction_rate) * TIME_INCREASE

                if numpy.random.uniform(0, 1) <= total_probability_of_event:

                    event = self.choose_event(speciation_rate, extinction_rate)

                    if event == "SP":

                        self._get_speciated(lineage, time_counter, multipliers=True)

                    elif event == "EX":

                        self._get_extinct(lineage, time_counter)

    def generate_tree_mode_4(self):

        speciation, extinction, periods = self.RM.mode_4()

        total_probability_of_event = (speciation + extinction) * TIME_INCREASE

        for time_counter in range(int(TOTAL_TIME / TIME_INCREASE)):

            if time_counter in periods:

                speciation = periods[time_counter][0]
                extinction = periods[time_counter][1]
                total_probability_of_event = (speciation + extinction) * TIME_INCREASE

            lineages_alive = [x for x in self.whole_species_tree.get_leaves() if x.is_alive == True]

            print(time_counter, speciation, extinction, len(lineages_alive))

            self.lineages_in_time[time_counter] = [x.name for x in lineages_alive]

            for lineage in lineages_alive:
                lineage.dist += TIME_INCREASE

                if numpy.random.uniform(0, 1) <= total_probability_of_event:
                    # An event takes place

                    event = self.choose_event(speciation, extinction)

                    if event == "SP":
                        self._get_speciated(lineage, time_counter)
                    elif event == "EX":
                        self._get_extinct(lineage, time_counter)

    def _get_speciated(self, lineage, time_counter, multipliers = False):

        c1 = lineage.add_child(dist=0)
        self.lineages_counter += 1
        c1.name = "n" + str(self.lineages_counter)
        c1.add_feature("is_alive", True)  # One indicates that the species is alive at present time
        c1.add_feature("active_gene_families", list())

        c2 = lineage.add_child(dist=0)
        self.lineages_counter += 1
        c2.name = "n" + str(self.lineages_counter)
        c2.add_feature("is_alive", True)  # One indicates that the species is alive at present time
        c2.add_feature("active_gene_families", list())

        if multipliers == True:
            c1.add_feature("speciation_multiplier", 1.0)
            c1.add_feature("extinction_multiplier", 1.0)
            c2.add_feature("speciation_multiplier", 1.0)
            c2.add_feature("extinction_multiplier", 1.0)

        if time_counter not in self.events:
            self.events[time_counter] = []
        self.events[time_counter].append(("SP", lineage.name, c1.name + "+" + c2.name))  # Store the event

        lineage.is_alive = False

    def _get_extinct(self, lineage, time_counter):

        if time_counter not in self.events:
            self.events[time_counter] = []
        self.events[time_counter].append(("EX", lineage.name, None))  # Store the event
        lineage.is_alive = False
        # lineage.name = lineage.name + "E"

    def get_extant_tree(self):

        leaves_alive = [x.name for x in self.whole_species_tree.get_leaves() if x.is_alive==True]

        if len(leaves_alive) < 3:
            return 0

        self.extant_species_tree = ete3.Tree(self.whole_species_tree.write(format=1), format=1)

        ### THIS IS A BIG BOTTLE NECK

        self.extant_species_tree.prune(leaves_alive, preserve_branch_length=True)

    def get_total_ages(self):

        # We are going to round till the same precision of the time ages

        precision = int(-1 * math.log10(TIME_INCREASE))

        whole_age = round(self.whole_species_tree.get_farthest_leaf()[1], precision)
        extant_age = round(self.extant_species_tree.get_farthest_leaf()[1], precision)

        return((whole_age, extant_age))

    def store_log(self, logfolder):

        logfile = os.path.join(logfolder,"ParametersLog.tsv")
        eventsfile = os.path.join(logfolder, "SpeciesTreeEvents.tsv")
        lineagesintime = os.path.join(logfolder, "LineagesInTime.tsv")
        wholetreefile = os.path.join(logfolder, "WholeTree")
        extanttreefile = os.path.join(logfolder, "ExtantTree")


        with open(logfile, "w") as f:

            for k, v in self.parameters.items():
                line = k + "\t" + str(v) + "\n"
                f.write(line)

        with open(wholetreefile,"w") as f:
            f.write(self.whole_species_tree.write(format=1))

        if self.parameters["PRUNING"] == "1":
            self.get_extant_tree()
            with open(extanttreefile, "w") as f:
                f.write(self.extant_species_tree.write(format=1))

        with open(eventsfile, "w") as f:

            f.write("### Evolutionary Events ###\n")

            for k,v in self.events.items():
                for event,ln,cld in v:

                    if ln == None:
                        ln = "None"

                    if cld == None:
                        cld = "None"

                    line = "\t".join([str(k),event,ln,cld])+"\n"
                    f.write(line)

        with open(lineagesintime, "w") as f:

            f.write("### Lineages alive in each unit of time ###\n")

            for k,v in self.lineages_in_time.items():
                all_lineages = ";".join(v)
                line = "\t".join([str(k),all_lineages])+"\n"
                f.write(line)

    def choose_event(self, speciation, extinction):

        if numpy.random.uniform(0,1) <= (speciation /(speciation + extinction)):
            return "SP"
        else:
            return "EX"
