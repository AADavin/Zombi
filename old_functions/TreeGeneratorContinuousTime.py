from old_functions.RatesManager import SpeciesEvolutionRates
import ete3
import numpy
import os
import random


class TreeGeneratorContinuousTime():

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

        c1 = tree.add_child(dist=0)
        self.lineages_counter += 1
        c1.name = "n" + str(self.lineages_counter)
        c1.add_feature("is_alive", True)  # One indicates that the species is alive at present time

        c2 = tree.add_child(dist=0)
        self.lineages_counter += 1
        c2.name = "n" + str(self.lineages_counter)
        c2.add_feature("is_alive", True)  # One indicates that the species is alive at present time

        tree.is_alive = False
        self.events[0] = []
        self.events[0].append(("SP", "Root", "n1+n2"))

    def new_tree_generator(self):

        speciation, extinction = self.RM.mode_0()

        time_counter = 0

        success = False

        while True:

            lineages_alive = [x for x in self.whole_species_tree.get_leaves() if x.is_alive == True]
            success = self.check_if_stop(lineages_alive, time_counter)

            if success == True:
                break

            time_to_next_event = self.get_time_to_next_event(len(lineages_alive), speciation, extinction)
            time_counter += time_to_next_event

            self.increase_distances(time_to_next_event, lineages_alive)

            event = self.choose_event(speciation, extinction)
            lineage = random.choice(lineages_alive)

            if event == "SP":
                self._get_speciated(lineage, time_counter)
            elif event == "EX":
                self._get_extinct(lineage, time_counter)

        # We add one more unit of time to avoid branches with length == 0

        lineages_alive = [x for x in self.whole_species_tree.get_leaves() if x.is_alive == True]
        for lineage in lineages_alive:
            lineage.dist += 1

        return success

    def increase_distances(self, time, lineages_alive):

        for lineage in lineages_alive:
            lineage.dist += time

    def check_if_stop(self, lineages_alive, time_counter):

        stopping_rule = int(self.parameters["STOPPING_RULE"])
        n_lineages = int(self.parameters["N_LINEAGES"])
        total_time = int(self.parameters["TOTAL_TIME"])
        all_lineages = len(self.whole_species_tree.get_leaves())
        dead_lineages = all_lineages - len(lineages_alive)

        success = False

        if stopping_rule == 0:
            if time_counter == total_time:
                if time_counter not in self.events:
                    self.events[time_counter] = []
                self.events[time_counter].append(("END", None, None))
                success = True


        elif stopping_rule == 1:
            if all_lineages >= n_lineages:
                if time_counter not in self.events:
                    self.events[time_counter] = []
                self.events[time_counter].append(("END", None, None))
                success = True


        elif stopping_rule == 2:
            if dead_lineages >= n_lineages:
                if time_counter not in self.events:
                    self.events[time_counter] = []
                self.events[time_counter].append(("END", None, None))
                success = True


        elif stopping_rule == 3:
            if len(lineages_alive) >= n_lineages:
                if time_counter not in self.events:
                    self.events[time_counter] = []
                self.events[time_counter].append(("END", None, None))
                success = True

        if len(lineages_alive) == 0:
            print("Simulation interrupted at time %s" % time_counter)
            print("All lineages are dead")
            if time_counter not in self.events:
                self.events[time_counter] = []
            self.events[time_counter].append(("END", None, None))
            success = False

        return success

    def get_time_to_next_event(self, n, s, e):

        total = 0.0
        for i in range(n):
            total += s
            total += e
        time = numpy.random.exponential(1/total)
        return time

    def _get_speciated(self, lineage, time_counter):

        c1 = lineage.add_child(dist=0)
        self.lineages_counter += 1
        c1.name = "n" + str(self.lineages_counter)
        c1.add_feature("is_alive", True)  # One indicates that the species is alive at present time

        c2 = lineage.add_child(dist=0)
        self.lineages_counter += 1
        c2.name = "n" + str(self.lineages_counter)
        c2.add_feature("is_alive", True)  # One indicates that the species is alive at present time

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
