import AuxiliarFunctions as af
import numpy
import ete3
import random

class EfficientSpeciesTreeGenerator():

    def __init__(self, parameters):

        self.whole_tree = ete3.Tree()

        self.lineages_counter = 0
        self.events = list()

        self.active_lineages = set()
        self.inactive_lineages = set()
        self.distances = dict()

        self.parameters = parameters


    def _start(self):

        self.active_lineages.add("Root")
        self.distances["Root"] = 0.0
        self._get_speciated("Root", 0.0)

    def run(self):

        self._start()

        speciation = self.parameters["SPECIATION"]
        extinction = self.parameters["EXTINCTION"]
        stopping_rule = self.parameters["STOPPING_RULE"]
        total_time = self.parameters["TOTAL_TIME"]
        total_lineages = self.parameters["TOTAL_LINEAGES"]
        max_lineages = self.parameters["MAX_LINEAGES"]

        time = 0
        success = False

        n_lineages_alive = 1

        while True:

            print("Time: %s ; Number of lineages alive: %s" % (str(time), str(n_lineages_alive)))

            time_to_next_event = self.get_time_to_next_event(n_lineages_alive, (speciation, extinction))

            if n_lineages_alive == 0:

                self.increase_distances(time_to_next_event)
                print("All dead")
                return success

            elif stopping_rule == 0 and time + time_to_next_event >= total_time:

                self.increase_distances(total_time - time)
                self.events.append((total_time, "F", "None"))
                success = True
                return success

            elif stopping_rule == 1 and n_lineages_alive >= total_lineages:

                self.increase_distances(time_to_next_event)
                self.events.append((time+time_to_next_event, "F", "None"))
                success = True
                return success

            elif n_lineages_alive >= max_lineages:

                print("Aborting. Max n of lineages attained")
                success = True
                return success

            else:
                # In this case we do the normal the computation

                time += time_to_next_event

                self.increase_distances(time_to_next_event)
                event = self.choose_event(speciation, extinction)
                lineage = random.sample(self.active_lineages, 1)[0]

                if event == "S":
                    self._get_speciated(lineage, time)
                    n_lineages_alive += 1

                elif event == "E":
                    self._get_extinct(lineage, time)
                    n_lineages_alive -= 1



    def run_1(self):

        self._start()

        speciation = self.parameters["SPECIATION"]
        extinction = self.parameters["EXTINCTION"]
        turnover = self.parameters["TURNOVER"]


        time_slices = self.parameters["POPULATION_SIZES"]
        total_time = time_slices[-1][0]
        current_time_slice = 0

        time = 0
        success = False
        action = 0

        n_lineages_alive = 1

        while True:

            print("Time: %s ; Number of lineages alive: %s" % (str(time), str(n_lineages_alive)))

            goal_N = time_slices[current_time_slice][1]  # The population size we have to attained

            if n_lineages_alive == goal_N:
                time_to_next_event = self.get_time_to_next_event(n_lineages_alive, [turnover])
                action = 0
            elif n_lineages_alive < goal_N:
                time_to_next_event = self.get_time_to_next_event(n_lineages_alive, [speciation])
                action = 1
            elif n_lineages_alive > goal_N:
                time_to_next_event = self.get_time_to_next_event(n_lineages_alive, [extinction])
                action = 2

            if n_lineages_alive == 0:

                self.increase_distances(time_to_next_event)
                print("All dead")
                success = False
                return success

            elif time + time_to_next_event >= total_time:

                self.increase_distances(total_time - time)
                self.events.append((total_time, "F", "None"))
                success = True
                return success

            else:

                # In this case we do the normal the computation
                time += time_to_next_event
                self.increase_distances(time_to_next_event)

                if time >= time_slices[current_time_slice][0]:
                    current_time_slice +=1

                if action == 0:

                    lineage1, lineage2 = random.sample(self.active_lineages, 2)
                    self._get_speciated(lineage1, time)
                    self._get_extinct(lineage2, time)

                elif action == 1:

                    lineage = random.sample(self.active_lineages, 1)[0]
                    self._get_speciated(lineage, time)
                    n_lineages_alive += 1

                elif action == 2:

                    lineage = random.sample(self.active_lineages, 1)[0]
                    self._get_extinct(lineage, time)
                    n_lineages_alive -= 1




    def increase_distances(self, time):

        for lineage in self.active_lineages:
            self.distances[lineage] += time

    def get_time_to_next_event(self, n, events):

        total = 0.0
        for i in range(n):
            for event in events:
                total += event
        time = numpy.random.exponential(1/total)
        return time

    def _get_speciated(self, lineage, time):

        self.lineages_counter += 1
        c1name = "n" + str(self.lineages_counter)
        self.active_lineages.add(c1name)

        self.lineages_counter += 1
        c2name = "n" + str(self.lineages_counter)
        self.active_lineages.add(c2name)

        self.distances[c1name] = 0.0
        self.distances[c2name] = 0.0

        self.inactive_lineages.add(";".join((lineage, c1name, c2name)))
        self.active_lineages.discard(lineage)

        self.events.append((time,"S", ";".join((lineage, c1name, c2name))))  # Store the event

    def _get_extinct(self, lineage, time):

        self.active_lineages.discard(lineage)
        self.inactive_lineages.add(lineage)
        self.events.append((time, "E", lineage))  # Store the event

    def choose_event(self, speciation, extinction):

        if numpy.random.uniform(0, 1) <= (speciation / (speciation + extinction)):
            return "S"
        else:
            return "E"


    def generate_whole_tree(self):

        # We create a hash of nodes for quicker access

        all_nodes = dict()

        root = self.whole_tree.get_tree_root()
        root.name = "Root"
        all_nodes[root.name] = root

        for time, event, nodes in self.events:
            if event == "S":
                p,c1name,c2name = nodes.split(";")

                node = all_nodes[p]

                c1 = node.add_child(dist=self.distances[c1name])
                c2 = node.add_child(dist=self.distances[c2name])

                c1.name = c1name
                c2.name = c2name

                all_nodes[c1name] = c1
                all_nodes[c2name] = c2


    def generate_extant_tree(self):

        leaves_alive = [x.name for x in self.whole_tree.get_leaves() if x.name in self.active_lineages]

        if len(leaves_alive) < 3:
            return 0

        self.extant_tree = ete3.Tree(self.whole_tree.write(format=1), format=1)
        self.extant_tree.prune(leaves_alive, preserve_branch_length=True)

    def write_events_file(self, events_file):

        with open(events_file, "w") as f:
            for item in self.events:
                line = "\t".join(map(str,item)) + "\n"
                f.write(line)

    def write_whole_tree(self, tree_file):

        with open(tree_file, "w") as f:
            f.write(self.whole_tree.write(format=1))

    def write_extant_tree(self, tree_file):

        with open(tree_file, "w") as f:
            f.write(self.extant_tree.write(format=1))


class SpeciesTreeGenerator():

    def __init__(self, parameters):

        self.whole_species_tree = ete3.Tree()

        self.lineages_counter = 0
        self.events = list()
        self.events.append(("0", "S", "Root;n1;n2"))
        self._start_tree(self.whole_species_tree)
        self.parameters = parameters

    def _start_tree(self, tree):

        tree.name = "Root"
        tree.add_feature("is_alive", True)

        c1 = tree.add_child(dist=0)
        self.lineages_counter += 1
        c1.name = "n" + str(self.lineages_counter)
        c1.add_feature("is_alive", True)  # One indicates that the species is alive at present time

        c2 = tree.add_child(dist=0)
        self.lineages_counter += 1
        c2.name = "n" + str(self.lineages_counter)
        c2.add_feature("is_alive", True)  # One indicates that the species is alive at present time

        tree.is_alive = False

    def generate_precomputed_tree(self, events):

        new_tree = ete3.Tree()
        self._start_tree(new_tree)

    def run(self):

        speciation = self.parameters["SPECIATION"]
        extinction = self.parameters["EXTINCTION"]
        stopping_rule = self.parameters["STOPPING_RULE"]
        total_time = self.parameters["TOTAL_TIME"]
        total_lineages = self.parameters["TOTAL_LINEAGES"]
        max_lineages = self.parameters["MAX_LINEAGES"]

        time = 0
        success = False

        while True:

            lineages_alive = [x for x in self.whole_species_tree.get_leaves() if x.is_alive == True]
            n_lineages_alive = len(lineages_alive)
            print("Time: %s ; Number of lineages alive: %s" % (str(time), str(n_lineages_alive)))

            time_to_next_event = self.get_time_to_next_event(n_lineages_alive, [speciation, extinction])

            if stopping_rule == 0 and time + time_to_next_event >= total_time:

                self.increase_distances(total_time - time, lineages_alive)
                self.events.append((total_time, "F", "None"))

                break

            elif stopping_rule == 1 and n_lineages_alive >= total_lineages:

                self.increase_distances(time_to_next_event, lineages_alive)
                self.events.append((time+time_to_next_event, "F", "None"))

                break


            elif n_lineages_alive == 0:

                self.increase_distances(time_to_next_event, lineages_alive)
                print("All dead")
                break

            elif n_lineages_alive >= max_lineages:

                print("Aborting. Max n of lineages attained")
                break

            else:
                # In this case we do the normal the computation

                time += time_to_next_event

                self.increase_distances(time_to_next_event, lineages_alive)

                event = self.choose_event(speciation, extinction)
                lineage = random.sample(lineages_alive, 1)[0]

                if event == "S":
                    self._get_speciated(lineage, time)

                elif event == "E":
                    self._get_extinct(lineage, time)

    def increase_distances(self, time, lineages_alive):

        for lineage in lineages_alive:
            lineage.dist += time

    def get_time_to_next_event(self, n, s, e):

        total = 0.0
        for i in range(n):
            total += s
            total += e
        time = numpy.random.exponential(1/total)
        return time

    def _get_speciated(self, lineage, time):

        c1 = lineage.add_child(dist=0)
        self.lineages_counter += 1
        c1.name = "n" + str(self.lineages_counter)
        c1.add_feature("is_alive", True)

        c2 = lineage.add_child(dist=0)
        self.lineages_counter += 1
        c2.name = "n" + str(self.lineages_counter)
        c2.add_feature("is_alive", True)

        self.events.append((time,"S", ";".join((lineage.name,c1.name,c2.name))))  # Store the event
        lineage.is_alive = False

    def _get_extinct(self, lineage, time):

        lineage.is_alive = False

        self.events.append((time, "E", lineage.name))  # Store the event

    def write_extant_tree(self, tree_file):

        leaves_alive = [x.name for x in self.whole_species_tree.get_leaves() if x.is_alive == True]

        if len(leaves_alive) < 3:
            return 0

        self.extant_species_tree = ete3.Tree(self.whole_species_tree.write(format=1), format=1)
        self.extant_species_tree.prune(leaves_alive, preserve_branch_length=True)

        with open(tree_file, "w") as f:
            f.write(self.extant_species_tree.write(format=1))

    def choose_event(self, speciation, extinction):

        if numpy.random.uniform(0, 1) <= (speciation / (speciation + extinction)):
            return "S"
        else:
            return "E"

    def write_events_file(self, events_file):

        with open(events_file, "w") as f:
            for item in self.events:
                line = "\t".join(map(str,item)) + "\n"
                f.write(line)

    def write_whole_tree(self, tree_file):

        with open(tree_file, "w") as f:
            f.write(self.whole_species_tree.write(format=1))
