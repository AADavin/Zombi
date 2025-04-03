import AuxiliarFunctions as af
import numpy
import ete3
import random


class SpeciesTreeGenerator():

    def __init__(self, parameters):

        self.parameters = parameters

        try:


            mseed = self.parameters["SEED"]
            if mseed != 0:
                random.seed(parameters["SEED"])
                numpy.random.seed(parameters["SEED"])
        except:
            pass

    def start(self):

        self.whole_tree = ete3.Tree()

        self.lineages_counter = 0
        self.events = list()

        self.active_lineages = set()
        self.inactive_lineages = set()
        self.distances = dict()

        self.active_lineages.add("Root")
        self.distances["Root"] = 0.0

    def run(self):

        self.start()

        speciation = af.obtain_value(self.parameters["SPECIATION"])
        extinction = af.obtain_value(self.parameters["EXTINCTION"])
        stopping_rule = self.parameters["STOPPING_RULE"]
        total_time = self.parameters["TOTAL_TIME"]
        total_lineages = self.parameters["TOTAL_LINEAGES"]
        max_lineages = self.parameters["MAX_LINEAGES"]

        time = 0

        n_lineages_alive = 1

        while True:

            if n_lineages_alive == 0:

                print("All dead")
                success = False
                return success

            if self.parameters["VERBOSE"] == 1:
                print("Time: %s ; Number of lineages alive: %s" % (str(time), str(n_lineages_alive)))

            time_to_next_event = self.get_time_to_next_event(n_lineages_alive, (speciation, extinction))

            if stopping_rule == 0 and time + time_to_next_event >= total_time:

                self.increase_distances(total_time - time)
                for lineage in self.active_lineages:
                    self.events.append((total_time, "F", lineage))
                success = True
                return success

            elif stopping_rule == 1 and n_lineages_alive == total_lineages:

                self.increase_distances(time_to_next_event)
                for lineage in self.active_lineages:
                    self.events.append((time + time_to_next_event, "F", lineage))
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
                lineage = random.sample(sorted(self.active_lineages), 1)[0]

                if event == "S":
                    self._get_speciated(lineage, time)
                    n_lineages_alive += 1

                elif event == "E":
                    self._get_extinct(lineage, time)
                    n_lineages_alive -= 1

    def run_b(self):

        # Speciation and extinction rates are branch-wise
        # Each time I create a new lineage, I have to generate new number for its rates
        # We create a dictionary to store the rates

        speciation = af.obtain_value(self.parameters["SPECIATION"])
        extinction = af.obtain_value(self.parameters["EXTINCTION"])

        self.branchwise_rates = dict()
        self.branchwise_rates["Root"] = (speciation, extinction)

        self.start()


        stopping_rule = self.parameters["STOPPING_RULE"]
        total_time = self.parameters["TOTAL_TIME"]
        total_lineages = self.parameters["TOTAL_LINEAGES"]
        max_lineages = self.parameters["MAX_LINEAGES"]

        time = 0

        n_lineages_alive = 1


        while True:

            if n_lineages_alive == 0:

                print("All dead")
                success = False
                return success

            if self.parameters["VERBOSE"] == 1:
                print("Time: %s ; Number of lineages alive: %s" % (str(time), str(n_lineages_alive)))

            time_to_next_event = self.get_time_to_next_event_advanced_modes()

            if stopping_rule == 0 and time + time_to_next_event >= total_time:

                self.increase_distances(total_time - time)
                for lineage in self.active_lineages:
                    self.events.append((total_time, "F", lineage))
                success = True
                return success

            elif stopping_rule == 1 and n_lineages_alive == total_lineages:

                self.increase_distances(time_to_next_event)
                for lineage in self.active_lineages:
                    self.events.append((time + time_to_next_event, "F", lineage))
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

                # Now we have to choose the lineage doing the event. This will be proportional to the value of the rates
                ###

                active_lineages = list(sorted(self.active_lineages))

                lineage = numpy.random.choice(active_lineages, 1, p=af.normalize(
                    [sum((self.branchwise_rates[x][0], self.branchwise_rates[x][1])) for x in active_lineages]))[0]

                myspeciation = self.branchwise_rates[lineage][0]
                myextinction = self.branchwise_rates[lineage][1]

                event = self.choose_event(myspeciation, myextinction)

                if event == "S":
                    c1, c2 = self._get_speciated(lineage, time)
                    n_lineages_alive += 1

                    self.branchwise_rates[c1] = (
                        af.obtain_value(self.parameters["SPECIATION"]),
                        af.obtain_value(self.parameters["EXTINCTION"]))
                    self.branchwise_rates[c2] = (
                        af.obtain_value(self.parameters["SPECIATION"]),
                        af.obtain_value(self.parameters["EXTINCTION"]))

                elif event == "E":
                    self._get_extinct(lineage, time)
                    n_lineages_alive -= 1


    def run_p(self):

        self.start()
        print("Computing tree with fine control of the number of lineages")

        speciation = af.obtain_value(self.parameters["SPECIATION"])
        extinction = af.obtain_value(self.parameters["EXTINCTION"])
        turnover = self.parameters["TURNOVER"]

        time_slices = self.parameters["LINEAGE_PROFILE"]
        total_time = time_slices[-1][0]
        current_time_slice = 0

        time = 0
        action = 0

        n_lineages_alive = 1

        while True:

            if self.parameters["VERBOSE"] == 1:
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
                for lineage in self.active_lineages:
                    self.events.append((total_time, "F", lineage))
                success = True
                return success

            else:

                # In this case we do the normal the computation
                time += time_to_next_event
                self.increase_distances(time_to_next_event)

                if time >= time_slices[current_time_slice][0]:
                    current_time_slice +=1

                if action == 0:

                    lineage1, lineage2 = random.sample(sorted(self.active_lineages), 2)
                    self._get_speciated(lineage1, time)
                    self._get_extinct(lineage2, time)

                elif action == 1:

                    lineage = random.sample(sorted(self.active_lineages), 1)[0]
                    self._get_speciated(lineage, time)
                    n_lineages_alive += 1

                elif action == 2:

                    lineage = random.sample(sorted(self.active_lineages), 1)[0]
                    self._get_extinct(lineage, time)
                    n_lineages_alive -= 1

    def run_m(self):

        self.start()

        speciation = af.obtain_value(self.parameters["SPECIATION"])
        extinction = af.obtain_value(self.parameters["EXTINCTION"])
        stopping_rule = self.parameters["STOPPING_RULE"]
        total_time = self.parameters["TOTAL_TIME"]
        total_lineages = self.parameters["TOTAL_LINEAGES"]
        max_lineages = self.parameters["MAX_LINEAGES"]

        handle = self.parameters["MASS_EXTINCTION"].split("-") # I should add the possibility of having several extinctions
        extinction_time, p_extinction = float(handle[0]), float(handle[1])
        time = 0

        n_lineages_alive = 1

        is_extinction_over = False

        while True:

            if n_lineages_alive == 0:

                print("All dead")
                success = False
                return success

            if self.parameters["VERBOSE"] == 1:
                print("Time: %s ; Number of lineages alive: %s" % (str(time), str(n_lineages_alive)))

            time_to_next_event = self.get_time_to_next_event(n_lineages_alive, (speciation, extinction))

            if time + time_to_next_event >= extinction_time and is_extinction_over == False:

                self.increase_distances(extinction_time - time)

                time = extinction_time
                extinctions = set()
                for lineage in self.active_lineages:
                    if numpy.random.uniform(0, 1) <= p_extinction:
                        extinctions.add(lineage)
                for lineage in extinctions:
                    self._get_extinct(lineage, time)
                    self.events.append((time, "E", lineage))
                    n_lineages_alive -= 1
                is_extinction_over = True

            elif stopping_rule == 0 and time + time_to_next_event >= total_time:

                self.increase_distances(total_time - time)
                for lineage in self.active_lineages:
                    self.events.append((total_time, "F", lineage))
                success = True
                return success

            elif stopping_rule == 1 and n_lineages_alive == total_lineages:

                self.increase_distances(time_to_next_event)
                for lineage in self.active_lineages:
                    self.events.append((time + time_to_next_event, "F", lineage))
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
                lineage = random.sample(sorted(self.active_lineages), 1)[0]

                if event == "S":
                    self._get_speciated(lineage, time)
                    n_lineages_alive += 1

                elif event == "E":
                    self._get_extinct(lineage, time)
                    n_lineages_alive -= 1

    def run_s(self):

        # Shift-birth-death model

        hspeciation = self.parameters["BASE_SPECIATION"]
        hextinction = self.parameters["BASE_EXTINCTION"]
        cat_speciation = int(self.parameters["NUM_SPECIATION_RATE_CATEGORIES"])
        cat_extinction = int(self.parameters["NUM_EXTINCTION_RATE_CATEGORIES"])
        s_speciation = af.obtain_value(self.parameters["SHIFT_SPECIATION_RATE_FREQUENCY"])
        s_extinction = af.obtain_value(self.parameters["SHIFT_EXTINCTION_RATE_FREQUENCY"])

        # We get the categories for speciations

        distribution, value = hspeciation.split(":")
        value = float(value)

        if distribution == "g":
            speciation_rates = af.normalize_middle(af.discretize(value,cat_speciation,"gamma"))
        elif distribution == "l":
            speciation_rates = af.normalize_middle(af.discretize(value,cat_speciation,"lognorm"))
        else:
            print("Unrecognized distribution. Please, use g or l")
            return False

        # We get the categories for extinctions

        distribution, value = hextinction.split(":")
        value = float(value)

        if distribution == "g":
            extinction_rates = af.normalize_middle(af.discretize(value,cat_extinction,"gamma"))
        elif distribution == "l":
            extinction_rates = af.normalize_middle(af.discretize(value,cat_extinction,"lognorm"))
        else:
            print("Unrecognized distribution. Please, use g or l")
            return False

        self.extinction_rates = extinction_rates
        self.speciation_rates = speciation_rates

        # We will write the shift events

        self.shift_events = list()

        # We put the initial category right in the middle

        self.category_position = dict()
        self.category_position["Root"] = (int(cat_speciation/2), int(cat_extinction/2))

        speciation = speciation_rates[self.category_position["Root"][0]]
        extinction = extinction_rates[self.category_position["Root"][1]]

        self.branchwise_rates = dict()
        self.branchwise_rates["Root"] = (speciation, extinction, s_speciation, s_extinction)

        self.start()

        stopping_rule = self.parameters["STOPPING_RULE"]
        total_time = self.parameters["TOTAL_TIME"]
        total_lineages = self.parameters["TOTAL_LINEAGES"]
        max_lineages = self.parameters["MAX_LINEAGES"]

        time = 0

        n_lineages_alive = 1

        while True:

            if n_lineages_alive == 0:

                print("All dead")
                success = False
                return success

            if self.parameters["VERBOSE"] == 1:
                print("Time: %s ; Number of lineages alive: %s" % (str(time), str(n_lineages_alive)))

            time_to_next_event = self.get_time_to_next_event_advanced_modes()

            if stopping_rule == 0 and time + time_to_next_event >= total_time:

                self.increase_distances(total_time - time)
                for lineage in self.active_lineages:
                    self.events.append((total_time, "F", lineage))
                success = True
                return success

            elif stopping_rule == 1 and n_lineages_alive == total_lineages:

                self.increase_distances(time_to_next_event)
                for lineage in self.active_lineages:
                    self.events.append((time + time_to_next_event, "F", lineage))
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

                # Now we have to choose the lineage doing the event. This will be proportional to the value of the rates
                ###

                active_lineages = list(sorted(self.active_lineages))

                lineage = numpy.random.choice(active_lineages, 1, p=af.normalize(
                    [sum(self.branchwise_rates[x])
                     for x in active_lineages]))[0]

                myspeciation = self.branchwise_rates[lineage][0]
                myextinction = self.branchwise_rates[lineage][1]
                myshiftspeciation = self.branchwise_rates[lineage][2]
                myshiftextinction = self.branchwise_rates[lineage][3]

                event = self.choose_event_s_mode(myspeciation, myextinction,
                                                myshiftspeciation, myshiftextinction)

                if event == "S":
                    c1, c2 = self._get_speciated(lineage, time)
                    n_lineages_alive += 1

                    # We inherit the values

                    self.branchwise_rates[c1] = (myspeciation, myextinction, myshiftspeciation, myshiftextinction)
                    self.branchwise_rates[c2] = (myspeciation, myextinction, myshiftspeciation, myshiftextinction)

                    self.category_position[c1] = self.category_position[lineage]
                    self.category_position[c2] = self.category_position[lineage]

                elif event == "E":
                    self._get_extinct(lineage, time)
                    n_lineages_alive -= 1

                elif event == "SS":

                    # Shift speciation

                    cat = self.category_position[lineage][0]
                    if cat == cat_speciation-1:
                        direction = numpy.random.choice([-1,0])
                    elif cat == 0:
                        direction = numpy.random.choice([0,1])
                    else:
                        direction = numpy.random.choice([-1,1], p = [0.5,0.5])

                    p_sp, p_ex = self.category_position[lineage]
                    self.category_position[lineage] = (p_sp + direction, p_ex)
                    new_speciation = speciation_rates[self.category_position[lineage][0]]
                    self.shift_events.append((time, "SS", lineage, self.branchwise_rates[lineage][0], new_speciation))
                    self.events.append((time, "SS",  lineage + ";" + str(self.branchwise_rates[lineage][0]) + "->" + str(new_speciation)))
                    self.branchwise_rates[lineage] = (new_speciation, myextinction, myshiftspeciation, myshiftextinction)


                elif event == "SE":

                    # Shift event
                    cat = self.category_position[lineage][1]
                    if cat == cat_extinction-1:
                        direction = numpy.random.choice([-1,0])
                    elif cat == 0:
                        direction = numpy.random.choice([0,1])
                    else:
                        direction = numpy.random.choice([-1,1], p = [0.5,0.5])
                    p_sp, p_ex = self.category_position[lineage]
                    self.category_position[lineage] = (p_sp, p_ex + direction)
                    new_extinction = extinction_rates[self.category_position[lineage][1]]
                    self.shift_events.append((time, "SE", lineage, self.branchwise_rates[lineage][1], new_extinction))
                    self.events.append((time, "SE", lineage + ";" + str(self.branchwise_rates[lineage][1]) + "->" + str(new_extinction)))
                    self.branchwise_rates[lineage] = (myspeciation, new_extinction, myshiftspeciation, myshiftextinction)



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

    def get_time_to_next_event_advanced_modes(self):
        # To obtain the time to next event in case that we have different rates per branch
        total = 0.0
        for lineage in self.active_lineages:
            total += sum(self.branchwise_rates[lineage])

        time = numpy.random.exponential(1 / total)
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

        return c1name, c2name

    def _get_extinct(self, lineage, time):

        self.active_lineages.discard(lineage)
        self.inactive_lineages.add(lineage)
        self.events.append((time, "E", lineage))  # Store the event

    def choose_event(self, speciation, extinction):

        if numpy.random.uniform(0, 1) <= (speciation / (speciation + extinction)):
            return "S"
        else:
            return "E"

    def choose_event_s_mode(self, speciation, extinction, shif_speciation, shif_extinction):

        draw = numpy.random.choice(["S", "E", "SS", "SE"], 1,
                                   p=af.normalize([speciation, extinction,
                                                   shif_speciation, shif_extinction]))
        return draw

    def generate_newick_trees(self):

        def find_descendant(surviving_nodes, node):

            found = 0
            mynode = surviving_nodes[node]["descendant"]
            collapsed_nodes = list()

            while found == 0:

                if surviving_nodes[mynode]["state"] == 1:
                    collapsed_nodes.append(mynode)
                    found = 1
                else:
                    collapsed_nodes.append(mynode)
                    mynode = surviving_nodes[mynode]["descendant"]

            return mynode, collapsed_nodes

        def get_extinct(surviving_nodes, node):

            extinct_nodes = list()

            if surviving_nodes[node]["extinct"] == "E":
                extinct_nodes.append(node)
            else:
                extinct_nodes.append(node)
                extinct_nodes += surviving_nodes[node]["extinct"].split(";")

            return ";".join(extinct_nodes)

        # Eric's algorithm

        # First we will iterate the events from the end

        events = self.events

        surviving_nodes = dict()
        times = dict()

        for current_time, event, nodes in events[::-1]:

            if event == "F":

                times[nodes] = float(current_time)
                surviving_nodes[nodes] = {"state": 1, "descendant": "None", "collapsed": "", "extinct": ""}

            elif event == "E":

                times[nodes] = float(current_time)
                surviving_nodes[nodes] = {"state": 0, "descendant": "None", "collapsed": "", "extinct": "E"}

            elif event == "S":

                p, c1, c2 = nodes.split(";")

                times[p] = float(current_time)

                if surviving_nodes[c1]["state"] == 1 and surviving_nodes[c2]["state"] == 1:
                    surviving_nodes[p] = {"state": 1, "descendant": c1 + ";" + c2, "collapsed": "", "extinct": ""}

                elif surviving_nodes[c1]["state"] == 0 and surviving_nodes[c2]["state"] == 0:

                    en1 = get_extinct(surviving_nodes, c1)
                    en2 = get_extinct(surviving_nodes, c2)

                    surviving_nodes[p] = {"state": 0, "descendant": "None", "collapsed": "", "extinct": en1 + ";" + en2}

                elif surviving_nodes[c1]["state"] == -1 and surviving_nodes[c2]["state"] == -1:

                    en1 = get_extinct(surviving_nodes, c1)
                    en2 = get_extinct(surviving_nodes, c2)

                    mynode1, cp_nodes1 = find_descendant(surviving_nodes, c1)
                    mynode2, cp_nodes2 = find_descendant(surviving_nodes, c2)
                    surviving_nodes[p] = {"state": 1, "descendant": mynode1 + ";" + mynode2, "collapsed": surviving_nodes[c1]["collapsed"] + "+" + surviving_nodes[c2]["collapsed"],
                                          "extinct": en1 + "+" + en2}

                elif surviving_nodes[c1]["state"] == 1 and surviving_nodes[c2]["state"] == 0:

                    en2 = get_extinct(surviving_nodes, c2)
                    surviving_nodes[p] = {"state": -1, "descendant": c1, "collapsed": p, "extinct": en2}

                elif surviving_nodes[c1]["state"] == 0 and surviving_nodes[c2]["state"] == 1:

                    en1 = get_extinct(surviving_nodes, c1)
                    surviving_nodes[p] = {"state": -1, "descendant": c2, "collapsed": p, "extinct": en1}

                elif surviving_nodes[c1]["state"] == 1 and surviving_nodes[c2]["state"] == -1:
                    mynode, cp_nodes = find_descendant(surviving_nodes, c2)
                    surviving_nodes[p] = {"state": 1, "descendant": c1 + ";" + mynode, "collapsed": "N+" + surviving_nodes[c2]["collapsed"],
                                          "extinct": ""}

                elif surviving_nodes[c1]["state"] == -1 and surviving_nodes[c2]["state"] == 1:
                    mynode, cp_nodes = find_descendant(surviving_nodes, c1)
                    surviving_nodes[p] = {"state": 1, "descendant": mynode + ";" + c2, "collapsed": surviving_nodes[c1]["collapsed"] + "+N",
                                          "extinct": ""}

                elif surviving_nodes[c1]["state"] == -1 and surviving_nodes[c2]["state"] == 0:
                    mynode, cp_nodes = find_descendant(surviving_nodes, c1)

                    en2 = get_extinct(surviving_nodes, c2)

                    surviving_nodes[p] = {"state": -1, "descendant": mynode, "collapsed": surviving_nodes[c1]["collapsed"] + ";" + p,
                                          "extinct": en2}

                elif surviving_nodes[c1]["state"] == 0 and surviving_nodes[c2]["state"] == -1:
                    mynode, cp_nodes = find_descendant(surviving_nodes, c2)

                    en1 = get_extinct(surviving_nodes, c2)

                    surviving_nodes[p] = {"state": -1, "descendant": mynode, "collapsed": surviving_nodes[c2]["collapsed"] + ";" + p,
                                          "extinct": en1}

        extanttree = ete3.Tree()
        wholetree = ete3.Tree()
        eroot = extanttree.get_tree_root()
        eroot.name = ""
        wroot = wholetree.get_tree_root()
        wroot.name = "Root"
        wroot.dist = float(events[0][0]) # We add the root distance

        t = (len(events))

        wquick_nodes = dict()
        equick_nodes = dict()

        wquick_nodes["Root"] = wroot

        # I create a dict for storing the collapsed nodes:

        map_collapsed = dict()
        map_extinct = dict()

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

                if state == -1:

                    extinct_nodes = surviving_nodes[p]["extinct"]

                    if extinct_nodes != "":
                        if "+" in extinct_nodes:
                            ep1, ep2 = extinct_nodes.split("+")
                            if ep1 != "N":
                                map_extinct[c1name] = ep1
                            if ep2 != "N":
                                map_extinct[c2name] = ep2

                            #print(p + " -> " +  " : " + ep1)
                            #print(p + " -> " +  " : " + ep2)

                        else:
                            pass
                            #print(p + " -> " +  " : " + extinct_nodes)

                if state == 1:  # Now the extant tree

                    c1name, c2name = surviving_nodes[p]["descendant"].split(";")

                    collapsed_nodes = surviving_nodes[p]["collapsed"]
                    if collapsed_nodes != "":
                        cp1, cp2 = collapsed_nodes.split("+")
                        if cp1 != "N":
                            map_collapsed[c1name] = cp1
                        if cp2 != "N":
                            map_collapsed[c2name] = cp2

                        #print(p + " -> " + c1name + " : " + cp1)
                        #print(p + " -> " + c2name + " : " + cp2)

                    extinct_nodes = surviving_nodes[p]["extinct"]

                    if extinct_nodes != "":
                        if "+" in extinct_nodes:
                            ep1, ep2 = extinct_nodes.split("+")
                            if ep1 != "N":
                                map_extinct[c1name] = ep1
                            if ep2 != "N":
                                map_extinct[c2name] = ep2

                            #print(p + " -> " + c1name + " : " + ep1)
                            #print(p + " -> " + c2name + " : " + ep2)

                        else:
                            pass

                            #print(p + " -> " + c1name + " : " + extinct_nodes)
                            #print(p + " -> " + c2name + " : " + extinct_nodes)

                    if eroot.name == "":
                        eroot.name = p
                        equick_nodes[p] = eroot
                        eroot.dist = events[0][0]

                    mynode = equick_nodes[p]
                    myc1 = mynode.add_child()
                    myc2 = mynode.add_child()
                    myc1.name = c1name
                    myc2.name = c2name
                    myc1.dist = times[c1name] - times[p]
                    myc2.dist = times[c2name] - times[p]
                    equick_nodes[c1name] = myc1
                    equick_nodes[c2name] = myc2

        # There is a more efficient way to do this

        if eroot.name != "Root":
            for time, event, nodes in events:
                if event == "S":
                    n0,n1,n2 = nodes.split(";")
                    if n0 == eroot.name:
                        eroot.dist = float(time)
                        break

        return wholetree.write(format=1, format_root_node = True), extanttree.write(format=1, format_root_node = True), map_collapsed

    def scale_trees(self, extant_tree, scaling):

        mextant_tree = ete3.Tree(extant_tree, format=1)
        eroot = mextant_tree.get_tree_root()
        eleaf = mextant_tree.get_leaves()[0]
        crown_length = eroot.get_distance(eleaf)

        mfactor = scaling / float(crown_length)

        for n in mextant_tree.traverse():
            n.dist *= mfactor

        scaled_events = list()

        beginning = ""

        for t, kind, nodes in self.events:
            if kind == "S" and eroot.name == nodes.split(";")[0]:
                beginning = float(t)
                break
        if beginning == "":
            beginning = 0


        for t, kind, nodes in self.events:
            scaled_events.append(( (float(t) - beginning) * mfactor, kind, nodes))

        #print(mfactor, beginning, eroot.name)

        return (mextant_tree.write(format=1), scaled_events)


    def quick_pruner():

        mextant_tree


    def write_events_file(self, events_file):

        header = ["TIME","EVENT","NODES"]
        header = "\t".join(map(str, header)) + "\n"

        with open(events_file, "w") as f:
            f.write(header)
            for item in self.events:
                line = "\t".join(map(str,item)) + "\n"
                f.write(line)

    def write_rates(self, rates_file):

        with open(rates_file, "w") as f:

            line = "\t".join(["lineage", "speciation", "extinction"]) + "\n"
            f.write(line)

            for lineage, values in self.branchwise_rates.items():

                speciation, extinction = values

                line = "\t".join(map(str,[lineage, speciation, extinction])) + "\n"
                f.write(line)

    def write_lengths(self, length_file, complete_tree, extant_tree):

        mcomplete_tree = ete3.Tree(complete_tree, format=1)
        mextant_tree = ete3.Tree(extant_tree, format=1)

        croot = mcomplete_tree.get_tree_root()
        cleaf = mcomplete_tree.get_leaves()[0]

        eroot = mextant_tree.get_tree_root()
        eleaf = mextant_tree.get_leaves()[0]

        crown_length = eroot.get_distance(eleaf)
        stem_length = eroot.dist
        root_length = croot.dist

        with open(length_file, "w") as f:

            f.write("### Extant Tree lengths \n")

            line = "\t".join(["Total_length", str(stem_length + crown_length)]) + "\n"
            f.write(line)

            line = "\t".join(["Root_length", str(root_length)]) + "\n"
            f.write(line)

            line = "\t".join(["Stem_length", str(stem_length)]) + "\n"
            f.write(line)

            line = "\t".join(["Crown_length", str(crown_length)]) + "\n"
            f.write(line)


            #line = "\t".join(["Total_ED", str(total_ed)]) + "\n"
            #f.write(line)

            #line = "\t".join(["Crown_ED", str(crown_ed)]) + "\n"
            #f.write(line)

    def write_scaled_files(self, scaled_tree, scaled_extant_tree_file, scaled_events, scaled_events_file):

        stree = ete3.Tree(scaled_tree, format=1)
        with open(scaled_extant_tree_file, "w") as f:
            f.write(stree.write(format=1, format_root_node = True))
        with open(scaled_events_file, "w") as f:

            header = ["SCALED_TIME","EVENT","NODES"]
            header = "\t".join(map(str, header)) + "\n"
            f.write(header)
            for tp in scaled_events:
                f.write("\t".join(list(map(str,tp)))+"\n")

    #def write_shifts(self, shift_file):

    #    with open(shift_file, "w") as f:

    #        header = ["NODE", "TIME_INTERVAL", "EFFECTIVE_RATE"]
    #        header = "\t".join(map(str, header)) + "\n"
    #        f.write(header)

    def write_categories(self, categories_file):

        with open(categories_file, "w") as f:

            f.write("SPECIATION_RATE_CATEGORIES\n")
            f.write("\t".join(map(str,self.speciation_rates))+"\n")
            f.write("EXTINCTION_RATE_CATEGORIES\n")
            f.write("\t".join(map(str,self.extinction_rates))+"\n")
