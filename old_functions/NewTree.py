class SpeciesTree(Tree):

    def start_tree(self):

        self.name = "Root"
        self.add_feature("is_alive", True)

    def generate(self, speciation, extinction):

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

        c1 = self.add_child(dist=0)
        self.lineages_counter += 1
        c1.name = "n" + str(self.lineages_counter)
        c1.add_feature("is_alive", True)  # One indicates that the species is alive at present time

        c2 = self.add_child(dist=0)
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




