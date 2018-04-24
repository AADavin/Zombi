import sys
import os
import ete3


def create_furry_tree(path_to_extant_tree, path_events_per_branch, mapping_file):

    # First we look at the branching place of all dead lineages

    map_branching = dict()

    with open(mapping_file) as f:
        for line in f:
            l1, l2 = line.strip().split("\t")
            map_branching[l1] = l2

    # Second, we look at the arriving transfers
    # We need to check that they arrive and that they arrive to a lineage that survives. This is not obvious


    events_per_branch = os.listdir(path_events_per_branch)

    transfer_events = list()

    for myfile in events_per_branch:
        with open(os.path.join(path_events_per_branch, myfile)) as f:
            f.readline()
            for line in f:
                t, e, nodes = line.strip().split("\t")
                if e == "AT":
                    transfer_events.append((t, nodes))

    # Third, we get the midpoints of the extant tree

    points = dict()

    with open(path_to_extant_tree) as f:
        extant_tree = ete3.Tree(f.readline(), format=1)
        root = extant_tree.get_tree_root()
        root.name = "Root"
        for node in extant_tree.traverse():
            if node.is_root():
                points["Root"] = 0
            else:
                points[node.name] = node.up.get_distance(root)

    # Now I am going to map all the dead lineages to the corresponding branches of the extant tree

    times_and_diversity = list()

    extinct_diversity = dict()

    for t, nodes in transfer_events:

        ln = nodes.split(";")[0]

        # We check that the leaving lineage is not among the lineages in the extant tree

        if ln in points:
            continue

        an = nodes.split(";")[4]

        # I need to recover the time of leaving the tree an the time to reentry the tree

        if ln in map_branching:
            mybranch = map_branching[ln]
        else:
            mybranch = ln
        try:
            if ln not in extinct_diversity:
                extinct_diversity[ln] = (points[mybranch], t, mybranch, ln, an)
            else:
                if t > extinct_diversity[ln][1]:
                    extinct_diversity[ln] = (points[mybranch], t, mybranch, ln, an)
        except:
            pass

    for k, v in extinct_diversity.items():
        print(k,v)











    # This is very similar to counting the diversity in the extant tree

    times_and_diversity = list()

    for node in extant_tree.traverse():
        if node.is_leaf():
            continue
        d = node.get_distance(root)
        times_and_diversity.append((d, "+"))

    times_and_diversity = sorted(times_and_diversity, key=lambda x: x[0])

    counter = 1

    for time, event in times_and_diversity:

        if event == "+":
            counter += 1

        elif event == "-":
            counter -=1

        line = "\t".join([str(time), str(counter)])

        #print(line)




if __name__ == "__main__":

    extanttree, events_per_branch, mapping = sys.argv[1:]
    create_furry_tree(extanttree, events_per_branch, mapping)

