import sys
import os
import ete3


def map_extinct_lineages(path_to_whole_tree, path_to_extant_tree):

    # This reads the whole tree

    whole_nodes_dict = dict()

    with open(path_to_whole_tree) as f:
        whole_tree = ete3.Tree(f.readline().strip(), format=1)
        r = whole_tree.get_tree_root()
        r.name = "Root"
        for node in whole_tree.traverse():
            whole_nodes_dict[node.name] = node

    # This reads the extant tree

    extant_nodes_dict = dict()

    with open(path_to_extant_tree) as f:
        extant_tree = ete3.Tree(f.readline().strip(), format=1)
        r = extant_tree.get_tree_root()
        r.name = "Root"
        for node in extant_tree.traverse():
            extant_nodes_dict[node.name] = node

    # Now I am going to map all the dead lineages to the corresponding branches of the extant tree

    extinct_leaves = [whole_nodes_dict[x].name for x in
                      {x.name for x in whole_tree.get_leaves()} - {x.name for x in extant_tree.get_leaves()}]

    missing_lineages = dict()

    middle_groups = dict()


    for leaf in extinct_leaves:
        extinct_lineages = list()
        mynode = leaf

        while mynode not in extant_nodes_dict:
            extinct_lineages.append(mynode)
            old_node = mynode
            mynode = whole_nodes_dict[mynode].up.name

        if mynode not in middle_groups:
            middle_groups[mynode] =

        for node in extinct_lineages:
            if node not in missing_lineages:
                missing_lineages[node] = mynode

    ## This maps the missing clades into the upper branch, but we want to map them into the lower branch
    ## For that, I check all the descendants of the missing lineages dict:

    good_mapping = dict()


    for extinct_lineage, parent_node in missing_lineages.items():
        ec1, ec2 = extant_nodes_dict[parent_node].get_children()
        wc1, wc2 = whole_nodes_dict[parent_node].get_children()
        dwc1 = [x.name for x in wc1.get_descendants()] + [wc1.name]
        dwc2 = [x.name for x in wc2.get_descendants()] + [wc2.name]
        if ec1.name in dwc1 and extinct_lineage in dwc1:
            good_mapping[extinct_lineage] = ec1.name
        elif ec2.name in dwc2 and extinct_lineage in dwc2:
            good_mapping[extinct_lineage] = ec2.name


    for k,v in good_mapping.items():
        print(k + "\t" + v)




if __name__ == "__main__":

    extanttree, wholetree = sys.argv[1:]
    map_extinct_lineages(wholetree, extanttree)

