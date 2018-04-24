import ete3
import sys

def count_lineages_in_time_extant(extant_tree_file, whole_tree_file):

    with open(whole_tree_file) as f:
        mytree = ete3.Tree(f.readline().strip(), format=1)

    root = mytree.get_tree_root()

    extant_leaves = {x.name for x in mytree.get_leaves()}
    lca = mytree.get_common_ancestor(extant_leaves)
    lca_t = lca.get_distance(root)

    with open(extant_tree_file) as f:
        mytree = ete3.Tree(f.readline().strip(), format=1)

    root = mytree.get_tree_root()

    nodes = list()

    for node in mytree.traverse():
        if node.is_leaf():
            continue
        d = node.get_distance(root)
        nodes.append((d + lca_t, node.name))

    times_and_nodes = sorted(nodes, key=lambda x: x[0])

    counter = 1

    for time, node in times_and_nodes:
        counter += 1
        line = "\t".join([str(time), str(counter)])
        print(line)


if __name__ == "__main__":


    try:
        extanttree, wholetree = sys.argv[1:]
        count_lineages_in_time_extant(extanttree, wholetree)
    except:

        print("count_lineages_in_time_extant(extanttree, wholetree)")