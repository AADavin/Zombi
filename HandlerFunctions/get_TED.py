import ete3
import sys

# Total Evolutionary Distance


def read_tree_from_file(tree_file):
    with open(tree_file) as f:
        tree = ete3.Tree(f.readline().strip(),format=1)
    return tree


def get_TED(whole_tree_file, sampled_tree_file):

    # This script only works if the names of the inner nodes are maintained between both trees
    # It outputs the name of the node in the first column and the TED in the second column

    whole_tree = read_tree_from_file(whole_tree_file)
    sampled_tree = read_tree_from_file(sampled_tree_file)

    # First thing we do is we get the nodes that are preserved between both datasets

    preserved_nodes = {x.name for x in whole_tree.traverse()}

    for snode in sampled_tree.iter_descendants(strategy="preorder"):
        ted = 0.0
        wnode = whole_tree&snode.name
        ted += wnode.dist
        nodeup = wnode.up
        while nodeup.name not in preserved_nodes and not nodeup.is_root():
            sisternode = wnode.get_sisters()[0]
            for sisternode in sisternode.traverse():
                ted += sisternode.dist
            ted += nodeup.dist
            wnode = nodeup
            nodeup = wnode.up

        print(snode.name + "\t" +  str(ted))

if __name__ == "__main__":

    args = sys.argv[1:]
    if len(args) < 2:
        print("python get_TED.py whole_tree sampled_tree")
    else:
        get_TED(args[0],args[1])

