from itertools import cycle, repeat



class Gene():
    def __init__(self):
        self.t = "G"
        self.len = 3
    def __str__(self):
        return "G_" + str(self.len)

class Intergene():
    def __init__(self):
        self.t = "I"
        self.len = 4
    def __str__(self):
        return "I_" + str(self.len)

my = list()

g0 = Gene()
g1 = Gene()
g2 = Gene()

g0.len = 3
g1.len = 4
g2.len = 3

i0 = Intergene()
i1 = Intergene()
i2 = Intergene()

i0.len = 3
i1.len = 4


my.append(g0)
my.append(i0)
my.append(g1)
my.append(i1)

intersections = list()


def start_cycling(position):

    p = 0

    for l in cycle(my):

        i = 0

        if l.len == 0:

            p+=1

            if p < position:
                continue

            yield "GI"

        else:

            for n in repeat(l.t, l.len):

                i += 1
                p += 1

                if p < position:
                    continue

                if i == 0 or i == l.len:
                    yield "GI"

                else:
                    yield l.t



#mm = start_cycling(4)

#for i in start_cycling(4):
#    print(mm)
#
#
import numpy

import ete3
import pyvolve

m = ete3.Tree("(n1_2:19.7901,((((n29_30:2.71742,n30_31:2.71742)n25_26:0.0797206,((n33_34:1.96966,(n39_40:0.516997,n40_41:0.516997)n34_35:1.45266)n27_28:0.749581,n28_29:2.71924)n26_27:0.0778982)n13_14:8.60907,n14_15:11.4062)n3_4:7.8954,((((((n41_42:0.167051,n42_43:0.167051)n35_36:1.59426,n36_37:1.76131)n23_24:1.59459,(n37_38:0.518771,n38_39:0.518771)n24_25:2.83713)n11_12:9.12547,(n15_16:9.86667,n16_17:9.86667)n12_13:2.61469)n9_10:1.94133,(n31_32:2.23383,n32_33:2.23383)n10_11:12.1889)n5_6:3.09489,((n19_20:6.61221,(n21_22:3.97597,n22_23:2.74445)n20_21:2.63624)n7_8:8.20344,(n17_18:9.60973,n18_19:9.60973)n8_9:5.20592)n6_7:2.70193)n4_5:1.78402)n2_3:16.8068);", format=1)
print(m.write(format=5))

tree = pyvolve.read_tree(tree=m.write(format=5))

pyvolvemap = dict()

def traverse(root):
    if root:
        if len(root.children) == 0:
            return None
        else:
            traverse(root.children[0])
            traverse(root.children[1])
            pyvolvemap[root.children[0].name +"+" + root.children[1].name] = root.name

traverse(tree)

good_mapping = dict()

for n in m.traverse(strategy="postorder"):
    if not n.is_leaf():
        c1, c2 = n.get_children()
        n1 = c1.name + "+" + c2.name
        n2 = c2.name + "+" + c1.name
        if n1 in pyvolvemap:
            good_mapping[pyvolvemap[n1]] = n.name
            n.name = pyvolvemap[n1]
        if n2 in pyvolvemap:
            good_mapping[pyvolvemap[n2]] = n.name
            n.name = pyvolvemap[n2]

for k,v in good_mapping.items():
    print(k,v)



