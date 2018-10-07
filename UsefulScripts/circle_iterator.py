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



mm = start_cycling(4)

for i in start_cycling(4):
    print(mm)





