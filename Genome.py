import AuxiliarFunctions as af
import numpy

class Genome():

    def __init__(self):

        self.genes = list()


    def start_genome(self, n):

        for i in range(1, n+1):

            if numpy.random.randint(2) == 0:
                sense = "+"
            else:
                sense = "-"

            name = "_".join([sense,str(i),"1"])
            self.genes.append(name)

    def select_random_position(self):

        return numpy.random.randint(len(self.genes))

    def select_random_length(self):

        return numpy.random.geometric(0.5)

    def obtain_affected_genes(self):

        position = self.select_random_position()
        length = self.select_random_length()

        total_length = len(self.genes)

        affected_genes = list()

        for i in range(position, position + length):
            if i >= total_length:
                affected_genes.append(i - total_length)
            else:
                affected_genes.append(i)

        return affected_genes

    def obtain_segment(self, affected_genes):

        segment = list()
        for x in affected_genes:
            segment.append(self.genes[x])
        return segment


    def duplicate_segment(self, homologous, time, affected_genes, adjacent = True):

        old_segment = [self.genes[x] for x in affected_genes]

        new_segment1 = list()
        new_segment2 = list()

        for i in affected_genes:

            sense, gf, cp = self.genes[i].split("_")

            homologous[gf]["Copies"] += 1
            name1 = sense + "_" + gf + "_" + str(homologous[gf]["Copies"])

            new_segment1.append(name1)

            homologous[gf]["Copies"] += 1
            name2 = sense + "_" + gf + "_" + str(homologous[gf]["Copies"])

            new_segment2.append(name2)

            homologous[gf]["Events"].append(("D", time, cp + "_" + name1.split("_")[2] + "_" + name2.split("_")[2]))

        # Now we have to delete the ancient segment that has been duplicated and insert in the new position
        # both segments

        elements_to_remove = [self.genes[x] for x in affected_genes]

        for element in elements_to_remove:

            self.genes.remove(element)

        position = affected_genes[0]

        for i, x in enumerate(new_segment1 + new_segment2):

            self.genes.insert(position + i + 1, x)


    def loss_segment(self, homologous, time, affected_genes):

        elements_to_remove = [self.genes[x] for x in affected_genes]

        for element in elements_to_remove:

            gf = element.split("_")[1]

            self.genes.remove(element)

            homologous[gf]["Events"].append(("L", time, element))

    def insert_segment(self, position, segment):

        for i, x in enumerate(segment):
            self.genes.insert(position + i + 1, x)

    def invert_segment(self, homologous,time, affected_genes):

        segment = [self.genes[x] for x in affected_genes]
        new_segment = list()

        for i, x in enumerate(segment):

            name = x.replace("-","A")
            name = name.replace("+","-")
            name = name.replace("A", "+")
            new_segment.append(name)

        new_segment = new_segment[::-1]

        for i,x in enumerate(new_segment):
           self.genes[affected_genes[i]] = x
           gf,element = x.split("_")[1:]
           homologous[gf]["Events"].append(("I", time, element))

    def translocate_segment(self, homologous, time, affected_genes): # Watch out, verify that it is outside

        # first we remove

        segment = [self.genes[x] for x in affected_genes]
        for element in segment:
            self.genes.remove(element)

        # second we insert

        position = self.select_random_position()

        for i, x in enumerate(segment):
            self.genes.insert(position + i + 1, x)
            gf, element = x.split("_")[1:]
            homologous[gf]["Events"].append(("C", time, element))

    def update_homologous(self, homologous):

        for i,gene in enumerate(self.genes):
            sense,gf,hml = gene.split("_")
            homologous[gf]["Copies"] += 1

            new_hml = homologous[gf]["Copies"]
            self.genes[i] = "_".join((sense,gf,str(new_hml)))

        return self.genes


    def write_genome(self, genome_file):

        with open(genome_file) as f:

            f.write("\t".join(("Position","Gene_family","Orientation","Id"))+"\t")

            for i,gene in enumerate(self.genes):

                sense, gf, id = gene.split("_")
                line = "\t".join((str(i),gf,sense,id)) + "\n"
                f.write(line)









'''

gnm = Genome()
gnm.start_genome(10)

duplication = 5
loss = 2
inversion = 5
traslocation = 3
transfer = 0

def choose_event(duplication, transfer, loss, inversion, traslocation):

    draw = numpy.random.choice(["D", "T", "L", "I", "C"], 1, p=af.normalize([duplication, transfer, loss, inversion, traslocation]))
    return draw

segment = ["+_100_1"]

for i in range(10):

    event = choose_event(duplication, transfer, loss, inversion, traslocation)
    a = gnm.obtain_affected_genes()

    if event == "D":
        gnm.duplicate_segment(i, a)

    elif event == "T":
        position = gnm.select_random_position()
        gnm.insert_segment(position, segment)

    elif event == "L":
        gnm.loss_segment(i, a)

    elif event == "I":
        gnm.invert_segment(a)

    elif event == "C":
        gnm.traslocate_segment(a)

    print(event, gnm.genes)

'''








