import AuxiliarFunctions as af
import argparse
import os
import random
import numpy
import ete3


class RatesCustomizer():

    def __init__(self, parameters, experiment_folder):

        self.parameters = parameters
        self.experiment_folder = experiment_folder

    def generate_transfers_file(self):

        rates = dict()
        alive_lineages = {"Root"}

        with open(os.path.join(experiment_folder,"T/Events.tsv")) as f:

            f.readline()
            for line in f:
                t, e, nodes = line.strip().split("\t")

                if e == "S":

                    sp, c1, c2 = nodes.split(";")

                    alive_lineages.remove(sp)
                    alive_lineages.add(c1)
                    alive_lineages.add(c2)

                    for node in alive_lineages:
                        if c1 == node:
                            continue
                        if c1 not in rates:
                            rates[c1] = dict()
                        if node not in rates:
                            rates[node] = dict()
                        rates[c1][node] = 1
                        rates[node][c1] = 1

                    for node in alive_lineages:
                        if c2 == node:
                            continue
                        if c2 not in rates:
                            rates[c2] = dict()
                        if node not in rates:
                            rates[node] = dict()
                        rates[c2][node] = 1
                        rates[node][c2] = 1

                if e == "E":
                    alive_lineages.remove(nodes)


        with open(os.path.join(experiment_folder,"CustomRates/Transfer_rates.tsv"),"w") as f:
            header = "\t".join(["DONOR","RECIPIENT","WEIGHT"]) +"\n"
            f.write(header)
            for k1,v1 in rates.items():
                for k2,v2 in v1.items():
                    line = "\t".join([k1, k2, str(v2)]) + "\n"
                    f.write(line)


    def generate_rates_file(self):

        tree_file = os.path.join(self.experiment_folder, "T/WholeTree.nwk")

        with open(tree_file) as f:
            mytree = ete3.Tree(f.readline().strip(), format=1)
            root = mytree.get_tree_root()
            root.name = "Root"

        branch_rates = dict()

        for node in mytree.traverse():

            d = af.obtain_value(self.parameters["DUPLICATION"])
            t = af.obtain_value(self.parameters["TRANSFER"])
            l = af.obtain_value(self.parameters["LOSS"])
            i = af.obtain_value(self.parameters["INVERSION"])
            c = af.obtain_value(self.parameters["TRANSLOCATION"])
            o = af.obtain_value(self.parameters["ORIGINATION"])

            branch_rates[node.name] = (d,t,l,i,c,o)

        with open(os.path.join(experiment_folder,"CustomRates/Event_rates.tsv"),"w") as f:

            line = "\t".join(["lineage","D","T","L","I","C","O"]) + "\n"
            f.write(line)

            for lineage, values in branch_rates.items():

                d,t,l,i,c,o = values

                line = "\t".join(map(str,[lineage, d,t,l,i,c,o])) + "\n"
                f.write(line)

    def generate_extension_file(self):

        tree_file = os.path.join(self.experiment_folder, "T/WholeTree.nwk")

        with open(tree_file) as f:
            mytree = ete3.Tree(f.readline().strip(), format=1)
            root = mytree.get_tree_root()
            root.name = "Root"

        branch_rates = dict()

        for node in mytree.traverse():

            d = af.obtain_value(self.parameters["DUPLICATION_EXTENSION"])
            t = af.obtain_value(self.parameters["TRANSFER_EXTENSION"])
            l = af.obtain_value(self.parameters["LOSS_EXTENSION"])
            i = af.obtain_value(self.parameters["INVERSION_EXTENSION"])
            c = af.obtain_value(self.parameters["TRANSLOCATION_EXTENSION"])

            branch_rates[node.name] = (d,t,l,i,c)

        with open(os.path.join(experiment_folder,"CustomRates/Extension_rates.tsv"),"w") as f:

            line = "\t".join(["lineage","D_E","T_E","L_E","I_E","C_E"]) + "\n"
            f.write(line)

            for lineage, values in branch_rates.items():

                 d,t,l,i,c = values
                 line = "\t".join(map(str,[lineage, d,t,l,i,c])) + "\n"
                 f.write(line)


if __name__ == "__main__":


    parser = argparse.ArgumentParser()

    parser.add_argument("parameters", type=str,  help="GenomeParameters.tsv")
    parser.add_argument("output", type=str, help="Name of the experiment folder")

    args = parser.parse_args()

    parameters_file, experiment_folder = args.parameters, args.output
    parameters = af.prepare_genome_parameters(af.read_parameters(parameters_file))

    rates_folder = os.path.join(experiment_folder,"CustomRates")

    if not os.path.isdir(rates_folder):
            os.mkdir(rates_folder)

    rc = RatesCustomizer(parameters, experiment_folder)

    rc.generate_transfers_file()
    rc.generate_rates_file()
    rc.generate_extension_file()



