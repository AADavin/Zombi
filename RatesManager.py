from globals import *
import ete3
import argparse
import numpy


class GeneEvolutionRates():

    def __init__(self, params):

        self.params = params

        if self.params["MODE"] == 0:
            self.mode_0() # Global rates
        elif self.params["MODE"] == 1:
            pass # Rates per family
        elif self.params["MODE"] == 2:
            pass # Aucorrelated rates per lineage
        elif self.params["MODE"] == 3:
            pass
        elif self.params["MODE"] == 4:
            pass
        elif self.params["MODE"] == 5:
            pass

    def get_value(self, distribution, p0, p1):

        r = 0.0

        if distribution == "fixed":
            r = p0
        elif distribution == "uniform":
            r = numpy.random.uniform(p0,p1)
        elif distribution == "normal":
            r = numpy.random.normal(p0,p1)

        return r

    def mode_0(self):

        # Global rates

        dr = self.get_value(self.params["DUPLICATION_D"], float(self.params["DUPLICATION_P0"]), float(self.params["DUPLICATION_P1"]))
        tr = self.get_value(self.params["TRANSFER_D"], float(self.params["TRANSFER_P0"]), float(self.params["TRANSFER_P1"]))
        lr = self.get_value(self.params["LOSS_D"], float(self.params["LOSS_P0"]), float(self.params["LOSS_P1"]))

        return dr, tr, lr

    def mode_1(self):
        # Rates per family
        return self.mode_0()

    def mode_2(self, tree_file):

        # Autocorrelated time and lineages

        lineages_rates = dict()

        with open(tree_file) as f:
            tree = ete3.Tree(f.readline().strip())

        dr = self.get_value(self.params["duplication_d"], self.params["duplication_p0"], self.params["duplication_p1"])
        tr = self.get_value(self.params["transfer_d"], self.params["transfer_p0"], self.params["transfer_p1"])
        lr = self.get_value(self.params["loss_d"], self.params["loss_p0"], self.params["loss_p1"])

        myroot = tree.get_tree_root()
        myroot.add_rates("rates", (dr, tr, lr))

        lineages_rates["Root"] = (dr,tr,lr)

        for node in tree.iter_descendants(strategy="preorder"):

            node.add_feature("rates")

            fdr, ftr, flr = node.up.rates

            dr = numpy.random.normal(fdr, node.dist)
            tr = numpy.random.normal(ftr, node.dist)
            lr = numpy.random.normal(flr, node.dist)

            node.rates = (dr, tr, lr)
            lineages_rates[node.name] = (dr, tr, lr)

        return lineages_rates


    def mode_3(self, tree_file):

        # Autocorrelated in lineages

        lineages_rates = dict()

        with open(tree_file) as f:
            tree = ete3.Tree(f.readline().strip())

        myroot = tree.get_tree_root()
        myroot.name = "Root"

        for node in tree.traverse():

            dr = self.get_value(self.params["duplication_d"], self.params["duplication_p0"],
                                self.params["duplication_p1"])
            tr = self.get_value(self.params["transfer_d"], self.params["transfer_p0"], self.params["transfer_p1"])
            lr = self.get_value(self.params["loss_d"], self.params["loss_p0"], self.params["loss_p1"])

            lineages_rates[node.name] = (dr, tr, lr)

        return lineages_rates

    def mode_4(self, rates_file):

        # User defined rates

        lineages_rates = dict()
        with open(rates_file) as f:
            for line in f:
                lineage,dr,tr,lr = line.strip().split("\t")
                lineages_rates[lineage] = (float(dr),float(tr),float(lr))


class SpeciesEvolutionRates():

    def __init__(self, params):

        self.caller_count = 0
        self.params = params

        if self.params["SPECIES_EVOLUTION_MODE"] == 0:
            self.SPECIES_EVOLUTION_MODE_0() # Global rates
        elif self.params["SPECIES_EVOLUTION_MODE"] == 1:
            pass # Lineage specific (time and lineage autocorrelated)
        elif self.params["SPECIES_EVOLUTION_MODE"] == 2:
            pass # Aucorrelated rates per lineage
        elif self.params["SPECIES_EVOLUTION_MODE"] == 3:
            pass # Uncorrelated
        elif self.params["SPECIES_EVOLUTION_MODE"] == 4:
            self.SPECIES_EVOLUTION_MODE_4() # User defined

    def mode_0(self):

        sr = self.get_value(self.params["SPECIATION_D"], float(self.params["SPECIATION_P0"]),
                            float(self.params["SPECIATION_P1"]))
        er = self.get_value(self.params["EXTINCTION_D"], float(self.params["EXTINCTION_P0"]),
                            float(self.params["EXTINCTION_P1"]))

        return sr, er

    def mode_1(self, speciation_mean = 0, extinction_mean = 0):

        if self.caller_count == 0:
            sr = self.get_value(self.params["SPECIATION_D"], float(self.params["SPECIATION_P0"]),
                                float(self.params["SPECIATION_P0"]))
            er = self.get_value(self.params["EXTINCTION_D"], float(self.params["EXTINCTION_P0"]),
                                float(self.params["EXTINCTION_P1"]))

            self.caller_count += 1

            return sr, er

        else:
            ms = self.get_value(self.params["MULTIPLIERS_SPECIATION_D"], speciation_mean, float(self.params["MULTIPLIERS_SPECIATION_P1"]))
            me = self.get_value(self.params["MULTIPLIERS_EXTINCTION_D"], extinction_mean, float(self.params["MULTIPLIERS_EXTINCTION_P1"]))

            return ms, me

    def mode_2(self):

        return self.mode_0()


    def mode_3(self):

         return self.mode_0()

    def mode_4(self):

        sr = self.get_value(self.params["SPECIATION_D"], float(self.params["SPECIATION_P0"]),
                            float(self.params["SPECIATION_P0"]))
        er = self.get_value(self.params["EXTINCTION_D"], float(self.params["EXTINCTION_P0"]),
                            float(self.params["EXTINCTION_P1"]))
        periods = dict()

        with open(self.parameters["USER_DEFINED_RATES"]) as f:

            for line in f:
                start, end, spec, ext = line.strip().split("\t")
                periods[int(float(start) / TIME_INCREASE)] = (sr * int(spec), er * int(ext))
                periods[int(float(end) / TIME_INCREASE)] = (float(self.params["SPECIATION_RATE"]),
                                                            float(self.params["EXTINCTION_RATE"]))
        return sr, er, periods

    def get_value(self, distribution, p0, p1):

        r = 0.0

        if distribution == "fixed":
            r = p0
        elif distribution == "uniform":
            r = abs(numpy.random.uniform(p0,p1))
        elif distribution == "normal":
            r = abs(numpy.random.normal(p0,p1))

        return r




#lvl0parser = argparse.ArgumentParser(add_help=False)
#lvl0parser.add_argument('file', type=open, action=LoadFromFile)
'''
lvl0parser.add_argument("--mode", type=int, default=0, help="Choose the mode of gene family evolution. Check up the manual for more information")
lvl0parser.add_argument("--duplication_d", type=str,  choices=["fixed", "uniform", "normal"], help="Duplication rate distribution")
lvl0parser.add_argument("--transfer_d", type=str, choices=["fixed", "uniform", "normal"], help="Transfer rate distribution")
lvl0parser.add_argument("--loss_d", type=str, choices=["fixed", "uniform", "normal"], help="Loss rate distribution")
lvl0parser.add_argument("--replacement_d", type=str, choices=["fixed", "uniform"], help="Replacement rate distribution")
lvl0parser.add_argument("--duplication_p0", type=float, help="Duplication rate parameter 0")
lvl0parser.add_argument("--duplication_p1", type=float, help="Duplication rate parameter 1")
lvl0parser.add_argument("--transfer_p0", type=float, help="Transfer rate parameter 0")
lvl0parser.add_argument("--transfer_p1", type=float, help="Transfer rate parameter 1")
lvl0parser.add_argument("--loss_p0", type=float, help="Loss rate parameter 0")
lvl0parser.add_argument("--loss_p1", type=float, help="Loss rate parameter 1")
lvl0parser.add_argument("--replacement_r", type=float, help="Replacement rate")
'''
#args = lvl0parser.parse_args()
#rm = GeneEvolutionRates(args)

