import ete3
import os
import AuxiliarFunctions as af
from RatesManager import GeneEvolutionRates
import numpy
import random
from globals import *

class SequenceSimulator():
    def __init__(self):
        pass

    def mode_0(self):
        pass

    def mode_1(self):
        pass

    def mode_2(self):
        pass

    def mode_3(self):
        pass



def substitution_heterogenity(tree,rates,mode):

    mytree = ete3.Tree(tree,format=1)

    multiplying_factors = dict()

    # First case, branch wise correlated
    #
    if mode == "1":
        for node in mytree.traverse():
            if node.is_root():
                continue
            multiplying_factors[node.name] = get_value(rates)
            node.dist = node.dist * multiplying_factors[node.name]

    ## Second case, completely uncorrelated

    if mode == "2":
        for node in mytree.traverse():
            if node.is_root():
                continue
            total_dist = 0.0
            for i in range(int(node.dist / TIME_INCREASE)):
                total_dist += get_value(rates)

            node.dist = total_dist * TIME_INCREASE

    ## Third case, time and lineage autocorrelated

    if mode == "3":
        for node in mytree.traverse(strategy="preorder"):
            if node.is_root():
                node.add_feature("substitution_rate", numpy.random.normal(0,1))
                continue

            v = abs(numpy.random.normal(node.up.substitution_rate,node.dist))
            node.add_feature("substitution_rate", v)
            node.dist = node.dist * node.substitution_rate

    return mytree.write(format=1)
