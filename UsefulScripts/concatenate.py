from fastareader import fastareader
import ete3
import os
import sys
import argparse

def concatenate(infile, extant_tree_path):

    with open(extant_tree_path) as f:
        extant_tree = ete3.Tree(f.readline().strip(), format=1)

    sps = extant_tree.get_leaf_names()

    concatenate = dict()

    with open(infile) as f:

        for line in f:

            gf = line.strip()

            gfs = dict()

            for h, s in fastareader(gf):
                name = h[1:].split("_")[0]
                gfs[name] = s

            for sp in sps:
                if sp not in gfs:
                    mys = "".join(["-" for x in range(100)])
                else:
                    gfs[sp] = mys

            for sp,seq in gfs.items():
                if sp not in concatenate:
                    concatenate[sp] = ""
                concatenate[sp] += seq





def fasta_reader(fasta_file):

    with open(fasta_file) as f:

        seq = ""
        for line in f:
            if ">" == line[0]:
                if seq != "":
                    yield header, seq
                    header = line.strip()
                    seq = ""
                else:
                    header = line.strip()
                    seq = ""
            else:
                seq += line.strip()

        yield header, seq








if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("e", type=str, help="ExtantTree")
    parser.add_argument("f", type=str, help="Folder with gene trees")

    args = parser.parse_args()
    infile, extant_tree,  =  args.f, args.e

    concatenate(infile, extant_tree)

