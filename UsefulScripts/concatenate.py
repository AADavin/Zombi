import ete3
import argparse
import sys
import os

def concatenate(inpath):

    concatenate = dict()

    fasta_files = [x for x in os.listdir(inpath) if "pruned" in x]
    
    for fasta_file in fasta_files:
        gf = fasta_file.split("_")[0]
        for h, s in fasta_reader(os.path.join(inpath, fasta_file)):
            name = h[1:].split("_")[0]
            if name not in concatenate:
                concatenate[name] = ""
            concatenate[name] += s

    for sp, seq in concatenate.items():
        print(">" + sp)
        print(seq)


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
    #parser.add_argument("e", type=str, help="ExtantTree")
    parser.add_argument("f", type=str, help="Folder with gene trees")

    args = parser.parse_args()
    inpath  =  args.f

    concatenate(inpath)

