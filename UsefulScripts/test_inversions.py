import ete3
import os

genome_path = "/Users/davin/Desktop/Zombi/TestInversions/G/Genomes/n2528_LENGTHS.tsv"
genome_path = "/Users/davin/Desktop/Zombi/TestInversions/lengths"

#all_genomes_files = os.listdir(genomes_path)

#for genome_file in all_genomes_files:
#    print(genome_file)

distribution = list()

with open(genome_path) as f:
    f.readline()
    for line in f:
        n,t,s = line.strip().split("\t")
        if t == "I":
            distribution.append(s)

for l in distribution

