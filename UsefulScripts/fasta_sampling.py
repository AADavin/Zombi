import random
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

mapping = dict()

with open("/Users/aadavin/Desktop/Last/GregTest/mapping_names") as f:
    for line in f:
        sp, nm = line.strip().split("\t")
        mapping[nm] = sp

myfile = "/Users/aadavin/Desktop/Last/GregTest/ZombiTest/S/concatenate.fasta"

with open("/Users/aadavin/Desktop/Last/GregTest/ZombiTest/S/concatenate_goodnames.fasta", "w") as f:
    n = 14500
    #sel = random.sample(range(n), 2000)
    for h, s in fasta_reader(myfile):
        f.write(">"+mapping[h[1:]] + "\n")
        #seq = "".join([s[x] for x in sel])
        #f.write(seq+"\n")
        f.write(s)