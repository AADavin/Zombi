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



myfile = "/Users/aadavin/Desktop/Last/GregTest/Tree/alignment_102_taxa.fst"
with open("/Users/aadavin/Desktop/Last/GregTest/Tree/alignment_sample2000.fst", "w") as f:
    n = 14745
    sel = random.sample(range(n), 2000)

    for h,s in fasta_reader(myfile):
        f.write(h + "\n")
        seq = "".join([s[x] for x in sel])
        f.write(seq+"\n")