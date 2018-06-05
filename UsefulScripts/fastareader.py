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