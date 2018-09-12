from Bio import AlignIO

input_handle = open("/Users/aadavin/Desktop/Last/GregTest/Tree/alignment_sample2000.fst", "rU")
output_handle = open("/Users/aadavin/Desktop/Last/GregTest/Tree/alignment_sample2000.phy", "w")

alignments = AlignIO.parse(input_handle, "fasta")
AlignIO.write(alignments, output_handle, "phylip-relaxed")

output_handle.close()
input_handle.close()