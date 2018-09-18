from Bio import AlignIO

input_handle = open("/Users/aadavin/Desktop/Last/GregTest/ZombiTest/S/concatenate2000.fasta", "rU")
output_handle = open("/Users/aadavin/Desktop/Last/GregTest/ZombiTest/PBanalysis/concatenate2000.phy", "w")

alignments = AlignIO.parse(input_handle, "fasta")
AlignIO.write(alignments, output_handle, "phylip-relaxed")

output_handle.close()
input_handle.close()