## This script is made for transorming several fasta files (containing aligned sequences) into one maf file where each chunks is associated to an artificial chromosome and artificial coordinates (each chunk must be separated by more than 100 bp)
## It takes as inout a list of fasta files
## It outputs one single maf file.

import  sys
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio import AlignIO
from Bio.AlignIO import MafIO
import math

if len(sys.argv)!=3:
	sys.exit("usage: python from_fasta_to_maf.py <list of fasta (50)> <branchID (for output file name>")

with open(sys.argv[1]) as f:
    fasta_list = f.read().splitlines() 

branchID=sys.argv[2]

for t in range(0,49,10):
	g=math.floor((t/10)+1)
	# Code for 1 maf output. Loop to get five!
	alignments=[]
	for j in range(t,t+10):
		alignments.append(AlignIO.read(fasta_list[j], "fasta"))
	for ali in alignments:
		for record in ali:
			record.id = "species" + record.id + ".chrm" + str(g)
		#print(record.id) 
	AlignIO.write(alignments, "my_exa.maf", "maf")

	# now change coordinates...
	alignmentsmaf = list(AlignIO.parse("my_exa.maf", "maf"))
	for i in range(0,10):
		#for multiple_alignment in AlignIO.parse("my_exa.maf", "maf"):
		starting_coord = i*2000000
		for seqrec in alignmentsmaf[i]:
			seqrec.annotations["start"] = starting_coord
	AlignIO.write(alignmentsmaf, 'alignment_{}_{}.maf'.format(g, branchID), "maf")

