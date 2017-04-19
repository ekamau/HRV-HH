#!/usr/local/bin/python

# Compares all sequences in aligned fasta file to first sequence in set
# .. determines number of sequences to compare, adjusts accordingly 
# Adopted from Dr. Matt Cotten, modified by Everlyn Kamau (18-04-2017)
# duplicate the last sequence in the alignment.

#1. import the required modules
import sys, os, csv, time
import os.path
import matplotlib.pyplot as plt
import numpy as np
sys.path.insert(0, "/Applications/biopython-1.62b")

from Bio import SeqIO
from Bio.Seq import Seqfrom Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import Entrez

#2. Ensure that the necessary input files are fed into the program

if len(sys.argv)!= 2:
	print("Usage: python HiLiter.py Sequences.fasta")
	sys.exit()
print "Running HiLiterB.py"

All_seqs= sys.argv[1]
outprefix = os.path.splitext(All_seqs)[0]

#3. Start an empty list (list comprenension), parse the sequence file into the program

seq_names =[]
sequences =[]
for record in SeqIO.parse(open(All_seqs, "rU"), "fasta"):
	seq_names.append(record.id)
	sequence_string = str(record.seq)	
	sequences.append(sequence_string)
total_number_sequences = len(sequences)
genome_length = len(sequences[0])

graph_names=[]
graph_names.append(seq_names[1])
graph_names.append(seq_names[2])

#4. Generate lists of differences between ref and test sequences
#5. function for comparing sequence a to b
def vergleichen(seqA, seqB, index):	
	for i in range (len(seqA)):
		if seqB[i]!=seqA[i] and seqB[i] == "A":
			color = "Green"
			line = index
			diff_list.append((i,line,color))
		elif seqB[i]!=seqA[i] and seqB[i] == "T":
			color = "Red"
			line = index
			diff_list.append((i,line,color))
		elif seqB[i]!=seqA[i] and seqB[i] == "G":
			color = "Black"
			line = index
			diff_list.append((i,line,color))
		elif seqB[i]!=seqA[i] and seqB[i] == "C":
			color = "DarkBlue"
			line = index
			diff_list.append((i,line,color))
		elif seqB[i]!=seqA[i] and seqB[i] == "N":
			color = "DarkGrey"
			line = index
			diff_list.append((i,line,color))
		elif seqB[i]!=seqA[i]:
			color = "LightGrey"	
			line = index
			diff_list.append((i,line,color))	
		else:
			continue

diff_list=[]

#6. Call vergleichen on each pair

for i in range(1, len(sequences)):
	vergleichen(sequences[0], sequences[i], i)

print diff_list

#7. Now generate graph of differences. 
ticks=[] 
for x in range(len(sequences)-1):
	ticks.append(x)
#ticks=["1","2","3"] 
fig = plt.figure()
#fig.subplots_adjust(left=0.1, bottom=None, right=None, wspace=None, top = None)
ax2 = fig.add_subplot(211)
ax2.set_ylim(0.9, len(sequences)-1) 
ax2.set_xlim(0,(genome_length))

ax2.set_yticks(ticks)
ax2.set_yticklabels(seq_names, size = 5)
plt.grid(color='silver', linewidth=0.4)

# Remove box
#ax2.set_xlabel("Genome Position of "+seq_names[0]+", also acting as reference", size = 8)
ax2.set_xlabel("Nucleotide Position of "+seq_names[0], size = 8)

ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)

# Remove_ticks
ax2.tick_params(axis='both', direction='out',size = 4)
ax2.get_xaxis().tick_bottom()   # remove unneeded ticks 
ax2.get_yaxis().tick_left()

# xtick label size
ax2.tick_params(axis='x', labelsize=8) 

for i in range(len(diff_list)):
	color = diff_list[i][2]
	ax2.broken_barh([(int(diff_list[i][0]), 1)] , (int(diff_list[i][1]), 0.75), facecolors=(color),edgecolors=('None')) 

#plt.savefig(outprefix+'_differences.png', dpi=250)

#Can save image as pdf
plt.savefig(outprefix+'_differences.pdf')

print "That's All Folks!"
