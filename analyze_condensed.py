__author__ = 'fsiegris'
# Date -- 18.09.2015
# Title -- Analysis functions -for- statistical analyses of tef scaffold matches
# 'condensed' from CoGe's SynMap on Sorghum chromosomes

import csv

from Bio import SeqIO, SeqFeature
import matplotlib.pyplot as plt

from tef_functions_FS2015 import *

lol = list(csv.reader(filter(lambda row: row[0]!='#', open(
    '../../i1sz/22790_24796.CDS-CDS.last.tdd10.cs0.filtered.dag.all.go_D20_g10_A3.aligncoords.gcoords.condensed',
    'r+'
)), delimiter='\t'))

comments = list(csv.reader(filter(lambda row: row[0]=='#', open(
    '../../i1sz/22790_24796.CDS-CDS.last.tdd10.cs0.filtered.dag.all.go_D20_g10_A3.aligncoords.gcoords.condensed',
    'r+'
)), delimiter='\t'))

comment = []

for m in comments:
    try:
        for n in range(0,int(m[5])):
            #print(m)
            comment.append(m)
    except:
        pass

# print(str(len(comment))+' '+str(len(comments))+' '+str(len(lol)))


D = {}
for q in range(1, len(lol)):
    """
    The 'r' entries in comment[q][4] should get extended in reverse order!
    """
    entry = [DagChain(
                      int(str(lol[q][5]).split('||')[0]), #chr
                      int(str(lol[q][5]).split('||')[1]), #start
                      int(str(lol[q][5]).split('||')[2]), #end
                      int(str(lol[q][1]).split('||')[4]), #ori
                      str(lol[q][1]).split('||')[0], # scaffold
                      int(str(lol[q][1]).split('||')[1]), # sstart
                      int(str(lol[q][5]).split('||')[2]), # send
                      int(str(lol[q][5]).split('||')[4]), # sori
                      int(comment[q][5]), # stretch
                      int(str(comment[q][0])[1]), # number
                      float(comment[q][1]) # score
                      )
            ]
    if str(lol[q][1]).split('||')[0] in D:
        D[(str(lol[q][1]).split('||')[0])].extend(entry)
    else:
        D[(str(lol[q][1]).split('||')[0])] = entry
print(D['scaffold105'])

# Here is a code-chunk to easely import fasta files
input_file='../../i1sz/GNYt98ter.41.closedgt1000.sorted'
output_file='../../i1sz/GNY98.pyout'
# Read in all E. tef scaffolds and echos scaffold name and nucleotide sequence
fasta_sequences = SeqIO.parse(open(input_file),'fasta')

# Search scaffolds on DAGchainer document and build up hash links to the fasta file nucleotide information
from numpy import histogram
missmatch = 0
a = []
# OK, I understand I should use english, but I have to violate this rule,
# from time to time.
print('\nTest phase over! Starting real thing: Раз Два Три')
input('\nPaused --- Press Enter to continue\n\n')
for fasta in fasta_sequences:
    try:
        a.append(find_synthenic_block(D[fasta.id],str(fasta.id)))
    except:
        missmatch=missmatch+1
    finally:
         pass
print('\nNumber of not matched scaffolds: '+str(missmatch))
print(histogram(a, bins=[0,1,2,3,4,5,10,20,100]))

# Plot a histogram of how many stretches have been found on different scaffolds
plt.hist(a)
plt.title("Histogram")
plt.xlabel("3+ genes stretches in scaffold")
plt.ylabel("Frequency")
plt.show()


