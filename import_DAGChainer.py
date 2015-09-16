__author__ = 'fsiegris'
# Date: 07.09.2015

import csv

from Bio import SeqIO, SeqFeature

from tef_functions_FS2015 import *

# Test importing the dag file and print out the first
# significant line.  To run on local and lab computer.

try:
    print('Starting test phase:')
    file = open(
        '../../i1sz/22790_24796.CDS-CDS.last.tdd10.cs0.filtered.dag.all.go',
        'r+'
    )
    print(file.readline())
    print(file.readline())
    file.close()
except NameError:
    print('File not readable')
except FileNotFoundError:
    print('File not located where it should')
except Exception:
    print('I give up on you!')
finally:
    print('Starting proof of principle:')


# Get out the different columns and calculate the difference
# of the scaffold and reference chromosome coordinates.

lol = list(csv.reader(open(
    '../../i1sz/22790_24796.CDS-CDS.last.tdd10.cs0.filtered.dag.all.go',
    'r+'
), delimiter='\t'))
for m in lol[1]:
    print(m)
print('All columns of first line\n\n')
for n in range(1, 10):
    print(lol[n][0]+' '+str(int(lol[n][2])-int(lol[n][6]))+' '+lol[n][4])
    print(str(lol[n][1]).split('||'))

for p in range(1, 10):
    print(str(lol[p][5]).split('||'))

input('\nPaused --- Press Enter to continue')

# Calculate for same Sorghum chromosome gene the distances
##############################################################################
# Line character length should be max 72 comments and 79 code: testing that
##############################################################################


# Build up a hash-table with scaffold name and ( chr, (start, end), orientation )

D = {}
for q in range(1, len(lol)):
    if str(lol[q][1]).split('||')[0] in D:
        D[(str(lol[q][1]).split('||')[0])].extend([(int(str(lol[q][5]).split('||')[0]),(int(str(lol[q][5]).split('||')[1]),int(str(lol[q][5]).split('||')[2])),int(str(lol[q][1]).split('||')[4]))])
    else:
        D[(str(lol[q][1]).split('||')[0])] = [(int(str(lol[q][5]).split('||')[0]),(int(str(lol[q][5]).split('||')[1]),int(str(lol[q][5]).split('||')[2])),int(str(lol[q][1]).split('||')[4]))]

print(D['scaffold14093'])
print( "Value : %s" %  D.keys())


input('\nPaused --- Press Enter to continue')


# OK, I understand I should use english, but I have to violate this rule,
# from time to time.
print('\n\nTest phase over! Starting real life: Раз Два Три\n\n')


# Here is a code-chunk to easely import fasta files

input_file='/windows/GNY98ter_41.closed'
output_file='/windows/GNY98.pyout'
fasta_sequences = SeqIO.parse(open(input_file),'fasta')
"""with open(output_file) as out_file:
    for fasta in fasta_sequences:
        #name, sequence = fasta.id, fasta.seq.tostring()
        #new_sequence = some_function(sequence)
        #write_fasta(out_file)
        print(fasta.id+' '+fasta.seq.tostring())"""

"""
Read in all E. tef scaffolds and echos scaffold name and nucleotide sequence

"""


# Search scaffolds on DAGchainer document and build up hash links to the fasta file nucleotide information
from numpy import histogram
missmatch = 0
a = []
for fasta in fasta_sequences:
    try:
        #print(' '+fasta.id[3:]+' '+str(D[fasta.id[3:]]))
        a.append(find_synthenic_block(D[fasta.id[3:]],str(fasta.id[3:])))
    except:
        missmatch=missmatch+1 #print(fasta.id[3:]+' '+'Not found.')
    finally:
         pass
print('\nNumber of not matched scaffolds: '+str(missmatch))
print(histogram(a))

# Now this is amazing how many synthenic stretches we find on a single scaffold,
# we have to find now the best one (longest-gene stretch on smallest space ?


