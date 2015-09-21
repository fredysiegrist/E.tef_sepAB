__author__ = 'fsiegris'
# Date: 07.09.2015

import csv, os, pwd

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
    print(file.readline())
    file.close()
except NameError:
    print('File not readable')
except FileNotFoundError:
    print('File not located where it should')
except Exception:
    print('I give up on you!')
finally:
    print('Starting integrity test of input data file:')


# Get out the different columns and calculate the difference
# of the scaffold and reference chromosome coordinates.

lol = list(csv.reader(filter(lambda row: row[0] != '#', open(
    '../../i1sz/22790_24796.CDS-CDS.last.tdd10.cs0.filtered.dag.all.go',
    'r+'
)), delimiter='\t'))
for m in lol[1]:
    print(m)
print('All columns of first line\n\n')
for n in range(1, 10):
    print(
        lol[n][0] + ' ' + str(int(lol[n][2]) - int(lol[n][6])) + ' ' + lol[n][
            4])
    print(str(lol[n][1]).split('||'))

for p in range(1, 10):
    print(str(lol[p][5]).split('||'))

input('\nPaused --- Press Enter to continue')

# Calculate for same Sorghum chromosome gene the distances
##############################################################################
# Line character length should be max 72 comments and 79 code: testing that
##############################################################################


# Build up a hash-table with scaffold name and (chr,(start,end),orientation)

D = {}
for q in range(1, len(lol)):
    entry = [SynMap(int(str(lol[q][5]).split('||')[0]),
                    int(str(lol[q][5]).split('||')[1]),
                    int(str(lol[q][5]).split('||')[2]),
                    int(str(lol[q][1]).split('||')[4]))]
    if str(lol[q][1]).split('||')[0] in D:
        D[(str(lol[q][1]).split('||')[0])].extend(entry)
    else:
        D[(str(lol[q][1]).split('||')[0])] = entry

# print(D['scaffold14093'])
print("Value : %s" % D.keys())

input('\nPaused --- Press Enter to continue')


# OK, I understand I should use english, but I have to violate this rule,
# from time to time.
print('\n\nTest phase over! Starting real life: One Two Tree\n\n')


# Here is a code-chunk to easely import fasta files
if pwd.getpwuid(os.getuid()).pw_gecos == 'fredi' or pwd.getpwuid(
        os.getuid()).pw_gecos == 'fsiegris':
    input_file = '/windows/GNYt98ter.41.closedgt1000.sorted'
elif pwd.getpwuid(os.getuid()).pw_gecos == 'fredy' or pwd.getpwuid(
        os.getuid()).pw_gecos == 'cannarozi':
    input_file = '../../i1sz/GNYt98ter.41.closedgt1000.sorted'
else:
    input_file = input(
        'Please enter directory and file of GNYt98ter.41.closedgt1000.sorted file')

output_file = '../../i1sz/GNY98.pyout'
fasta_sequences = SeqIO.parse(open(input_file), 'fasta')

# Read in all E. tef scaffolds and echos scaffold name and nucleotide sequence


# Search scaffolds on DAGchainer document and build up hash links to the fasta file nucleotide information
from numpy import histogram

missmatch = 0
a = []
for fasta in fasta_sequences:
    try:
        # print(' '+fasta.id[3:])#+' '+str(D[fasta.id[3:]]))
        a.append(find_synthenic_block(D[fasta.id], fasta))
    except:
        missmatch = missmatch + 1  # print(fasta.id[3:]+' '+'Not found.')
    finally:
        pass
print('\nNumber of not matched scaffolds: ' + str(missmatch))
print(histogram(a, bins=[0, 1, 2, 3, 4, 5, 10, 20, 100]))

import matplotlib.pyplot as plt

plt.hist(a)
plt.title("Histogram")
plt.xlabel("3+ genes stretches in scaffold")
plt.ylabel("Frequency")
plt.show()


# Now this is amazing how many synthenic stretches we find on a single scaffold,
# we have to find now the best one (longest-gene stretch on smallest space ?
