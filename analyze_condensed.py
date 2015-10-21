__author__ = 'fsiegris'
# Date -- 18.09.2015
# Title -- Analysis functions -for- statistical analyses of tef scaffold
# matches 'condensed' from CoGe's SynMap on Sorghum chromosomes

import csv

from Bio import SeqIO
from numpy import histogram
import matplotlib.pyplot as plt

from tef_functions_FS2015 import *

# Setting Debug parameters
printLastSingleton = False
plotHist = True

# Loading Data
# Check for console arguments, then UserIDs then ask for input.
(dagfile, output_file, input_fileS, input_fileC) = documentNames()

# Read-in CoGe SynMap condensed file as list of lists.
lol = list(csv.reader(filter(lambda row: row[0] != '#', open(
    dagfile,
    'r'
)), delimiter='\t'))

# Read in the comment lines from condensed CoGe file
comments = list(csv.reader(filter(lambda row: row[0] == '#', open(
    dagfile,
    'r'
)), delimiter='\t'))

# Duplicate each comment row for every entry.
comment = []
block = 1
for m in comments:
    block = block + 1
    m.append(block)
    try:
        for n in range(0, int(m[5])):
            comment.append(m)
    except IndentationError:
        pass  # empty comment first line
    except IndexError:
        pass

# Search scaffolds on DAGchainer document and build up hash links
# to the fasta file nucleotide information.
D = {}
for q in range(0, len(lol)):
    """
    The 'r' entries are extended in reverse order.
    """
    if str(comment[q][4]) == 'f':
        entry = dagList(q, lol, comment)
    elif str(comment[q][4]) == 'r':
        # Read the block from the last entry to the first
        if q == 0 or str(comment[q]) != str(comment[q - 1]):  # new block
            s = int(comment[q][5]) - 1  # start with last
        else:
            s = s - 2  # get to the entry before
        r = q + s  # step back one per one
        entry = dagList(r, lol, comment)
    if str(lol[q][1]).split('||')[0] in D:
        D[(str(lol[q][1]).split('||')[0])].extend(entry)
    else:
        D[(str(lol[q][1]).split('||')[0])] = entry

# Read in all E. tef scaffolds and echos scaffold name and nucleotide sequence
fasta_sequences = SeqIO.parse(open(input_fileS), 'fasta')

# Read in all Sorghum chromosomes and sort chromosomes by number
chr_sequences = SeqIO.parse(open(input_fileC), 'fasta')

# Sort the cromosomes numerically.
chr_sequences = list(sorted(chr_sequences, key=lambda x: int(x.id[3:])))

# Setup:
missmatch = 0
a = []
stretches = []
scaf_fasta = {}

# Loading fasta sequences in SynMap file to dictionary
# and find synthenic Stretch and decide on best one.

print('\nLoading over! Starting Stretch search')

for fasta in fasta_sequences:
    try:
        # Try to find scaffold name in dictionary of CoGe DagChain and
        # build up list of stretches found
        singleton = find_synthenic_block(D[fasta.id], str(fasta.id))
        stretches.append(singleton[0])
        a.append(singleton[1])
        # build up fasta dictionary with known scaffold ids
        scaf_fasta[fasta.id] = fasta
    except KeyError:  # Except for missing keys in CoGe file and count !
        missmatch = missmatch + 1

if printLastSingleton:  # Check class module for sequence information.
    print('Chr ' + singleton[0][0].cfa(chr_sequences))
    print('Sca ' + singleton[0][0].sfa(scaf_fasta[singleton[0][0].scaffold]))

# Report histogram data
print('\nNumber of not matched scaffolds: ' + str(missmatch))
print(histogram(a, bins=[0, 1, 2, 3, 4, 5, 10, 20, 100]))

# Save Stretch-list in a file for R import.
try:
    with open(output_file, 'w') as csvfile:
        stretchwriter = csv.writer(csvfile, delimiter='\t',
                                   quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for stretch in stretches:
            for entry in stretch:
                stretchwriter.writerow(entry.aslist())
except FileNotFoundError:
    print('File not found')

# Plot a histogram of how many stretches have been found on different scaffolds
if plotHist:
    plt.hist(a)
    plt.title("Histogram")
    plt.xlabel("3+ genes stretches in scaffold")
    plt.ylabel("Frequency")
    plt.show()
