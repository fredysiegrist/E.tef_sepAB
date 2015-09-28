__author__ = 'fsiegris'
# Date -- 18.09.2015
# Title -- Analysis functions -for- statistical analyses of tef scaffold
# matches 'condensed' from CoGe's SynMap on Sorghum chromosomes

import csv, pwd, os, sys

from Bio import SeqIO, SeqFeature
import matplotlib.pyplot as plt

from tef_functions_FS2015 import *

try:
    dagfile = sys.argv[1]
except:
    if os.getuid() == 1000 or pwd.getpwuid(
            os.getuid()).pw_gecos == 'fsiegris':
        dagfile = '../../i1sz/22790_24796.CDS-CDS.last.tdd10.cs0.filtered.dag.all.go_D20_g10_A3.aligncoords.gcoords.condensed'
    elif pwd.getpwuid(
            os.getuid()).pw_gecos == 'Fredy Siegrist,,,' or pwd.getpwuid(
        os.getuid()).pw_gecos == 'Gina Cannarozzi,,,':
        dagfile = '/home/fredy/i1sz/22790_24796.CDS-CDS.last.tdd10.cs0.filtered.dag.all.go_D20_g10_A3.aligncoords.gcoords.condensed'
    else:
        dagfile = input(
            'Please enter directory and file of \
            *.dag.all.go_D20_g10_A3.aligncoords.gcoords.condensed file')
lol = list(csv.reader(filter(lambda row: row[0] != '#', open(
    dagfile,
    'r'
)), delimiter='\t'))

comments = list(csv.reader(filter(lambda row: row[0] == '#', open(
    dagfile,
    'r'
)), delimiter='\t'))

comment = []
block = 1
for m in comments:
    block = block + 1
    m.append(block)
    try:
        for n in range(0, int(m[5])):
            comment.append(m)
    except IndentationError:
        pass    # empty comment first line
    except IndexError:
        pass

# print(str(len(comment))+' '+str(len(comments))+' '+str(len(lol)))


D = {}
for q in range(0, len(lol)):
    """
    The 'r' entries are extended in reverse order.
    """
    if str(comment[q][4]) == 'f':
        entry = [DagChain(
            int(str(lol[q][5]).split('||')[0]),  # chr
            int(str(lol[q][5]).split('||')[1]),  # start
            int(str(lol[q][5]).split('||')[2]),  # end
            int(str(lol[q][5]).split('||')[4]),  # ori
            str(lol[q][1]).split('||')[0],  # scaffold
            int(str(lol[q][1]).split('||')[1]),  # sstart
            int(str(lol[q][1]).split('||')[2]),  # send
            int(str(lol[q][1]).split('||')[4]),  # sori
            int(comment[q][5]),  # stretch
            int(str(comment[q][0])[1]),  # number
            float(comment[q][1]),  # score
            str(comment[q][4]),  # reversion
            int(comment[q][6]) # block
        )
        ]
    elif str(comment[q][4]) == 'r':
        if q == 0 or str(comment[q]) != str(comment[q - 1]):  # new block
            s = int(comment[q][5]) - 1
        else:
            s = s - 2
        # print(str(s)+'\t'+'\t'.join(comment[q]))
        r = q + s
        entry = [DagChain(
            int(str(lol[r][5]).split('||')[0]),  # chr
            int(str(lol[r][5]).split('||')[1]),  # start
            int(str(lol[r][5]).split('||')[2]),  # end
            int(str(lol[r][5]).split('||')[4]),  # ori
            str(lol[r][1]).split('||')[0],  # scaffold
            int(str(lol[r][1]).split('||')[1]),  # sstart
            int(str(lol[r][1]).split('||')[2]),  # send
            int(str(lol[r][1]).split('||')[4]),  # sori
            int(comment[q][5]),  # stretch
            int(str(comment[q][0])[1]),  # number
            float(comment[q][1]),  # score
            str(comment[q][4]),  # reversion
            int(comment[q][6]) # block
        )
        ]
    if str(lol[q][1]).split('||')[0] in D:
        D[(str(lol[q][1]).split('||')[0])].extend(entry)
    else:
        D[(str(lol[q][1]).split('||')[0])] = entry
# print(D['scaffold105'])
# input()
# print(s.encode(pwd.getpwuid(os.getuid()).pw_gecos))

# Here is a code-chunk to easely import fasta files dependent on current user
try:
    input_file = sys.argv[2]
except:
    if os.getuid() == 1000 or pwd.getpwuid(
            os.getuid()).pw_gecos == 'fsiegris':
        input_file = '/windows/GNYt98ter.41.closedgt1000.sorted'
    elif pwd.getpwuid(
            os.getuid()).pw_gecos == 'Fredy Siegrist,,,' or pwd.getpwuid(
        os.getuid()).pw_gecos == 'Gina Cannarozzi,,,':
        input_file = '/home/fredy/i1sz/GNYt98ter.41.closedgt1000.sorted'
    else:
        input_file = input(
            'Please enter directory and file of GNYt98ter.41.closedgt1000.sorted file')


# Read in all E. tef scaffolds and echos scaffold name and nucleotide sequence
fasta_sequences = SeqIO.parse(open(input_file), 'fasta')

# Read in all Sorghum chromosomes and sort chromosomes by number
try:
    input_file = sys.argv[3]
except:
    if os.getuid() == 1000 or pwd.getpwuid(
            os.getuid()).pw_gecos == 'fsiegris':
        input_file = '../../i1sz/Sorghum_bicolor.faa'
    elif pwd.getpwuid(
            os.getuid()).pw_gecos == 'Fredy Siegrist,,,' or pwd.getpwuid(
        os.getuid()).pw_gecos == 'Gina Cannarozzi,,,':
        input_file = '/home/fredy/i1sz/Sorghum_bicolor.faa'
    else:
        input_file = input(
            'Please enter directory and file of Sorghum_bicolor.faa file')
chr_sequences = SeqIO.parse(open(input_file), 'fasta')
# Sort the cromosomes numerically.
chr_sequences = list(sorted(chr_sequences, key=lambda x: int(x.id[3:])))


# Search scaffolds on DAGchainer document and build up hash links
# to the fasta file nucleotide information
from numpy import histogram

missmatch = 0
a = []
stretches = []

print('\nTest phase over! Starting real thing: One two three')

for fasta in fasta_sequences:
    try:
        singleton = find_synthenic_block(D[fasta.id], fasta, chr_sequences)
        stretches.append(singleton[0])
        a.append(singleton[1])
    except KeyError:  # Only except for missing keys!
        missmatch = missmatch + 1

print('\nNumber of not matched scaffolds: ' + str(missmatch))
print(histogram(a, bins=[0, 1, 2, 3, 4, 5, 10, 20, 100]))

# Save Stretch-list in a file for R import.
output_file = '../../i1sz/stretches_condensed.csv'
try:
    with open(output_file, 'w') as csvfile:
        stretchwriter = csv.writer(csvfile, delimiter='\t',
                                   quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for stretch in stretches:
            for entry in stretch:
                stret = [str(entry.chr), str(entry.start), str(entry.end),
                         str(entry.ori),
                         str(entry.scaffold), str(entry.sstart),
                         str(entry.send),
                         str(entry.sori), str(entry.stretch),
                         str(entry.number),
                         str(entry.score), str(entry.reversion),
                         str(entry.genes_in_row), str(entry.block), str(entry.position)]#,
                         #str(entry.cfa), str(entry.sfa)]
                stretchwriter.writerow(stret)
                # stretchwriter.writerow("\# "+str(entry[0].scafname))
except FileNotFoundError:
    print('File not found')

# Plot a histogram of how many stretches have been found on different scaffolds

plt.hist(a)
plt.title("Histogram")
plt.xlabel("3+ genes stretches in scaffold")
plt.ylabel("Frequency")
plt.show()
