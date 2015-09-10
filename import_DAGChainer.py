__author__ = 'fsiegris'
# Date: 07.09.2015

# from recipe-577444-1 import *

import csv

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

# OK, I understand I should use english, but I have to violate this rule,
# from time to time.
print('\n\nTest phase over! Starting real life: Раз Два Три\n\n')

# Calculate for same Sorghum chromosome gene the distances
##############################################################################
# Line character length should be max 72 comments and 79 code: testing that
##############################################################################


"""
# Here is a code-chunk to easely import fasta files
from Bio import SeqIO

fasta_sequences = SeqIO.parse(open(input_file),'fasta')
with open(output_file) as out_file:
    for fasta in fasta_sequences:
        name, sequence = fasta.id, fasta.seq.tostring()
        new_sequence = some_function(sequence)
        write_fasta(out_file)
"""