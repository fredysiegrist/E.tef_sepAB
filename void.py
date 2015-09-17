__author__ = 'fsiegris'
class Etef:
    """ Object representing a E. tef scaffold """

    def __init__(self, width, height,
                 color='black', emphasis=None, highlight=0):
        """Initiates the class Etef as example,
        though tef is a plant and not a class.
        """
        __version__ = "$Revision$"
        # $Source$
        if (width == 0 and height == 0 and
                    color == 'red' and emphasis == 'strong' or
                    highlight > 100):
            raise ValueError("sorry, you lose")
        if width == 0 and height == 0 and (color == 'red' or
                                                   emphasis is None):
            raise ValueError("I don't think so -- values are %s, %s" %
                             (width, height))
        Etef.__init__(self, width, height,
                      color, emphasis, highlight)


# Here are some classes from "https://www.biostars.org/p/710/" as examples

class Dna:
    ''' Object representing a FASTA record. '''

    def __init__(self, header, sequence):
        self.head = header
        self.seq = sequence

    def __repr__(self):
        return '[HTML]' % (self.head)

    def __str__(self, separator=''):
        return '>%s\n%s' % (self.head, separator.join(self.seq))

    def __len__(self):
        return len(''.join(self.seq))

    @property
    def sequence(self, separator=''):
        return separator.join(self.seq)


class Fasta:
    ''' A FASTA iterator/generates DNA objects. '''

    def __init__(self, handle):
        self.handle = handle

    def __repr__(self):
        return '[HTML]' % handle

    def __iter__(self):
        header, sequence = '', []
        for line in self.handle:
            if line[0] == '>':
                if sequence: yield Dna(header, sequence)
                header = line[1:-1]
                sequence = []
            else:
                sequence.append(line.strip())
        yield Dna(header, sequence)


class FastaSeq:
    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence


def get_seqs(file):
    items = []
    index = 0
    for line in file:
        if line.startswith(">"):
            if index >= 1:
                items.append(aninstance)
            index += 1
            name = line[:-1]
            seq = ''
            aninstance = FastaSeq(name, seq)
        else:
            seq += line[:-1]
            aninstance = FastaSeq(name, seq)
    items.append(aninstance)

    return items



"""with open(output_file) as out_file:
    for fasta in fasta_sequences:
        #name, sequence = fasta.id, fasta.seq.tostring()
        #new_sequence = some_function(sequence)
        #write_fasta(out_file)
        print(fasta.id+' '+fasta.seq.tostring())"""
import csv

from Bio import SeqIO, SeqFeature

from tef_functions_FS2015 import *

lol = list(csv.reader(filter(lambda row: row[0]!='#', open(
    '../../i1sz/22790_24796.CDS-CDS.last.tdd10.cs0.filtered.dag.all.go_D20_g10_A3.aligncoords.gcoords.condensed',
    'r+'
)), delimiter='\t'))

D = {}
for q in range(1, len(lol)):
    entry = [SynMap(int(str(lol[q][5]).split('||')[0]),int(str(lol[q][5]).split('||')[1]),int(str(lol[q][5]).split('||')[2]),int(str(lol[q][1]).split('||')[4]))]
    if str(lol[q][1]).split('||')[0] in D:
        D[(str(lol[q][1]).split('||')[0])].extend(entry)
    else:
        D[(str(lol[q][1]).split('||')[0])] = entry
print(D['scaffold100'][0])
from operator import itemgetter, attrgetter, methodcaller

print(sorted(D['scaffold105'], key=attrgetter('chr','start','end')))


input('\nPaused --- Press Enter to continue')

input_file='../../i1sz/GNYt98ter.41.closedgt1000.sorted'
output_file='../../i1sz/GNY98.pyout'
fasta_sequences = SeqIO.parse(open(input_file),'fasta')



for fasta in fasta_sequences:
    find_synthenic_block(D[fasta.id],str(fasta.id))
