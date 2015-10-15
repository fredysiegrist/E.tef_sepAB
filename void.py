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







                genes_in_row = genes_in_row + 1
                if (genes_in_row == 2):
                    if (entry_old.start == 0):
                        cordstart = entry.start
                        if typecheck:
                            scafstart = entry.sstart
                    else:
                        cordstart = entry_old.start
                        if typecheck:
                            scafstart = entry_old.sstart
                        else:
                            scafstart = 0





            # if ((entry.start - entry_old.end)<0):
            # You will kill me for the next statement :-)
            if False:  # typecheck:
                orient = (
                    entry.ori == entry.sori and entry_old.ori == entry.ori or
                    entry.ori != entry.sori and entry_old.ori == entry.ori
                )
            else:


"""
     chr    start      end ori     scaffold sstart   send sori stretch number
1732   2  7423978 71237703   1 scaffold6518  26669 132614   -1       3      1
4045   1 27532747 68908460   1 scaffold4243   9731  29813   -1       3      1
     score reversion gir block position
1732   144         r   3  3775        2
4045   150         r   3  2637        1
"""



findmaxpos <- function(x) {
    if (max(cond[cond[,5] %in% slimmaster[x,2],15])>actualpos[x]) {
        actualpos[x]=actualpos[x]+1
    }
    else {
        actualpos[x]=1
    }
    outlier[x] <- outlierpos[outlierpos[,15]==actualpos[x],]
}



betterlist[match(outlier$scaffold,  betterlist$scaffold),] <- outlier
betterlist <- betterlist[order((betterlist[,1]*1e10)+betterlist[,2]+betterlist[,3]),]


betterremaster <- betterlist[,c(1, 5, 12)]
newcon <- match(slimmaster$CHR2, betterremaster$scaffold, nomatch=NULL)


# Damn thing won't work: which are  the connec values of updated scaffolds, or the newcon positio of them
updated <- (1:length(connec))[abs(newcon - connec) > 100]  # match(outlier$scaffold, cond$scaffold)
plot(1:length(connec), (1:length(connec))-newcon, pch=16, xlab="ordinal # in better master", ylab="position difference", cex=1, main="ordering in master and 'better' master file", col=colors()[as.numeric(betterremaster$chr)+98])
legend(0, 3, unique(betterremaster$chr), col=colors()[99:108], pch=16, ncol=5)
points(updated, updated-newcon[updated], pch="°", col="green", cex=2)
