__author__ = 'fsiegris'


# date: 09.09.2015
# function source file for tef practical

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


def find_synthenic_block(coordlist, scafname, D=120000, mingen=3, plot=0):
    import matplotlib.pyplot as plt
    # from time import sleep
    # Maximum distance between two matches (-D): plant-default 120000 bp
    #D = 120000
    entry_old = (0, (0, 0), 0)
    genes_in_row = 0
    synthenic_hits = 0
    cordstart = 0
    found_stretch = []
    for entry in sorted(coordlist):
        if (entry_old[0] == 0 or entry[0] == entry_old[0]):  # check if on same chromosome
            # print(str(entry[1][0])+'-'+str(entry_old[1][1])+'*'+str(entry[2])+'*'+str(entry_old[2])+' '+str(cordstart)+' Decision: '+str(abs((entry[1][0]-entry_old[1][1])*entry[2]*entry_old[2])<=D))
            if (abs((entry[1][0] - entry_old[1][1]) * entry[2] * entry_old[2]) <= D):
                genes_in_row = genes_in_row + 1
                if (genes_in_row == 1):
                    if (entry_old[1][0] == 0):
                        cordstart = entry[1][0]
                    else:
                        cordstart = entry_old[1][0]
            else:
                synthenic_hits = __found_synthenic(synthenic_hits, genes_in_row, entry, entry_old, cordstart, scafname, mingen)
                found_stretch.append(genes_in_row)
                genes_in_row = 0
        else:
            synthenic_hits = __found_synthenic(synthenic_hits, genes_in_row, entry, entry_old, cordstart, scafname, mingen)
            found_stretch.append(genes_in_row)
            genes_in_row = 0
            # print('chromosome jump '+str(entry_old[0])+' '+str(entry[0]))
            cordstart = 0
        entry_old = entry
    if plot:
        plt.hist(found_stretch)
        plt.title("Histogram of "+scafname)
        plt.xlabel("genes in row on 3+stretches")
        plt.ylabel("Frequency")
        plt.show()
    return (synthenic_hits)


def __found_synthenic(synthenic_hits, genes_in_row, entry, entry_old, cordstart, scafname, mingen):
    if (genes_in_row >= mingen):
        print('got ' + str(genes_in_row) + '-genes stretch of synthenic genes on chromosome ' + str(
            entry_old[0]) + ' coordinates start: ' + str(cordstart) + ' end: ' + str(
            entry_old[1][1]) + ' on ' + scafname)
        # sleep(1)
        synthenic_hits = synthenic_hits + 1
    return (synthenic_hits)
