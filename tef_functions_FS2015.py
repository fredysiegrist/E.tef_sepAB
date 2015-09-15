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
            index+=1
            name = line[:-1]
            seq = ''
            aninstance = FastaSeq(name, seq)
        else:
            seq += line[:-1]
            aninstance = FastaSeq(name, seq)
    items.append(aninstance)

    return items