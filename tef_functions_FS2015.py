__author__ = 'fsiegris'


# Date: 09.09.2015
# Functions for tef practical.

class DagChain:
    """
        Class is designed to take additional information as coordinates on
        chromosome and scaffold and block identity information and can be
        supplemented with nucleotide sequence extracted from fasta files for
        calculation of evolutionary distance with PAML or RSD (on proteins?)
    """

    def __init__(self, chr, start, end, ori, scaffold, sstart, send, sori,
                 stretch, number, score, reversion, block):
        self.chr = chr
        self.start = start
        self.end = end
        self.ori = ori
        self.coord = (start, end)
        self.scaffold = scaffold
        self.sstart = sstart
        self.send = send
        self.scoord = (sstart, send)
        self.sori = sori
        self.stretch = stretch
        self.number = number
        self.score = score
        self.reversion = reversion
        self.block = block

    def chr(self):
        return self.chr

    def start(self):
        return self.start

    def end(self):
        return self.end

    def ori(self):
        return self.ori

    def length(self):
        return int((self.end - self.start) * self.ori)

    def scaffold(self):
        return self.scaffold

    def sstart(self):
        return self.sstart

    def send(self):
        return self.send

    def sori(self):
        return self.sori

    def slength(self):
        return int((self.send - self.sstart) * self.sori)

    def stretch(self):
        return self.stretch

    def number(self):
        return self.number

    def score(self):
        return self.score

    def reversion(self):
        return self.reversion

    def block(self):
        return self.block

    def all(self):
        return (
            self.chr, self.start, self.end, self.ori, self.scaffold,
            self.sstart,
            self.send, self.scoord, self.sori, self.stretch, self.number,
            self.score, self.block)

    def __repr__(self):
        return repr((self.chr, (self.start, self.end), self.ori))

    def __lt__(self, other):
        (self.chr, self.start, self.end, self.sstart, self.send) < (
            other.chr, other.start, other.end, other.sstart, self.send)


class Stretch(DagChain):
    """
    A class for returning the values found for the best (and additional)
    nth best stretches of genes found on both scaffolds and reference genome
    (chromosomes).
    """

    def __init__(self, dag, genes_in_row, newblock, position):
        DagChain.__init__(self, dag.chr, dag.start, dag.end, dag.ori,
                          dag.scaffold, dag.sstart, dag.send,
                          dag.sori, dag.stretch, dag.number,
                          dag.score, dag.reversion, dag.block)
        self.genes_in_row = genes_in_row
        self.newblock = newblock
        self.position = position
        self.len = self.end - self.start
        self.slen = self.send - self.sstart

    def gir(self):
        return self.genes_in_row

    def genes_in_row(self):
        return self.genes_in_row

    def newblock(self):
        return self.newblock

    def position(self):
        return self.position

    def cfa(self, chrfastafile):
        """
        :param chrfastafile: SeqIO object with fasta information of reference genome.
        """
        return str(
            chrfastafile[int(self.chr - 1)].seq[self.start - 1:self.end - 1])

    def sfa(self, scafastafile):
        """
        :param scafastafile: SeqIO object with fasta information of scaffold.
        """
        return (str(scafastafile.seq[self.sstart - 1:self.send - 1]))

    def all(self):
        return (
            self.chr, self.start, self.end, self.ori, self.scaffold,
            self.sstart,
            self.send, self.sori, self.stretch, self.number, self.score,
            self.genes_in_row, self.newblock, self.position, self.cfa,
            self.sfa)

    def __repr__(self):
        return repr((self.position, self.block, self.genes_in_row, self.score,
                     self.len,
                     self.slen, self.chr, self.coord, self.scaffold,
                     self.scoord))
    def aslist(self):
        return [str(self.chr), str(self.start), str(self.end),
                      str(self.ori),
                      str(self.scaffold), str(self.sstart),
                      str(self.send),
                      str(self.sori), str(self.stretch),
                      str(self.number),
                      str(self.score), str(self.reversion),
                      str(self.genes_in_row), str(self.block),
                      str(self.position)]

    def str(self):
        return ('got ' + str(
            self.genes_in_row) +
                '-genes stretch of synthenic genes on chromosome ' + str(
            self.chr) + ' coordinates start: ' + str(
            self.start) + ' end: ' + str(
            self.end) + ' on ' + self.scaffold + ' in ' + str(self.sori) +
                str(self.ori) + ' orientation:\n' + str(
            self.cfa) + '   ' + str(self.sfa))

    def __lt__(self, other):
        (self.genes_in_row, self.score, self.end - self.start,
         self.send - self.sstart, self.coord, self.scoord) < (
            other.genes_in_row, other.score, other.end - other.start,
            other.send - other.sstart, other.coord, other.scoord)


def find_synthenic_block(coordlist, scafname, D=12000000, mingen=3,
                         plot=False):
    """
    Checkes the nucleotide distance of the blocks and applies minimal genes
    in a syngenic stretch cutoff to a list of DagChain objects.
    Returns the sorted (by DAG score OR genes_in_row) list of Stretch objects.
    Maximum distance between two matches (-D): plant-default 120000 bp
    100XD = 12 000 000 applied to accept alternative by gene distance
    criteria used in generation of condensed list.
    :param coordlist: List of all DagChain objects for one scaffold.
    :param D: Integer for maximal distance between two genes for cut-off.
    :param mingen: Integer for minimal genes to form a stretch.
    :param plot: Boolean if histogram of number of genes_in_row should be
                 plotted for every entry - to be used for subsets.
    :return: Tuple of list of sorted Stretch objects and number of geq mingen
             stretches for statistical reporting.
    """

    if plot:
        import matplotlib.pyplot as plt

    # Initiation of looping variables for new Stretch attributes.
    entry_old = DagChain(0, 0, 0, 0, '', 0, 0, 0, 0, 0, 0, 'f', 0)
    block = 1
    genes_in_row = 1
    synthenic_hits = []
    cordstart = 0
    scafstart = 0
    found_stretch = []
    ecount = 0
    for entry in coordlist:
        ecount = ecount + 1
        chcheck = entry.chr == entry_old.chr and entry_old.number == entry.number and entry_old.score == entry.score and entry_old.block == entry.block
        if (  # Checks if on same chromosome or block.
              (entry_old.chr == 0 or (chcheck))):
            if (  # Checks for distance between two mapping entries.
                  abs(entry.start - entry_old.end) <= D):
                genes_in_row = genes_in_row + 1
                if (genes_in_row == 2):
                    # Sets the start coordinates for the Stretch.
                    if (entry_old.start == 0):
                        cordstart = entry.start
                        if entry.reversion == 'f':
                            scafstart = entry.sstart
                        else:
                            scafstart = entry.send
                    else:
                        cordstart = entry_old.start
                        if entry.reversion == 'f':
                            scafstart = entry_old.sstart
                        else:
                            scafstart = entry_old.send
            elif (  # Checks if there are enough genes in row for a good
                    # stretch if distance criteria is no longer met.
                    entry != entry_old or entry.block != entry_old.block):
                synthenic_hits = __found_synthenic(synthenic_hits,
                                                   genes_in_row, entry,
                                                   cordstart,
                                                   scafname, mingen,
                                                   scafstart, block)
                # print("GATE2 " + str(entry.sori))
                if genes_in_row >= mingen:
                    found_stretch.append(genes_in_row)
                genes_in_row = 1
                cordstart = 0
                scafstart = 0
        else:  # Checks if there were enough genes in row if chromosome or
            # block has changed between entries.
            synthenic_hits = __found_synthenic(synthenic_hits, genes_in_row,
                                               entry_old, cordstart,
                                               scafname, mingen, scafstart,
                                               block)
            # print("GATE3 " + str(entry.sori))
            if genes_in_row >= mingen:
                found_stretch.append(genes_in_row)
            genes_in_row = 1
            # Sets back the start coordinates for new block or chromosome.
            cordstart = 0
            scafstart = 0
            block = block + 1
        if ecount == len(coordlist):  # Final check at the end of the list.
            synthenic_hits = __found_synthenic(synthenic_hits, genes_in_row,
                                               entry, cordstart,
                                               scafname, mingen, scafstart,
                                               block)
            if genes_in_row >= mingen:
                found_stretch.append(genes_in_row)
            # print("GATE1 " + str(entry.sori))
        entry_old = entry
    if plot:
        plt.hist(found_stretch)
        plt.title("Histogram of " + scafname)
        plt.xlabel("genes in row on 3+stretches")
        plt.ylabel("Frequency")
        plt.show()
    return (__decide_best_stretch(synthenic_hits), len(found_stretch))


def __switch_starts_ends(entry, cordstart, scafstart):
    """
    Sets strings for easy readability of orientations and switch start and
    end coordinates in reversed block occurence.
    :param entry: Last DagChain entry in successful check.
    :param cordstart: Recorded start of the stretch.
    :param scafstart: Recorded start of the scaffold.
    :return: Tuple of ori strings and start/end coordinates
    """
    oridict = {1: '+', -1: '-'}
    ori = oridict[entry.ori]
    sori = oridict[entry.sori]
    if cordstart > entry.start:
        end = cordstart
        cordstart = entry.end
    else:
        end = entry.end
    if scafstart > entry.sstart:
        send = scafstart
        scafstart = entry.send
    else:
        send = entry.send
    return ((ori, sori, cordstart, end, scafstart, send))


def __found_synthenic(synthenic_hits, genes_in_row, entry,
                      cordstart, scafname, mingen, scafstart, block):
    """
    Determines the start and end coordinates on the scaffold dependent on
    whether the block entrie was reversed or not and checkes if the
    minimal genes in rows criteria is met, in this case the found good stretch
    is appended to the first input argument and returned.
    :param synthenic_hits: List of Stretch objects that met criteria so far.
    :param genes_in_row: Integer counting number of genes in actual stretch.
    :param entry: DagChain object of last maping element.
    :param cordstart: Integer with start position of stretch on reference.
    :param scafname: String with current tef scaffold name.
    :param mingen: Integer of minimal genes in row criteria.
    :param scafstart: Integer with start or end coordinate on scaffold.
    :param block: Integer counting the current block the DagChain elements are
                  coming from
    :return: Appended list of Stretch objects (1. attribute).
    """
    if (genes_in_row >= mingen):
        (ori, sori, cordstart, end, scafstart, send) = __switch_starts_ends(
            entry, cordstart, scafstart)
        goodstretch = Stretch(
            DagChain(entry.chr, cordstart, end,
                     ori, scafname, scafstart,
                     send, sori,
                     entry.stretch, entry.number,
                     entry.score, entry.reversion, entry.block),
            genes_in_row, block, 0)
        synthenic_hits.append(goodstretch)
    return synthenic_hits


def __decide_best_stretch(listofstretches):
    """
    Sorts Stretches by DAG score OR genes_in_row.
    :param listofstretches: Unsorted list of Stretch objects.
    :return: List of Stretch objects sorted on score and genes_in_row, etc.
    """
    from operator import attrgetter
    p = 0
    sortedstretchlist = []
    for entry in sorted(listofstretches,
                        key=attrgetter('score', 'genes_in_row', 'len',
                                       'slen', 'coord', 'scoord', 'end'),
                        reverse=False):
        p = p + 1
        entry.position = p
        sortedstretchlist.append(entry)
    return (sortedstretchlist)


def documentNames():
    """
    Sets the directories and file names for input and output
    :return: Tuple of strings: dagfile, output_file, input_fileS, input_fileC
    """
    import pwd, os, sys

    try:
        dagFile = sys.argv[1]
    except IndexError:
        if os.getuid() == 1000 or pwd.getpwuid(
                os.getuid()).pw_gecos == 'fsiegris':
            dagFile = '../../i1sz/22790_24796.CDS-CDS.last.tdd10.cs0.filtered.dag.all.go_D20_g10_A3.aligncoords.gcoords.condensed'
        elif pwd.getpwuid(
                os.getuid()).pw_gecos == 'Fredy Siegrist,,,' or pwd.getpwuid(
            os.getuid()).pw_gecos == 'Gina Cannarozzi,,,':
            dagFile = '/home/fredy/i1sz/22790_24796.CDS-CDS.last.tdd10.cs0.filtered.dag.all.go_D20_g10_A3.aligncoords.gcoords.condensed'
        else:
            dagFile = input(
                'Please enter directory and file of \
                *.dag.all.go_D20_g10_A3.aligncoords.gcoords.condensed file')

    # Here is a code-chunk to easely import fasta files dependent on current user
    try:
        inputFileS = sys.argv[3]
    except:
        if os.getuid() == 1000 or pwd.getpwuid(
                os.getuid()).pw_gecos == 'fsiegris':
            inputFileS = '/windows/GNYt98ter.41.closedgt1000.sorted'
        elif pwd.getpwuid(
                os.getuid()).pw_gecos == 'Fredy Siegrist,,,' or pwd.getpwuid(
            os.getuid()).pw_gecos == 'Gina Cannarozzi,,,':
            inputFileS = '/home/fredy/i1sz/GNYt98ter.41.closedgt1000.sorted'
        else:
            inputFileS = input(
                'Please enter directory and file of GNYt98ter.41.closedgt1000.sorted file')
    try:
        inputFileC = sys.argv[4]
    except IndexError:
        if os.getuid() == 1000 or pwd.getpwuid(
                os.getuid()).pw_gecos == 'fsiegris':
            inputFileC = '/debian/Sorghum_bicolor.faa'
        elif pwd.getpwuid(
                os.getuid()).pw_gecos == 'Fredy Siegrist,,,' or pwd.getpwuid(
            os.getuid()).pw_gecos == 'Gina Cannarozzi,,,':
            inputFileC = '/home/fredy/i1sz/Sorghum_bicolor.faa'
        else:
            inputFileC = input(
                'Please enter directory and file of Sorghum_bicolor.faa file')

    try:
        outputFile = sys.argv[2]
    except IndexError:
        if os.getuid() == 1000 or pwd.getpwuid(
                os.getuid()).pw_gecos == 'fsiegris':
            outputFile = '../../i1sz/stretches_condensed_ac.csv'
        elif pwd.getpwuid(
                os.getuid()).pw_gecos == 'Fredy Siegrist,,,' or pwd.getpwuid(
            os.getuid()).pw_gecos == 'Gina Cannarozzi,,,':
            outputFile = '/home/fredy/i1sz/stretches_condensed_ac.csv'
        else:
            outputFile = input(
                'Please enter directory and file for output')

    return (dagFile, outputFile, inputFileS, inputFileC)

def dagList(x, lol, comment):
    """
    Extract list element and returns single element list of DagChain Object
    :param x: line counter
    :param lol: List of lists from condensed file
    :param comment: Multiplied Comment list from condensed file
    :return: List (1 element) of DagChain
    """
    return [DagChain(
            int(str(lol[x][5]).split('||')[0]),  # chr
            int(str(lol[x][5]).split('||')[1]),  # start
            int(str(lol[x][5]).split('||')[2]),  # end
            int(str(lol[x][5]).split('||')[4]),  # ori
            str(lol[x][1]).split('||')[0],  # scaffold
            int(str(lol[x][1]).split('||')[1]),  # sstart
            int(str(lol[x][1]).split('||')[2]),  # send
            int(str(lol[x][1]).split('||')[4]),  # sori
            int(comment[x][5]),  # stretch
            int(str(comment[x][0])[1]),  # number
            1 / float(comment[x][1]),  # score
            str(comment[x][4]),  # reversion
            int(comment[x][6])  # block
        )
    ]