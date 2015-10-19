__author__ = 'fsiegris'


# Date: 09.09.2015
# Functions for tef practical.

class SynMap:
    """
        Simple class used for raw mapping inputs without block grouping of
        condensed output.
    """

    def __init__(self, chr, start, end, ori):
        self.chr = chr
        self.start = start
        self.end = end
        self.ori = ori
        self.coord = (start, end)

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

    def __repr__(self):
        return repr((self.chr, (self.start, self.end), self.ori))

    def __lt__(self, other):
        (self.chr, self.start, self.end, self.ori) < (
            other.chr, other.start, other.end, other.ori)

    def __le__(self, other):
        (self.chr, self.start, self.end, self.ori) <= (
            other.chr, other.start, other.end, other.ori)

    def __ge__(self, other):
        (self.chr, self.start, self.end, self.ori) >= (
            other.chr, other.start, other.end, other.ori)

    def __gt__(self, other):
        (self.chr, self.start, self.end, self.ori) > (
            other.chr, other.start, other.end, other.ori)

    def __eq__(self, other):
        (self.chr, self.start, self.end, self.ori) == (
            other.chr, other.start, other.end, other.ori)

    def __ne__(self, other):
        (self.chr, self.start, self.end) != (other.chr, other.start, other.end)


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

    def __init__(self, dag, genes_in_row, newblock, position, cfa, sfa):
        DagChain.__init__(self, dag.chr, dag.start, dag.end, dag.ori,
                          dag.scaffold, dag.sstart, dag.send,
                          dag.sori, dag.stretch, dag.number,
                          dag.score, dag.reversion, dag.block)
        self.genes_in_row = genes_in_row
        self.newblock = newblock
        self.position = position
        self.cfa = cfa
        self.sfa = sfa
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


def find_synthenic_block(coordlist, fa, chfa, D=12000000, mingen=3,
                         plot=False):
    """
    Checkes the nucleotide distance of the blocks and applies minimal genes
    in a syngenic stretch cutoff to a list of DagChain objects.
    Returns the sorted (by DAG score OR genes_in_row) list of Stretch objects.
    Maximum distance between two matches (-D): plant-default 120000 bp
    100XD = 12 000 000 applied to accept alternative by gene distance
    criteria used in generation of condensed list.
    :param coordlist: List of all DagChain objects for one scaffold.
    :param fa: SeqIO Fasta object of tef scaffold.
    :param chfa: SeqIO Fasta object of reference genomes with ordered chr.
    :param D: Integer for maximal distance between two genes for cut-off.
    :param mingen: Integer for minimal genes to form a stretch.
    :param plot: Boolean if histogram of number of genes_in_row should be
                 plotted for every entry - to be used for subsets.
    :return: Tuple of list of sorted Stretch objects and number of geq mingen
             stretches for statistical reporting.
    """

    if plot:
        import matplotlib.pyplot as plt

    from operator import attrgetter

    # Setup of Stretch attributes.
    scafname = str(fa.id)

    # Integrity check for new version of input files (condensed).
    typecheck = type(coordlist[0]) == type(
        DagChain(0, 0, 0, 0, '', 0, 0, 0, 0, 0, 0, 'f', 0))
    if typecheck:
        entry_old = DagChain(0, 0, 0, 0, '', 0, 0, 0, 0, 0, 0, 'f', 0)
    elif type(coordlist[0]) == type(SynMap(0, 0, 0, 0)):
        entry_old = SynMap(0, 0, 0, 0)
    else:
        print('Wrong list type!')

    # Initiation of looping variables for new Stretch attributes.
    block = 1
    genes_in_row = 1
    synthenic_hits = []
    cordstart = 0
    scafstart = 0
    found_stretch = []
    ecount = 0
    # if not typecheck: # was used to sort mapped elements in old version
    # coordlist = sorted(coordlist, key=attrgetter('chr', 'start', 'end'))
    for entry in coordlist:
        # print(entry)
        ecount = ecount + 1
        if typecheck:
            # print(str(entry.chr == entry_old.chr) + str(entry_old.number ==
            # entry.number) + str(entry_old.score == entry.score))
            chcheck = entry.chr == entry_old.chr and entry_old.number == entry.number and entry_old.score == entry.score and entry_old.block == entry.block

        else:
            chcheck = entry.chr == entry_old.chr
        if (  # Checks if on same chromosome or block.
              (entry_old.chr == 0 or (chcheck))):
            orient = True
            # print(' '+str(entry.chr)+' '+str(entry_old.chr)+' \ '+
            # str(entry.start - entry_old.end)+' '+str(entry.coord)+'-'+
            # str(entry_old.coord)+'*'+str(entry.ori)+'*'+str(entry_old.ori)
            # +' '+str(cordstart)+' Decision: '+str((entry.start -
            # entry_old.end) <= D and entry.start > entry_old.start))
            if (  # Checks for distance between two mapping entries.
                  abs(entry.start - entry_old.end) <= D and
                          orient == True):
                genes_in_row = genes_in_row + 1
                if (genes_in_row == 2):
                    # Sets the start coordinates for the Stretch.
                    if (entry_old.start == 0):
                        cordstart = entry.start
                        if typecheck:
                            if entry.reversion == 'f':
                                scafstart = entry.sstart
                            else:
                                scafstart = entry.send
                    else:
                        cordstart = entry_old.start
                        if typecheck:
                            scafstart = entry_old.sstart
                        else:
                            scafstart = 0
            elif (  # Checks if there are enough genes in row for a good
                    # stretch if distance criteria is no longer met.
                    entry != entry_old or entry.block != entry_old.block):
                synthenic_hits = __found_synthenic(synthenic_hits,
                                                   genes_in_row, entry,
                                                   cordstart,
                                                   scafname, mingen, fa,
                                                   scafstart, block, chfa)
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
                                               scafname, mingen, fa, scafstart,
                                               block, chfa)
            # print("GATE3 " + str(entry.sori))
            if genes_in_row >= mingen:
                found_stretch.append(genes_in_row)
            genes_in_row = 1
            # print('chromosome jump '+str(entry_old[0])+' '+str(entry[0]))
            # Sets back the start coordinates for new block or chromosome.
            cordstart = 0
            scafstart = 0
            block = block + 1
        if ecount == len(coordlist):  # Final check at the end of the list.
            if genes_in_row >= mingen:
                found_stretch.append(genes_in_row)
            synthenic_hits = __found_synthenic(synthenic_hits, genes_in_row,
                                               entry, cordstart,
                                               scafname, mingen, fa, scafstart,
                                               block, chfa)
            # print("GATE1 " + str(entry.sori))
        entry_old = entry
    if plot:
        plt.hist(found_stretch)
        plt.title("Histogram of " + scafname)
        plt.xlabel("genes in row on 3+stretches")
        plt.ylabel("Frequency")
        plt.show()
    return (__decide_best_stretch(synthenic_hits), len(found_stretch))


def __found_synthenic(synthenic_hits, genes_in_row, entry,
                      cordstart, scafname, mingen, fa, scafstart, block, chfa):
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
    :param fa: SeqIO object with fasta information of scaffold.
    :param scafstart: Integer with start or end coordinate on scaffold.
    :param block: Integer counting the current block the DagChain elements are
                  coming from
    :param chfa: SeqIO object with fasta information of reference genome.
    :return:
    """
    try:
        if entry.ori == 1:
            ori = '+'
        else:
            ori = '-'
        if entry.sori == 1:
            sori = '+'
        else:
            sori = '-'
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
        sequ = __get_scaffold_fa(fa, scafstart, send)
        chseq = __get_chr_fa(chfa, entry.chr, cordstart, end)
        # sequence = sequ[:50] + ' ... ' + sequ[-50:]
        # seqcoord = ' ' + str(scafstart) + ' ' + str(send)
    except:
        if (genes_in_row >= mingen):
            print('Orientation not checked for!')
        chseq = 'NA'
    finally:
        if (genes_in_row >= mingen):
            """
            print('got ' + str(
                genes_in_row) +
                  '-genes stretch of synthenic genes on chromosome ' + str(
                entry.chr) + ' coordinates start: ' + str(
                cordstart) + ' end: ' + str(
                entry.end) + ' on ' + scafname + ' in ' + sori + ori +
                ' orientation:\n' + sequence + seqcoord + ' '+ str(block))
            """
            goodstretch = Stretch(
                DagChain(entry.chr, cordstart, end,
                         entry.ori, scafname, scafstart,
                         send, entry.sori,
                         entry.stretch, entry.number,
                         entry.score, entry.reversion, entry.block),
                genes_in_row, block, 0,
                chseq, sequ)
            synthenic_hits.append(goodstretch)
    return synthenic_hits


def __get_scaffold_fa(fasta, start, end):
    """
    Returns a nucleotide string (fasta) between the longest
    syntenic stretch between first and last mapped gene on tef scaffold.
    :param fasta: SeqIO fasta object for given stretch.
    :param start: String of lower Fasta nucleotide coordinate on scaffold.
    :param end: String of higher Fasta nucleotide coordinate on scaffold.
    :return: String of nucleotide sequence for stretch between mapped elements
             on scaffold.
    """
    from Bio import SeqIO, SeqFeature
    return (str(fasta.seq[start - 1:end - 1]))


def __get_chr_fa(fasta, chr, start, end):
    """
    Returns a nucleotide string (fasta) between the longest syntenic stretch
    between first and last mapped gene on reference genome.
    :param fasta: SeqIO fasta object of reference genome for given stretch.
    :param chr: String or integer of chromosome number.
    :param start: String of higher Fasta nucleotide coordinate on chromosome.
    :param end: String of lower Fasta nucleotide coordinate on chromosome.
    :return: String of nucleotide sequence for stretch between mapped elements
             on reference genome.
    """
    return str(fasta[int(chr - 1)].seq[start - 1:end - 1])


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
                        key=attrgetter('len', 'score', 'genes_in_row', 
                                       'slen', 'coord', 'scoord', 'end'),
                        reverse=False):
        p = p + 1
        entry.position = p
        sortedstretchlist.append(entry)
    return (sortedstretchlist)
