__author__ = 'fsiegris'


# date: 09.09.2015
# function source file for tef practical

class SynMap:
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
    We need additional information on nucleotide sequence of chromosome and scaffold
    for calculation of evolutionary distance with PAML or RSD (on proteins?)
    """

    def __init__(self, chr, start, end, ori, scaffold, sstart, send, sori,
                 stretch, number, score):
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

    def all(self):
        return (
            self.chr, self.start, self.end, self.ori, self.scaffold,
            self.sstart,
            self.send, self.scoord, self.sori, self.stretch, self.number,
            self.score)

    def __repr__(self):
        return repr((self.chr, (self.start, self.end), self.ori))

    def __lt__(self, other):
        (self.chr, self.start, self.end, self.sstart, self.send) < (
            other.chr, other.start, other.end, other.sstart, self.send)


class Stretch(DagChain):
    """
    A class for returning the values found for the best (or additional)
    stretches of genes found on both scaffolds and reference genome
    (chromosomes).
    """

    def __init__(self, dag, genes_in_row, block, position, cfa, sfa):
        super(Stretch, self).__init__(dag.chr, dag.start, dag.end, dag.ori,
                                      dag.scaffold, dag.sstart, dag.send,
                                      dag.sori, dag.stretch, dag.number,
                                      dag.score)
        self.genes_in_row = genes_in_row
        self.block = block
        self.position = position
        self.cfa = cfa
        self.sfa = sfa
        self.len = self.end - self.start
        self.slen = self.send - self.sstart

    def gir(self):
        return self.genes_in_row

    def block(self):
        return self.block

    def position(self):
        return self.position

    def all(self):
        return (
            self.chr, self.start, self.end, self.ori, self.scaffold,
            self.sstart,
            self.send, self.sori, self.stretch, self.number, self.score,
            self.genes_in_row, self.block, self.position, self.cfa, self.sfa)

    def __repr__(self):
        return repr((self.position, self.block, self.genes_in_row, self.score, self.len,
                    self.slen, self.chr, self.coord, self.scaffold, self.scoord))

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


def find_synthenic_block(coordlist, fa, chfa, D=120000, mingen=3, plot=False):
    if plot:
        import matplotlib.pyplot as plt

    from operator import itemgetter, attrgetter, methodcaller

    scafname = str(fa.id)
    # Maximum distance between two matches (-D): plant-default 120000 bp
    # D = 120000
    typecheck = type(coordlist[0]) == type(
        DagChain(0, 0, 0, 0, '', 0, 0, 0, 0, 0, 0))
    if typecheck:
        entry_old = DagChain(0, 0, 0, 0, '', 0, 0, 0, 0, 0, 0)
    elif type(coordlist[0]) == type(SynMap(0, 0, 0, 0)):
        entry_old = SynMap(0, 0, 0, 0)
    else:
        print('Wrong list type!')
    block = 0
    genes_in_row = 1
    synthenic_hits = []
    cordstart = 0
    scafstart = 0
    found_stretch = []
    ecount = 0
    if not typecheck:
        coordlist = sorted(coordlist, key=attrgetter('chr', 'start', 'end'))
    for entry in coordlist:
        ecount = ecount + 1
        if typecheck:
            # print(str(entry.chr == entry_old.chr) + str(entry_old.number ==
            # entry.number) + str(entry_old.score == entry.score))
            chcheck = entry.chr == entry_old.chr and entry_old.number == entry.number and entry_old.score == entry.score
        else:
            chcheck = entry.chr == entry_old.chr
        if (  # check if on same chromosome or block
              (entry_old.chr == 0 or (chcheck))):
            # !!! Will not work if number and score are the same for
            # different blocks;
            # if ((entry.start - entry_old.end)<0):
            block = block + 1
            # You will kill me for the next statement :-)
            if False:  # typecheck:
                orient = (
                    entry.ori == entry.sori and entry_old.ori == entry.ori or
                    entry.ori != entry.sori and entry_old.ori == entry.ori
                )
            else:
                orient = True
            # print(' '+str(entry.chr)+' '+str(entry_old.chr)+' \ '+
            # str(entry.start - entry_old.end)+' '+str(entry.coord)+'-'+
            # str(entry_old.coord)+'*'+str(entry.ori)+'*'+str(entry_old.ori)
            # +' '+str(cordstart)+' Decision: '+str((entry.start -
            # entry_old.end) <= D and entry.start > entry_old.start))
            if (
                            abs(entry.start - entry_old.end) <= D and
                    # entry.start > entry_old.start and
                            orient == True):
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
            elif (entry != entry_old):  # for identical entries in 1st vers.
                synthenic_hits = __found_synthenic(synthenic_hits,
                                                   genes_in_row, entry,
                                                   entry_old, cordstart,
                                                   scafname, mingen, fa,
                                                   scafstart, block, chfa)
                if genes_in_row >= mingen:
                    found_stretch.append(genes_in_row)
                genes_in_row = 1
                cordstart = 0
                scafstart = 0
        else:
            synthenic_hits = __found_synthenic(synthenic_hits, genes_in_row,
                                               entry_old, entry_ancient,
                                               cordstart,
                                               scafname, mingen, fa, scafstart,
                                               block, chfa)
            if genes_in_row >= mingen:
                found_stretch.append(genes_in_row)
            genes_in_row = 1
            # print('chromosome jump '+str(entry_old[0])+' '+str(entry[0]))
            cordstart = 0
            scafstart = 0
        if ecount == len(coordlist):
            synthenic_hits = __found_synthenic(synthenic_hits, genes_in_row,
                                               entry_old, entry_ancient,
                                               cordstart,
                                               scafname, mingen, fa, scafstart,
                                               block, chfa)
            if genes_in_row >= mingen:
                found_stretch.append(genes_in_row)
        entry_ancient = entry_old
        entry_old = entry
    print(__decide_best_stretch(synthenic_hits))
    if plot:
        plt.hist(found_stretch)
        plt.title("Histogram of " + scafname)
        plt.xlabel("genes in row on 3+stretches")
        plt.ylabel("Frequency")
        plt.show()
    return (len(found_stretch))


def __found_synthenic(synthenic_hits, genes_in_row, entry, entry_old,
                      cordstart, scafname, mingen, fa, scafstart, block, chfa):
    try:
        if entry_old.ori == 1:
            ori = '+'
        else:
            ori = '-'
        if entry_old.sori == 1:
            sori = '+'
        else:
            sori = '-'
        sequ = __get_scaffold_fa(fa, scafstart, entry_old.send)
        chseq = __get_chr_fa(chfa, entry_old.chr, entry_old.start,
                             entry_old.end)
        sequence = sequ[:50] + ' ... ' + sequ[-50:]
        seqcoord = ' ' + str(scafstart) + ' ' + str(entry.send)
    except:
        if (genes_in_row >= mingen):
            print('Orientation not checked for!')
        ori = str(entry_old.ori)
        sori = str(entry.ori)
        chseq = 'NA'
        sequence = 'NA'
        seqcoord = ' / NA'
    finally:
        if (genes_in_row >= mingen):
            """print('got ' + str(
                genes_in_row) +
                  '-genes stretch of synthenic genes on chromosome ' + str(
                entry_old.chr) + ' coordinates start: ' + str(
                cordstart) + ' end: ' + str(
                entry_old.end) + ' on ' + scafname + ' in ' + sori +
                  ori + ' orientation:\n' + sequence + seqcoord + ' '+ str(block))"""
            goodstretch = Stretch(
                DagChain(entry_old.chr, cordstart, entry_old.end,
                         entry_old.ori, scafname, scafstart,
                         entry_old.send, entry_old.sori,
                         entry_old.stretch, entry_old.number,
                         entry_old.score), genes_in_row, block, 0,
                chseq, sequ)
            synthenic_hits.append(goodstretch)
    return synthenic_hits


def __get_scaffold_fa(fasta, start, end):
    from Bio import SeqIO, SeqFeature
    return (str(fasta.seq[start - 1:end - 1]))


def __get_chr_fa(fasta, chr, start, end):
    return str(fasta[int(chr - 1)].seq[start - 1:end - 1])


def __decide_best_stretch(listofstretches):
    from operator import attrgetter
    p = 0
    sortedstretchlist = []
    for entry in sorted(listofstretches, key=attrgetter('genes_in_row', 'score', 'len',
         'slen', 'coord', 'scoord', 'end'), reverse=True):
        p = p + 1
        entry.position = p
        sortedstretchlist.append(entry)
    return (sortedstretchlist)
