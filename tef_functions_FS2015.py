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

    def __repr__(self):
        return repr((self.chr, (self.start, self.end), self.ori))

    def __lt__(self, other):
        (self.chr, self.start, self.end, self.sstart, self.ssend) < (
        other.chr, other.start, other.end, other.sstart, self.send)


def find_synthenic_block(coordlist, fa, D=120000, mingen=3, plot=0):
    # from time import sleep
    import matplotlib.pyplot as plt

    from operator import itemgetter, attrgetter, methodcaller

    scafname=str(fa.id)
    # Maximum distance between two matches (-D): plant-default 120000 bp
    # D = 120000
    typecheck = type(coordlist[0])==type(DagChain(0, 0, 0, 0, '', 0, 0, 0, 0, 0, 0))
    if typecheck:
        entry_old = DagChain(0, 0, 0, 0, '', 0, 0, 0, 0, 0, 0)
    elif type(coordlist[0])==type(SynMap(0, 0, 0, 0)):
        entry_old = SynMap(0, 0, 0, 0)
    else:
        print('Wrong list type!')
    genes_in_row = 1
    synthenic_hits = 0
    cordstart = 0
    scafstart = 0
    found_stretch = []
    for entry in sorted(coordlist, key=attrgetter('chr', 'start', 'end')):
        if (                    # check if on same chromosome
            (entry_old.chr == 0 or (entry.chr == entry_old.chr and
             entry_old.number == entry.number and entry_old.score == entry.score))):
            # !!! Will not work if number and score are the same for different blocks
            # if ((entry.start - entry_old.end)<0):
            if typecheck:
                orient=(
                    entry.ori==entry.sori and entry_old.ori==entry.ori or
                    entry.ori!=entry.sori and entry_old.ori==entry.ori
                )
            else:
                orient=True
            #    print(' '+str(entry.chr)+' '+str(entry_old.chr)+' \ '+str(entry.start - entry_old.end)+' '+str(entry.coord)+'-'+str(entry_old.coord)+'*'+str(entry.ori)+'*'+str(entry_old.ori)+' '+str(cordstart)+' Decision: '+str((entry.start - entry_old.end) <= D and entry.start > entry_old.start))
            if ((
                    entry.start - entry_old.end) <= D and
                    entry.start > entry_old.start and
                    orient==True):
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
            elif (entry != entry_old):  # for identical entires in 1st vers.
                synthenic_hits = __found_synthenic(synthenic_hits,
                                                   genes_in_row, entry,
                                                   entry_old, cordstart,
                                                   scafname, mingen, fa, scafstart)
                found_stretch.append(genes_in_row)
                genes_in_row = 1
        else:
            synthenic_hits = __found_synthenic(synthenic_hits, genes_in_row,
                                               entry_old, entry_ancient, cordstart,
                                               scafname, mingen, fa, scafstart)
            found_stretch.append(genes_in_row)
            genes_in_row = 1
            print('chromosome jump '+str(entry_old[0])+' '+str(entry[0]))
            cordstart = 0
            scafstart = 0
        entry_ancient = entry_old
        entry_old = entry
    if plot:
        plt.hist(found_stretch)
        plt.title("Histogram of " + scafname)
        plt.xlabel("genes in row on 3+stretches")
        plt.ylabel("Frequency")
        plt.show()
    return (synthenic_hits)


def __found_synthenic(synthenic_hits, genes_in_row, entry, entry_old,
                      cordstart, scafname, mingen, fa, scafstart):
    try:
        if entry_old.ori==1:
            ori='+'
        else:
            ori='-'
        if entry_old.sori==1:
            sori='+'
        else:
            sori='-'
        sequ=__get_scaffold_fa(fa, scafstart, entry.send)
        sequence=sequ[:50]+' ... '+sequ[-50:]
    except:
        if (genes_in_row >= mingen):
            print('Orientation not checked for!')
        ori=str(entry_old.ori)
        sori=str(entry.ori)
        sequence=('NA')
    finally:
        if (genes_in_row >= mingen):
            print('got ' + str(
                genes_in_row) + '-genes stretch of synthenic genes on chromosome ' + str(
                entry_old.chr) + ' coordinates start: ' + str(
                cordstart) + ' end: ' + str(
                entry_old.end) + ' on ' + scafname + ' in ' + sori +
                ori + ' orientation:\n' + sequence  + str(scafstart) + ' ' + str(entry.send))
            # sleep(1)
            synthenic_hits = synthenic_hits + 1
    return (synthenic_hits)

def __get_scaffold_fa(fasta, start, end):
    from Bio import SeqIO, SeqFeature
    return(str(fasta.seq[start:end]))

def __get_chr_fa(chr, start, end):
    from Bio import SeqIO, SeqFeature
    return('')
