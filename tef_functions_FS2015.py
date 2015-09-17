__author__ = 'fsiegris'


# date: 09.09.2015
# function source file for tef practical

class SynMap:
    def __init__(self, chr, start, end, ori):
        self.chr = chr
        self.start = start
        self.end = end
        self.ori = ori
        self.coord = (start,end)
    def chr(self):
        return self.chr
    def start(self):
        return self.start
    def end(self):
        return self.end
    def ori(self):
        return self.ori
    def length(self):
        return int((self.end-self.start)*self.ori)
    def __repr__(self):
        return repr((self.chr, (self.start, self.end), self.ori))
    def __lt__(self, other):
        (self.chr, self.start, self.end, self.ori) < (other.chr, other.start, other.end, other.ori)
    def __le__(self, other):
        (self.chr, self.start, self.end, self.ori) <= (other.chr, other.start, other.end, other.ori)
    def __ge__(self, other):
        (self.chr, self.start, self.end, self.ori) >= (other.chr, other.start, other.end, other.ori)
    def __gt__(self, other):
        (self.chr, self.start, self.end, self.ori) > (other.chr, other.start, other.end, other.ori)
    def __eq__(self, other):
        (self.chr, self.start, self.end, self.ori) == (other.chr, other.start, other.end, other.ori)
    def __ne__(self, other):
        (self.chr, self.start, self.end) != (other.chr, other.start, other.end)


def find_synthenic_block(coordlist, scafname, D=120000, mingen=3, plot=0):
    # from time import sleep
    import matplotlib.pyplot as plt
    from operator import itemgetter, attrgetter, methodcaller
    # Maximum distance between two matches (-D): plant-default 120000 bp
    #D = 120000
    entry_old = SynMap(0, 0, 0, 0)
    genes_in_row = 1
    synthenic_hits = 0
    cordstart = 0
    found_stretch = []
    for entry in sorted(coordlist, key=attrgetter('chr','start','end')):
        if (entry_old.chr == 0 or entry.chr == entry_old.chr):  # check if on same chromosome
            #if ((entry.start - entry_old.end)<0):
            #    print(' '+str(entry.chr)+' '+str(entry_old.chr)+' \ '+str(entry.start - entry_old.end)+' '+str(entry.coord)+'-'+str(entry_old.coord)+'*'+str(entry.ori)+'*'+str(entry_old.ori)+' '+str(cordstart)+' Decision: '+str((entry.start - entry_old.end) <= D and entry.start > entry_old.start))
            if ((entry.start - entry_old.end) <= D and entry.start > entry_old.start):
                genes_in_row = genes_in_row + 1
                if (genes_in_row == 2):
                    if (entry_old.start == 0):
                        cordstart = entry.start
                    else:
                        cordstart = entry_old.start
            elif (entry!=entry_old):
                synthenic_hits = __found_synthenic(synthenic_hits, genes_in_row, entry, entry_old, cordstart, scafname, mingen)
                found_stretch.append(genes_in_row)
                genes_in_row = 1
        else:
            synthenic_hits = __found_synthenic(synthenic_hits, genes_in_row, entry, entry_old, cordstart, scafname, mingen)
            found_stretch.append(genes_in_row)
            genes_in_row = 1
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
            entry_old.chr) + ' coordinates start: ' + str(cordstart) + ' end: ' + str(
            entry_old.end) + ' on ' + scafname)
        # sleep(1)
        synthenic_hits = synthenic_hits + 1
    return (synthenic_hits)
