from bagatelle import openBam, Baminfo, dhsbam, mhsbam
import pyBigWig
import os

"""
Note: pyBigWig only allow float as values
"""

def dhstobw(bamfile, bwfile, library='Duke'):

    bamfor = Baminfo.Baminfo(bamfile)

    bw = pyBigWig.open(bwfile, "w")

    pass


def midtobw(bamfile, bwfile, maxinsert, mininsert, paired=False):

    bamfor = Baminfo.Baminfo(bamfile)

    bw = pyBigWig.open(bwfile, "w")

    bw.addHeader(list(bamfor.chrlen.items()))

    for chromosome in bamfor.chrlen:

        end = bamfor.chrlen[chromosome]

        mhsmidcount = mhsbam.mhsmidcount(bamfile=bamfile, chromosome=chromosome, start=1,
                                         end=end, maxinsert=maxinsert, mininsert=mininsert, paired=paired)

        if mhsmidcount:

            starts = list()

            values = list()

            for start in sorted(mhsmidcount):

                starts.append(start)

                values.append(float(mhsmidcount[start]))

            bw.addEntries(chromosome, starts=starts, values=values,
                          span=1, step=1)

    bw.close()


def mhsmidtobw(bamfile, bwfile, maxinsert=80, mininsert=1, paired=False):

    midtobw(bamfile=bamfile, bwfile=bwfile, maxinsert=maxinsert, mininsert=mininsert, paired=paired)


def chipmidtobw(bamfile, bwfile, maxinsert=180, mininsert=130, paired=False):

    midtobw(bamfile=bamfile, bwfile=bwfile, maxinsert=maxinsert, mininsert=mininsert, paired=paired)

