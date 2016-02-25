from bagatelle import openBam, Baminfo, dhsbam, mhsbam
import pyBigWig
import os

"""
Note: pyBigWig only allow float as values
"""

def dhstobw(bamfile, bwfile, library='Duke'):

    """

    :param bamfile:
    :param bwfile:
    :param library:Duke or Washington

        Duke: |=====>
                        <=====|

        Washington: |===========|

        Out put cutting site '|'
    :return:
    """

    bamfor = Baminfo.Baminfo(bamfile)

    bw = pyBigWig.open(bwfile, "w")

    bw.addHeader(list(bamfor.chrlen.items()))

    for chromosome in bamfor.chrlen:

        end = bamfor.chrlen[chromosome]

        dhscut = dhsbam.dhcutcount(bamfile=bamfile, chromosome=chromosome, start=1,
                                         end=end, library=library)

        if dhscut:

            starts = list()

            values = list()

            for start in sorted(dhscut):

                starts.append(start)

                values.append(float(dhscut[start]))

            bw.addEntries(chromosome, starts=starts, values=values,
                          span=1, step=1)

    bw.close()


def dhsstrandtobw(bamfile,rbwfile, fbwfile, library='Duke', nagtiverev=True):

    """

    :param bamfile:
    :param bwfile:
    :param library:Duke or Washington

        Duke: |=====>
                        <=====|

        Washington: |===========|

        Out put cutting site '|'
    :param nagtiverev: - Strand use nagative score
    :return:
    """

    bamfor = Baminfo.Baminfo(bamfile)

    bwf = pyBigWig.open(fbwfile, "w")

    bwr = pyBigWig.open(rbwfile, "w")

    bwf.addHeader(list(bamfor.chrlen.items()))

    bwr.addHeader(list(bamfor.chrlen.items()))

    for chromosome in bamfor.chrlen:

        end = bamfor.chrlen[chromosome]

        dhscut = dhsbam.dhstrandcutcount(bamfile=bamfile, chromosome=chromosome, start=1,
                                         end=end, library=library)

        if dhscut['+']:

            starts = list()

            values = list()

            for start in sorted(dhscut):

                starts.append(start)

                values.append(float(dhscut[start]))

            bwf.addEntries(chromosome, starts=starts, values=values,
                          span=1, step=1)

        if dhscut['-']:

            starts = list()

            values = list()

            for start in sorted(dhscut):

                starts.append(start)
                if nagtiverev:
                    values.append(0-float(dhscut[start]))
                else:
                    values.append(float(dhscut[start]))

            bwr.addEntries(chromosome, starts=starts, values=values,
                          span=1, step=1)

    bwf.close()

    bwr.close()


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

