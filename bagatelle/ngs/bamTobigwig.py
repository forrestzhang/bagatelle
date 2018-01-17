import pyBigWig

from bagatelle.ngs import Baminfo, mhsbam
from bagatelle.ngs import dhsbam
from bagatelle.ngs import chipbam
from bagatelle.ngs import kernel
import numpy as np

"""
Note: pyBigWig only allow float as values
"""


def dhstobw(bamfile, bwfile, library='Duke'):
    # Washington is under processing
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


def dhsstrandtobw(bamfile, rbwfile, fbwfile, library='Duke', nagtiverev=True):
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

            for start in sorted(dhscut['+']):
                starts.append(start)

                values.append(float(dhscut['+'][start]))

            bwf.addEntries(chromosome, starts=starts, values=values,
                           span=1, step=1)

        if dhscut['-']:

            starts = list()

            values = list()

            for start in sorted(dhscut['-']):

                starts.append(start)
                if nagtiverev:
                    values.append(0 - float(dhscut['-'][start]))
                else:
                    values.append(float(dhscut['-'][start]))

            bwr.addEntries(chromosome, starts=starts, values=values,
                           span=1, step=1)

    bwf.close()

    bwr.close()


def dhextendtobw(bamfile, bwfile, extendsize=20, library='Duke'):
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

    bw = pyBigWig.open(bwfile, "w")

    bw.addHeader(list(bamfor.chrlen.items()))

    for chromosome in bamfor.chrlen:

        end = bamfor.chrlen[chromosome]

        dhscut = dhsbam.dhstrandcutcount(bamfile=bamfile, chromosome=chromosome, start=1,
                                         end=end, library=library)

        dhsext = dict()

        if dhscut['+']:

            # starts = list()
            #
            # values = list()

            for start in sorted(dhscut['+']):

                # starts.append(start)
                #
                # values.append(float(dhscut['+'][start]))
                for i in range(0, extendsize):

                    nowsite = start + i

                    if 0 < nowsite < end:

                        if nowsite in dhsext:
                            dhsext[nowsite] += float(dhscut['+'][start])
                        else:
                            dhsext[nowsite] = float(dhscut['+'][start])

        if dhscut['-']:

            starts = list()

            values = list()

            for start in sorted(dhscut['-']):

                for i in range(0 - extendsize, 0):
                    nowsite = start + i

                    if 0 < nowsite < end:

                        if nowsite in dhsext:

                            dhsext[nowsite] += float(dhscut['-'][start])
                        else:
                            dhsext[nowsite] = float(dhscut['-'][start])

        starts = list()

        values = list()

        for nowsite in sorted(dhsext):
            starts.append(nowsite)

            values.append(dhsext[nowsite])

        bw.addEntries(chromosome, starts=starts, values=values,
                      span=1, step=1)

    bw.close()


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


def midextendtobw(bamfile, bwfile, maxinsert, mininsert, paired=False, extend=80):
    bamfor = Baminfo.Baminfo(bamfile)

    bw = pyBigWig.open(bwfile, "w")

    bw.addHeader(list(bamfor.chrlen.items()))

    for chromosome in bamfor.chrlen:

        end = bamfor.chrlen[chromosome]

        mhsmidcount = mhsbam.mhsmidcount(bamfile=bamfile, chromosome=chromosome, start=1,
                                         end=end, maxinsert=maxinsert, mininsert=mininsert, paired=paired)

        if mhsmidcount:

            # starts = list()
            #
            # values = list()
            #
            # for start in sorted(mhsmidcount):
            #     starts.append(start)
            #
            #     values.append(float(mhsmidcount[start]))

            extended = dict()

            for midnow in mhsmidcount:

                halfextend = int(extend/2)

                for i in range(0-halfextend, halfextend+1):

                    nowsite = i + midnow

                    if nowsite < 1:
                        continue
                    if nowsite > end:
                        continue
                    if i in extended:
                        extended[i] += 1
                    else:
                        extended[i] = 1

            starts = list()

            values = list()

            for start in sorted(extended):
                starts.append(start)

                values.append(float(extended[start]))


            bw.addEntries(chromosome, starts=starts, values=values,
                          span=1, step=1)

    bw.close()


def coveragetobw(bamfile, bwfile, maxinsert, mininsert, paired=False):
    bamfor = Baminfo.Baminfo(bamfile)

    bw = pyBigWig.open(bwfile, "w")

    bw.addHeader(list(bamfor.chrlen.items()))

    for chromosome in bamfor.chrlen:

        end = bamfor.chrlen[chromosome]

        coveragecount = mhsbam.coveragecount(bamfile=bamfile, chromosome=chromosome, start=1,
                                             end=end, maxinsert=maxinsert, mininsert=mininsert, paired=paired)

        if coveragecount:

            starts = list()

            values = list()

            for start in sorted(coveragecount):
                starts.append(start)

                values.append(float(coveragecount[start]))

            bw.addEntries(chromosome, starts=starts, values=values,
                          span=1, step=1)

    bw.close()


def kernelsmooth(scorecount, regionstart, regionend, chr_length, kernelsize):
    """

    :param scorecount: dict
    :param regionstart:
    :param regionend:
    :param sitecount:
    :param chr_length:
    :param kernelsize:
    :return:
    """

    startsite = regionstart

    endsite = regionend

    regionstart = regionstart - kernelsize * 2

    regionend = regionend + kernelsize * 2

    if regionstart < 1:
        regionstart = 1

    if regionend > chr_length:
        regionend = chr_length

    renewlength = regionend - regionstart + 1

    kernelnow = kernel.smooth_kernel(kernelsize)

    readcount = list()

    kernel_score = list()

    for w in sorted(kernelnow):
        kernel_score.append(kernelnow[w])

    for n in range(regionstart, regionend + 1):

        nowcount = 0

        if n in scorecount:
            nowcount = scorecount[n]

        readcount.append(nowcount)

    nowsmoothed = np.correlate(np.array(readcount), kernel_score, "same")

    outputscore = dict()

    for j in range(0, renewlength):

        nowsite = j + regionstart

        nowscore = nowsmoothed[j]

        if (startsite <= nowsite <= endsite):
            outputscore[nowsite] = nowscore

    return outputscore


def mhsmidtobw(bamfile, bwfile, maxinsert=80, mininsert=1, paired=False):
    midtobw(bamfile=bamfile, bwfile=bwfile, maxinsert=maxinsert, mininsert=mininsert, paired=paired)


def chipmidtobw(bamfile, bwfile, maxinsert=180, mininsert=130, paired=False):
    midtobw(bamfile=bamfile, bwfile=bwfile, maxinsert=maxinsert, mininsert=mininsert, paired=paired)


def chipcvtobw(bamfile, bwfile, maxinsert=180, mininsert=130, paired=False, extend=150):
    """
    ChIP coverage to bw file
    :param bamfile:
    :param bwfile:
    :param maxinsert:
    :param mininsert:
    :param paired:
    :param extend:
    :return:
    """
    chipbam.chipcoveragetobw(bamfile=bamfile, bwfile=bwfile, maxinsert=maxinsert,
                             mininsert=mininsert, paired=paired, extend=extend)


def mhsmidkernelsmooth(bamfile, bwfile, maxinsert=80, mininsert=1, paired=False, kernelsize=30):
    bamfor = Baminfo.Baminfo(bamfile)

    bw = pyBigWig.open(bwfile, "w")

    bw.addHeader(list(bamfor.chrlen.items()))

    for chromosome in bamfor.chrlen:

        end = bamfor.chrlen[chromosome]

        mhsmidcount = mhsbam.mhsmidcount(bamfile=bamfile, chromosome=chromosome, start=1,
                                         end=end, maxinsert=maxinsert, mininsert=mininsert, paired=paired)

        mhsmidsmoothed = kernelsmooth(mhsmidcount, 1, end, end, kernelsize)

        if mhsmidsmoothed:

            starts = list()

            values = list()

            for start in sorted(mhsmidsmoothed):
                starts.append(start)

                values.append(float(mhsmidsmoothed[start]))

            bw.addEntries(chromosome, starts=starts, values=values,
                          span=1, step=1)

    bw.close()


def dhscutkernelsmooth(bamfile, bwfile, library='Duke', kernelsize=200):
    bamfor = Baminfo.Baminfo(bamfile)

    bw = pyBigWig.open(bwfile, "w")

    bw.addHeader(list(bamfor.chrlen.items()))

    for chromosome in bamfor.chrlen:

        end = bamfor.chrlen[chromosome]

        dhscut = dhsbam.dhcutcount(bamfile=bamfile, chromosome=chromosome, start=1,
                                   end=end, library=library)

        dhscutsmoothed = kernelsmooth(dhscut, 1, end, end, kernelsize)

        if dhscutsmoothed:

            starts = list()

            values = list()

            for start in sorted(dhscutsmoothed):
                starts.append(start)

                values.append(float(dhscutsmoothed[start]))

            bw.addEntries(chromosome, starts=starts, values=values,
                          span=1, step=1)

    bw.close()
