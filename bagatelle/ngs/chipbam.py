import pyBigWig

from bagatelle.ngs import openBam
from bagatelle.ngs import Baminfo

def chipcount(bamfile, bwfile, maxinsert=180, mininsert=130, paired=False, extend=150):
    """
    :param bamfile:
    :param chromosome:
    :param start:
    :param end:
    :param paired:
    :return:
    """

    samfile = openBam.openBam(bamfile)

    readscount = dict()

    if paired:

        pass

    else:

        pass

    return readscount



def chipcoverage(bamfile, chromosome, start, end, maxinsert=180, mininsert=130, paired=False, extend=150):
    """
    :param bamfile:
    :param chromosome:
    :param start:
    :param end:
    :param paired:
    :return:
    """

    samfile = openBam.openBam(bamfile)

    readscount = dict()


    if paired:

        for aligened_read in samfile.fetch(reference=str(chromosome), start=start, end=end):

            if aligened_read.is_proper_pair:

                if not aligened_read.is_reverse:

                    if mininsert <= aligened_read.isize <= maxinsert:

                        pair_start = aligened_read.pos

                        for i in range(pair_start, pair_start+aligened_read.isize):

                            if start <= i <= end:

                                if i in readscount:

                                    readscount[i] += 1

                                else:

                                    readscount[i] = 1

    else:

        for alignend_read in samfile.fetch(reference=str(chromosome), start=start, end=end):

            for i in range(alignend_read.pos, alignend_read.pos+extend):

                if start <= i <= end:

                    if i in readscount:

                        readscount[i] += 1

                    else:

                        readscount[i] = 1

    return readscount


def chipcoveragetobw(bamfile, bwfile, maxinsert,
                         mininsert, paired, extend):

    bamfor = Baminfo.Baminfo(bamfile)

    bw = pyBigWig.open(bwfile, "w")

    bw.addHeader(list(bamfor.chrlen.items()))

    for chromosome in bamfor.chrlen:

        end = bamfor.chrlen[chromosome]

        chipcount = chipcoverage(bamfile=bamfile, chromosome=chromosome, start=1,
                                 end=end, maxinsert=maxinsert, mininsert=mininsert,
                                 paired=paired, extend=extend)

        if chipcount:

            starts = list()

            values = list()

            for start in sorted(chipcount):

                starts.append(start)

                values.append(float(chipcount[start]))

            bw.addEntries(chromosome, starts=starts, values=values,
                          span=1, step=1)

    bw.close()