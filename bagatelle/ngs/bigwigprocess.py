from bagatelle.ngs import regioncount
from bagatelle.ngs import sorttools
from os import path
import pyBigWig
import numpy as np

def getbwbybed(bedfile, inbwfile, outbwfile):

    """
    :param bedfile:
    :param inbwfile:
    :param outbwfile:
    :return:
    """

    outbw = pyBigWig.open(outbwfile, 'w')

    inbw = pyBigWig.open(inbwfile, 'r')

    chromslist = list()

    for chrom in sorted(inbw.chroms(), key=sorttools.sortChr):
        chromslist.append((chrom, inbw.chroms(chrom)))

    bwheader = chromslist

    outbw.addHeader(bwheader)

    inbw.close()

    with open(bedfile) as bedio:

        for bedlin in bedio:

            infor = bedlin.rstrip().split('\t')

            chrom = infor[0]

            start = int(infor[1])

            end = int(infor[2])

            values = regioncount.getbwbyregion(chrom, start, end, inbwfile)

            outbw.addEntries(chrom, starts=list(range(start+1, end+1)), values=values, span=1, step=1)

    outbw.close()


def bwadjust(inbwfile, outbwfile, ratio):

    inbw = pyBigWig.open(inbwfile, 'r')

    outbw = pyBigWig.open(outbwfile, 'w')

    chromslist = list()

    for chrom in sorted(inbw.chroms(), key=sorttools.sortChr):
        chromslist.append((chrom, inbw.chroms(chrom)))

    bwheader = chromslist

    outbw.addHeader(bwheader)

    for chrom in sorted(inbw.chroms(), key=sorttools.sortChr):
        # sort chromosome by using sorttools.sort
        #
        values = inbw.values(chrom, 0, inbw.chroms(chrom), numpy=True) * ratio

        outbw.addEntries(chrom, starts=list(range(1, inbw.chroms(chrom)+1)), values=values.tolist(), span=1, step=1)

    outbw.close()

    inbw.close()


def bwcompare(treatbwfile, controlbwfile, outbwfile):

    treatbw = pyBigWig.open(treatbwfile, 'r')

    controlbw = pyBigWig.open(controlbwfile, 'r')

    outbw = pyBigWig.open(outbwfile, 'w')

    chromslist = list()

    for chrom in sorted(treatbw.chroms(), key=sorttools.sortChr):
        chromslist.append((chrom, treatbw.chroms(chrom)))

    bwheader = chromslist

    outbw.addHeader(bwheader)

    for chrom in sorted(treatbw.chroms(), key=sorttools.sortChr):

        treatvalue = treatbw.values(chrom, 0, treatbw.chroms(chrom), numpy=True)

        controlvalue = controlbw.values(chrom, 0, controlbw.chroms(chrom), numpy=True)

        outvalue = treatvalue - controlvalue

        outbw.addEntries(chrom, starts=list(range(1, treatvalue.chroms(chrom) + 1)), values=outvalue.tolist(), span=1, step=1)

    treatbw.close()

    controlbw.close()

    outbw.close()