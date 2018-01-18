from bagatelle.ngs import regioncount
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

    bwheader = list(inbw.chroms().items())

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

    bwheader = list(inbw.chroms().items())

    outbw.addHeader(bwheader)

    for chrom in inbw.chroms():

        values = inbw.values(chrom, 0, inbw.chroms(chrom), numpy=True) * ratio

        outbw.addEntries(chrom, starts=list(range(1, inbw.chroms(chrom)+1)), values=values.tolist(), span=1, step=1)

    outbw.close()

    inbw.close()

