from bagatelle.ngs import regioncount
from os import path
import pyBigWig

def getbwbybed(bedfile, inbwfile, outbwfile):

    """
    :param bedfile:
    :param inbwfile:
    :param outbwfile:
    :return:
    """

    outbw = pyBigWig.open(outbwfile, 'w')

    inbw = pyBigWig.open(inbwfile, 'r')

    bwheader = inbw.header()

    outbw.addHeader(bwheader)

    inbw.close()

    with open(bedfile) as bedio:

        for bedlin in bedio:

            infor = bedlin.rstrip().split('\t')

            chrom = infor[0]

            start = int(infor[1]) + 1

            end = int(infor[2]) + 1

            values = regioncount.getbwbyregion(chrom, start, end, inbwfile)

            outbw.addEntries(chrom, start=start, values=values, span=1, step=1)

    outbw.close()
