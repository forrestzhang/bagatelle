import sys
import pysam
from os import path
import argparse


def main():

    args = check_options(get_options())

    inbamfle = args.inbamfile

    outbamfile = args.outbamfile

    maxinsert = args.maxinsert

    mininsert = args.mininsert

    print("start size selection, input=%s output=%s maxinsert=%s mininsert=%s"
          % (inbamfle, outbamfile, maxinsert, mininsert))

    insersize(inbamfle, outbamfile, maxinsert, mininsert)

    print("start size selection, input=%s output=%s maxinsert=%s mininsert=%s finished!"
          % (inbamfle, outbamfile, maxinsert, mininsert))



def check_options(parser):

    args = parser.parse_args()

    if not path.exists(args.inbamfile):

        print("Can not located bamfile, please input full path of bamfile")

        parser.print_help()

        sys.exit(1)

    return args

def get_options():

    parser = argparse.ArgumentParser(description="bagatelle package for insert size selection", prog='bagatelle')

    parser.add_argument('-i', dest='inbamfile', required=True, type=str, help='input bamfile')

    parser.add_argument('-o', dest='outbamfile', required=True, type=str, help='output bamfile')

    parser.add_argument('-l', dest='mininsert', required=True, type=int, default=1, help='min insert size')

    parser.add_argument('-m', dest='maxinsert', required=True, type=int, default=200, help='max insert size')

    return parser

def insersize(inbamfile, outbamfile, maxinsert, mininsert, paired=True):
    """

    :param bamfile:
    :param chromosome:
    :param start:
    :param end:
    :param paired:
    :return: dictionary of mhs middle site count
    """

    samfile = pysam.AlignmentFile(inbamfile)

    outbam = pysam.AlignmentFile(outbamfile, "wb", template=samfile)

    if paired:

        for aligened_read in samfile.fetch():

            if aligened_read.is_proper_pair:

                insize = abs(aligened_read.isize)

                if mininsert <= insize <= maxinsert:

                        outbam.write(aligened_read)

    samfile.close()
    outbam.close()


if __name__ == '__main__':

    try:

        main()

    except KeyboardInterrupt:

        sys.stderr.write("User interrupt\n")

        sys.exit(0)