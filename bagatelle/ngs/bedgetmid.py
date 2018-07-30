import sys
import pysam
from os import path
import argparse



def main():

    args = check_options(get_options())

    inbedfle = args.inbedfile

    outbedfile = args.outbedfile

    getmid(inbedfle, outbedfile)


def check_options(parser):

    args = parser.parse_args()

    if not path.exists(args.inbedfile):

        print("Can not located bed file, please input full path of bed file")

        parser.print_help()

        sys.exit(1)

    return args


def get_options():

    parser = argparse.ArgumentParser(description="bagatelle package for get middle site of bed file", prog='bagatelle')

    parser.add_argument('-i', dest='inbedfile', required=True, type=str, help='input bedfile')

    parser.add_argument('-o', dest='outbedfile', required=True, type=str, help='output bedfile')

    return parser


def getmid(inbedfle, outbedfile):

    outio = open(outbedfile, 'w')

    with open(inbedfle) as inio:
        for lin in inio:

            infor = lin.rstrip().split('\t')

            chrom = infor[0]

            start = int(infor[1])

            end = int(infor[2])

            mid1 = int((start+end)/2)

            mid2 = int((start+end)/2 + 1)


            print(chrom, mid1, mid2, sep='\t', file=outio)

    outio.close()

if __name__ == '__main__':

    try:

        main()

    except KeyboardInterrupt:

        sys.stderr.write("User interrupt\n")

        sys.exit(0)