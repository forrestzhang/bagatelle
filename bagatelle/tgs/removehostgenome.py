import re
import argparse
from os import path
import sys

def removehost(samfile):


    # samfile generat by minialign
    # example
    # minialign -t 40  -x pacbio K12.fa m54143_180307_064412.subreads.bam.fastq > m54143_180307_064412_K12.sam

    filename = samfile

    cfile = filename.replace('.sam', '_nohost.fa')

    cio = open(cfile, 'w')

    with open(filename) as inio:
        for i in inio:
            if '@' in i:
                continue

            infor = i.rstrip().split('\t')

            if infor[2] == '*':

                outst = '>' + infor[0] + '\n' + infor[9] + '\n'

                cio.write(outst)

    cio.close()



def main():

    args = check_options(get_options())

    removehost(args.input)

def get_options():

    parser = argparse.ArgumentParser(description="Remove vector in PacBio or Nanopore", prog='RemoveVector')

    parser.add_argument('-i', '--input', dest='input', help='input file generate by minialign',
                        required=True, type=str)

    return(parser)


def check_options(parser):

    args = parser.parse_args()

    if args.input:

        if not path.exists(args.input):

            print("Can not locate input file, please input input file.\n")

            parser.print_help()

            sys.exit(1)

    return args


if __name__ == "__main__":

    try:

        main()

    except KeyboardInterrupt:

        sys.stderr.write("User interrupt\n")

        sys.exit(0)