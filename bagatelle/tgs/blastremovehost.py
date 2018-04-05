from pyfasta import Fasta
import argparse
from os import path
import sys


# blastn -db K12DH10B.fasta  -query removevector_minscore45.1--1.fa  -outfmt 6 -num_threads 40 | cut -f 1 > removevector_minscore45.1--1.bls


def removehost(fasta, bed):
    removeregion = dict()

    with open(bed) as bedin:
        for i in bedin:




                removeregion[i.rstrip()[0]] = 1


    fa = Fasta(fasta)

    outfile = 'removehost_' + fasta

    outio = open(outfile, 'w')

    for seqname in fa.keys():

        if seqname in removeregion:

            continue

        else:

            outst = '>' + seqname + '\n' + str(fa[seqname]) + '\n'

            outio.write(outst)

    outio.close()


def main():
    args = check_options(get_options())

    removehost(args.fasta, args.bed)


def get_options():
    parser = argparse.ArgumentParser(description="Remove vector in PacBio or Nanopore", prog='RemoveVector')

    parser.add_argument('-f', '--fasta', dest='fasta', help='input fasta',
                        required=True, type=str)
    parser.add_argument('-b', '--bls', dest='bed', help='blast reslut',
                        required=True, type=str)


    return (parser)


def check_options(parser):
    args = parser.parse_args()

    if args.fasta:

        if not path.exists(args.fasta):
            print("Can not locate input file, please input input file.\n")

            parser.print_help()

            sys.exit(1)

    if args.bed:

        if not path.exists(args.bed):
            print("Can not locate input file, please input bed file.\n")

            parser.print_help()

            sys.exit(1)

    return args


if __name__ == "__main__":

    try:

        main()

    except KeyboardInterrupt:

        sys.stderr.write("User interrupt\n")

        sys.exit(0)