from pyfasta import Fasta
import argparse
from os import path
import sys

# blastn -db vector.fa  -query minscore45.1--1.fa -outfmt 6 -num_threads 40 | awk '{print $1"\t"$7"\t"$8}'  > minscore45.1--1.bls
# bedtools sort -i minscore45.1--1.bls | bedtools merge -i - -d 7000 > minscore45.1--1.bed

def removehost(fasta, bed, minilength):

    removeregion = dict()

    with open(bed) as bedin:
        for i in bedin:
            (seqname, seqstart, seqend) = i.rstrip().split('\t')[0:3]

            if seqname in removeregion:

                removeregion[seqname].append(i.rstrip())
            else:
                removeregion[seqname] = list()
                removeregion[seqname].append(i.rstrip())

    fa = Fasta(fasta)

    outfile = 'removevector_'+ fasta

    outio = open(outfile,'w')

    for seqname in fa.keys():



        if seqname in removeregion:

            if len(removeregion[seqname]) > 1:

                continue

            else:

                (seqname, seqstart, seqend) = removeregion[seqname].split('\t')[0:3]

                seqlen = len(removeregion[seqname])

                mod = list()

                # sl = skip left, sr = skip right

                if seqend - seqstart < minilength:

                    continue

                else:

                    lseq = fa[seqname][:seqstart]
                    rseq = fa[seqname][seqend:]
                    outst = ''
                    if len(lseq) > minilength:

                        lname = seqname+'_L'
                        outst = outst + '>' + lname + '\n' + lseq + '\n'

                    if len(rseq) > minilength:
                        rname = seqname + '_R'
                        outst = outst + '>' + rname + '\n' + rseq + '\n'



        else:

                outst = '>' + seqname + '\n' + str(fa[seqname]) + '\n'

                outio.write(outst)

    outio.close()



def main():

    args = check_options(get_options())

    removehost(args.fasta, args.bed, args.minilength)

def get_options():

    parser = argparse.ArgumentParser(description="Remove vector in PacBio or Nanopore", prog='RemoveVector')

    parser.add_argument('-f', '--fasta', dest='fasta', help='input fasta',
                        required=True, type=str)
    parser.add_argument('-b', '--bed', dest='bed', help='bedfile',
                        required=True, type=str)
    parser.add_argument('-l', '--length', dest='minilength', help='bedfile', default=1000,
                        required=True, type=int)

    return(parser)


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