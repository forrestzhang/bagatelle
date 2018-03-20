import re
import argparse
from os import path
import sys

def removevector(samfile):


    # samfile generat by minialign
    # example
    # minialign -t 40  -x pacbio vector.fa m54143_180307_064412.subreads.bam.fastq > m54143_180307_064412.sam

    filename = samfile

    cfile = filename.replace('.sam', '_clean.fa')

    lfile = filename.replace('.sam', '_L.fa')

    rfile = filename.replace('.sam', '_R.fa')

    rmfile = filename.replace('.sam', '_removevec.fa')

    cio = open(cfile, 'w')

    lio = open(lfile, 'w')

    rio = open(rfile, 'w')

    allio = open(rmfile, 'w')

    searchl = re.compile("^(\d*)[H,S]")
    searchr = re.compile("(\d*)[H,S]$")

    with open(filename) as inio:
        for i in inio:
            if '@' in i:
                continue

            infor = i.rstrip().split('\t')

            if infor[2] == '*':

                #             print(infor[0], 'clean')

                outst = '>' + infor[0] + '\n' + infor[9] + '\n'

                cio.write(outst)

                allio.write(outst)

            else:

                #             print(infor[5])
                #             print(infor[5],infor[9])
                ml = searchl.search(infor[5])
                mr = searchr.search(infor[5])

                if ml:
                    # get left length
                    ll = int(ml.group(1))

                    # print('l',ml.group(1))
                    seqname = infor[0] + '_L'
                    seq = (infor[9][:ll])
                    outst = '>' + seqname + '\n' + seq + '\n'
                    lio.write(outst)
                    allio.write(outst)
                if mr:
                    # get right length

                    rl = 0 - int(mr.group(1))

                    # print('r',mr.group(1))
                    seqname = infor[0] + '_R'
                    seq = (infor[9][rl:])
                    outst = '>' + seqname + '\n' + seq + '\n'
                    rio.write(outst)
                    allio.write(outst)
    cio.close()
    lio.close()
    rio.close()
    allio.close()


def main():

    args = check_options(get_options())

    removevector(args.input)

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