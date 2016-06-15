from Bio import SeqIO
import glob
import os.path


def abitofa(abifile, fafile):

    inio = open(abifile, 'rb')

    outio = open(fafile, 'w')

    for record in SeqIO.parse(inio, 'abi'):

        SeqIO.write(record, outio, "fasta")

    inio.close()

    outio.close()


def abidirtofa(abidir, fafile, ender="ab1"):

    abidir = os.path.dirname(abidir)

    searcher = abidir+"/*"+ender

    abifiles = glob.glob(searcher)

    print(abifiles)

    outio = open(fafile, 'w')

    for abifile in abifiles:

        inio = open(abifile, 'rb')

        for record in SeqIO.parse(inio, 'abi'):

            SeqIO.write(record, outio, "fasta")

        inio.close()

    outio.close()
