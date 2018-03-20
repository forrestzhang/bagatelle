from os import path
import re
import gzip

def fix_fq_pair_inmem(infq1, infq2, outfq1='fix1.fq', outfq2='fix2.fq', gzip=False):

    print(path.abspath(infq1), path.abspath(infq2))

    infq1 = path.abspath(infq1)

    infq2 = path.abspath(infq2)

    paired = dict()

    namere = re.compile('^@')

    checklist = ['f1name','f1seq','f1+','f1q', 'f2name', 'f2seq', 'f2+', 'f2q']

    ch1 = ['f1name', 'f1seq', 'f1+', 'f1q']

    ch2 = ['f2name', 'f2seq', 'f2+', 'f2q']

    if gzip:

        outio1 = gzip.open(path.abspath(outfq1), 'wb')

    else:
        outio1 = open(path.abspath(outfq1), 'w')

        outio2 = open(path.abspath(outfq2), 'w')

        with open(infq1) as inio1, open(infq2) as inio2:

            for f1, f2 in zip(inio1, inio2):

                f1 = f1.rstrip()

                f2 = f2.rstrip()

                if namere.match(f1):
                    pairname = f1.split(' ')[0]
                    fl1 = 0
                    if pairname in paired:
                        paired[pairname]['f1name'] = f1
                    else:
                        paired[pairname] = dict()
                        paired[pairname]['f1name'] = f1
                    fl1 +=1
                else:
                    if fl1 == 1:
                        paired[pairname]['f1seq'] = f1

                    if fl1 == 2:
                        paired[pairname]['f1+'] = f1

                    if fl1 == 3:
                        paired[pairname]['f1q'] = f1

                    fl1 += 1


                if namere.match(f2):
                    pairname = f2.split(' ')[0]
                    fl2 = 0
                    if pairname in paired:
                        paired[pairname]['f2name'] = f2
                    else:
                        paired[pairname] = dict()
                        paired[pairname]['f2name'] = f2
                    fl2 +=1
                else:
                    if fl2 == 1:
                        paired[pairname]['f2seq'] = f2

                    if fl2 == 2:
                        paired[pairname]['f2+'] = f2

                    if fl2 == 3:
                        paired[pairname]['f2q'] = f2

                    fl2 += 1


                dellist = list()

                for pairname in paired:

                    printornot = True

                    for check in checklist:

                        if check not in paired[pairname]:

                            printornot = False

                    if printornot:

                        for c1 in ch1:
                            print(c1)
                            print(paired[pairname][c1], file=outio1)

                        for c2 in ch2:
                            print(paired[pairname][c2], file=outio2)

                        #del paired[pairname]
                        dellist.append(pairname)

                for pairname in dellist:

                    del paired[pairname]

        outio1.close()

        outio2.close()

if __name__ == '__main__':

    fix_fq_pair_inmem(infq1='../tests/inputfile/infq1.fq', infq2='../tests/inputfile/infq2.fq',
                      outfq1='../tests/inputfile/outfq1.fq',  outfq2='../tests/inputfile/outfq2.fq',)