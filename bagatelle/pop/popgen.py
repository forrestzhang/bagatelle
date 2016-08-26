import pandas as pd
import math


## Start Tajima's D

def _count_diff(list1, list2):
    diff_count = 0

    for i in range(0, len(list1)):

        if list1[i] != list2[i]:
            diff_count += 1
    return diff_count


def _pi(datafram):
    totaldiff = 0

    pair = 0

    for col1 in datafram.columns:
        # print(d3[col])
        for col2 in datafram.columns:

            if col1 == col2:

                continue

            else:

                diff = _count_diff(datafram[col1].tolist(), datafram[col2].tolist())

                #             print(col1, col2, diff)
                totaldiff += diff

                pair += 1

    return (totaldiff / pair)


def _segsite(datafram):
    segsite = 0

    for col in datafram.index:

        vc = datafram.ix[col].value_counts()

        if len(vc) > 1:
            segsite += 1

    return segsite


def _tajimad(piscore, segsite, n):
    a1 = 0

    for i in range(1, n):
        a1 += 1 / i

    a2 = 0

    for i in range(1, n):
        a2 += 1 / (i * i)

    b1 = (n + 1) / (3 * (n - 1))

    b2 = (2 * (n * n + n + 3)) / (9 * n * (n - 1))

    c1 = b1 - 1 / a1

    c2 = b2 - (n + 2) / (a1 * n) + (a2 / (a1 * a1))

    e1 = c1 / a1

    e2 = c2 / ((a1 * a1) + a2)

    d = (piscore - (segsite / a1)) / (math.sqrt((e1 * segsite) + (e2 * segsite * (segsite - 1))))
    #     print(a1,a2,b1,b2,'\n',c1, c2,e1,e2)
    return d


def tajimaD(datafram):
    """
    :param datafram: Pandas dataframe

     +------+----+----+----+----+
     |index |sp1 |sp2 |sp3 |sp4 |
     +======+====+====+====+====+
     |  1   |A   |A   |A   |T   |
     +------+----+----+----+----+
     | 2    |A   |C   |A   |A   |
     +------+----+----+----+----+
     |  3   |A   |C   |C   |G   |
     +------+----+----+----+----+
     |  4   |A   |A   |A   |A   |
     +------+----+----+----+----+
    :type datafram: pandas dataframe
    :return: Tajima's D
    :rtype: float
    """
    try:

        if len(datafram) == 0:

            tajimad = 0

        else:

            n = len(datafram.columns)

            segsite = _segsite(datafram)

            pi = _pi(datafram)

            if segsite == 0:

                tajimad = 0

                print('segsite=', segsite)

            else:

                tajimad = _tajimad(pi, segsite, n)

    except:

        print("Error:", pi, segsite, n)

        return '-'
    else:

        return tajimad


def tajimaworker(par):
    scare = par['scare']

    tmpma = par['tmpma']

    tjimad = tajimaD(tmpma)
    # print(scare,tjimad)
    return (scare, tjimad)


## End Tajima's D


def Fst_region(datafram, subpopulationlist):
    """
    :param datafram: Pandas dataframe
    :param subpopulationlist: colnames list of subpopulation in panda dataframe ['sp1','sp2']

     +------+----+----+----+----+
     |index |sp1 |sp2 |sp3 |sp4 |
     +======+====+====+====+====+
     |  1   |A   |A   |A   |T   |
     +------+----+----+----+----+
     | 2    |A   |C   |A   |A   |
     +------+----+----+----+----+
     |  3   |A   |C   |C   |G   |
     +------+----+----+----+----+
     |  4   |A   |A   |A   |A   |
     +------+----+----+----+----+
    :type subpopulationlist: list
    :type datafram: Pandas dataframe
    :return: Fst
    :rtype: float
    """
    pis = _pi(datafram[subpopulationlist])

    pit = _pi(datafram)

    fst = (pit - pis) / pit

    if fst < 0:
        fst = 0

    return fst
