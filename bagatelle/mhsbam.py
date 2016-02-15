from bagatelle import openBam


def mhsmidcount(bamfile, chromosome, start, end, maxinsert=80, mininsert=1, paired=False):
    """

    :param bamfile:
    :param chromosome:
    :param start:
    :param end:
    :param paired:
    :return: dictionary of mhs middle site count
    """

    samfile = openBam.openBam(bamfile)

    readscount = dict()

    if paired:

        for aligened_read in samfile.fetch(reference=str(chromosome), start=start, end=end):

            if aligened_read.is_proper_pair:

                if not aligened_read.is_reverse:

                    if mininsert <= aligened_read.isize <= maxinsert:

                        pair_start = aligened_read.pos

                        if aligened_read.isize % 2 == 0:

                            middle1 = pair_start + aligened_read.isize / 2

                            middle2 = pair_start + aligened_read.isize / 2 - 1

                        else:

                            middle1 = pair_start + int(aligened_read.isize / 2)

                            middle2 = pair_start + int(aligened_read.isize / 2)

                        middleint1 = int(middle1)

                        middleint2 = int(middle2)

                        if start<= middleint1 <=end:

                            if middleint1 in readscount:

                                readscount[middleint1] = readscount[middleint1] + 1

                            else:

                                readscount[middleint1] = 1

                        if start<= middleint2 <=end:

                            if middleint2 in readscount:

                                readscount[middleint2] = readscount[middleint2] + 1

                            else:

                                readscount[middleint2] = 1

    else:

        for alignend_read in samfile.fetch(reference=str(chromosome), start=start, end=end):

            if mininsert <= alignend_read.alen <= maxinsert:

                if alignend_read.alen % 2 == 0:

                    middle1 = alignend_read.pos + alignend_read.alen / 2

                    middle2 = alignend_read.pos + alignend_read.alen / 2 + 1

                else:

                    middle1 = alignend_read.pos + int(alignend_read.alen / 2)

                    middle2 = alignend_read.pos + int(alignend_read.alen / 2)


                middleint1 = int(middle1)

                middleint2 = int(middle2)

                if (start <= middleint1 <=end):

                    if middleint1 in readscount:

                        readscount[middleint1] = readscount[middleint1] + 1

                    else:

                         readscount[middleint1] = 1

                if (start <= middleint2 <=end):

                    if middleint2 in readscount:

                        readscount[middleint2] = readscount[middleint2] + 1

                    else:

                        readscount[middleint2] = 1

    return readscount
