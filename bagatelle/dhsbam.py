from bagatelle import openBam


def dhcutcount(bamfile, chromosome, start, end, library='Duke'):
    """

    :param bamfile: bamfile
    :param chromosome:
    :param start:
    :param end:
    :param library: Duke or Washington

        Duke: |=====>
                        <=====|

        Washington: |===========|

        Out put cutting site '|'

    :return: dictionary of cutting site count
    """

    samfile = openBam.openBam(bamfile)

    readscount = dict()

    if library == 'Duke':

        for aligned_read in samfile.fetch(reference=str(chromosome), start=start, end=end):

            if aligned_read.is_reverse:

                site = aligned_read.aend

            else:

                site = aligned_read.pos

            site = site + 1

            if site in readscount:

                readscount[site] = readscount[site] + 1

            else:

                readscount[site] = 1

    elif library == 'Washington':

        pass

    else:

        pass

    return readscount



def dhstrandcutcount(bamfile, chromosome, start, end, library='Duke'):
    """

    :param bamfile: bamfile
    :param chromosome:
    :param start:
    :param end:
    :param library: Duke or Washington

        Duke: |=====>
                        <=====|

        Washington: |===========|

        Out put cutting site '|'

    :return: dictionary of cutting site count
    """

    samfile = openBam.openBam(bamfile)

    readscount = dict()

    readscount['+'] = dict()

    readscount['-'] = dict()

    if library == 'Duke':

        for aligned_read in samfile.fetch(reference=str(chromosome), start=start, end=end):

            if aligned_read.is_reverse:

                site = aligned_read.aend

                site = site + 1

                if site in readscount['-']:

                    readscount['-'][site] = readscount['-'][site] + 1

                else:

                    readscount['-'][site] = 1

            else:

                site = aligned_read.pos

                site = site + 1

                if site in readscount['+']:

                    readscount['+'][site] = readscount['+'][site] + 1

                else:

                    readscount['+'][site] = 1



    elif library == 'Washington':

        pass

    else:

        pass


    return readscount