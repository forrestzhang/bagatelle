from bagatelle.ngs import openBam


def chipcount(bamfile, chromosome, start, end, paired=False):
    """
    :param bamfile:
    :param chromosome:
    :param start:
    :param end:
    :param paired:
    :return:
    """

    samfile = openBam.openBam(bamfile)

    readscount = dict()

    if paired:

        pass

    else:

        pass

    return readscount
