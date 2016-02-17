import pyBigWig

# count score up and down stream of middlesite or region


def midcount(bwfile, bedfile, up=1000, down=1000, winsize=50,
             maxfilter=9999999999, minfilter=0, stranded=False,
             nantozero=True, outfile="midcount.txt"):

    """
    :param bwfile:
    :param bedfile:
    :param up:
    :param down:
    :param winsize: windowsize
    :param stranded:
    :param maxfilter: if score higher than maxfilter, use maxfilter
    :param minfilter: if score lower than minfilter, use minfilter
    :return:

    Given bed file and bw file count score in flanking region of bed middle site(s)


    """

    bw = pyBigWig.open(bwfile, "r")

    bedio = open(bedfile, 'r')

    socorefile = open(outfile, 'w')

    for bed in bedio.readlines():

        bed = bed.rstrip('\n')

        bedinfor = bed.split("\t")

        chromosome = bedinfor[0]

        start = int(bedinfor[1])

        end = int(bedinfor[2])

        midsite = int((start + end)/2)

        nbins = int((up+down)/winsize+1)

        ctregionstart = midsite - up

        ctregionend = midsite + end

        if ctregionstart < 0 :

            print("%d is too close to chromosome start" % start)

            continue

        if ctregionend > bw.chroms(chromosome):

            print("%d is too close to chromosome end" % end)

            continue

        socorelist1 = bw.stats(chromosome, ctregionstart, ctregionend, nBins=nbins)

        scorelist = list()

        for i in socorelist1:

            if not i:

                if nantozero:

                    i = 0

                scorelist.append(str(i))

                continue

            if i < minfilter:

                i = minfilter

            if i > maxfilter:

                i = maxfilter

            scorelist.append(str(i))



        if stranded:
            # Reverse - strand score
            if bedinfor[5] == '-':

                scorelist = scorelist.reverse()

        header = list()

        for j in range(0-up, down+1, winsize):

            header.append(str(j))

        print('chromosome', 'start', 'end', "\t".join(header), sep='\t', file=socorefile)

        print(chromosome, start, end, "\t".join(scorelist), sep='\t', file=socorefile)

    bedio.close()

    socorefile.close()

    bw.close()