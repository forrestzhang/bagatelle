import pyBigWig
from multiprocessing import Pool
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

    header = list()

    for j in range(0-up, down+1, winsize):

        header.append(str(j))

    print('chromosome', 'start', 'end', "\t".join(header), sep='\t', file=socorefile)

    chrlen = bw.chroms()

    # print(chrlen)
    nbins = int((up+down)/winsize+1)

    for bed in bedio.readlines():

        bed = bed.rstrip('\n')

        bedinfor = bed.split("\t")

        chromosome = bedinfor[0]

        start = int(bedinfor[1])

        end = int(bedinfor[2])

        midsite = int((start + end)/2)



        ctregionstart = midsite - up

        ctregionend = midsite + down

        if ctregionstart < 0 :

            print("%d is too close to %s start %d" % (start, chromosome, chrlen[chromosome]))

            continue

        if ctregionend > chrlen[chromosome]:

            print("%d is too close to %s end %d" % (end, chromosome, chrlen[chromosome]))
            # print(chromosome, start, end,  ctregionstart, ctregionend, )
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


        print(chromosome, start, end, "\t".join(scorelist), sep='\t', file=socorefile)

    bedio.close()

    socorefile.close()

    bw.close()


def midcountmp(bwfile, bedfile, up=1000, down=1000, winsize=50,
             maxfilter=9999999999, minfilter=0, stranded=False,
             nantozero=True, outfile="midcount.txt",threads=2):

    bw = pyBigWig.open(bwfile, "r")

    bedio = open(bedfile, 'r')

    socorefile = open(outfile, 'w')

    header = list()

    for j in range(0-up, down+1, winsize):

        header.append(str(j))

    print('chromosome', 'start', 'end', "\t".join(header), sep='\t', file=socorefile)

    workerlist = list()

    nbins = int((up+down)/winsize+1)

    chrlen = bw.chroms()

    for worker in range(0, threads):

        workerinfor = dict()

        workerinfor['reglist'] = list()

        workerinfor['bwfile'] = bwfile

        workerinfor['nbins'] = nbins

        workerinfor['minfilter'] = minfilter

        workerinfor['maxfilter'] = maxfilter

        workerinfor['stranded'] = stranded

        workerinfor['nantozero'] = nantozero

        workerlist.append(workerinfor)

    nowworker = 0

    for bed in bedio.readlines():

        bed = bed.rstrip('\n')

        bedinfor = bed.split("\t")

        chromosome = bedinfor[0]

        start = int(bedinfor[1])

        end = int(bedinfor[2])

        midsite = int((start + end)/2)

        ctregionstart = midsite - up

        ctregionend = midsite + down

        if ctregionstart < 0 :

            print("%d is too close to %s start %d" % (start, chromosome, chrlen[chromosome]))

            continue

        if ctregionend > chrlen[chromosome]:

            print("%d is too close to %s end %d" % (end, chromosome, chrlen[chromosome]))
            # print(chromosome, start, end,  ctregionstart, ctregionend, )
            continue

        reg = chromosome+"_"+str(ctregionstart) + "_" + str(ctregionend)

        if stranded:

            reg = reg+'_'+bedinfor[5]

        nowth = nowworker % threads

        workerlist[nowth]['reglist'].append(reg)

        nowworker += 1

    pool = Pool(threads)

    for res in pool.imap_unordered(midcountworker, workerlist):

        for resultstring in res:

            print(resultstring, file=socorefile)

    socorefile.close()


def midcountworker(workerinfor):

    nbins = workerinfor['nbins']

    minfilter =workerinfor['minfilter']

    maxfilter = workerinfor['maxfilter']

    stranded = workerinfor['stranded']

    nantozero = workerinfor['nantozero']

    bw = pyBigWig.open(workerinfor['bwfile'])

    result = list()

    for region in workerinfor['reglist']:

        regioninfo = region.split("_")

        socorelist1 = bw.stats(regioninfo[0], int(regioninfo[1]), int(regioninfo[2]), nBins=nbins)

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
            if regioninfo[3] == '-':

                scorelist = scorelist.reverse()


        scorestring = region.replace('_','\t') + '\t' +'\t'.join(scorelist)

        result.append(scorestring)

    return result