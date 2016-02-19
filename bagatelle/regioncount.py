import pyBigWig
from multiprocessing import Pool
import gzip
import shutil
import os
# count score up and down stream of middlesite or region


def midcount(bwfile, bedfile, up=1000, down=1000, winsize=50,
             maxfilter=9999999999, minfilter=0, stranded=False,
             nantozero=True, outfile="midcount.txt", gziped=True):

    """
    :param bwfile: bigwig file
    :param bedfile: bed file of peak/hotspots
    :param up:
    :param down:
    :param winsize: windowsize
    :param stranded:
    :param maxfilter: if score higher than maxfilter, use maxfilter
    :param minfilter: if score lower than minfilter, use minfilter
    :return:

    Given bed file and bw file count score in flanking region of bed middle site(s)

              up                mid                   down
    ==========*=================|======================*===========

    """

    bw = pyBigWig.open(bwfile, "r")

    bedio = open(bedfile, 'r')

    if gziped:

        gzfile = outfile+'.gz'

        socorefile = gzip.open(gzfile, 'w')

    else:
        socorefile = open(outfile, 'w')

    header = list()

    for j in range(0-up, down+1, winsize):

        header.append(str(j))

    headerstring = 'chromosome' + "\t"+'start'+"\t"+'end'+"\t"+"\t".join(header)+'\n'

    if gziped:
        socorefile.write(headerstring.encode('ascii'))
    else:
        socorefile.write(headerstring)

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

                else:

                    i = 'Nan'

            else:
                if i < minfilter:

                    i = minfilter

                if i > maxfilter:

                    i = maxfilter

            scorelist.append(str(i))



        if stranded:
            # Reverse - strand score
            if bedinfor[5] == '-':

                scorelist = reversed(scorelist)


        # print(chromosome, start, end, "\t".join(scorelist), sep='\t', file=socorefile)

        resultstring = chromosome+ '\t' + str(start) + '\t' + str(end) + '\t' +'\t'.join(scorelist)+"\n"

        if gziped:
                socorefile.write(resultstring.encode('ascii'))
        else:
                socorefile.write(resultstring)


    bedio.close()

    socorefile.close()

    bw.close()


def midcountmp(bwfile, bedfile, up=1000, down=1000, winsize=50,
             maxfilter=9999999999, minfilter=0, stranded=False,
             nantozero=True, outfile="midcount.txt",threads=2, gziped=True):

    """

    :param bwfile:
    :param bedfile:
    :param up:
    :param down:
    :param winsize:
    :param maxfilter:
    :param minfilter:
    :param stranded:
    :param nantozero:
    :param outfile:
    :param threads:
    :param gziped:
    :return:

    Generate matrix file
    """

    bw = pyBigWig.open(bwfile, "r")

    bedio = open(bedfile, 'r')

    if gziped:

        gzfile = outfile+'.gz'

        socorefile = gzip.open(gzfile, 'w')

    else:
        socorefile = open(outfile, 'w')

    header = list()

    for j in range(0-up, down+1, winsize):

        header.append(str(j))

    headerstring = 'chromosome' + "\t"+'start'+"\t"+'end'+"\t"+"\t".join(header)+'\n'

    if gziped:
        socorefile.write(headerstring.encode('ascii'))
    else:
        socorefile.write(headerstring)

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

    bedio.close()

    pool = Pool(threads)

    for res in pool.imap_unordered(midcountworker, workerlist):

        for resultstring in res:

            # print(resultstring, file=socorefile)
            resultstring = resultstring + "\n"

            if gziped:
                socorefile.write(resultstring.encode('ascii'))
            else:
                socorefile.write(resultstring)

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

        if (int(regioninfo[2]) - int(regioninfo[1]) + 1 ) == nbins:

            socorelist1 = bw.values(regioninfo[0], int(regioninfo[1]), int(regioninfo[2]))

        else:
            socorelist1 = bw.stats(regioninfo[0], int(regioninfo[1]), int(regioninfo[2]), nBins=nbins)

        scorelist = list()

        for i in socorelist1:

            if not i:

                if nantozero:

                    i = 0

                else:

                    i = 'Nan'

            else:
                if i < minfilter:

                    i = minfilter

                if i > maxfilter:

                    i = maxfilter

            scorelist.append(str(i))

        if stranded:
            # Reverse - strand score
            if regioninfo[3] == '-':

                scorelist = reversed(scorelist)

        scorestring = '\t'.join(regioninfo[0:3]) + '\t' +'\t'.join(scorelist)

        result.append(scorestring)

    return result
