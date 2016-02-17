

from optparse import OptionParser
import pybedtools
import logging
import sys
import os
import pprint
from BCBio.GFF import GFFExaminer
from collections import defaultdict

def main():

    opt = opt_check(get_optparser())

    gfffile = opt.gfffile

    bedfile = opt.bedfile

    gffprefn = pybedtools.BedTool(gfffile).remove_invalid()

    gffd = dict()

    featurecluster = dict()

    featurenoparent = dict()

    examiner = GFFExaminer()

    outio = open(opt.outfile, 'w')

    openBED = open(bedfile, 'r')

    summitbedoutfile = "summit"+bedfile

    bedoutio = open(summitbedoutfile,'w')

    for line in openBED:

        if line.startswith("#"):

            continue

        strings = line.strip().split("\t")

        peakchromosome = strings[0]

        peakstart = int(strings[1])

        peakend = int(strings[2])

        mid1 = int((peakstart + peakend)/2)

        mid2 = int((peakstart + peakend)/2)

        print (peakchromosome,mid1,mid2, sep="\t", file=bedoutio)

    bedoutio.close()

    derivesfeature = dict()

    for gffinf in gffprefn:

        # a = gffinf.fields

        attrs = gffinf.attrs

        #    continue
        featuretype = gffinf[2]

        if "Parent" not in attrs:

            featurenoparent[featuretype] = 1

        if featuretype in gffd:

            gffd[featuretype] += 1

        else:

            gffd[featuretype] = 1

        if "Derives_from" in attrs:

            derivesfeature[featuretype] = 1
            # print (featuretype)
    #
    # pprint.pprint(featurenoparent)

    updown = list()

    overlap = list()

    skip = list()

    other = list()

    if opt.profile:

        profile = open(opt.profile, 'r')

        for lin in profile.readlines():

            lin = lin.rstrip('\n')

            (typenow,inf) = lin.split(':')

            if typenow == 'updown':

                updown = inf.split(',')

            if typenow == 'overlap':

                overlap = inf.split(',')

            if typenow == 'skip':

                skip = inf.split(',')

            if typenow == 'other':

                other = inf.split(',')




    else:

        for featuretype in gffd:

            while True:

                print ("#"*36)

                print ("Find ",featuretype, gffd[featuretype],"in genome")

                print ("please choose model: \n"
                       "1) calculate up and downstream, \n"
                       "2) overlap, \n"
                       "3) skip, \n"
                       "4) count this type as other")

                if featuretype == 'chromosome':

                    choose = eval(input("suggest 3: ")) or 3

                elif featuretype in derivesfeature:

                    choose = eval(input("suggest 3: ")) or 3

                elif featuretype in featurenoparent:

                    choose = eval(input("suggest 1: ")) or 1

                elif featuretype in ['exon',
                              'CDS',
                              'intron',
                              'five_prime_UTR',
                              'three_prime_UTR']:

                    choose = eval(input("suggest 2: ")) or 2

                else:

                    choose = eval(input("suggest 2: ")) or 2

                # choose = input()

                choose = int(choose)

                if choose == 1:

                    updown.append(featuretype)

                    print (featuretype, "calculate up and downstream")

                    print ()

                    break

                elif choose == 2:

                    overlap.append(featuretype)

                    print (featuretype, "overlap")

                    print ()

                    break

                elif choose == 3:

                    skip.append(featuretype)

                    print (featuretype, "skip this featuretype")

                    print ()

                    break

                elif choose == 4:

                    other.append(featuretype)

                    print (featuretype, " count this featuretype as other")

                    print ()

                    break

                else:

                    print ("Please input 1,2,3,4 model")

                    print ()

    mkintron = False

    if 'intron' not in gffd:

        while True:

            print ("Do not find intron annotation, suggest make intron annotation. y(es) or n(o)")

            intronyes = eval(input("suggest yes: ")) or 'yes'

            if intronyes=='yes' or intronyes=='y':

                mkintron = True

                break

            elif intronyes=='no' or intronyes=='n':

                mkintron = False

                break

            else:

                continue

    updonwfile = gfffile+"updown"


    updownio = open(updonwfile,'w')

    gffinio = open(gfffile,'r')

    for line in gffinio:

        if line.startswith("#"):

            continue

        line = line.rstrip('\n')

        linecontain = line.split("\t")

        if linecontain[2] in updown:

            print (line, file=updownio)

    updownio.close()

    if mkintron:

        gffinio1 = open(gfffile,'r')

        print ("make intron file")

        nointronfile = gfffile+"notinron"

        genefile = gfffile+"gene"

        nointronio = open(nointronfile,'w')

        geneio = open(genefile,'w')

        for line in gffinio1:

            if line.startswith("#"):

                continue

            line = line.rstrip('\n')

            linecontain = line.split("\t")

            if linecontain[2] in ['exon',
                          'CDS',
                          'five_prime_UTR',
                          'three_prime_UTR']:

                print (line, file=nointronio)

            if linecontain[2] in ['gene']:

                print (line, file=geneio)

        genefn = pybedtools.BedTool(genefile)

        nointronfn = pybedtools.BedTool(nointronfile)

        intronfn = genefn.subtract(nointronfn)

        intronfn.saveas('tmp_intron.gff')

        # intronfile = gfffile+"intron"

        intronin = open('tmp_intron.gff','r')

        # intronout = open (intronfile, 'w')



        gffinio2 = open(gfffile,'r')

        overlapfile = gfffile+"overlap"

        overlapio = open(overlapfile, 'w')

        for line in gffinio2:

            if line.startswith("#"):

                continue

            line = line.rstrip('\n')

            linecontain = line.split("\t")

            if linecontain[2] in overlap:

                print (line,file=overlapio)

            if linecontain[2] in updown:

                print (line, file=overlapio)

            if linecontain[2] in other:

                print (line, file=overlapio)

        for line in intronin:

            line = line.rstrip('\n')

            b = line.replace('\tgene\t', '\tintron\t')

            print (b, file=overlapio)

        overlapio.close()

        overlap.append('intron')

        # os.remove('tmp_intron.gff')

        os.remove(genefile)

        # os.remove(nointronfile)

    else:


        gffinio2 = open(gfffile,'r')

        overlapfile = gfffile+"overlap"

        overlapio = open(overlapfile, 'w')

        for line in gffinio2:

            if line.startswith("#"):

                continue

            line = line.rstrip('\n')

            linecontain = line.split("\t")

            if linecontain[2] in overlap:

                print (line, file=overlapio)

            if linecontain[2] in updown:

                print (line, file=overlapio)

            if linecontain[2] in other:

                print (line, file=overlapio)

        overlapio.close()

    print ("updown", updown, file=outio)

    print ("overlap", overlap, file=outio)

    print ("skip", skip, file=outio)

    print ("other", other, file=outio)

    overlapfn = pybedtools.BedTool(overlapfile).sort()

    updownfn = pybedtools.BedTool(updonwfile).sort()

    summitfn = pybedtools.BedTool(summitbedoutfile).sort()

    intergenic_summitfn = summitfn.subtract(overlapfn)

    intergenic_summitfn.saveas("intergenic_summitfn.txt")

    nearbyfn =intergenic_summitfn.closest(updownfn, d=True, stream=True)

    # nearbyfn.saveas("nearby.txt")

    d = defaultdict(set)

    bedfields = summitfn.field_count()

    type_idx = bedfields + 2

    bedintersectgff = summitfn.intersect(overlapfn, wao=True)

    for feature in bedintersectgff:

        featuretype = feature[type_idx]

        key = '\t'.join(feature[:bedfields])

        if featuretype in overlap:

            d[key].update([featuretype])
            # print ("overlap")

        elif featuretype in updown:

            d[key].update([featuretype])
            # print ("updown")

        elif featuretype in other:

            d[key].update(['other'])
            # print ("other")

        elif featuretype in skip:
            print (featuretype,"skip")
            # d[key].update(['.'])
            continue

        else:

            continue

            # d[key].update(['.'])

    npeaks = float(len(d))

    count_d = defaultdict(int)

    for peak, featuretypes in list(d.items()):

        if featuretypes == set('.'):

            featuretype = 'unannotated'

            continue

        else:

            featuretype = labelfilter(featuretypes)

        count_d[featuretype] += 1

    results = list(count_d.items())

    # results.sort(key=lambda x: x[1])
    results = sorted(results)
    labels, counts = list(zip(*results))

    labels = []
    counts_to_use = []

    nearpeakd = defaultdict(set)

    for nearpeak in nearbyfn:
        #Chr,peakstart, peakend, genechr,genestart,geneend, genestrand
        #print (nearpeak[0],nearpeak[1], nearpeak[2],nearpeak[bedfields], nearpeak[bedfields+3],nearpeak[bedfields+4],nearpeak[bedfields+6])

        peakkey = '\t'.join(nearpeak[:bedfields])

        if peakkey in d:

            continue

        genestrand = nearpeak[bedfields+6]

        distance = int(nearpeak[-1])

        typenow = 'error'

        if distance == 0:

            continue

        if int(nearpeak[bedfields+3]) <= int(nearpeak[1]) <= int(nearpeak[2]) <= int(nearpeak[bedfields+4]):

            print ("error")

            print (nearpeak[0],nearpeak[1], nearpeak[2],nearpeak[bedfields], nearpeak[bedfields+3],nearpeak[bedfields+4],nearpeak[bedfields+6])


        if genestrand == '+':

            if int(nearpeak[1]) >= int(nearpeak[bedfields+4]):

                #typenow = 'downstrand'
                if distance <= 1000:

                    typenow = nearpeak[bedfields+2]+"_"+'TTS_1000'

                elif distance <= 3000:

                    typenow = nearpeak[bedfields+2]+"_"+'TTS_3000'

                else:

                    typenow = 'intergentic'

            elif int(nearpeak[2]) <= int(nearpeak[bedfields+3]):

                if distance <= 1000:

                    typenow = nearpeak[bedfields+2]+"_"+'TSS_1000'

                elif distance <= 3000:

                    typenow = nearpeak[bedfields+2]+"_"+'TSS_3000'

                else:

                    typenow = 'intergentic'

            else:

                print ("error", nearpeak[0],nearpeak[1], nearpeak[2],nearpeak[bedfields], nearpeak[bedfields+3],nearpeak[bedfields+4],nearpeak[bedfields+6])

        elif genestrand == '-':

            if int(nearpeak[1]) >= int(nearpeak[bedfields+4]):

                #typenow = 'downstrand'
                if distance <= 1000:

                    typenow = nearpeak[bedfields+2]+"_"+'TSS_1000'

                elif distance <= 3000:

                    typenow = nearpeak[bedfields+2]+"_"+'TSS_3000'

                else:

                    typenow = 'intergentic'

            elif int(nearpeak[2]) <= int(nearpeak[bedfields+3]):

                if distance <= 1000:

                    typenow = nearpeak[bedfields+2]+"_"+'TTS_1000'

                elif distance <= 3000:

                    typenow = nearpeak[bedfields+2]+"_"+'TTS_3000'

                else:

                    typenow = 'intergentic'
            else:

                print ("error", nearpeak[0],nearpeak[1], nearpeak[2],nearpeak[bedfields], nearpeak[bedfields+3],nearpeak[bedfields+4],nearpeak[bedfields+6])
        else:

            print ("error", nearpeak[0],nearpeak[1], nearpeak[2],nearpeak[bedfields], nearpeak[bedfields+3],nearpeak[bedfields+4],nearpeak[bedfields+6])

        nearpeakd[peakkey].update([typenow])

    for peakid in nearpeakd:

        if peakid in d:

            print ("error peakid in nearpeakd", peakid, nearpeakd[peakid], d[peakid])

    for peakid in d:

        if peakid in nearpeakd:

            print ("error peakid in d", peakid, nearpeakd[peakid], d[peakid])

    discount = defaultdict(int)

    for peak, distypes in list(nearpeakd.items()):

        distype = labelfilter(distypes)

        discount[distype] += 1

    disres = list(discount.items())



    for label, count in results:

        print(label, count, file=outio)


    for label, count in disres:

        print(label, count, file=outio)

    outio.close()


def labelmaker(x):

    x.difference_update('.')

    label = []

    for i in list(x):

        if i == 'three_prime_UTR':

            i = "3'UTR"

        if i == 'five_prime_UTR':

            i = "5'UTR"

        label.append(i)

    return ','.join(sorted(label))


def labelfilter(x):
    x.difference_update('.')
    label = 'error'

    if 'other' in list(x):

        label = 'other'

    elif 'three_prime_UTR' in list(x):

        label = "3'UTR"

    elif 'five_prime_UTR' in list(x):

        label = "5'UTR"

    elif 'CDS' in list(x):

        label = 'CDS'

    elif 'intron' in list(x):

        label = 'intron'
    #
    # elif 'TSS_1000' in list(x):
    #
    #     label = 'TSS_1000'
    #
    # elif 'TSS_3000' in list(x):
    #
    #     label = 'TSS_3000'

    else:

        tmplab = []

        for i in list(x):

            tmplab.append(i)

        label = ','.join(sorted(tmplab))

    return label

    # return ', '.join(sorted(label))

def get_optparser():

    usage = """usage: %prog  <-g gffile> [-n outfile]
    Example %prog -g hg18.gff3 -o concerto_config.txt
    """

    overtureopt = OptionParser(version="%prog 0.1", description="", usage=usage, add_help_option=False)

    overtureopt.add_option("-h", "--help", action="help", help="show this help message and exit.")

    overtureopt.add_option("-g", "--gfffile", dest="gfffile", type="string", help="genome annotation file, gff3 format")

    overtureopt.add_option("-o", "--output", dest="outfile", type="string", help="output file, default=concerto_config.txt",
                            default="concerto_config.txt")

    overtureopt.add_option("-b", "--bedfile", dest="bedfile", type="string", help="peak file, bed format file")

    overtureopt.add_option("-p", "--profile", dest="profile", type="string", help="profile of annotation")

    return overtureopt

def opt_check(overtureopt):

    (opt, args) = overtureopt.parse_args()

    if not opt.gfffile:

        logging.error("you need input a genome annotation gff3 format file, '-g genome.gff'")

        overtureopt.print_help()

        sys.exit(1)

    if not os.path.isfile(opt.gfffile):

        logging.error("No such file: %s" % opt.gfffile)

        sys.exit(1)

    if not opt.bedfile:

        logging.error("you need input a peak bed format file, '-b test.bed'")

        overtureopt.print_help()

        sys.exit(1)

    if not os.path.isfile(opt.bedfile):

        logging.error("No such file: %s" % opt.bedfile)

        sys.exit(1)

    if opt.profile:

        if not os.path.isfile(opt.profile):

            logging.error("No such file: %s" % opt.profile)

            sys.exit(1)


    return opt

if __name__ == "__main__":

    # try:

        main()

    # except KeyboardInterrupt:
    #
    #     sys.stderr.write("User interrupt\n")
    #
    #     sys.exit(0)
