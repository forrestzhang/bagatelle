import unittest
import sys
from os import path

from bagatelle import regioncount



class TestbamTobigwig(unittest.TestCase):

    def setUp(self):

        self.inputpath = path.join(path.abspath('./'), 'inputfile')

        self.nostrandbed = path.join(self.inputpath, 'testmidcount.bed')

        self.strandbed = path.join(self.inputpath, 'testmidcountstrand.bed')

        self.npsbw = path.join(self.inputpath, 'npstest.bw')

        self.outfile = path.join(self.inputpath, 'npsscore.txt')

        self.outfilenth1 = path.join(self.inputpath, 'npsscore1th.txt')

        self.outfilest = path.join(self.inputpath, 'npsscorest.txt')

    def testmidcountnostrand(self):

        regioncount.midcount(bwfile=self.npsbw, bedfile=self.nostrandbed, up=50, down=50,
                             winsize=5, outfile=self.outfilenth1)

    def testmidcountnostrandthread(self):

        regioncount.midcountmp(bwfile=self.npsbw, bedfile=self.nostrandbed, up=50, down=50,
                             winsize=1, outfile=self.outfile,threads=3, gziped=False)

    def testmidcountnostrandthreadgzip(self):

        regioncount.midcountmp(bwfile=self.npsbw, bedfile=self.nostrandbed, up=50, down=50,
                             winsize=5, outfile=self.outfile,threads=3)

    def testmidcountstrand(self):

        regioncount.midcount(bwfile=self.npsbw, bedfile=self.strandbed, up=50, down=50,
                             winsize=5, outfile=self.outfilest, stranded=True)

    def testmidcountstrandmp(self):

        regioncount.midcountmp(bwfile=self.npsbw, bedfile=self.strandbed, up=50, down=50,
                             winsize=1, outfile=self.outfilest, stranded=True)

if __name__ == '__main__':

    unittest.main(exit=False)