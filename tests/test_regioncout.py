import unittest
import sys
from os import path

from bagatelle import regioncount



class TestbamTobigwig(unittest.TestCase):

    def setUp(self):

        self.inputpath = path.join(path.abspath('./'), 'inputfile')

        self.nostrandbed = path.join(self.inputpath, 'testmidcount.bed')

        self.npsbw = path.join(self.inputpath, 'npstest.bw')

        self.outfile = path.join(self.inputpath, 'npsscore.txt')

    def testmidcountnostrand(self):

        regioncount.midcount(bwfile=self.npsbw, bedfile=self.nostrandbed, up=50, down=50,
                             winsize=5, outfile=self.outfile)

    def testmidcountnostrandthread(self):

        regioncount.midcountmp(bwfile=self.npsbw, bedfile=self.nostrandbed, up=50, down=50,
                             winsize=5, outfile=self.outfile,threads=3)

if __name__ == '__main__':

    unittest.main(exit=False)