import unittest
import sys
from os import path

from bagatelle import bamTobigwig


class TestbamTobigwig(unittest.TestCase):

    def setUp(self):

        self.inputpath = path.join(path.abspath('./'), 'inputfile')

        self.npsbam = path.join(self.inputpath, 'npstest.bam')

        self.npsbw = path.join(self.inputpath, 'npstest.bw')

    def test_npsLong(self):

        bamTobigwig.chipmidtobw(bamfile=self.npsbam, bwfile=self.npsbw, paired=True)


if __name__ == '__main__':

    unittest.main(exit=False)