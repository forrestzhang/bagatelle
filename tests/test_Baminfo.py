import unittest
import sys
from os import path

from bagatelle import Baminfo
from bagatelle import mhsbam


class TestBaminfo(unittest.TestCase):

    def setUp(self):

        self.inputpath = path.join(path.abspath('./'), 'inputfile')

        self.pairedbamnps = path.join(self.inputpath, 'npstest.bam')

    def test_getchrlen(self):

        baminfo = Baminfo.Baminfo(self.pairedbamnps)

        self.assertEqual(baminfo.chrlen['Chr1'], 30427671)

    def test_openmhs_pairedshort(self):

        res = mhsbam.mhsmidcount(self.pairedbamnps, chromosome='Chr1', start=1, end=300, paired=True)

        self.assertEqual(res[115], 2)

    def test_openmhs_pairedlong(self):

        res = mhsbam.mhsmidcount(self.pairedbamnps, chromosome='Chr1', start=1, end=300, maxinsert=180, paired=True)

        self.assertEqual(res[251], 3)

    def test_openmhs_pairedlong_empty(self):

        res = mhsbam.mhsmidcount(self.pairedbamnps, chromosome='Chr2', start=1, end=300, maxinsert=180, paired=True)

        self.assertFalse(res)
        # self.assertTrue(res)

if __name__ == '__main__':

    unittest.main(exit=False)