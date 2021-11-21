import unittest as test
import tools.hashbox as hb

class TestStringMethods(test.TestCase):

    def testHash(self):
        self.assertEqual(hb.hash("ACCAT"), 788)
        self.assertEqual(hb.hash("ACGT"), 228)
        self.assertEqual(hb.hash("AAAAAAAAA"), 0)
        self.assertEqual(hb.hash(""), 0)
        self.assertEqual(hb.hash("GCAATTTA"), 16134)
        self.assertEqual(hb.hash("GCGT"), 230)
        self.assertEqual(hb.dehash(788, 5), "ACCAT")
        self.assertEqual(hb.dehash(228, 4), "ACGT")
        self.assertEqual(hb.dehash(0, 9), "AAAAAAAAA")
        self.assertEqual(hb.dehash(0, 0), "")
        self.assertEqual(hb.dehash(16134, 8), "GCAATTTA")
        self.assertEqual(hb.dehash(230, 4), "GCGT")
        self.assertEqual(hb.dehash(hb.hash("ACCAT"), 5), "ACCAT")
        self.assertEqual(hb.dehash(hb.hash("ACGT"), 4), "ACGT")
        self.assertEqual(hb.dehash(hb.hash("AAAAAAAAA"), 9), "AAAAAAAAA")
        self.assertEqual(hb.dehash(hb.hash(""), 0), "")
        self.assertEqual(hb.dehash(hb.hash("GCAATTTA"), 8), "GCAATTTA")
        self.assertEqual(hb.dehash(hb.hash("GCGT"), 4), "GCGT")

    def testHashCreation(self):
        self.assertEqual(hb.createHashList("GCAATTTA", 4, False), [hb.hash("GCAA"), hb.hash("CAAT"), hb.hash("AATT"),
                                                                   hb.hash("ATTT"), hb.hash("TTTA")])
        self.assertEqual(hb.createHashList("GCAATTTA", 4, True), [hb.hash("GCAA"), hb.hash("TTTA"), hb.hash("CAAT"),
                                                                  hb.hash("AATT"), hb.hash("ATTT")])

    def testAlignAgnosticScore(self):
        self.assertEqual(hb.getSimilarity("CCCAAACCC", "CCC", 3), 2)
        self.assertEqual(hb.getSimilarity("CCCAAACCC", "CCCAAA", 3), 5)
        self.assertEqual(hb.getSimilarity("CCCAAACCC", "CCCAAA", 6), 1)
        self.assertEqual(hb.getSimilarity("ACGTACGTACGTACGT", "TTACTTACTT", 4), 0)
        self.assertEqual(hb.getSimilarity("", "", 1), 0)
        self.assertEqual(hb.getSimilarity("CCCAAACCC", "CCC", 1), 18)
        self.assertEqual(hb.getSimilarity("AATTAATCGTGTACTGA", "AATTATACTGTAATCG", 5), 4)

    def testConvertPos(self):
        self.assertEqual(hb.convertPosToAlignedPos(1, "A--A"), 3)
        self.assertEqual(hb.convertPosToAlignedPos(0, "A--A"), 0)
        self.assertEqual(hb.convertPosToAlignedPos(2, "CCC------"), 2)
        self.assertEqual(hb.convertPosToAlignedPos(3, "---AA-A-C---CC"), 8)
        self.assertEqual(hb.convertPosToAlignedPos(4, "---AA-A-C---CC"), 12)

    def testAlignConsciousScore(self):
        self.assertEqual(hb.getSimilarityExcludeAlignment("CCCAAACCC", "CCC", 3, "CCCAAACCC", "CCC------"), 1)
        self.assertEqual(hb.getSimilarityExcludeAlignment("CCCAAACCC", "CCC", 3, "CCCAAACCC", "------CCC"), 1)
        self.assertEqual(hb.getSimilarityExcludeAlignment("CCCAAACCC", "CCCAAA", 3, "CCCAAACCC", "CCCAAA---"), 1)
        self.assertEqual(hb.getSimilarityExcludeAlignment("CCCAAACCC", "CCCAAA", 3, "CCCAAACCC", "---AAACCC"), 1)
        self.assertEqual(hb.getSimilarityExcludeAlignment("CCCAAACCC", "CCCAAA", 6, "CCCAAACCC", "CCCAAA---"), 0)
        #self.assertEqual(hb.getSimilarityExcludeAlignment("CCCAAACCC", "AAACCC", 6, "CCCAAACCC", "---AAACCC"), 0)
        self.assertEqual(hb.getSimilarityExcludeAlignment("CCCAAACCC", "CCC", 1, "CCCAAACCC", "CCC------"), 15)
