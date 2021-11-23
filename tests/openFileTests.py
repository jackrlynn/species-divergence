import os
import unittest as test
import dataNavigation as nav

class TestStringMethods(test.TestCase):

    def testAllSpeciesIncludedTraining(self):
        directory = '/data/hfe/training/'
        for file in os.listdir(directory):
            file_name = file.replace('.txt', '')
            self.run(nav.getSpeciesDivergence(file_name, 'Artibeus jamaicensis'))

    def testAllSpeciesIncludedTesting(self):
        directory = '/data/hfe/testing/'
        for file in os.listdir(directory):
            f1 = file.replace('.txt', '')
            for file2 in os.listdir(directory):
                f2 = file2.replace('.txt', '')
                self.assertTrue(nav.getSpeciesDivergence(f1, f2) > -1)