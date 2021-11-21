import os
import unittest as test
import dataNavigation as nav

class TestStringMethods(test.TestCase):

    def allSpeciesIncludedTest(self):
        directory = 'data/hfe/training/'
        for file in os.listdir(directory):
            file_name = file.replace('.txt', '')
            self.run(nav.getSpeciesDivergence(file_name, 'Artibeus jamaicensis'))