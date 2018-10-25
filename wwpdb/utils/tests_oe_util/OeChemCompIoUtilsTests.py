##
# File:    OeChemCompIoUtilsTests.py
# Author:  J. Westbrook
# Date:    23-Jan-2012
# Version: 0.001
#
# Updated:
#     6-Jun-2016 jdw general cleanup
#
#
##
"""
Test cases for persistent storage of serialized OE molecule objects.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.01"

import sys
import unittest
import traceback
import time
import os

try:
    from wwpdb.utils.oe_util.build.OeChemCompIoUtils import OeChemCompIoUtils
    skiptests = False
except ImportError:
    skiptests = True
    

@unittest.skipIf(skiptests, 'Could not import openeye')
class OeChemCompIoUtilsTests(unittest.TestCase):

    def setUp(self):
        self.__lfh = sys.stdout
        self.__verbose = True
        self.__topCachePath = "../../../../../reference/components/ligand-dict-v3"
        self.__pathList = ['../data/ATP.cif', '../data/GTP.cif', '../data/ARG.cif']
        self.__idList = ['MSE', 'GTP', 'TRP']

    def tearDown(self):
        pass

    def testMakeFromFiles(self):
        """Test case -  create OE molecules from the input chem comp definition path list.
        """
        self.__lfh.write("\nStarting %s %s\n" % (self.__class__.__name__,
                                                 sys._getframe().f_code.co_name))
        try:
            oemList = []
            oeU = OeChemCompIoUtils(topCachePath=self.__topCachePath, verbose=self.__verbose, log=self.__lfh)
            oemList = oeU.getFromPathList(self.__pathList, use3D=True, coordType='model')
            for oem in oemList:
                self.__lfh.write("Title              = %s\n" % oem.getTitle())
                self.__lfh.write("SMILES (canonical) = %s\n" % oem.getCanSMILES())
                self.__lfh.write("SMILES (isomeric)  = %s\n" % oem.getIsoSMILES())
        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testMakeFromIds(self):
        """Test case -  create OE molecules from the input chem comp definition path list.
        """
        self.__lfh.write("\nStarting %s %s\n" % (self.__class__.__name__,
                                                 sys._getframe().f_code.co_name))
        try:
            oemList = []
            oeU = OeChemCompIoUtils(topCachePath=self.__topCachePath, verbose=self.__verbose, log=self.__lfh)
            oemList = oeU.getFromIdList(self.__idList, use3D=True, coordType='model')
            for oem in oemList:
                self.__lfh.write("Title              = %s\n" % oem.getTitle())
                self.__lfh.write("SMILES (canonical) = %s\n" % oem.getCanSMILES())
                self.__lfh.write("SMILES (isomeric)  = %s\n" % oem.getIsoSMILES())
        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()


def suite1():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeChemCompIoUtilsTests("testMakeFromFiles"))
    suiteSelect.addTest(OeChemCompIoUtilsTests("testMakeFromIds"))
    return suiteSelect

if __name__ == '__main__':
    #
    mySuite = suite1()
    unittest.TextTestRunner(verbosity=2).run(mySuite)