##
# File:    OeAlignDepictTests.py
# Author:  jdw
# Date:    2-Oct-2011
# Version: 0.001
#
# Updates:
#  11-Mar-2012 jdw update to revised calling interface
#   1-Nov-2014 jdw test for SDF alignment -
#   6-Jun-2016 jdw general cleanup
#
##
"""
A collection of tests for the OEAlignDepict and related classes which perform
MCSS comparison and depiction.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.01"


import sys
import unittest
import traceback
import os

try:
    from wwpdb.utils.oe_util.oedepict.OeAlignDepict import OeDepictMCSAlign
    skiptests = False
except ImportError:
    skiptests = True
    

@unittest.skipIf(skiptests, 'Could not import openeye')
class OeAlignDepictTests(unittest.TestCase):

    def setUp(self):
        self.__lfh = sys.stderr
        self.__verbose = True
        #
        # Chemical component repository path -
        self.__topCachePath = "../../../../../reference/components/ligand-dict-v3"
        self.__rnaPairFile = './examples/rna-linking-components.txt'
        #
        self.__refId = 'C'
        #
        self.__idList = ['cg1', 'atp', 'gtp', 'A', 'C', 'G', 'DG']
        self.__pairIdList = [('c', 'cg1'), ('c', 'atp'), ('c', 'gtp'), ('c', 'A'), ('c', 'C'), ('c', 'G'), ('c', 'DG')]

    def tearDown(self):
        pass

    def __readPairList(self, fn='./examples/rna-linking-components.txt'):
        pairList = []
        ifh = open(fn, 'r')
        for line in ifh:
            fields = line.split()
            pairList.append((fields[1], fields[0]))
        ifh.close()
        return pairList

    def testMCSAlignPairDepict(self):
        """Test case -  Simple pairwise MCSS alignment  -  Each aligned pair output to a separate image file
        """
        self.__lfh.write("\nStarting %s %s\n" % (self.__class__.__name__,
                                                 sys._getframe().f_code.co_name))
        try:
            oed = OeDepictMCSAlign(verbose=self.__verbose, log=self.__lfh)
            oed.setRefId(self.__refId, cachePath=self.__topCachePath)
            for fitId in self.__idList:
                oed.setFitId(fitId, cachePath=self.__topCachePath)
                fName = "ref-" + self.__refId + "-trg-" + fitId + ".png"
                aML = oed.alignPair(imagePath=fName)
                if len(aML) > 0:
                    for (rCC, rAt, tCC, tAt) in aML:
                        self.__lfh.write("%5s %-5s %5s %-5s\n" % (rCC, rAt, tCC, tAt))
        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testMCSAlignListDepict(self):
        """Test case -  List view of pairwise MCS alignment - multipage output
        """
        self.__lfh.write("\nStarting %s %s\n" % (self.__class__.__name__,
                                                 sys._getframe().f_code.co_name))
        try:
            oed = OeDepictMCSAlign(verbose=self.__verbose, log=self.__lfh)
            oed.setRefId(self.__refId, cachePath=self.__topCachePath)
            oed.setFitIdList(self.__idList, cachePath=self.__topCachePath)
            imageFile = "list-example-mcs-alignment.pdf"
            aML = oed.alignPairList(imagePath=imageFile)
            if len(aML) > 0:
                for (rCC, rAt, tCC, tAt) in aML:
                    self.__lfh.write("%5s %-5s %5s %-5s\n" % (rCC, rAt, tCC, tAt))

        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testMCSAlignPairListDepict(self):
        """Test case -  List view of  pairwise MCS alignment using pair id list input
        """
        self.__lfh.write("\nStarting %s %s\n" % (self.__class__.__name__,
                                                 sys._getframe().f_code.co_name))
        try:
            oed = OeDepictMCSAlign(verbose=self.__verbose, log=self.__lfh)
            oed.setPairIdList(self.__pairIdList, cachePath=self.__topCachePath)
            imageFile = "pair-list-example-mcs-alignment.pdf"
            aML = oed.alignPairList(imagePath=imageFile)
            if len(aML) > 0:
                for (rCC, rAt, tCC, tAt) in aML:
                    self.__lfh.write("%5s %-5s %5s %-5s\n" % (rCC, rAt, tCC, tAt))

        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testMCSAlignRnaPairListDepict(self):
        """Test case -  Modified RNA nucleotide alignment with parent nucleotied using pair list input
        """
        self.__lfh.write("\nStarting %s %s\n" % (self.__class__.__name__,
                                                 sys._getframe().f_code.co_name))
        try:
            pairIdList = self.__readPairList(fn=self.__rnaPairFile)
            oed = OeDepictMCSAlign(verbose=self.__verbose, log=self.__lfh)
            oed.setPairIdList(pairIdList, cachePath=self.__topCachePath)
            imageFile = "rna-modified-pair-alignment.pdf"
            aML = oed.alignPairList(imagePath=imageFile)
            if len(aML) > 0:
                for (rCC, rAt, tCC, tAt) in aML:
                    self.__lfh.write("%5s %-5s %5s %-5s\n" % (rCC, rAt, tCC, tAt))

        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testMCSAlignAtomMap(self):
        """Test case -  match test with return of atom maps
        """
        self.__lfh.write("\nStarting %s %s\n" % (self.__class__.__name__,
                                                 sys._getframe().f_code.co_name))
        try:
            pairIdList = self.__readPairList(fn=self.__rnaPairFile)
            oed = OeDepictMCSAlign(verbose=self.__verbose, log=self.__lfh)
            for refId, fitId in pairIdList:
                oed.setRefId(refId, cachePath=self.__topCachePath)
                oed.setFitId(fitId, cachePath=self.__topCachePath)
                aML = oed.testAlign()
                if len(aML) > 0:
                    self.__lfh.write("Match suceeded for: %s %s\n" % (refId, fitId))
                else:
                    self.__lfh.write("Match failed for: %s %s\n" % (refId, fitId))
        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testSdfMCSAlignAtomMap(self):
        """Test case -  match test with return of atom maps for foreign SDF with a CC definition
        """
        self.__lfh.write("\nStarting %s %s\n" % (self.__class__.__name__, sys._getframe().f_code.co_name))
        try:
            # self.__extPairTup=('../data/ATP.sdf','ATP')
            self.__extPairTup = ('../data/ATP.sdf', '../data/ATP.cif')
            refPath = self.__extPairTup[0]
            fitId = self.__extPairTup[1]
            fitPath = self.__extPairTup[1]
            oed = OeDepictMCSAlign(verbose=self.__verbose, log=self.__lfh)
            # oed.setSearchType(sType='exact')
            oed.setSearchType(sType='exact')
            #
            oed.setRefPath(refPath, type="SDF")
            # oed.setFitId(fitId,cachePath=self.__topCachePath)
            oed.setFitPath(fitPath, title="ATP")
            aML = oed.testAlign(suppressHydrogens=False, unique=False, minFrac=0.9)
            if len(aML) > 0:
                self.__lfh.write("Match suceeded for: %s %s\n" % (refPath, fitId))
                for ii, aM in enumerate(aML):
                    self.__lfh.write("Mapping: %4d %s\n" % (ii, aM))
            else:
                self.__lfh.write("Match failed for: %s %s\n" % (refPath, fitId))
        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testMCSRelaxAlignPairDepict(self):
        """ Test case -  Relaxed pairwise MCSS alignment  -
        """
        self.__lfh.write("\nStarting %s %s\n" % (self.__class__.__name__,
                                                 sys._getframe().f_code.co_name))
        try:
            oed = OeDepictMCSAlign(verbose=self.__verbose, log=self.__lfh)
            oed.setSearchType(sType='relaxed')
            oed.setRefPath(ccPath='./examples/PRDCC_000225.cif', title="PRD_000225", suppressHydrogens=False)
            oed.setFitPath(ccPath='./examples/L_LDI_990_.comp.cif', title='L_LDI_990', suppressHydrogens=False)
            fName = "relaxed-fit.png"
            aML = oed.alignPair(imagePath=fName)
            if len(aML) > 0:
                for (rCC, rAt, tCC, tAt) in aML:
                    self.__lfh.write("%5s %-5s %5s %-5s\n" % (rCC, rAt, tCC, tAt))
        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()


def suiteAlignPairRelax():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeAlignDepictTests("testMCSRelaxAlignPairDepict"))
    return suiteSelect


def suiteAlignExtPairRelax():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeAlignDepictTests("testSdfMCSAlignAtomMap"))
    return suiteSelect


def suiteAlignPair():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeAlignDepictTests("testMCSAlignPairDepict"))
    suiteSelect.addTest(OeAlignDepictTests("testMCSAlignListDepict"))
    suiteSelect.addTest(OeAlignDepictTests("testMCSAlignPairListDepict"))
    suiteSelect.addTest(OeAlignDepictTests("testMCSAlignRnaPairListDepict"))
    suiteSelect.addTest(OeAlignDepictTests("testMCSAlignAtomMap"))
    return suiteSelect


def suite():
    return unittest.makeSuite(OeAlignDepictTests, 'test')

if __name__ == '__main__':
    # unittest.main()
    #
    if (False):
        mySuite1 = suiteAlignPair()
        unittest.TextTestRunner(verbosity=2).run(mySuite1)
        #
        mySuite1 = suiteAlignPairRelax()
        unittest.TextTestRunner(verbosity=2).run(mySuite1)
        #

    mySuite1 = suiteAlignExtPairRelax()
    unittest.TextTestRunner(verbosity=2).run(mySuite1)