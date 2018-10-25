##
# File:    OeAlignDepictUtilsTests.py
# Author:  jdw
# Date:    30-Jun-2013
# Version: 0.001
#
# Updates:
#     4-Nov-2013 jdw add example for specifying input definition paths and output
#                    image paths.
#
##
"""
A collection of tests for the OEAlignDepictUtils and related classes which perform
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
import sys
import time
import os
import os.path
import shutil

try:
    from wwpdb.utils.oe_util.oedepict.OeAlignDepictUtils import OeDepictMCSAlign, OeDepictMCSAlignMulti, OeDepictMCSAlignSingle, OeTestMCSAlign
    skiptests = False
except ImportError:
    skiptests = True


@unittest.skipIf(skiptests, "Could not import openeye")
class OeAlignDepictUtilsTests(unittest.TestCase):

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
        #
        self.__refPathTup = ('A', '../data/A.cif', './A-ref.svg')
        self.__fitPathTupList = [('A', '../data/A.cif', './A-fit.svg'), ('T', '../data/T.cif', './T-fit.svg'),
                                 ('ATP', '../data/ATP.cif', './ATP-fit.svg'), ('GTP', '../data/GTP.cif', './GTP-fit.svg')]
        #

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
        """ Test case -  Simple pairwise MCSS alignment  -  Each aligned pair output to a separate image file
        """
        self.__lfh.write("\nStarting %s %s\n" % (self.__class__.__name__,
                                                 sys._getframe().f_code.co_name))
        try:
            oed = OeDepictMCSAlign(verbose=self.__verbose, log=self.__lfh)
            oed.setDisplayOptions(labelAtomName=True, labelAtomCIPStereo=True,
                                  labelAtomIndex=False, labelBondIndex=False,
                                  highlightStyleFit='ballAndStickInverse',
                                  bondDisplayWidth=0.5)

            oed.setRefId(self.__refId, cachePath=self.__topCachePath)
            for fitId in self.__idList:
                oed.setFitId(fitId, cachePath=self.__topCachePath)
                fName = "ref-" + self.__refId + "-trg-" + fitId + ".svg"
                aML = oed.alignPair(imagePath=fName)
                if len(aML) > 0:
                    for (rCC, rAt, tCC, tAt) in aML:
                        self.__lfh.write("%5s %-5s %5s %-5s\n" % (rCC, rAt, tCC, tAt))
        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testMCSAlignFitListDepictMulti(self):
        """Test case -  List view of pairwise MCS alignment - multipage output
        """
        self.__lfh.write("\nStarting %s %s\n" % (self.__class__.__name__,
                                                 sys._getframe().f_code.co_name))
        try:
            oed = OeDepictMCSAlignMulti(verbose=self.__verbose, log=self.__lfh)
            oed.setRefId(self.__refId, cachePath=self.__topCachePath)
            oed.setFitIdList(self.__idList, cachePath=self.__topCachePath)
            oed.setDisplayOptions(labelAtomName=True, labelAtomCIPStereo=True,
                                  labelAtomIndex=False, labelBondIndex=False,
                                  highlightStyleFit='ballAndStickInverse',
                                  gridRows=3, gridCols=3, bondDisplayWidth=0.5)
            imageFile = "list-example-mcs-alignment.pdf"
            aML = oed.alignOneWithListMulti(imagePath=imageFile)
            if len(aML) > 0:
                for (rCC, rAt, tCC, tAt) in aML:
                    self.__lfh.write("%5s %-5s %5s %-5s\n" % (rCC, rAt, tCC, tAt))

        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testMCSAlignPairListDepictPortrait(self):
        """Test case -  List view MCS alignment using pair id list input
        """
        self.__lfh.write("\nStarting %s %s\n" % (self.__class__.__name__,
                                                 sys._getframe().f_code.co_name))
        try:
            oed = OeDepictMCSAlignMulti(verbose=self.__verbose, log=self.__lfh)
            oed.setPairIdList(self.__pairIdList, cachePath=self.__topCachePath)
            oed.setDisplayOptions(labelAtomName=True, labelAtomCIPStereo=True,
                                  labelAtomIndex=False, labelBondIndex=False,
                                  highlightStyleFit='ballAndStickInverse',
                                  pageOrientation='portrait',
                                  gridRows=4, bondDisplayWidth=0.5)

            imageFile = "pair-list-example-mcs-alignment-portrait.pdf"
            aML = oed.alignPairListMulti(imagePath=imageFile)
            if len(aML) > 0:
                for (rCC, rAt, tCC, tAt) in aML:
                    self.__lfh.write("%5s %-5s %5s %-5s\n" % (rCC, rAt, tCC, tAt))
        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testMCSAlignRnaPairListDepict(self):
        """Test case -  Modified RNA nucleotide alignment with parent nucleotied using pair list input
        """
        self.__lfh.write("\nStarting %s %s\n" % (self.__class__.__name__, sys._getframe().f_code.co_name))
        try:
            pairIdList = self.__readPairList(fn=self.__rnaPairFile)
            oed = OeDepictMCSAlignMulti(verbose=self.__verbose, log=self.__lfh)
            oed.setPairIdList(pairIdList, cachePath=self.__topCachePath)
            oed.setDisplayOptions(labelAtomName=True, labelAtomCIPStereo=True,
                                  labelAtomIndex=False, labelBondIndex=False,
                                  highlightStyleFit='ballAndStickInverse',
                                  bondDisplayWidth=0.5)

            imageFile = "rna-modified-pair-alignment.pdf"
            aML = oed.alignPairListMulti(imagePath=imageFile)
            if len(aML) > 0:
                for (rCC, rAt, tCC, tAt) in aML:
                    self.__lfh.write("%5s %-5s %5s %-5s\n" % (rCC, rAt, tCC, tAt))
        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()

    #
    def testMCSAlignAtomMap(self):
        """Test case -  match test with return of atom maps
        """
        self.__lfh.write("\nStarting %s %s\n" % (self.__class__.__name__, sys._getframe().f_code.co_name))
        try:
            pairIdList = self.__readPairList(fn=self.__rnaPairFile)
            oed = OeTestMCSAlign(verbose=self.__verbose, log=self.__lfh)
            for refId, fitId in pairIdList:
                oed.setRefId(refId, cachePath=self.__topCachePath)
                oed.setFitId(fitId, cachePath=self.__topCachePath)
                aML = oed.doAlign()
                if len(aML) > 0:
                    self.__lfh.write("Match length %3d for: %s %s\n" % (len(aML), refId, fitId))
                else:
                    self.__lfh.write("Match failed for: %s %s\n" % (refId, fitId))
        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testMCSRelaxAlignPairDepict(self):
        """Test case -  Relaxed pairwise MCSS alignment  -
        """
        self.__lfh.write("\nStarting %s %s\n" % (self.__class__.__name__,
                                                 sys._getframe().f_code.co_name))
        try:
            oed = OeDepictMCSAlign(verbose=self.__verbose, log=self.__lfh)
            oed.setSearchType(sType='relaxed')
            oed.setRefPath(refId='PRD_000225', ccPath='./examples/PRDCC_000225.cif', title="PRD_000225", suppressHydrogens=False)
            oed.setFitPath(fitId='LD_LDI_990', ccPath='./examples/L_LDI_990_.comp.cif', title='L_LDI_990', suppressHydrogens=False)
            fName = "relaxed-fit.svg"
            oed.setDisplayOptions(imageSizeX=2000, imageSizeY=1000, labelAtomName=True, labelAtomCIPStereo=True,
                                  labelAtomIndex=False, labelBondIndex=False,
                                  highlightStyleFit='ballAndStickInverse',
                                  bondDisplayWidth=1.0)
            aML = oed.alignPair(imagePath=fName)
            if len(aML) > 0:
                for (rCC, rAt, tCC, tAt) in aML:
                    self.__lfh.write("%5s %-5s %5s %-5s\n" % (rCC, rAt, tCC, tAt))
        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testMCSAlignListMultiDepict(self):
        """Test case -  List view of MCS alignment --- on multi-pages --
        """
        self.__lfh.write("\nStarting %s %s\n" % (self.__class__.__name__,
                                                 sys._getframe().f_code.co_name))
        try:
            oed = OeDepictMCSAlignMulti(verbose=self.__verbose, log=self.__lfh)
            oed.setDisplayOptions(labelAtomName=True, labelAtomCIPStereo=True,
                                  labelAtomIndex=False, labelBondIndex=False,
                                  highlightStyleFit='ballAndStickInverse',
                                  gridRows=3, gridCols=3, bondDisplayWidth=0.5)
            oed.setPairIdList(self.__pairIdList, cachePath=self.__topCachePath)
            imageFile = "mcs-align-with-list-multi.pdf"
            aML = oed.alignOneWithListMulti(imagePath=imageFile)
            if len(aML) > 0:
                for (rCC, rAt, tCC, tAt) in aML:
                    self.__lfh.write("%5s %-5s %5s %-5s\n" % (rCC, rAt, tCC, tAt))

        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testMCSAlignListDepict(self):
        """Test case -  List view of MCS alignment on single image -
        """
        self.__lfh.write("\nStarting %s %s\n" % (self.__class__.__name__,
                                                 sys._getframe().f_code.co_name))
        try:
            oed = OeDepictMCSAlign(verbose=self.__verbose, log=self.__lfh)
            oed.setDisplayOptions(imageSizeX=1500, imageSizeY=1500, labelAtomName=True, labelAtomCIPStereo=True,
                                  labelAtomIndex=False, labelBondIndex=False,
                                  highlightStyleFit='ballAndStickInverse',
                                  gridRows=3, gridCols=3, bondDisplayWidth=1.0)
            oed.setPairIdList(self.__pairIdList, cachePath=self.__topCachePath)
            imageFile = "mcs-align-with-list-single.svg"
            aML = oed.alignOneWithList(imagePath=imageFile)
            if len(aML) > 0:
                for (rCC, rAt, tCC, tAt) in aML:
                    self.__lfh.write("%5s %-5s %5s %-5s\n" % (rCC, rAt, tCC, tAt))

        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testMCSAlignListDepictSingle(self):
        """Test case -  View MCS alignment in single image files.

                        Input details specified as chemical component Id pair lists -
        """
        self.__lfh.write("\nStarting %s %s\n" % (self.__class__.__name__,
                                                 sys._getframe().f_code.co_name))
        try:
            oed = OeDepictMCSAlignSingle(verbose=self.__verbose, log=self.__lfh)
            oed.setDisplayOptions(imageSizeX=500, imageSizeY=500, labelAtomName=True, labelAtomCIPStereo=True,
                                  labelAtomIndex=False, labelBondIndex=False,
                                  highlightStyleFit='ballAndStickInverse',
                                  bondDisplayWidth=1.0)
            oed.setPairIdList(self.__pairIdList, cachePath=self.__topCachePath)
            aML = oed.alignOneWithList()
            if len(aML) > 0:
                for (rCC, rAt, tCC, tAt) in aML:
                    self.__lfh.write("%5s %-5s %5s %-5s\n" % (rCC, rAt, tCC, tAt))

        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testMCSAlignPathListDepictSingle(self):
        """Test case -  View MCS alignments in single image files.

                        Input details specified chemical component and image file paths.

                        Image file paths must end with a recognized image format (e.g. svg, png, jpg)
        """
        self.__lfh.write("\nStarting %s %s\n" % (self.__class__.__name__,
                                                 sys._getframe().f_code.co_name))
        try:
            # Use inverse  highlighting of matching/non-matching atoms/bonds -
            hOpt = 'inverse'
            #
            oed = OeDepictMCSAlignSingle(verbose=self.__verbose, log=self.__lfh)
            if hOpt == 'inverse':
                oed.setDisplayOptions(imageSizeX=500, imageSizeY=500, labelAtomName=True, labelAtomCIPStereo=True,
                                      labelAtomIndex=False, labelBondIndex=False,
                                      highlightStyleFit='ballAndStickInverse',
                                      bondDisplayWidth=1.0)
            elif hOpt == 'match':
                oed.setDisplayOptions(imageSizeX=500, imageSizeY=500, labelAtomName=True, labelAtomCIPStereo=True,
                                      labelAtomIndex=False, labelBondIndex=False,
                                      highlightStyleFit='ballAndStick',
                                      bondDisplayWidth=1.0)

            else:
                oed.setDisplayOptions(imageSizeX=500, imageSizeY=500, labelAtomName=True, labelAtomCIPStereo=True,
                                      labelAtomIndex=False, labelBondIndex=False,
                                      highlightStyleFit='None',
                                      bondDisplayWidth=1.0)

            oed.setRefPath(self.__refPathTup[0], self.__refPathTup[1], title=self.__refPathTup[0], imagePath=self.__refPathTup[2])
            for fitPathTup in self.__fitPathTupList:
                oed.addFitPath(fitPathTup[0], fitPathTup[1], title=fitPathTup[0], imagePath=fitPathTup[2])

            aML = oed.alignOneWithList()
            #
            # Write out no atom correspondences --
            #
            if len(aML) > 0:
                for (rCC, rAt, tCC, tAt) in aML:
                    self.__lfh.write("%5s %-5s %5s %-5s\n" % (rCC, rAt, tCC, tAt))

        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()


def suiteAlignTypeTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeAlignDepictUtilsTests("testMCSAlignAtomMap"))
    suiteSelect.addTest(OeAlignDepictUtilsTests("testMCSRelaxAlignPairDepict"))
    return suiteSelect


def suiteAlignPair():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeAlignDepictUtilsTests("testMCSAlignPairDepict"))
    suiteSelect.addTest(OeAlignDepictUtilsTests("testMCSAlignFitListDepictMulti"))
    suiteSelect.addTest(OeAlignDepictUtilsTests("testMCSAlignPairListDepictPortrait"))
    suiteSelect.addTest(OeAlignDepictUtilsTests("testMCSAlignRnaPairListDepict"))
    return suiteSelect


def suiteAlignWithList():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeAlignDepictUtilsTests("testMCSAlignListMultiDepict"))
    suiteSelect.addTest(OeAlignDepictUtilsTests("testMCSAlignListDepict"))
    return suiteSelect


def suiteAlignWithListSingle():
    suiteSelect = unittest.TestSuite()
    # suiteSelect.addTest(OeAlignDepictUtilsTests("testMCSAlignListDepictSingle"))
    suiteSelect.addTest(OeAlignDepictUtilsTests("testMCSAlignPathListDepictSingle"))
    return suiteSelect


if __name__ == '__main__':
    #
    if False:
        mySuite1 = suiteAlignPair()
        unittest.TextTestRunner(verbosity=2).run(mySuite1)
        #
        mySuite1 = suiteAlignWithList()
        unittest.TextTestRunner(verbosity=2).run(mySuite1)

        mySuite1 = suiteAlignTypeTests()
        unittest.TextTestRunner(verbosity=2).run(mySuite1)
        #
        mySuite1 = suiteAlignWithListSingle()
        unittest.TextTestRunner(verbosity=2).run(mySuite1)

    mySuite1 = suiteAlignWithListSingle()
    unittest.TextTestRunner(verbosity=2).run(mySuite1)