##
# File:    OeAlignDepictTests.py
# Author:  jdw
# Date:    2-Oct-2011
# Version: 0.001
#
# Updates:
#  11-Mar-2012 jdw update to revised calling interface
#   1-Nov-2014 jdw test for SDF alignment -
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


import unittest
import traceback
import sys
import platform
import os
import os.path

try:
    from openeye.oechem import OEFloatArray  # noqa: F401 pylint: disable=unused-import

    skiptests = False
except ImportError:
    skiptests = True

if not skiptests:
    from wwpdb.utils.oe_util.oedepict.OeAlignDepict import OeDepictMCSAlign


@unittest.skipIf(skiptests, "Could not import openeye")
class OeAlignDepictTests(unittest.TestCase):
    def setUp(self):
        self.__lfh = sys.stderr
        self.__verbose = False
        #

        self.__here = os.path.abspath(os.path.dirname(__file__))
        self.__examples = os.path.join(self.__here, "examples")
        self.__testoutput = os.path.join(self.__here, "test-output", platform.python_version())
        self.__datadir = os.path.join(self.__here, "data")
        if not os.path.exists(self.__testoutput):
            os.makedirs(self.__testoutput)

        # Chemical component repository path -
        self.__topCachePath = os.path.join(self.__here, "ligand-dict-v3")
        self.__rnaPairFile = os.path.join(self.__examples, "rna-linking-components.txt")
        #
        self.__refId = "C"
        #
        self.__idList = ["cg1", "atp", "gtp", "A", "C", "G", "DG"]
        self.__pairIdList = [("c", "cg1"), ("c", "atp"), ("c", "gtp"), ("c", "A"), ("c", "C"), ("c", "G"), ("c", "DG")]

    def tearDown(self):
        pass

    def __readPairList(self, fn="./examples/rna-linking-components.txt"):
        pairList = []
        ifh = open(fn, "r")
        for line in ifh:
            fields = line.split()
            pairList.append((fields[1], fields[0]))
        ifh.close()
        return pairList

    def testMCSAlignPairDepict(self):
        """Test case -  Simple pairwise MCSS alignment  -  Each aligned pair output to a separate image file"""
        self.__lfh.write("\nStarting OeAlignDepictTests tstMCSalignPairDepict\n")
        try:
            oed = OeDepictMCSAlign(verbose=self.__verbose, log=self.__lfh)
            oed.setRefId(self.__refId, cachePath=self.__topCachePath)
            for fitId in self.__idList:
                oed.setFitId(fitId, cachePath=self.__topCachePath)
                fName = os.path.join(self.__testoutput, "ref-" + self.__refId + "-trg-" + fitId + ".png")
                aML = oed.alignPair(imagePath=fName)
                if len(aML) > 0:
                    for (rCC, rAt, tCC, tAt) in aML:
                        self.__lfh.write("%5s %-5s %5s %-5s\n" % (rCC, rAt, tCC, tAt))
        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testMCSAlignListDepict(self):
        """Test case -  List view of pairwise MCS alignment - multipage output"""
        self.__lfh.write("\nStarting OeAlignDepictTests testMCSAlignListDepict\n")
        try:
            oed = OeDepictMCSAlign(verbose=self.__verbose, log=self.__lfh)
            oed.setRefId(self.__refId, cachePath=self.__topCachePath)
            oed.setFitIdList(self.__idList, cachePath=self.__topCachePath)
            imageFile = os.path.join(self.__testoutput, "list-example-mcs-alignment.pdf")
            aML = oed.alignPairList(imagePath=imageFile)
            if len(aML) > 0:
                for (rCC, rAt, tCC, tAt) in aML:
                    self.__lfh.write("%5s %-5s %5s %-5s\n" % (rCC, rAt, tCC, tAt))

        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testMCSAlignPairListDepict(self):
        """Test case -  List view of  pairwise MCS alignment using pair id list input"""
        self.__lfh.write("\nStarting OeAlignDepictTests testMCSAlignPairListDepict\n")
        try:
            oed = OeDepictMCSAlign(verbose=self.__verbose, log=self.__lfh)
            oed.setPairIdList(self.__pairIdList, cachePath=self.__topCachePath)
            imageFile = os.path.join(self.__testoutput, "pair-list-example-mcs-alignment.pdf")
            aML = oed.alignPairList(imagePath=imageFile)
            if len(aML) > 0:
                for (rCC, rAt, tCC, tAt) in aML:
                    self.__lfh.write("%5s %-5s %5s %-5s\n" % (rCC, rAt, tCC, tAt))

        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testMCSAlignRnaPairListDepict(self):
        """Test case -  Modified RNA nucleotide alignment with parent nucleotied using pair list input"""
        self.__lfh.write("\nStarting OeAlignDepictTests testMCSAlignRnaPairListDepict\n")
        try:
            pairIdList = self.__readPairList(fn=self.__rnaPairFile)
            oed = OeDepictMCSAlign(verbose=self.__verbose, log=self.__lfh)
            oed.setPairIdList(pairIdList, cachePath=self.__topCachePath)
            imageFile = os.path.join(self.__testoutput, "rna-modified-pair-alignment.pdf")
            aML = oed.alignPairList(imagePath=imageFile)
            if len(aML) > 0:
                for (rCC, rAt, tCC, tAt) in aML:
                    self.__lfh.write("%5s %-5s %5s %-5s\n" % (rCC, rAt, tCC, tAt))

        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testMCSAlignAtomMap(self):
        """Test case -  match test with return of atom maps"""
        self.__lfh.write("\nStarting OeAlignDepictTests testMCSAlignAtomMap\n")
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
        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testSdfMCSAlignAtomMap(self):
        """Test case -  match test with return of atom maps for foreign SDF with a CC definition"""
        self.__lfh.write("\nStarting OeAlignDepictTests testSdfMCSAlignAtomMap\n")
        try:
            extPairTup = (os.path.join(self.__datadir, "ATP.sdf"), os.path.join(self.__datadir, "ATP.cif"))
            refPath = extPairTup[0]
            fitId = extPairTup[1]
            fitPath = extPairTup[1]
            oed = OeDepictMCSAlign(verbose=self.__verbose, log=self.__lfh)
            # oed.setSearchType(sType='exact')
            oed.setSearchType(sType="exact")
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
        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testMCSRelaxAlignPairDepict(self):
        """Test case -  Relaxed pairwise MCSS alignment  -"""
        self.__lfh.write("\nStarting OeAlignDepictTests testMCSRelaceAlignPairDepict\n")
        try:
            oed = OeDepictMCSAlign(verbose=self.__verbose, log=self.__lfh)
            oed.setSearchType(sType="relaxed")
            oed.setRefPath(ccPath=os.path.join(self.__examples, "PRDCC_000225.cif"), title="PRD_000225", suppressHydrogens=False)
            oed.setFitPath(ccPath=os.path.join(self.__examples, "L_LDI_990_.comp.cif"), title="L_LDI_990", suppressHydrogens=False)
            fName = os.path.join(self.__testoutput, "relaxed-fit.png")
            aML = oed.alignPair(imagePath=fName)
            if len(aML) > 0:
                for (rCC, rAt, tCC, tAt) in aML:
                    self.__lfh.write("%5s %-5s %5s %-5s\n" % (rCC, rAt, tCC, tAt))
        except Exception as e:  # noqa: F841 pylint: disable=unused-variable
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
    return unittest.makeSuite(OeAlignDepictTests, "test")


if __name__ == "__main__":
    # unittest.main()
    #
    mySuite1 = suiteAlignPair()
    unittest.TextTestRunner(verbosity=2).run(mySuite1)
    #
    mySuite1 = suiteAlignPairRelax()
    unittest.TextTestRunner(verbosity=2).run(mySuite1)
    #

    mySuite1 = suiteAlignExtPairRelax()
    unittest.TextTestRunner(verbosity=2).run(mySuite1)
