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


import builtins
import os
import sys
import tempfile
import unittest
from typing import Any, List, Sequence, Tuple

# if sys.version_info[0] < 3:
#    from io import open as open

try:
    from openeye.oechem import OEFloatArray  # noqa: F401 pylint: disable=unused-import

    skiptests = False
except ImportError:
    skiptests = True

if not skiptests:
    from wwpdb.utils.oe_util.oedepict.OeAlignDepict import OeDepictMCSAlign


@unittest.skipIf(skiptests, "Could not import openeye")
class OeAlignDepictTests(unittest.TestCase):
    def setUp(self) -> None:
        self.__lfh = sys.stderr
        self.__verbose = True
        self.__here = os.path.abspath(os.path.dirname(__file__))
        self.__examples = os.path.join(self.__here, "examples")
        self.__tmpDir = tempfile.TemporaryDirectory()
        self.__testoutput = self.__tmpDir.name
        self.__datadir = os.path.join(self.__here, "data")
        # Chemical component repository path -
        self.__topCachePath = os.path.join(self.__here, "ligand-dict-v3")
        self.__rnaPairFile = os.path.join(self.__examples, "rna-linking-components.txt")
        self.__refId = "C"
        self.__idList: List[str] = ["cg1", "atp", "gtp", "A", "C", "G", "DG"]
        self.__pairIdList: List[Tuple[str, str]] = [
            ("c", "cg1"),
            ("c", "atp"),
            ("c", "gtp"),
            ("c", "A"),
            ("c", "C"),
            ("c", "G"),
            ("c", "DG"),
        ]

    def tearDown(self) -> None:
        self.__tmpDir.cleanup()

    def __readPairList(self, fn: str = "./examples/rna-linking-components.txt") -> List[Tuple[str, str]]:
        pairList: List[Tuple[str, str]] = []
        with builtins.open(fn, encoding="utf-8") as ifh:
            for line in ifh:
                fields = line.split()
                pairList.append((fields[1], fields[0]))
        return pairList

    def __assertImage(self, imagePath: str) -> None:
        """Check that a non-empty image file was written."""
        self.assertTrue(os.path.isfile(imagePath), "missing image %s" % imagePath)
        self.assertGreater(os.path.getsize(imagePath), 0)

    def __checkAtomMap(self, aML: Sequence[Tuple[Any, ...]], refIds: Sequence[str], fitIds: Sequence[str]) -> None:
        """Check (refId, refAtom, fitId, fitAtom) alignment tuples against the expected identifiers."""
        self.assertTrue(len(aML) > 0, "empty atom map")
        refU = {rId.upper() for rId in refIds}
        fitU = {fId.upper() for fId in fitIds}
        for rCC, rAt, tCC, tAt in aML:
            self.__lfh.write("%5s %-5s %5s %-5s\n" % (rCC, rAt, tCC, tAt))
            self.assertIn(str(rCC).upper(), refU)
            self.assertIn(str(tCC).upper(), fitU)
            self.assertTrue(len(rAt) > 0)
            self.assertTrue(len(tAt) > 0)

    def testMCSAlignPairDepict(self) -> None:
        """Test case -  Simple pairwise MCSS alignment  -  Each aligned pair output to a separate image file"""
        self.__lfh.write("\nStarting OeAlignDepictTests testMCSAlignPairDepict\n")
        oed = OeDepictMCSAlign(verbose=self.__verbose, log=self.__lfh)  # pylint: disable=possibly-used-before-assignment

        oed.setRefId(self.__refId, cachePath=self.__topCachePath)
        for fitId in self.__idList:
            oed.setFitId(fitId, cachePath=self.__topCachePath)
            fName = os.path.join(self.__testoutput, "ref-" + self.__refId + "-trg-" + fitId + ".png")
            aML = oed.alignPair(imagePath=fName)
            self.__checkAtomMap(aML, [self.__refId], [fitId])
            self.__assertImage(fName)

    def testMCSAlignListDepict(self) -> None:
        """Test case -  List view of pairwise MCS alignment - multipage output"""
        self.__lfh.write("\nStarting OeAlignDepictTests testMCSAlignListDepict\n")
        oed = OeDepictMCSAlign(verbose=self.__verbose, log=self.__lfh)
        oed.setRefId(self.__refId, cachePath=self.__topCachePath)
        oed.setFitIdList(self.__idList, cachePath=self.__topCachePath)
        imageFile = os.path.join(self.__testoutput, "list-example-mcs-alignment.pdf")
        aML = oed.alignPairList(imagePath=imageFile)
        self.__checkAtomMap(aML, [self.__refId], self.__idList)
        # every fit molecule shares a substructure with the reference
        self.assertEqual({str(tup[2]).upper() for tup in aML}, {fId.upper() for fId in self.__idList})
        self.__assertImage(imageFile)

    def testMCSAlignPairListDepict(self) -> None:
        """Test case -  List view of  pairwise MCS alignment using pair id list input"""
        self.__lfh.write("\nStarting OeAlignDepictTests testMCSAlignPairListDepict\n")
        oed = OeDepictMCSAlign(verbose=self.__verbose, log=self.__lfh)
        oed.setPairIdList(self.__pairIdList, cachePath=self.__topCachePath)
        imageFile = os.path.join(self.__testoutput, "pair-list-example-mcs-alignment.pdf")
        aML = oed.alignPairList(imagePath=imageFile)
        self.__checkAtomMap(aML, [tup[0] for tup in self.__pairIdList], [tup[1] for tup in self.__pairIdList])
        self.assertEqual({str(tup[2]).upper() for tup in aML}, {tup[1].upper() for tup in self.__pairIdList})
        self.__assertImage(imageFile)

    def testMCSAlignRnaPairListDepict(self) -> None:
        """Test case -  Modified RNA nucleotide alignment with parent nucleotied using pair list input"""
        self.__lfh.write("\nStarting OeAlignDepictTests testMCSAlignRnaPairListDepict\n")
        pairIdList = self.__readPairList(fn=self.__rnaPairFile)
        self.assertTrue(len(pairIdList) > 0)
        oed = OeDepictMCSAlign(verbose=self.__verbose, log=self.__lfh)
        oed.setPairIdList(pairIdList, cachePath=self.__topCachePath)
        imageFile = os.path.join(self.__testoutput, "rna-modified-pair-alignment.pdf")
        aML = oed.alignPairList(imagePath=imageFile)
        self.__checkAtomMap(aML, [tup[0] for tup in pairIdList], [tup[1] for tup in pairIdList])
        self.__assertImage(imageFile)

    def testMCSAlignAtomMap(self) -> None:
        """Test case -  match test with return of atom maps"""
        self.__lfh.write("\nStarting OeAlignDepictTests testMCSAlignAtomMap\n")
        pairIdList = self.__readPairList(fn=self.__rnaPairFile)
        self.assertTrue(len(pairIdList) > 0)
        oed = OeDepictMCSAlign(verbose=self.__verbose, log=self.__lfh)
        matchedPairs: List[Tuple[str, str]] = []
        for refId, fitId in pairIdList:
            oed.setRefId(refId, cachePath=self.__topCachePath)
            oed.setFitId(fitId, cachePath=self.__topCachePath)
            aML = oed.testAlign()
            if len(aML) > 0:
                self.__lfh.write("Match suceeded for: %s %s\n" % (refId, fitId))
                matchedPairs.append((refId, fitId))
                # (refId, refIdx, refType, refName, fitId, fitIdx, fitName)
                for aM in aML:
                    self.assertEqual(len(aM), 7)
                    self.assertEqual(aM[0], refId)
                    self.assertEqual(str(aM[4]).upper(), fitId.upper())
            else:
                self.__lfh.write("Match failed for: %s %s\n" % (refId, fitId))
        self.assertTrue(len(matchedPairs) > 0)
        # Parent nucleotides are complete substructures of these simple modified forms
        for expected in [("C", "5MC"), ("G", "GTP"), ("A", "MA6")]:
            self.assertIn(expected, matchedPairs)

    def testSdfMCSAlignAtomMap(self) -> None:
        """Test case -  match test with return of atom maps for foreign SDF with a CC definition"""
        self.__lfh.write("\nStarting OeAlignDepictTests testSdfMCSAlignAtomMap\n")
        # extPairTup=('../data/ATP.sdf','ATP')
        extPairTup = (os.path.join(self.__datadir, "ATP.sdf"), os.path.join(self.__datadir, "ATP.cif"))
        refPath = extPairTup[0]
        fitId = extPairTup[1]
        fitPath = extPairTup[1]
        oed = OeDepictMCSAlign(verbose=self.__verbose, log=self.__lfh)
        # oed.setSearchType(sType='exact')
        self.assertEqual(oed.setSearchType(sType="exact"), "exact")
        oed.setRefPath(refPath, type="SDF")
        # oed.setFitId(fitId,cachePath=self.__topCachePath)
        oed.setFitPath(fitPath, title="ATP")
        aML = oed.testAlign(suppressHydrogens=False, unique=False, minFrac=0.9)
        self.assertTrue(len(aML) > 0, "Match failed for: %s %s" % (refPath, fitId))
        self.__lfh.write("Match suceeded for: %s %s\n" % (refPath, fitId))
        for ii, aM in enumerate(aML):
            self.__lfh.write("Mapping: %4d %s\n" % (ii, aM))
            self.assertEqual(aM[0], "ATP")
            self.assertEqual(aM[4], "ATP")
        # All heavy atoms of ATP (31 atoms) are mapped by the exact search
        self.assertGreaterEqual(len(aML), 31)

    def testMCSRelaxAlignPairDepict(self) -> None:
        """Test case -  Relaxed pairwise MCSS alignment  -"""
        self.__lfh.write("\nStarting OeAlignDepictTests testMCSRelaxAlignPairDepict\n")
        oed = OeDepictMCSAlign(verbose=self.__verbose, log=self.__lfh)
        self.assertEqual(oed.setSearchType(sType="relaxed"), "relaxed")
        oed.setRefPath(
            ccPath=os.path.join(self.__examples, "PRDCC_000225.cif"), title="PRD_000225", suppressHydrogens=False
        )
        oed.setFitPath(
            ccPath=os.path.join(self.__examples, "L_LDI_990_.comp.cif"), title="L_LDI_990", suppressHydrogens=False
        )
        fName = os.path.join(self.__testoutput, "relaxed-fit.png")
        aML = oed.alignPair(imagePath=fName)
        self.assertTrue(len(aML) > 0)
        for rCC, rAt, tCC, tAt in aML:
            self.__lfh.write("%5s %-5s %5s %-5s\n" % (rCC, rAt, tCC, tAt))
        self.__assertImage(fName)


def suiteAlignPairRelax() -> unittest.TestSuite:
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeAlignDepictTests("testMCSRelaxAlignPairDepict"))
    return suiteSelect


def suiteAlignExtPairRelax() -> unittest.TestSuite:
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeAlignDepictTests("testSdfMCSAlignAtomMap"))
    return suiteSelect


def suiteAlignPair() -> unittest.TestSuite:
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeAlignDepictTests("testMCSAlignPairDepict"))
    suiteSelect.addTest(OeAlignDepictTests("testMCSAlignListDepict"))
    suiteSelect.addTest(OeAlignDepictTests("testMCSAlignPairListDepict"))
    suiteSelect.addTest(OeAlignDepictTests("testMCSAlignRnaPairListDepict"))
    suiteSelect.addTest(OeAlignDepictTests("testMCSAlignAtomMap"))
    return suiteSelect


def suite() -> unittest.TestSuite:
    return unittest.TestLoader().loadTestsFromTestCase(OeAlignDepictTests)


if __name__ == "__main__":
    # unittest.main()
    #
    mySuite1 = suiteAlignPair()
    unittest.TextTestRunner(verbosity=2).run(mySuite1)
    mySuite1 = suiteAlignPairRelax()
    unittest.TextTestRunner(verbosity=2).run(mySuite1)
    mySuite1 = suiteAlignExtPairRelax()
    unittest.TextTestRunner(verbosity=2).run(mySuite1)
