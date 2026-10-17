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


import os
import os.path
import sys
import tempfile
import unittest
from typing import Any, Dict, List, Sequence, Tuple

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
        self.__verbose = False

        self.__here = os.path.abspath(os.path.dirname(__file__))
        self.__examples = os.path.join(self.__here, "examples")
        self.__datadir = os.path.join(self.__here, "data")
        self.__tmpDir = tempfile.TemporaryDirectory()  # pylint: disable=consider-using-with
        self.__testoutput = self.__tmpDir.name

        # Chemical component repository path -
        self.__topCachePath = os.path.join(self.__here, "ligand-dict-v3")
        self.__rnaPairFile = os.path.join(self.__examples, "rna-linking-components.txt")
        self.__refId = "C"
        self.__idList = ["cg1", "atp", "gtp", "A", "C", "G", "DG"]
        self.__pairIdList = [("c", "cg1"), ("c", "atp"), ("c", "gtp"), ("c", "A"), ("c", "C"), ("c", "G"), ("c", "DG")]
        # Number of matched atoms for each reference C / fit pair (heavy atoms only)
        self.__expectedPairMatchCount: Dict[str, int] = {
            "cg1": 19,
            "atp": 20,
            "gtp": 19,
            "A": 20,
            "C": 21,
            "G": 19,
            "DG": 18,
        }

    def tearDown(self) -> None:
        self.__tmpDir.cleanup()

    def __readPairList(self, fn: str) -> List[Tuple[str, str]]:
        pairList: List[Tuple[str, str]] = []
        with open(fn, encoding="utf-8") as ifh:
            for line in ifh:
                fields = line.split()
                pairList.append((fields[1], fields[0]))
        return pairList

    def __assertImage(self, path: str) -> None:
        self.assertTrue(os.path.isfile(path), path)
        self.assertGreater(os.path.getsize(path), 0, path)

    def __writeMap(self, aML: Sequence[Tuple[Any, ...]]) -> None:
        for rCC, rAt, tCC, tAt in aML:
            self.__lfh.write("%5s %-5s %5s %-5s\n" % (rCC, rAt, tCC, tAt))

    def testMCSAlignPairDepict(self) -> None:
        """Test case -  Simple pairwise MCSS alignment  -  Each aligned pair output to a separate image file"""
        self.__lfh.write("\nStarting OeAlignDepictTests tstMCSalignPairDepict\n")
        oed = OeDepictMCSAlign(verbose=self.__verbose, log=self.__lfh)  # pylint: disable=possibly-used-before-assignment
        oed.setRefId(self.__refId, cachePath=self.__topCachePath)
        for fitId in self.__idList:
            oed.setFitId(fitId, cachePath=self.__topCachePath)
            fName = os.path.join(self.__testoutput, "ref-" + self.__refId + "-trg-" + fitId + ".png")
            aML = oed.alignPair(imagePath=fName)
            self.__writeMap(aML)
            self.assertEqual(len(aML), self.__expectedPairMatchCount[fitId], fitId)
            self.assertEqual({tup[0] for tup in aML}, {self.__refId})
            self.assertEqual({tup[2] for tup in aML}, {fitId.upper()})
            self.__assertImage(fName)
            if fitId == self.__refId:
                # Self alignment maps every atom onto the atom of the same name
                for _rCC, rAt, _tCC, tAt in aML:
                    self.assertEqual(rAt, tAt)

    def testMCSAlignListDepict(self) -> None:
        """Test case -  List view of pairwise MCS alignment - multipage output"""
        self.__lfh.write("\nStarting OeAlignDepictTests testMCSAlignListDepict\n")
        oed = OeDepictMCSAlign(verbose=self.__verbose, log=self.__lfh)
        oed.setRefId(self.__refId, cachePath=self.__topCachePath)
        oed.setFitIdList(self.__idList, cachePath=self.__topCachePath)
        imageFile = os.path.join(self.__testoutput, "list-example-mcs-alignment.pdf")
        aML = oed.alignPairList(imagePath=imageFile)
        self.__writeMap(aML)
        self.assertEqual(len(aML), sum(self.__expectedPairMatchCount.values()))
        self.assertEqual({tup[2] for tup in aML}, set(self.__idList))
        self.__assertImage(imageFile)

    def testMCSAlignPairListDepict(self) -> None:
        """Test case -  List view of  pairwise MCS alignment using pair id list input"""
        self.__lfh.write("\nStarting OeAlignDepictTests testMCSAlignPairListDepict\n")
        oed = OeDepictMCSAlign(verbose=self.__verbose, log=self.__lfh)
        oed.setPairIdList(self.__pairIdList, cachePath=self.__topCachePath)
        imageFile = os.path.join(self.__testoutput, "pair-list-example-mcs-alignment.pdf")
        aML = oed.alignPairList(imagePath=imageFile)
        self.__writeMap(aML)
        self.assertEqual(len(aML), sum(self.__expectedPairMatchCount.values()))
        self.assertEqual({tup[0] for tup in aML}, {"C"})
        self.assertEqual({tup[2] for tup in aML}, {fitId for _refId, fitId in self.__pairIdList})
        self.__assertImage(imageFile)

    def testMCSAlignRnaPairListDepict(self) -> None:
        """Test case -  Modified RNA nucleotide alignment with parent nucleotied using pair list input"""
        self.__lfh.write("\nStarting OeAlignDepictTests testMCSAlignRnaPairListDepict\n")
        pairIdList = self.__readPairList(fn=self.__rnaPairFile)
        self.assertEqual(len(pairIdList), 128)
        oed = OeDepictMCSAlign(verbose=self.__verbose, log=self.__lfh)
        oed.setPairIdList(pairIdList, cachePath=self.__topCachePath)
        imageFile = os.path.join(self.__testoutput, "rna-modified-pair-alignment.pdf")
        aML = oed.alignPairList(imagePath=imageFile)
        self.__writeMap(aML)
        self.assertGreater(len(aML), 0)
        # All modified nucleotides except URD (vs U) satisfy the full minimum match fraction
        matchedFitIds = {tup[2] for tup in aML}
        self.assertEqual(matchedFitIds, {fitId for _refId, fitId in pairIdList} - {"URD"})
        self.__assertImage(imageFile)

    def testMCSAlignAtomMap(self) -> None:
        """Test case -  match test with return of atom maps"""
        self.__lfh.write("\nStarting OeAlignDepictTests testMCSAlignAtomMap\n")
        pairIdList = self.__readPairList(fn=self.__rnaPairFile)
        oed = OeDepictMCSAlign(verbose=self.__verbose, log=self.__lfh)
        nMatched = 0
        for refId, fitId in pairIdList:
            oed.setRefId(refId, cachePath=self.__topCachePath)
            oed.setFitId(fitId, cachePath=self.__topCachePath)
            aML = oed.testAlign()
            if len(aML) > 0:
                nMatched += 1
                self.__lfh.write("Match suceeded for: %s %s\n" % (refId, fitId))
                # (refId, refIdx, refType, refName, fitId, fitIdx, fitName)
                self.assertEqual({tup[0] for tup in aML}, {refId})
                self.assertEqual({tup[4] for tup in aML}, {fitId})
            else:
                self.__lfh.write("Match failed for: %s %s\n" % (refId, fitId))
        # testAlign requires a complete match of the smaller molecule (minFrac=1.0)
        self.assertEqual(nMatched, 71)

    def testSdfMCSAlignAtomMap(self) -> None:
        """Test case -  match test with return of atom maps for foreign SDF with a CC definition"""
        self.__lfh.write("\nStarting OeAlignDepictTests testSdfMCSAlignAtomMap\n")
        refPath = os.path.join(self.__datadir, "ATP.sdf")
        fitPath = os.path.join(self.__datadir, "ATP.cif")
        oed = OeDepictMCSAlign(verbose=self.__verbose, log=self.__lfh)
        self.assertEqual(oed.setSearchType(sType="exact"), "exact")
        oed.setRefPath(refPath, type="SDF")
        oed.setFitPath(fitPath, title="ATP")
        aML = oed.testAlign(suppressHydrogens=False, unique=False, minFrac=0.9)
        for ii, aM in enumerate(aML):
            self.__lfh.write("Mapping: %4d %s\n" % (ii, aM))
        # ATP with all hydrogens - 31 heavy atoms + 16 hydrogens
        self.assertEqual(len(aML), 47)
        self.assertEqual({tup[0] for tup in aML}, {"ATP"})
        self.assertEqual({tup[4] for tup in aML}, {"ATP"})

    def testMCSRelaxAlignPairDepict(self) -> None:
        """Test case -  Relaxed pairwise MCSS alignment  -"""
        self.__lfh.write("\nStarting OeAlignDepictTests testMCSRelaceAlignPairDepict\n")
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
        self.__writeMap(aML)
        self.assertEqual(len(aML), 72)
        self.assertEqual({tup[0] for tup in aML}, {"PRD_000225"})
        self.assertEqual({tup[2] for tup in aML}, {"LDI"})
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
