# pylint: disable=possibly-used-before-assignment
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
    from wwpdb.utils.oe_util.oedepict.OeAlignDepictUtils import (
        OeDepictMCSAlign,
        OeDepictMCSAlignMulti,
        OeDepictMCSAlignSingle,
        OeTestMCSAlign,
    )


@unittest.skipIf(skiptests, "Could not import openeye")
class OeAlignDepictUtilsTests(unittest.TestCase):
    def setUp(self) -> None:
        self.__lfh = sys.stderr
        self.__verbose = True
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
        self.__refPathTup = ("A", os.path.join(self.__datadir, "A.cif"), os.path.join(self.__testoutput, "A-ref.svg"))
        self.__fitPathTupList = [
            ("A", os.path.join(self.__datadir, "A.cif"), os.path.join(self.__testoutput, "A-fit.svg")),
            ("T", os.path.join(self.__datadir, "T.cif"), os.path.join(self.__testoutput, "T-fit.svg")),
            ("ATP", os.path.join(self.__datadir, "ATP.cif"), os.path.join(self.__testoutput, "ATP-fit.svg")),
            ("GTP", os.path.join(self.__datadir, "GTP.cif"), os.path.join(self.__testoutput, "GTP-fit.svg")),
        ]

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

    def __setPairDisplayOptions(self, oed: Any) -> None:
        oed.setDisplayOptions(
            labelAtomName=True,
            labelAtomCIPStereo=True,
            labelAtomIndex=False,
            labelBondIndex=False,
            highlightStyleFit="ballAndStickInverse",
            bondDisplayWidth=0.5,
        )

    def testMCSAlignPairDepict(self) -> None:
        """Test case -  Simple pairwise MCSS alignment  -  Each aligned pair output to a separate image file"""
        self.__lfh.write("\nStarting OeAlignDepictUtilsTests testMCSAlignPairDepict\n")
        for fitId in self.__idList:
            # A new aligner is used for each pair - alignPair() only depicts the first pair
            # it is given when an instance is reused (see testMCSAlignPairDepictReuse).
            oed = OeDepictMCSAlign(verbose=self.__verbose, log=self.__lfh)
            self.__setPairDisplayOptions(oed)
            self.assertTrue(oed.setRefId(self.__refId, cachePath=self.__topCachePath))
            self.assertTrue(oed.setFitId(fitId, cachePath=self.__topCachePath))
            fName = os.path.join(self.__testoutput, "ref-" + self.__refId + "-trg-" + fitId + ".svg")
            aML = oed.alignPair(imagePath=fName)
            self.__writeMap(aML)
            self.assertEqual(len(aML), self.__expectedPairMatchCount[fitId], fitId)
            self.assertEqual({tup[0] for tup in aML}, {self.__refId})
            self.assertEqual({tup[2] for tup in aML}, {fitId.upper()})
            self.__assertImage(fName)

    # alignPair() appends the current ref/fit pair to the internal pair list only when that
    # list is empty (OeAlignDepictUtils.py OeDepictMCSAlign.alignPair), so a reused instance
    # keeps re-depicting the first fit molecule after setFitId() is called again.
    @unittest.expectedFailure
    def testMCSAlignPairDepictReuse(self) -> None:
        """Test case -  Reuse one aligner for successive pairwise alignments"""
        oed = OeDepictMCSAlign(verbose=self.__verbose, log=self.__lfh)
        self.__setPairDisplayOptions(oed)
        self.assertTrue(oed.setRefId(self.__refId, cachePath=self.__topCachePath))
        for fitId in ["atp", "DG"]:
            self.assertTrue(oed.setFitId(fitId, cachePath=self.__topCachePath))
            aML = oed.alignPair(imagePath=os.path.join(self.__testoutput, "reuse-" + fitId + ".svg"))
            self.assertEqual({tup[2] for tup in aML}, {fitId.upper()})

    def testMCSAlignFitListDepictMulti(self) -> None:
        """Test case -  List view of pairwise MCS alignment - multipage output"""
        self.__lfh.write("\nStarting OeAlignDepictUtilsTests testMCSAlignFitListDepictMulti\n")
        oed = OeDepictMCSAlignMulti(verbose=self.__verbose, log=self.__lfh)
        self.assertTrue(oed.setRefId(self.__refId, cachePath=self.__topCachePath))
        self.assertTrue(oed.setFitIdList(self.__idList, cachePath=self.__topCachePath))
        oed.setDisplayOptions(
            labelAtomName=True,
            labelAtomCIPStereo=True,
            labelAtomIndex=False,
            labelBondIndex=False,
            highlightStyleFit="ballAndStickInverse",
            gridRows=3,
            gridCols=3,
            bondDisplayWidth=0.5,
        )
        imageFile = os.path.join(self.__testoutput, "list-example-mcs-alignment.pdf")
        aML = oed.alignOneWithListMulti(imagePath=imageFile)
        self.__writeMap(aML)
        self.assertEqual(len(aML), sum(self.__expectedPairMatchCount.values()))
        self.assertEqual({tup[2] for tup in aML}, {fitId.upper() for fitId in self.__idList})
        self.__assertImage(imageFile)

    def testMCSAlignPairListDepictPortrait(self) -> None:
        """Test case -  List view MCS alignment using pair id list input"""
        self.__lfh.write("\nStarting OeAlignDepictUtilsTests testMCSAlignPairListDepictPortrait\n")
        oed = OeDepictMCSAlignMulti(verbose=self.__verbose, log=self.__lfh)
        self.assertTrue(oed.setPairIdList(self.__pairIdList, cachePath=self.__topCachePath))
        oed.setDisplayOptions(
            labelAtomName=True,
            labelAtomCIPStereo=True,
            labelAtomIndex=False,
            labelBondIndex=False,
            highlightStyleFit="ballAndStickInverse",
            pageOrientation="portrait",
            gridRows=4,
            bondDisplayWidth=0.5,
        )

        imageFile = os.path.join(self.__testoutput, "pair-list-example-mcs-alignment-portrait.pdf")
        aML = oed.alignPairListMulti(imagePath=imageFile)
        self.__writeMap(aML)
        self.assertEqual(len(aML), sum(self.__expectedPairMatchCount.values()))
        self.assertEqual({tup[0] for tup in aML}, {"c"})
        self.assertEqual({tup[2] for tup in aML}, set(self.__idList))
        self.__assertImage(imageFile)

    def testMCSAlignRnaPairListDepict(self) -> None:
        """Test case -  Modified RNA nucleotide alignment with parent nucleotied using pair list input"""
        self.__lfh.write("\nStarting OeAlignDepictUtilsTests testCSAlignRnaPairListDepict\n")
        pairIdList = self.__readPairList(fn=self.__rnaPairFile)
        self.assertEqual(len(pairIdList), 128)
        oed = OeDepictMCSAlignMulti(verbose=self.__verbose, log=self.__lfh)
        self.assertTrue(oed.setPairIdList(pairIdList, cachePath=self.__topCachePath))
        self.__setPairDisplayOptions(oed)

        imageFile = os.path.join(self.__testoutput, "rna-modified-pair-alignment.pdf")
        aML = oed.alignPairListMulti(imagePath=imageFile)
        self.__writeMap(aML)
        # Every modified nucleotide aligns with its parent
        self.assertEqual({tup[2] for tup in aML}, {fitId for _refId, fitId in pairIdList})
        self.assertGreater(len(aML), len(pairIdList))
        self.__assertImage(imageFile)

    def testMCSAlignAtomMap(self) -> None:
        """Test case -  match test with return of atom maps"""
        self.__lfh.write("\nStarting OeAlignDepictUtilsTests testMCSAlignAtomMap\n")
        pairIdList = self.__readPairList(fn=self.__rnaPairFile)
        oed = OeTestMCSAlign(verbose=self.__verbose, log=self.__lfh)
        failedL: List[Tuple[str, str]] = []
        for refId, fitId in pairIdList:
            self.assertTrue(oed.setRefId(refId, cachePath=self.__topCachePath))
            self.assertTrue(oed.setFitId(fitId, cachePath=self.__topCachePath))
            aML = oed.doAlign()
            if len(aML) > 0:
                self.__lfh.write("Match length %3d for: %s %s\n" % (len(aML), refId, fitId))
                self.assertEqual({tup[0] for tup in aML}, {refId})
                self.assertEqual({tup[2] for tup in aML}, {fitId})
            else:
                self.__lfh.write("Match failed for: %s %s\n" % (refId, fitId))
                failedL.append((refId, fitId))
        self.assertEqual(failedL, [])

    def testMCSRelaxAlignPairDepict(self) -> None:
        """Test case -  Relaxed pairwise MCSS alignment  -"""
        self.__lfh.write("\nStarting OeAlignDepictUtilsTests testMCSRelaxAlignPairDepict\n")
        oed = OeDepictMCSAlign(verbose=self.__verbose, log=self.__lfh)
        self.assertEqual(oed.setSearchType(sType="relaxed"), "relaxed")
        self.assertTrue(
            oed.setRefPath(
                refId="PRD_000225",
                ccPath=os.path.join(self.__examples, "PRDCC_000225.cif"),
                title="PRD_000225",
                suppressHydrogens=False,
            )
        )
        self.assertTrue(
            oed.setFitPath(
                fitId="LD_LDI_990",
                ccPath=os.path.join(self.__examples, "L_LDI_990_.comp.cif"),
                title="L_LDI_990",
                suppressHydrogens=False,
            )
        )
        fName = os.path.join(self.__testoutput, "relaxed-fit.svg")
        oed.setDisplayOptions(
            imageSizeX=2000,
            imageSizeY=1000,
            labelAtomName=True,
            labelAtomCIPStereo=True,
            labelAtomIndex=False,
            labelBondIndex=False,
            highlightStyleFit="ballAndStickInverse",
            bondDisplayWidth=1.0,
        )
        aML = oed.alignPair(imagePath=fName)
        self.__writeMap(aML)
        self.assertEqual(len(aML), 72)
        self.assertEqual({tup[0] for tup in aML}, {"PRD_000225"})
        self.__assertImage(fName)

    def testMCSAlignListMultiDepict(self) -> None:
        """Test case -  List view of MCS alignment --- on multi-pages --"""
        self.__lfh.write("\nStarting OeAlignDepictUtilsTests testMCSAlignListMultiDepict\n")
        oed = OeDepictMCSAlignMulti(verbose=self.__verbose, log=self.__lfh)
        oed.setDisplayOptions(
            labelAtomName=True,
            labelAtomCIPStereo=True,
            labelAtomIndex=False,
            labelBondIndex=False,
            highlightStyleFit="ballAndStickInverse",
            gridRows=3,
            gridCols=3,
            bondDisplayWidth=0.5,
        )
        self.assertTrue(oed.setPairIdList(self.__pairIdList, cachePath=self.__topCachePath))
        imageFile = os.path.join(self.__testoutput, "mcs-align-with-list-multi.pdf")
        aML = oed.alignOneWithListMulti(imagePath=imageFile)
        self.__writeMap(aML)
        self.assertEqual(len(aML), sum(self.__expectedPairMatchCount.values()))
        self.assertEqual({tup[2] for tup in aML}, set(self.__idList))
        self.__assertImage(imageFile)

    def testMCSAlignListDepict(self) -> None:
        """Test case -  List view of MCS alignment on single image -"""
        self.__lfh.write("\nStarting OeAlignDepictUtilsTests testMCSAlignListDepict\n")
        oed = OeDepictMCSAlign(verbose=self.__verbose, log=self.__lfh)
        oed.setDisplayOptions(
            imageSizeX=1500,
            imageSizeY=1500,
            labelAtomName=True,
            labelAtomCIPStereo=True,
            labelAtomIndex=False,
            labelBondIndex=False,
            highlightStyleFit="ballAndStickInverse",
            gridRows=3,
            gridCols=3,
            bondDisplayWidth=1.0,
        )
        self.assertTrue(oed.setPairIdList(self.__pairIdList, cachePath=self.__topCachePath))
        imageFile = os.path.join(self.__testoutput, "mcs-align-with-list-single.svg")
        aML = oed.alignOneWithList(imagePath=imageFile)
        self.__writeMap(aML)
        self.assertEqual(len(aML), sum(self.__expectedPairMatchCount.values()))
        self.assertEqual({tup[2] for tup in aML}, set(self.__idList))
        self.__assertImage(imageFile)

    def testMCSAlignListDepictSingle(self) -> None:
        """Test case -  View MCS alignment in single image files.

        Input details specified as chemical component Id pair lists -
        """
        self.__lfh.write("\nStarting OeAlignDepictUtilsTests testMCSAlignListDep\n")
        oed = OeDepictMCSAlignSingle(verbose=self.__verbose, log=self.__lfh)
        oed.setDisplayOptions(
            imageSizeX=500,
            imageSizeY=500,
            labelAtomName=True,
            labelAtomCIPStereo=True,
            labelAtomIndex=False,
            labelBondIndex=False,
            highlightStyleFit="ballAndStickInverse",
            bondDisplayWidth=1.0,
        )
        imageDirPath = os.path.join(self.__testoutput, "single")
        os.makedirs(imageDirPath)
        self.assertTrue(oed.setPairIdList(self.__pairIdList, cachePath=self.__topCachePath, imageDirPath=imageDirPath))
        aML = oed.alignOneWithList()
        self.__writeMap(aML)
        self.assertEqual(len(aML), sum(self.__expectedPairMatchCount.values()))
        self.assertEqual({tup[2] for tup in aML}, set(self.__idList))
        # One reference image plus one image per fit molecule (reference C is also a fit)
        expectedFiles = sorted(fitId.upper() + ".svg" for fitId in self.__idList)
        self.assertEqual(sorted(os.listdir(imageDirPath)), expectedFiles)
        for fileName in expectedFiles:
            self.__assertImage(os.path.join(imageDirPath, fileName))

    def testMCSAlignPathListDepictSingle(self) -> None:
        """Test case -  View MCS alignments in single image files.

        Input details specified chemical component and image file paths.

        Image file paths must end with a recognized image format (e.g. svg, png, jpg)
        """
        self.__lfh.write("\nStarting OeAlignDepictUtilsTests testMCSAlignPathListDepictSingle\n")
        # Use inverse  highlighting of matching/non-matching atoms/bonds -
        oed = OeDepictMCSAlignSingle(verbose=self.__verbose, log=self.__lfh)
        oed.setDisplayOptions(
            imageSizeX=500,
            imageSizeY=500,
            labelAtomName=True,
            labelAtomCIPStereo=True,
            labelAtomIndex=False,
            labelBondIndex=False,
            highlightStyleFit="ballAndStickInverse",
            bondDisplayWidth=1.0,
        )

        self.assertTrue(
            oed.setRefPath(
                self.__refPathTup[0], self.__refPathTup[1], title=self.__refPathTup[0], imagePath=self.__refPathTup[2]
            )
        )
        for fitPathTup in self.__fitPathTupList:
            self.assertTrue(oed.addFitPath(fitPathTup[0], fitPathTup[1], title=fitPathTup[0], imagePath=fitPathTup[2]))

        aML = oed.alignOneWithList()
        #
        # Write out no atom correspondences --
        #
        self.__writeMap(aML)
        self.assertGreater(len(aML), 0)
        self.assertEqual({tup[0] for tup in aML}, {self.__refPathTup[0]})
        self.assertEqual({tup[2] for tup in aML}, {fitPathTup[0] for fitPathTup in self.__fitPathTupList})
        self.__assertImage(self.__refPathTup[2])
        for fitPathTup in self.__fitPathTupList:
            self.__assertImage(fitPathTup[2])


def suiteAlignTypeTests() -> unittest.TestSuite:
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeAlignDepictUtilsTests("testMCSAlignAtomMap"))
    suiteSelect.addTest(OeAlignDepictUtilsTests("testMCSRelaxAlignPairDepict"))
    return suiteSelect


def suiteAlignPair() -> unittest.TestSuite:
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeAlignDepictUtilsTests("testMCSAlignPairDepict"))
    suiteSelect.addTest(OeAlignDepictUtilsTests("testMCSAlignFitListDepictMulti"))
    suiteSelect.addTest(OeAlignDepictUtilsTests("testMCSAlignPairListDepictPortrait"))
    suiteSelect.addTest(OeAlignDepictUtilsTests("testMCSAlignRnaPairListDepict"))
    return suiteSelect


def suiteAlignWithList() -> unittest.TestSuite:
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeAlignDepictUtilsTests("testMCSAlignListMultiDepict"))
    suiteSelect.addTest(OeAlignDepictUtilsTests("testMCSAlignListDepict"))
    return suiteSelect


def suiteAlignWithListSingle() -> unittest.TestSuite:
    suiteSelect = unittest.TestSuite()
    # suiteSelect.addTest(OeAlignDepictUtilsTests("testMCSAlignListDepictSingle"))
    suiteSelect.addTest(OeAlignDepictUtilsTests("testMCSAlignPathListDepictSingle"))
    return suiteSelect


if __name__ == "__main__":
    mySuite1 = suiteAlignPair()
    unittest.TextTestRunner(verbosity=2).run(mySuite1)
    mySuite1 = suiteAlignWithList()
    unittest.TextTestRunner(verbosity=2).run(mySuite1)

    mySuite1 = suiteAlignTypeTests()
    unittest.TextTestRunner(verbosity=2).run(mySuite1)
    mySuite1 = suiteAlignWithListSingle()
    unittest.TextTestRunner(verbosity=2).run(mySuite1)

    mySuite1 = suiteAlignWithListSingle()
    unittest.TextTestRunner(verbosity=2).run(mySuite1)
