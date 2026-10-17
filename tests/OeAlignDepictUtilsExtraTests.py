# pylint: disable=possibly-used-before-assignment
##
# File:    OeAlignDepictUtilsExtraTests.py
# Date:    13-Oct-2026
#
##
"""
Additional tests for OeAlignDepictUtils covering setters, alternate display and
search options, grid overflow handling and failure branches.
"""

import io
import os
import os.path
import tempfile
import unittest
from unittest import mock

try:
    from openeye.oechem import OEFloatArray  # noqa: F401 pylint: disable=unused-import

    skiptests = False
except ImportError:
    skiptests = True

if not skiptests:
    from openeye.oechem import OEGraphMol, OESmilesToMol
    from wwpdb.utils.cc_dict_util.timeout.TimeoutMultiProc import TimeoutException

    from wwpdb.utils.oe_util.oedepict.OeAlignDepictUtils import (
        OeDepictMCSAlign,
        OeDepictMCSAlignMulti,
        OeDepictMCSAlignSingle,
        OeTestMCSAlign,
    )

HERE = os.path.abspath(os.path.dirname(__file__))
CACHE_PATH = os.path.join(HERE, "ligand-dict-v3")
DATA_DIR = os.path.join(HERE, "data")
EXAMPLES_DIR = os.path.join(HERE, "examples")
MODULE = "wwpdb.utils.oe_util.oedepict.OeAlignDepictUtils"


@unittest.skipIf(skiptests, "Could not import openeye")
class OeAlignDepictUtilsExtraTests(unittest.TestCase):
    def setUp(self) -> None:
        self.__tmpDir = tempfile.TemporaryDirectory()  # pylint: disable=consider-using-with
        self.__outDir = self.__tmpDir.name
        self.__log = io.StringIO()

    def tearDown(self) -> None:
        self.__tmpDir.cleanup()

    def __out(self, fileName: str) -> str:
        return os.path.join(self.__outDir, fileName)

    def __assertImage(self, path: str) -> None:
        self.assertTrue(os.path.exists(path), path)
        self.assertGreater(os.path.getsize(path), 0)

    def __smilesMol(self, smiles: str) -> "OEGraphMol":
        mol = OEGraphMol()
        self.assertTrue(OESmilesToMol(mol, smiles))
        return mol

    # ------------------------------------------------------------------
    # Setter success paths
    # ------------------------------------------------------------------
    def testSetRefAndFitMol(self) -> None:
        """setRefMol/setFitMol with and without titles, then a test alignment"""
        oed = OeTestMCSAlign(verbose=True, log=self.__log)
        self.assertTrue(oed.setRefMol(self.__smilesMol("c1ccccc1O"), "PHN", title="phenol"))
        self.assertTrue(oed.setFitMol(self.__smilesMol("c1ccccc1N"), "ANL", title="aniline"))
        self.assertTrue(oed.setRefMol(self.__smilesMol("c1ccccc1O"), "PHN"))
        self.assertTrue(oed.setFitMol(self.__smilesMol("c1ccccc1N"), "ANL"))
        aML = oed.doAlign()
        self.assertGreaterEqual(len(aML), 6)
        for rId, _rAt, fId, _fAt in aML:
            self.assertEqual(rId, "PHN")
            self.assertEqual(fId, "ANL")

    def testSetMolFailures(self) -> None:
        """setRefMol/setFitMol return False when the molecule object is invalid"""
        oed = OeTestMCSAlign(verbose=True, log=self.__log)
        self.assertFalse(oed.setRefMol(None, "XXX"))
        self.assertFalse(oed.setFitMol(None, "XXX"))
        self.assertIn("Traceback", self.__log.getvalue())

    def testSetIdAndPathTitleVariants(self) -> None:
        """setRefId/setFitId with titles and setRefPath/setFitPath without titles"""
        oed = OeTestMCSAlign(verbose=True, log=self.__log)
        self.assertTrue(oed.setRefId("A", title="Adenosine", cachePath=CACHE_PATH, suppressHydrogens=True))
        self.assertTrue(oed.setFitId("ATP", title="ATP title", cachePath=CACHE_PATH, suppressHydrogens=True))
        aML = oed.doAlign()
        self.assertGreater(len(aML), 0)

        self.assertTrue(oed.setRefPath("A", os.path.join(DATA_DIR, "A.cif")))
        self.assertTrue(oed.setFitPath("T", os.path.join(DATA_DIR, "T.cif")))
        aML = oed.doAlign()
        self.assertGreater(len(aML), 0)
        self.assertEqual(aML[0][0], "A")

    def testSetIdAndPathFailures(self) -> None:
        """Setters return False when the definition cannot be read"""
        badPath = os.path.join(self.__outDir, "missing.cif")
        oed = OeTestMCSAlign(verbose=True, log=self.__log)
        self.assertFalse(oed.setRefId("ZZZZ", cachePath=self.__outDir))
        self.assertFalse(oed.setFitId("ZZZZ", cachePath=self.__outDir))
        self.assertFalse(oed.setRefPath("ZZZZ", badPath))
        self.assertFalse(oed.setFitPath("ZZZZ", badPath))
        self.assertIn("Traceback", self.__log.getvalue())

    def testFitIdListFailureWithoutRef(self) -> None:
        """setFitIdList fails when no reference definition path has been set"""
        oed = OeTestMCSAlign(verbose=True, log=self.__log)
        self.assertFalse(oed.setFitIdList(["A"], cachePath=CACHE_PATH))

    def testPairIdListFailure(self) -> None:
        """setPairIdList fails on malformed identifiers"""
        oed = OeTestMCSAlign(verbose=True, log=self.__log)
        self.assertFalse(oed.setPairIdList([(None, "A")], cachePath=CACHE_PATH))  # type: ignore[list-item]  # deliberately malformed id

    def testAddFitPathVariants(self) -> None:
        """addFitPath with no reference, without a title, and with an unreadable path"""
        oed = OeTestMCSAlign(verbose=True, log=self.__log)
        self.assertFalse(oed.addFitPath("A", os.path.join(DATA_DIR, "A.cif")))
        self.assertTrue(oed.setRefPath("A", os.path.join(DATA_DIR, "A.cif")))
        self.assertTrue(oed.addFitPath("t", os.path.join(DATA_DIR, "T.cif")))
        self.assertFalse(oed.addFitPath("bad", os.path.join(self.__outDir, "missing.cif")))
        aML = oed.doAlignList()
        self.assertGreater(len(aML), 0)
        self.assertEqual({fId for _r, _ra, fId, _fa in aML}, {"T"})

    def testSetFitPathList(self) -> None:
        """setFitPathList with and without titles followed by doAlignList"""
        oed = OeTestMCSAlign(verbose=True, log=self.__log)
        fitTupList = [
            ("a", os.path.join(DATA_DIR, "A.cif"), None),
            ("gtp", os.path.join(DATA_DIR, "GTP.cif"), "GTP title"),
        ]
        self.assertFalse(oed.setFitPathList(fitTupList))
        self.assertTrue(oed.setRefPath("ATP", os.path.join(DATA_DIR, "ATP.cif")))
        self.assertTrue(oed.setFitPathList(fitTupList))
        aML = oed.doAlignList()
        fitIds = {fId for _r, _ra, fId, _fa in aML}
        self.assertEqual(fitIds, {"A", "GTP"})
        # Unreadable entry -
        self.assertFalse(oed.setFitPathList([("bad", os.path.join(self.__outDir, "missing.cif"), None)]))

    # ------------------------------------------------------------------
    # OeTestMCSAlign alternate search methods
    # ------------------------------------------------------------------
    def testSearchTypes(self) -> None:
        """doAlign with exact and unrecognized search types"""
        oed = OeTestMCSAlign(verbose=False, log=self.__log)
        for sType in ["exact", "unknown", "relaxed", "default"]:
            self.assertEqual(oed.setSearchType(sType=sType), sType)
            self.assertTrue(oed.setRefId("ATP", cachePath=CACHE_PATH))
            self.assertTrue(oed.setFitId("GTP", cachePath=CACHE_PATH))
            aML = oed.doAlign()
            self.assertGreater(len(aML), 0, sType)

    def testDoTestAlign(self) -> None:
        """doTestAlign on an identical pair maps every atom"""
        oed = OeTestMCSAlign(verbose=False, log=self.__log)
        self.assertTrue(oed.setRefPath("ATP", os.path.join(DATA_DIR, "ATP.cif")))
        self.assertTrue(oed.setFitPath("ATP", os.path.join(DATA_DIR, "ATP.cif")))
        aML = oed.doTestAlign()
        self.assertGreater(len(aML), 0)
        self.assertIn("doTestAlign", self.__log.getvalue())

    def testDoAlignWithAnal(self) -> None:
        """doAlignWithAnal reports unmatched atoms and charge differences"""
        oed = OeTestMCSAlign(verbose=False, log=self.__log, maxMatches=64)
        self.assertTrue(oed.setRefId("GTP", cachePath=CACHE_PATH))
        self.assertTrue(oed.setFitId("GTP", cachePath=CACHE_PATH))
        rL = oed.doAlignWithAnal()
        self.assertEqual(len(rL), 7)
        atomMap, uRefAtomList, uRefNList, chgDifRefD, uFitAtomList, uFitNList, chgDifFitD = rL
        self.assertGreater(len(atomMap), 0)
        self.assertIsInstance(uRefAtomList, list)
        self.assertIsInstance(uRefNList, list)
        self.assertEqual(chgDifRefD, {})
        self.assertIsInstance(uFitAtomList, list)
        self.assertIsInstance(uFitNList, list)
        self.assertEqual(chgDifFitD, {})

        # Partial match (fit has extra phosphates) - reference atoms all matched,
        # fit has unmatched atoms with neighbor lists
        self.assertTrue(oed.setRefId("ATP", cachePath=CACHE_PATH))
        self.assertTrue(oed.setFitId("A", cachePath=CACHE_PATH))
        rL = oed.doAlignWithAnal()
        atomMap, uRefAtomList, uRefNList = rL[0], rL[1], rL[2]
        self.assertIn("doAlignWithAnal", self.__log.getvalue())
        self.assertGreater(len(atomMap), 0)
        self.assertGreater(len(uRefAtomList), 0)
        self.assertGreater(len(uRefNList), 0)

    def testDoAlignList(self) -> None:
        """doAlignList over a pair id list"""
        oed = OeTestMCSAlign(verbose=False, log=self.__log)
        self.assertTrue(oed.setPairIdList([("a", "atp"), ("g", "gtp")], cachePath=CACHE_PATH))
        aML = oed.doAlignList()
        refIds = {rId for rId, _ra, _f, _fa in aML}
        self.assertEqual(refIds, {"a", "g"})

    # ------------------------------------------------------------------
    # Single image grid depiction
    # ------------------------------------------------------------------
    def testHighlightStyles(self) -> None:
        """alignPair with each supported highlight style and color option"""
        styleList = ["ballAndStick", "stick", "stickInverse", "ballAndStickInverse", "none"]
        colorList = ["blue", "green", "pink", "other"]
        oed = OeDepictMCSAlign(verbose=False, log=self.__log)
        self.assertTrue(oed.setRefId("C", cachePath=CACHE_PATH))
        self.assertTrue(oed.setFitId("DG", cachePath=CACHE_PATH))
        for ii, style in enumerate(styleList):
            color = colorList[ii % len(colorList)]
            notColor = colorList[(ii + 1) % len(colorList)]
            oed.setDisplayOptions(
                highlightStyleRef=style,
                highlightStyleFit=style,
                highLightMatchColorRef=color,
                highLightMatchColorFit=color,
                highLightNotMatchColorRef=notColor,
                highLightNotMatchColorFit=notColor,
            )
            fName = self.__out("style-%s.svg" % style)
            aML = oed.alignPair(imagePath=fName)
            self.assertGreater(len(aML), 0)
            self.__assertImage(fName)

    def testAlignPairListGridOverflow(self) -> None:
        """alignPairList stops when the grid has no more free cells"""
        oed = OeDepictMCSAlign(verbose=False, log=self.__log)
        oed.setDisplayOptions(gridRows=1, cellBorders=False)
        self.assertTrue(oed.setPairIdList([("c", "dg"), ("a", "atp")], cachePath=CACHE_PATH))
        fName = self.__out("pair-list-overflow.svg")
        aML = oed.alignPairList(imagePath=fName)
        self.__assertImage(fName)
        # Only the first pair fits in a 1x2 grid
        self.assertEqual({rId for rId, _ra, _f, _fa in aML}, {"c"})

    def testAlignOneWithListGridOverflow(self) -> None:
        """alignOneWithList stops when the fit molecules overflow the grid"""
        oed = OeDepictMCSAlign(verbose=False, log=self.__log)
        oed.setDisplayOptions(gridRows=1, gridCols=2)
        self.assertTrue(oed.setRefId("C", cachePath=CACHE_PATH))
        self.assertTrue(oed.setFitIdList(["dg", "cg1"], cachePath=CACHE_PATH, imageDirPath=self.__outDir))
        fName = self.__out("one-with-list-overflow.svg")
        aML = oed.alignOneWithList(imagePath=fName)
        self.__assertImage(fName)
        self.assertEqual({fId for _r, _ra, fId, _fa in aML}, {"DG"})

    # ------------------------------------------------------------------
    # Multi-page depiction
    # ------------------------------------------------------------------
    def testMultiNewPageForFit(self) -> None:
        """alignOneWithListMulti starts a new page when fit molecules fill the grid"""
        oed = OeDepictMCSAlignMulti(verbose=False, log=self.__log)
        oed.setDisplayOptions(gridRows=1, gridCols=2, cellBorders=False)
        self.assertTrue(oed.setRefId("C", cachePath=CACHE_PATH))
        self.assertTrue(oed.setFitIdList(["dg", "cg1", "c"], cachePath=CACHE_PATH))
        fName = self.__out("multi-new-page.pdf")
        aML = oed.alignOneWithListMulti(imagePath=fName)
        self.__assertImage(fName)
        self.assertEqual({fId for _r, _ra, fId, _fa in aML}, {"DG", "CG1", "C"})

    def testMultiPairListTimeout(self) -> None:
        """alignPairListMulti handles a timeout exception and returns an empty map"""
        oed = OeDepictMCSAlignMulti(verbose=True, log=self.__log)
        self.assertTrue(oed.setPairIdList([("c", "dg")], cachePath=CACHE_PATH))
        fName = self.__out("timeout.pdf")
        with mock.patch(MODULE + ".OEWriteMultiPageImage", side_effect=TimeoutException("timeout")):
            aML = oed.alignPairListMulti(imagePath=fName)
        self.assertEqual(aML, [])
        self.assertIn("alignPairListMulti timeout exception", self.__log.getvalue())

    def testMultiPairListException(self) -> None:
        """alignPairListMulti handles a general exception and returns an empty map"""
        oed = OeDepictMCSAlignMulti(verbose=True, log=self.__log)
        self.assertTrue(oed.setPairIdList([("c", "dg")], cachePath=CACHE_PATH))
        with mock.patch(MODULE + ".OEWriteMultiPageImage", side_effect=RuntimeError("boom")):
            aML = oed.alignPairListMulti(imagePath=self.__out("err.pdf"))
        self.assertEqual(aML, [])
        self.assertIn("alignPairListMulti general exception boom", self.__log.getvalue())

    def testMultiOneWithListTimeout(self) -> None:
        """alignOneWithListMulti handles a timeout exception and returns an empty map"""
        oed = OeDepictMCSAlignMulti(verbose=True, log=self.__log)
        self.assertTrue(oed.setRefId("C", cachePath=CACHE_PATH))
        self.assertTrue(oed.setFitIdList(["dg"], cachePath=CACHE_PATH))
        with mock.patch(MODULE + ".OEWriteMultiPageImage", side_effect=TimeoutException("timeout")):
            aML = oed.alignOneWithListMulti(imagePath=self.__out("timeout.pdf"))
        self.assertEqual(aML, [])
        self.assertIn("alignOneWithListMulti timeout exception", self.__log.getvalue())

    def testMultiOneWithListException(self) -> None:
        """alignOneWithListMulti handles a general exception and returns an empty map"""
        oed = OeDepictMCSAlignMulti(verbose=True, log=self.__log)
        self.assertTrue(oed.setRefId("C", cachePath=CACHE_PATH))
        self.assertTrue(oed.setFitIdList(["dg"], cachePath=CACHE_PATH))
        with mock.patch(MODULE + ".OEWriteMultiPageImage", side_effect=RuntimeError("boom")):
            aML = oed.alignOneWithListMulti(imagePath=self.__out("err.pdf"))
        self.assertEqual(aML, [])
        self.assertIn("alignOneWithListMulti general exception here boom", self.__log.getvalue())

    # ------------------------------------------------------------------
    # Single-file depiction
    # ------------------------------------------------------------------
    def testSingleAlignPair(self) -> None:
        """OeDepictMCSAlignSingle.alignPair writes separate ref and fit images"""
        refImage = self.__out("single-ref.svg")
        fitImage = self.__out("single-fit.svg")
        oed = OeDepictMCSAlignSingle(verbose=False, log=self.__log)
        oed.setDisplayOptions(highlightStyleFit="ballAndStick")
        self.assertTrue(oed.setRefId("A", cachePath=CACHE_PATH, imagePath=refImage))
        self.assertTrue(oed.setFitId("ATP", cachePath=CACHE_PATH, imagePath=fitImage))
        aML = oed.alignPair()
        self.assertGreater(len(aML), 0)
        self.__assertImage(refImage)
        self.__assertImage(fitImage)

    def testSingleAlignPairList(self) -> None:
        """OeDepictMCSAlignSingle.alignPairList writes images to the requested directory"""
        oed = OeDepictMCSAlignSingle(verbose=False, log=self.__log)
        self.assertTrue(oed.setPairIdList([("a", "atp"), ("a", "t")], cachePath=CACHE_PATH, imageDirPath=self.__outDir))
        aML = oed.alignPairList()
        self.assertGreater(len(aML), 0)
        self.__assertImage(self.__out("A.svg"))
        self.__assertImage(self.__out("ATP.svg"))


if __name__ == "__main__":
    unittest.main()
