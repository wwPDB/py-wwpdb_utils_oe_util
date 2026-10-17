##
# File:    OeDepictExtraTests.py
#
# Additional behaviour tests for OeDepict and OeDepictMultiPage
##
"""
Additional tests for 2D depictions with various display options.
"""

import os
import tempfile
import unittest
from typing import Any, List, Tuple

try:
    from openeye.oechem import OEFloatArray  # noqa: F401 pylint: disable=unused-import

    skiptests = False
except ImportError:
    skiptests = True

if not skiptests:
    from wwpdb.utils.oe_util.build.OeChemCompIoUtils import OeChemCompIoUtils
    from wwpdb.utils.oe_util.oedepict.OeDepict import OeDepict, OeDepictMultiPage


@unittest.skipIf(skiptests, "Cannot import openeye.oechem for tests")
class OeDepictExtraTests(unittest.TestCase):
    def setUp(self) -> None:
        here = os.path.abspath(os.path.dirname(__file__))
        self.__topCachePath = os.path.join(here, "ligand-dict-v3")
        self.__tmpDir = tempfile.TemporaryDirectory()  # pylint: disable=consider-using-with
        self.__outDir = self.__tmpDir.name

    def tearDown(self) -> None:
        self.__tmpDir.cleanup()

    def __getMolTitleList(self, idList: List[str]) -> List[Tuple[str, Any, str]]:
        """Fresh molecules for each call (hydrogen suppression modifies the stored molecule)"""
        oeU = OeChemCompIoUtils(topCachePath=self.__topCachePath, verbose=False)  # pylint: disable=possibly-used-before-assignment
        oemList = oeU.getFromIdList(idList, use3D=True, coordType="ideal")
        self.assertEqual(len(oemList), len(idList))
        molTitleList: List[Tuple[str, Any, str]] = []
        for oem in oemList:
            ccId = oem.getCcId()
            assert ccId is not None
            molTitleList.append((ccId, oem, "Title-" + ccId))
        return molTitleList

    def __depictSingle(self, fileName: str, useTitle: bool = True, idList: Any = None, **kwargs: Any) -> str:
        ids: List[str] = ["ATP"] if idList is None else idList
        oed = OeDepict(verbose=False, useTitle=useTitle)  # pylint: disable=possibly-used-before-assignment
        oed.setMolTitleList(self.__getMolTitleList(ids))
        oed.setDisplayOptions(**kwargs)
        oed.setGridOptions(rows=1, cols=len(ids))
        oed.prepare()
        imagePath = os.path.join(self.__outDir, fileName)
        oed.write(imagePath)
        self.assertTrue(os.path.isfile(imagePath))
        self.assertGreater(os.path.getsize(imagePath), 0)
        return imagePath

    def __depictMulti(self, fileName: str, useTitle: bool = True, nCopies: int = 1, **kwargs: Any) -> str:
        oed = OeDepictMultiPage(verbose=False, useTitle=useTitle)  # pylint: disable=possibly-used-before-assignment
        molTitleList = self.__getMolTitleList(["ATP", "GTP"] * nCopies)
        oed.setMolTitleList(molTitleList)
        oed.setDisplayOptions(**kwargs)
        oed.setGridOptions(rows=1, cols=1)
        oed.prepare()
        imagePath = os.path.join(self.__outDir, fileName)
        oed.write(imagePath)
        self.assertTrue(os.path.isfile(imagePath))
        self.assertGreater(os.path.getsize(imagePath), 0)
        return imagePath

    @staticmethod
    def __read(path: str) -> str:
        with open(path, encoding="utf-8") as ifh:
            return ifh.read()

    def testTitleShownAndHidden(self) -> None:
        """useTitle controls whether the title appears in the image"""
        withTitle = self.__read(self.__depictSingle("t.svg", useTitle=True))
        noTitle = self.__read(self.__depictSingle("nt.svg", useTitle=False))
        self.assertIn("Title-ATP", withTitle)
        self.assertNotIn("Title-ATP", noTitle)

    def testLabelAtomIndex(self) -> None:
        """Atom index labels change the depiction"""
        base = self.__read(self.__depictSingle("base.svg"))
        lab = self.__read(self.__depictSingle("aidx.svg", labelAtomIndex=True))
        self.assertNotEqual(base, lab)
        self.assertGreater(lab.count("<text"), base.count("<text"))

    def testLabelBondIndex(self) -> None:
        """Bond index labels change the depiction"""
        base = self.__read(self.__depictSingle("base.svg"))
        lab = self.__read(self.__depictSingle("bidx.svg", labelBondIndex=True))
        self.assertNotEqual(base, lab)
        self.assertGreater(lab.count("<text"), base.count("<text"))

    def testLabelAtomName(self) -> None:
        """Atom name labels include the chemical component atom names"""
        lab = self.__read(self.__depictSingle("aname.svg", labelAtomName=True))
        self.assertIn("PG", lab)
        self.assertIn("O1G", lab)

    @staticmethod
    def __stripDate(svg: str) -> str:
        return "\n".join(line for line in svg.splitlines() if "Creation Date" not in line)

    def testSuppressHydrogens(self) -> None:
        """Suppressing hydrogens renders the same drawing as the default

        OEPrepareDepiction() already suppresses explicit hydrogens by default, so the
        suppressHydrogens option does not change the rendered image.
        """
        base = self.__read(self.__depictSingle("base.svg"))
        noH = self.__read(self.__depictSingle("noh.svg", suppressHydrogens=True))
        self.assertEqual(self.__stripDate(base), self.__stripDate(noH))

    def testAllOptionsPng(self) -> None:
        """All options together render to a PNG with valid signature"""
        path = self.__depictSingle(
            "all.png",
            idList=["ATP", "GTP"],
            labelAtomIndex=True,
            labelBondIndex=True,
            labelAtomCIPStereo=True,
            suppressHydrogens=True,
            cellBorders=False,
            bondDisplayWidth=1.0,
            imageSizeX=400,
            imageSizeY=200,
        )
        with open(path, "rb") as ifh:
            self.assertEqual(ifh.read(8), b"\x89PNG\r\n\x1a\n")

    @unittest.expectedFailure
    def testDefaultGridMoreCellsThanMols(self) -> None:
        """Default 2x2 grid with a single molecule

        Known issue: OeDepict.prepare() iterates over all grid cells and indexes the
        molecule list, raising IndexError when there are fewer molecules than cells.
        """
        oed = OeDepict(verbose=False)
        oed.setMolTitleList(self.__getMolTitleList(["ATP"]))
        oed.prepare()
        imagePath = os.path.join(self.__outDir, "grid.svg")
        oed.write(imagePath)
        self.assertGreater(os.path.getsize(imagePath), 0)

    def testMultiPageOptions(self) -> None:
        """Multi-page depiction with hydrogen suppression, labels and no titles"""
        path = self.__depictMulti(
            "multi.pdf",
            useTitle=False,
            suppressHydrogens=True,
            labelAtomIndex=True,
            labelBondIndex=True,
            pageOrientation="portrait",
        )
        with open(path, "rb") as ifh:
            self.assertEqual(ifh.read(4), b"%PDF")

    def testMultiPageTitles(self) -> None:
        """Multi-page depiction with titles over several pages"""
        path = self.__depictMulti("multiT.pdf", useTitle=True, nCopies=2, labelAtomName=True)
        with open(path, "rb") as ifh:
            self.assertEqual(ifh.read(4), b"%PDF")


def suiteExtra() -> unittest.TestSuite:
    return unittest.TestLoader().loadTestsFromTestCase(OeDepictExtraTests)


if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suiteExtra())
