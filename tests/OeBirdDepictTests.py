# pylint: disable=possibly-used-before-assignment
##
#
# File:    OeBirdDepictTests.py
# Author:  jdw
# Date:    5-May-2013
# Version: 0.001
#
# Updates:
#        2-July-2013 jdw  -
#
##
"""
A collection of tests for the OEDepictAlignUtils class with FAMILY PRD data -

The full BIRD family/PRD repository is not part of the test data, so the
PdbxBirdIndex is replaced by a mock describing a small synthetic family index
built on the PRD chemical component definitions in tests/prdcc-v3.
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
import time
import unittest
from typing import Dict, List, Optional, Tuple
from unittest import mock

try:
    from openeye.oechem import OEFloatArray  # noqa: F401 pylint: disable=unused-import

    skiptests = False
except ImportError:
    skiptests = True

if not skiptests:
    from mmcif_utils.bird.PdbxBirdIndex import PdbxBirdIndex

    from wwpdb.utils.oe_util.build.OeChemCompIoUtils import OeChemCompIoUtils
    from wwpdb.utils.oe_util.oedepict.OeAlignDepictUtils import OeDepictMCSAlignMulti
    from wwpdb.utils.oe_util.oedepict.OeDepict import OeDepictMultiPage

# (familyId, prdId, ccId, ccPath)
FamilyMember = Tuple[str, str, str, str]
FamilyDict = Dict[str, List[FamilyMember]]

# Synthetic family index - familyId -> PRD id list
FAMILY_PRD_D: Dict[str, List[str]] = {
    "FAM_000001": ["PRD_000001", "PRD_000010", "PRD_000011"],  # actinomycins
    "FAM_000002": ["PRD_000180"],  # antipain (single member)
    "FAM_000003": ["PRD_999999"],  # no chemical component definition
}


@unittest.skipIf(skiptests, "Could not import openeye")
class OeBirdDepictTests(unittest.TestCase):
    def setUp(self) -> None:
        self.__lfh = sys.stderr
        self.__verbose = True
        self.__here = os.path.abspath(os.path.dirname(__file__))
        self.__prdCcPath = os.path.join(self.__here, "prdcc-v3")
        self.__tmpDir = tempfile.TemporaryDirectory()  # pylint: disable=consider-using-with
        self.__testoutput = self.__tmpDir.name

    def tearDown(self) -> None:
        self.__tmpDir.cleanup()

    def __chemCompId(self, prdId: str) -> Optional[str]:
        ccId = prdId.replace("PRD_", "PRDCC_")
        return ccId if self.__chemCompPath(prdId) else None

    def __chemCompPath(self, prdId: str) -> Optional[str]:
        ccId = prdId.replace("PRD_", "PRDCC_")
        pth = os.path.join(self.__prdCcPath, ccId[-1], ccId + ".cif")
        return pth if os.access(pth, os.R_OK) else None

    def __makeMockIndex(self) -> mock.MagicMock:
        """A PdbxBirdIndex instance mock serving the synthetic family index"""
        bI: mock.MagicMock = mock.create_autospec(PdbxBirdIndex, instance=True)
        bI.getFamilyList.return_value = sorted(FAMILY_PRD_D.keys())
        bI.getPrdIdList.side_effect = lambda familyId: FAMILY_PRD_D.get(familyId, [])
        bI.getPrdPathList.return_value = []
        bI.getChemCompId.side_effect = self.__chemCompId
        bI.getChemCompPath.side_effect = self.__chemCompPath
        return bI

    def __buildBirdIndex(self) -> FamilyDict:
        """Build index of family identifier correspondences -

        Returns -  a dictionary by family_id with valid id correpsondences.
        """
        self.__lfh.write("\nStarting OeBirdDepictTests __buildBirdIndex\n")
        fD: FamilyDict = {}
        indexPath = os.path.join(self.__testoutput, "bird-index.pic")
        with mock.patch("%s.PdbxBirdIndex" % __name__, return_value=self.__makeMockIndex()) as mockCls:
            bI = PdbxBirdIndex(indexPath=indexPath, verbose=self.__verbose, log=self.__lfh)
        mockCls.assert_called_once_with(indexPath=indexPath, verbose=self.__verbose, log=self.__lfh)
        familyIdL = bI.getFamilyList()
        for familyId in familyIdL:
            prdIdList = bI.getPrdIdList(familyId)
            bI.getPrdPathList(familyId)
            for prdId in prdIdList:
                ccId = bI.getChemCompId(prdId)
                ccPath = bI.getChemCompPath(prdId)
                self.__lfh.write("Family %r prdId %r ccId %r ccPath %r\n" % (familyId, prdId, ccId, ccPath))
                if familyId not in fD:
                    fD[familyId] = []
                if (prdId is not None) and (ccId is not None) and (ccPath is not None):
                    fD[familyId].append((familyId, prdId, ccId, ccPath))

        return fD

    def __depictFamilies(self, fD: FamilyDict) -> Dict[str, str]:
        """Depict each family - returns dictionary of familyId -> image path"""
        imageD: Dict[str, str] = {}
        for fId in sorted(fD.keys()):
            fmList = fD[fId]
            if len(fmList) < 1:
                continue

            (familyId, refPrdId, _refCcId, refCcPath) = fmList[0]
            imageFileName = os.path.join(self.__testoutput, familyId + "-members.pdf")
            if len(fmList) > 1:
                oed = OeDepictMCSAlignMulti(verbose=self.__verbose, log=self.__lfh)
                oed.setDisplayOptions(
                    labelAtomName=False,
                    labelAtomCIPStereo=True,
                    labelAtomIndex=False,
                    labelBondIndex=False,
                    highlightStyleFit="ballAndStickInverse",
                    gridRows=3,
                    gridCols=3,
                    bondDisplayWidth=0.5,
                )

                oed.setRefPath(refId=refPrdId, ccPath=refCcPath, title=refPrdId, suppressHydrogens=True)
                for fm in fmList[1:]:
                    (familyId, prdId, _ccId, ccPath) = fm
                    oed.addFitPath(fitId=prdId, ccPath=ccPath, title=prdId, suppressHydrogens=True)

                aML = oed.alignOneWithListMulti(imagePath=imageFileName)
                self.assertGreater(len(aML), 0, familyId)
                fitIdS = set()
                for rCC, rAt, tCC, tAt in aML:
                    self.__lfh.write("%5s %-5s %5s %-5s\n" % (rCC, rAt, tCC, tAt))
                    self.assertEqual(rCC, refPrdId)
                    fitIdS.add(tCC)
                # Not every member need have a MCS match with the reference
                # (e.g. PRD_000010 has a different chromophore tautomer than PRD_000001)
                self.assertTrue(fitIdS)
                self.assertTrue(fitIdS.issubset({fm[1] for fm in fmList[1:]}), fitIdS)
            else:
                oeU = OeChemCompIoUtils(verbose=self.__verbose, log=self.__lfh)
                oemList = oeU.getFromPathList([refCcPath], use3D=False)
                self.assertEqual(len(oemList), 1)
                oedM = OeDepictMultiPage(verbose=self.__verbose, log=self.__lfh)
                oedM.setMolTitleList([(refPrdId, oemList[0], refPrdId)])
                oedM.setDisplayOptions(
                    labelAtomName=False,
                    labelAtomCIPStereo=True,
                    labelAtomIndex=False,
                    labelBondIndex=False,
                    gridRows=3,
                    gridCols=3,
                    bondDisplayWidth=0.5,
                )
                oedM.prepare()
                oedM.write(imageFileName)
            self.assertTrue(os.path.isfile(imageFileName), imageFileName)
            self.assertGreater(os.path.getsize(imageFileName), 0)
            with open(imageFileName, "rb") as ifh:
                self.assertEqual(ifh.read(4), b"%PDF")
            imageD[familyId] = imageFileName
        return imageD

    def testFamilyDepiction(self) -> None:
        """Test case -  aligned family members --"""
        self.__lfh.write("\nStarting OeBirdDepictTests testFamilyDepiction\n")
        fD = self.__buildBirdIndex()
        self.assertEqual(sorted(fD.keys()), sorted(FAMILY_PRD_D.keys()))
        self.assertEqual([fm[1] for fm in fD["FAM_000001"]], FAMILY_PRD_D["FAM_000001"])
        self.assertEqual(fD["FAM_000003"], [])
        imageD = self.__depictFamilies(fD)
        self.assertEqual(sorted(imageD.keys()), ["FAM_000001", "FAM_000002"])

    def testFamilyDepictionHTMLIndex(self) -> None:
        """Test case -  aligned family members --"""
        self.__lfh.write("\nStarting OeBirdDepictTests testFamilyDepictionHTMLIndex\n")
        fD = self.__buildBirdIndex()
        # Produce the family images this index refers to
        self.__depictFamilies(fD)
        tS = time.strftime("%Y %m %d %H:%M:%S", time.localtime())
        htmlPath = os.path.join(self.__testoutput, "index.html")
        nLinks = 0
        with open(htmlPath, "w", encoding="utf-8") as ofh:
            ofh.write("<html>\n")
            ofh.write("<body>\n")
            ofh.write("<h4>Index of family chemical diagrams produced on: %s</h4>\n" % tS)
            ofh.write("<ul>\n")
            for fId in sorted(fD.keys()):
                fmList = fD[fId]
                if len(fmList) < 1:
                    continue
                (familyId, _refPrdId, _refCcId, _refCcPath) = fmList[0]
                imageFileName = familyId + "-members.pdf"
                if os.access(os.path.join(self.__testoutput, imageFileName), os.R_OK):
                    nLinks += 1
                    ofh.write(
                        '<li> <a href="%s">%s</a> with %2d members.</li>\n' % (imageFileName, familyId, len(fmList))
                    )
            ofh.write("</ul>\n")
            ofh.write("</body>\n")
            ofh.write("</html>\n")

        self.assertEqual(nLinks, 2)
        with open(htmlPath, encoding="utf-8") as ifh:
            html = ifh.read()
        self.assertIn('<a href="FAM_000001-members.pdf">FAM_000001</a> with  3 members.', html)
        self.assertIn('<a href="FAM_000002-members.pdf">FAM_000002</a> with  1 members.', html)
        self.assertNotIn("FAM_000003", html)


def suiteDepictFamily() -> unittest.TestSuite:
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeBirdDepictTests("testFamilyDepiction"))
    return suiteSelect


def suiteHTMLIndexFamily() -> unittest.TestSuite:
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeBirdDepictTests("testFamilyDepictionHTMLIndex"))
    return suiteSelect


if __name__ == "__main__":
    mySuite1 = suiteDepictFamily()
    unittest.TextTestRunner(verbosity=2).run(mySuite1)
    mySuite1 = suiteHTMLIndexFamily()
    unittest.TextTestRunner(verbosity=2).run(mySuite1)
