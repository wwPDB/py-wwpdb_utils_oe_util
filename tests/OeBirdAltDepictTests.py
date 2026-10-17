# pylint: disable=possibly-used-before-assignment
##
#
# File:    OeBirdAltDepictTests.py
# Author:  jdw
# Date:    4-July-2013
# Updates:
#
##
"""
A collection of tests for the OEDepictAlignUtils class with FAMILY PRD data -

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
from typing import Dict, List, Tuple

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


@unittest.skipIf(skiptests, "Could not import openeye")
@unittest.skip("Input files not available - family-members-alternative")
class OeBirdDepictTests(unittest.TestCase):
    def setUp(self) -> None:
        self.__lfh = sys.stderr
        self.__verbose = True
        self.__here = os.path.abspath(os.path.dirname(__file__))
        self.__altFamilyIndexPath = os.path.join(self.__here, "examples", "family-members-alternate.dat")
        self.__tmpDir = tempfile.TemporaryDirectory()  # pylint: disable=consider-using-with
        self.__testoutput = self.__tmpDir.name

    def tearDown(self) -> None:
        self.__tmpDir.cleanup()

    def __buildBirdIndex(self) -> FamilyDict:
        """Build index of family identifier correspondences -

        Returns -  a dictionary by alternate family_id with valid id correpsondences.
        """
        self.__lfh.write("\nStarting OeBirdDepictTests __buildBirdIndex\n")
        fD: FamilyDict = {}
        afD: FamilyDict = {}
        pD: Dict[str, Tuple[str, str]] = {}
        bI = PdbxBirdIndex(
            indexPath=os.path.join(self.__testoutput, "bird-index.pic"), verbose=self.__verbose, log=self.__lfh
        )
        familyIdL = bI.getFamilyList()
        self.assertGreater(len(familyIdL), 0)
        for familyId in familyIdL:
            prdIdList = bI.getPrdIdList(familyId)
            for prdId in prdIdList:
                ccId = bI.getChemCompId(prdId)
                ccPath = bI.getChemCompPath(prdId)
                self.__lfh.write("Family %r prdId %r ccId %r ccPath %r\n" % (familyId, prdId, ccId, ccPath))
                if familyId not in fD:
                    fD[familyId] = []
                if (prdId is not None) and (ccId is not None) and (ccPath is not None):
                    fD[familyId].append((familyId, prdId, ccId, ccPath))
                    pD[prdId] = (ccId, ccPath)

        with open(self.__altFamilyIndexPath, encoding="utf-8") as ifh:
            for line in ifh:
                fields = line.split()
                fId = fields[0].strip()
                pIdList = [f.strip() for f in fields[1:]]
                afD[fId] = []
                for pId in pIdList:
                    self.assertIn(pId, pD)
                    afD[fId].append((fId, pId, pD[pId][0], pD[pId][1]))

        self.assertGreater(len(afD), 0)
        return afD

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
                for rCC, rAt, tCC, tAt in aML:
                    self.__lfh.write("%5s %-5s %5s %-5s\n" % (rCC, rAt, tCC, tAt))
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
            imageD[familyId] = imageFileName
        return imageD

    def testFamilyDepiction(self) -> None:
        """Test case -  aligned family members --"""
        self.__lfh.write("\nStarting OeBirdDepictTests testFamilyDepiction\n")
        fD = self.__buildBirdIndex()
        imageD = self.__depictFamilies(fD)
        self.assertEqual(sorted(imageD.keys()), sorted(fId for fId, fmList in fD.items() if fmList))

    def testFamilyDepictionHTMLIndex(self) -> None:
        """Test case -  aligned family members --"""
        self.__lfh.write("\nStarting OeBirdDepictTests testFamilyDepictionHTMLIndex\n")
        fD = self.__buildBirdIndex()
        imageD = self.__depictFamilies(fD)
        tS = time.strftime("%Y %m %d %H:%M:%S", time.localtime())
        htmlPath = os.path.join(self.__testoutput, "index.html")
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
                    ofh.write(
                        '<li> <a href="%s">%s</a> with %2d members.</li>\n' % (imageFileName, familyId, len(fmList))
                    )
            ofh.write("</ul>\n")
            ofh.write("</body>\n")
            ofh.write("</html>\n")

        with open(htmlPath, encoding="utf-8") as ifh:
            html = ifh.read()
        for familyId in imageD:
            self.assertIn('<a href="%s-members.pdf">' % familyId, html)


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
