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

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.01"


import sys
import unittest
import traceback
import platform
import time
import os
import os.path

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


@unittest.skipIf(skiptests, "Could not import openeye")
class OeBirdDepictTests(unittest.TestCase):
    def setUp(self):
        self.__lfh = sys.stderr
        self.__verbose = True
        self.__here = os.path.abspath(os.path.dirname(__file__))
        self.__testoutput = os.path.join(self.__here, "test-output", platform.python_version())

    def tearDown(self):
        pass

    def __testBuildBirdIndex(self):
        """Test case -  build index of family identifier correspondences -

        Returns -  a dictionary by family_id with valid id correpsondences.
        """
        self.__lfh.write("\nStarting OeBirdDepictTests __testBuildBirdIndex\n")
        fD = {}
        try:
            bI = PdbxBirdIndex(indexPath=os.path.join(self.__testoutput, "bird-index.pic"), verbose=self.__verbose, log=self.__lfh)
            familyIdL = bI.getFamilyList()
            for familyId in familyIdL:
                prdIdList = bI.getPrdIdList(familyId)
                prdPathList = bI.getPrdPathList(familyId)  # noqa: F841 pylint: disable=unused-variable
                for prdId in prdIdList:
                    ccId = bI.getChemCompId(prdId)
                    ccPath = bI.getChemCompPath(prdId)
                    self.__lfh.write("Family %r prdId %r ccId %r ccPath %r\n" % (familyId, prdId, ccId, ccPath))
                    if familyId not in fD:
                        fD[familyId] = []
                    if (prdId is not None) and (ccId is not None) and (ccPath is not None):
                        fD[familyId].append((familyId, prdId, ccId, ccPath))
        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()

        return fD

    def testFamilyDepiction(self):
        """Test case -  aligned family members --"""
        self.__lfh.write("\nStarting OeBirdDepictTests testFamilyDepiction\n")
        try:
            fD = self.__testBuildBirdIndex()
            famList = sorted(fD.keys())
            for fId in famList:
                fmList = fD[fId]
                if len(fmList) < 1:
                    continue

                (familyId, refPrdId, _refCcId, refCcPath) = fmList[0]
                imageFileName = familyId + "-members.pdf"
                if os.access(imageFileName, os.F_OK):
                    continue
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
                    if len(aML) > 0:
                        for (rCC, rAt, tCC, tAt) in aML:
                            self.__lfh.write("%5s %-5s %5s %-5s\n" % (rCC, rAt, tCC, tAt))
                else:
                    oeU = OeChemCompIoUtils(verbose=self.__verbose, log=self.__lfh)
                    oemList = oeU.getFromPathList([refCcPath], use3D=False)
                    oed = OeDepictMultiPage(verbose=self.__verbose, log=self.__lfh)
                    oed.setMolTitleList([(refPrdId, oemList[0], refPrdId)])
                    oed.setDisplayOptions(labelAtomName=False, labelAtomCIPStereo=True, labelAtomIndex=False, labelBondIndex=False, gridRows=3, gridCols=3, bondDisplayWidth=0.5)
                    oed.prepare()
                    oed.write(imageFileName)
        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testFamilyDepictionHTMLIndex(self):
        """Test case -  aligned family members --"""
        self.__lfh.write("\nStarting OeBirdDepictTests testFamilyDepictionHTMLIndex\n")
        tS = time.strftime("%Y %m %d %H:%M:%S", time.localtime())
        ofh = open(os.path.join(self.__testoutput, "index.html"), "w")
        ofh.write("<html>\n")
        ofh.write("<body>\n")
        ofh.write("<h4>Index of family chemical diagrams produced on: %s</h4>\n" % tS)
        ofh.write("<ul>\n")
        try:
            fD = self.__testBuildBirdIndex()
            famList = sorted(fD.keys())
            for fId in famList:
                fmList = fD[fId]
                if len(fmList) < 1:
                    continue
                (familyId, _refPrdId, _refCcId, _refCcPath) = fmList[0]
                imageFileName = familyId + "-members.pdf"
                if os.access(imageFileName, os.R_OK):
                    ofh.write('<li> <a href="%s">%s</a> with %2d members.</li>\n' % (imageFileName, familyId, len(fmList)))
            ofh.write("</ul>\n")
            ofh.write("</body>\n")
            ofh.write("</html>\n")
            ofh.close()

        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()


def suiteDepictFamily():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeBirdDepictTests("testFamilyDepiction"))
    return suiteSelect


def suiteHTMLIndexFamily():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeBirdDepictTests("testFamilyDepictionHTMLIndex"))
    return suiteSelect


if __name__ == "__main__":
    mySuite1 = suiteDepictFamily()
    unittest.TextTestRunner(verbosity=2).run(mySuite1)
    mySuite1 = suiteHTMLIndexFamily()
    unittest.TextTestRunner(verbosity=2).run(mySuite1)
