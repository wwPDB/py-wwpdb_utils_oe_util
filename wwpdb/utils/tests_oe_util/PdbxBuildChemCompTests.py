##
# File:    PdbxBuildChemCompTests.py
# Author:  J. Westbrook
# Date:    9-Jun-2016
# Version: 0.001
#
# Updated:
#
#
#
##
"""
Test cases for OE -> CC building tools -

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.01"

import sys
import unittest
import traceback
import os

try:
    from openeye.oechem import OEFloatArray  # noqa: F401 pylint: disable=unused-import

    skiptests = False
except ImportError:
    skiptests = True

if not skiptests:
    from wwpdb.utils.oe_util.build.PdbxBuildChemComp import PdbxBuildChemComp
    from wwpdb.utils.oe_util.build.OeChemCompIoUtils import OeChemCompIoUtils


@unittest.skipIf(skiptests, "Requires openeye library")
class PdbxBuildChemCompTests(unittest.TestCase):
    def setUp(self):
        self.__lfh = sys.stdout
        self.__verbose = True
        self.__here = os.path.abspath(os.path.dirname(__file__))
        self.__topCachePath = os.path.join(self.__here, "ligand-dict-v3")

        self.__pathList = ["../data/ATP.cif", "../data/GTP.cif", "../data/ARG.cif"]
        self.__idList = ["MSE", "GTP", "TRP"]

    def tearDown(self):
        pass

    def __molListPrint(self, oeMolList):
        inKyD = {}
        for oem in oeMolList:
            self.__lfh.write("\nTitle              %s = %s\n" % (oem.getCcId(), oem.getTitle()))
            self.__lfh.write("Formula            %s = %s\n" % (oem.getCcId(), oem.getFormula()))
            self.__lfh.write("SMILES (canonical) %s = %s\n" % (oem.getCcId(), oem.getCanSMILES()))
            self.__lfh.write("SMILES (isomeric)  %s = %s\n" % (oem.getCcId(), oem.getIsoSMILES()))
            inKy = oem.getInChIKey()
            self.__lfh.write("InChIKey (std)     %s = %s\n" % (oem.getCcId(), inKy))
            inKy = oem.getInChIKey()
            if inKy not in inKyD:
                inKyD[inKy] = []
            inKyD[inKy].append(inKy)
            #
        self.__lfh.write("Unique InChiKeys =  %d\n" % (len(inKyD)))

    def testBuildFromFiles(self):
        """Test case -  create OE molecules from the input chem comp definition path list."""
        self.__lfh.write("\nStarting PdbxBuildChemCompTests testBuildFromFiles\n")
        try:
            oemList = []

            oeU = OeChemCompIoUtils(topCachePath=self.__topCachePath, verbose=self.__verbose, log=self.__lfh)
            oemList = oeU.getFromPathList(self.__pathList, use3D=True, coordType="model")
            self.__molListPrint(oemList)
            for oem in oemList:
                ccId = oem.getCcId()
                fp = "FF_" + ccId + ".cif"
                ccB = PdbxBuildChemComp(verbose=self.__verbose, log=self.__lfh)
                ccB.setOeMol(oem.getMol(), ccId, name=ccId)
                ccB.write(filePath=fp)
        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testBuildFromIds(self):
        """Test case -  create OE molecules from the input chem comp definition path list."""
        self.__lfh.write("\nStarting PdbxBuildChemCompTests testBuildFroIds\n")
        try:
            oemList = []
            oeU = OeChemCompIoUtils(topCachePath=self.__topCachePath, verbose=self.__verbose, log=self.__lfh)
            oemList = oeU.getFromIdList(self.__idList, use3D=True, coordType="model")
            self.__molListPrint(oemList)
            for oem in oemList:
                ccId = oem.getCcId()
                fp = "ID_" + ccId + ".cif"
                ccB = PdbxBuildChemComp(verbose=self.__verbose, log=self.__lfh)
                ccB.setOeMol(oem.getMol(), ccId, name=ccId)
                ccB.write(filePath=fp)
        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()


def suiteBuildTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(PdbxBuildChemCompTests("testBuildFromIds"))
    suiteSelect.addTest(PdbxBuildChemCompTests("testBuildFromFiles"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = suiteBuildTests()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
