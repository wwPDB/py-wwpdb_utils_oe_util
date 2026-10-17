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


import io
import os
import tempfile
import unittest
from typing import TYPE_CHECKING, Dict, List

try:
    from openeye.oechem import OEFloatArray  # noqa: F401 pylint: disable=unused-import

    skiptests = False
except ImportError:
    skiptests = True

if not skiptests:
    from mmcif_utils.persist.PdbxPyIoAdapter import PdbxPyIoAdapter as PdbxIoAdapter

    if TYPE_CHECKING:
        from wwpdb.utils.oe_util.build.OeBuildMol import OeBuildMol  # noqa: F401 pylint: disable=unused-import
    from wwpdb.utils.oe_util.build.OeChemCompIoUtils import OeChemCompIoUtils
    from wwpdb.utils.oe_util.build.PdbxBuildChemComp import PdbxBuildChemComp

# Expected formula and standard InChIKey for each component used below
EXPECTED: Dict[str, Dict[str, str]] = {
    "ATP": {"formula": "C10H16N5O13P3", "inchikey": "ZKHQWZAMYRWXGA-KQYNXXCUSA-N"},
    "GTP": {"formula": "C10H16N5O14P3", "inchikey": "XKMLYUALXHKNFT-UUOKFMHZSA-N"},
    "ARG": {"formula": "C6H15N4O2+", "inchikey": "ODKSFYDXXFIFQN-BYPYZUCNSA-O"},
    "MSE": {"formula": "C5H11NO2Se", "inchikey": "RJFAYQIBOAGBLC-BYPYZUCNSA-N"},
    "TRP": {"formula": "C11H12N2O2", "inchikey": "QIVBCDIJIAJPQS-VIFPVBQESA-N"},
}


@unittest.skipIf(skiptests, "Requires openeye library")
class PdbxBuildChemCompTests(unittest.TestCase):
    def setUp(self) -> None:
        self.__lfh = io.StringIO()
        self.__verbose = True
        self.__here = os.path.abspath(os.path.dirname(__file__))
        self.__topCachePath = os.path.join(self.__here, "ligand-dict-v3")
        self.__datadir = os.path.join(self.__here, "data")
        self.__tmpDir = tempfile.TemporaryDirectory()  # pylint: disable=consider-using-with
        self.__testoutput = self.__tmpDir.name

        self.__pathList: List[str] = [
            os.path.join(self.__datadir, "ATP.cif"),
            os.path.join(self.__datadir, "GTP.cif"),
            os.path.join(self.__datadir, "ARG.cif"),
        ]
        self.__idList: List[str] = ["MSE", "GTP", "TRP"]

    def tearDown(self) -> None:
        self.__tmpDir.cleanup()

    def __molListCheck(self, oeMolList: List["OeBuildMol"], ccIdList: List[str]) -> None:
        self.assertEqual([oem.getCcId() for oem in oeMolList], ccIdList)
        inKyD: Dict[str, List[str]] = {}
        for oem in oeMolList:
            ccId = oem.getCcId()
            assert ccId is not None
            self.assertEqual(oem.getTitle(), ccId)
            self.assertEqual(oem.getFormula(), EXPECTED[ccId]["formula"])
            self.assertGreater(len(oem.getCanSMILES()), 0)
            self.assertGreater(len(oem.getIsoSMILES()), 0)
            inKy = oem.getInChIKey()
            self.assertEqual(inKy, EXPECTED[ccId]["inchikey"])
            inKyD.setdefault(inKy, []).append(ccId)
        self.assertEqual(len(inKyD), len(ccIdList))

    def __writeAndVerify(self, oeMolList: List["OeBuildMol"], prefix: str) -> None:
        for oem in oeMolList:
            ccId = oem.getCcId()
            assert ccId is not None
            fp = os.path.join(self.__testoutput, prefix + ccId + ".cif")
            ccB = PdbxBuildChemComp(verbose=self.__verbose, log=self.__lfh)  # pylint: disable=possibly-used-before-assignment
            ccB.setOeMol(oem.getMol(), ccId, name=ccId)
            self.assertTrue(ccB.write(filePath=fp))
            self.assertTrue(os.path.isfile(fp))
            self.assertGreater(os.path.getsize(fp), 0)
            # Read back the written definition and compare with the source molecule
            myReader = PdbxIoAdapter(self.__verbose, self.__lfh)  # pylint: disable=possibly-used-before-assignment
            self.assertTrue(myReader.read(pdbxFilePath=fp))
            containerList = myReader.getContainerList()
            self.assertEqual(len(containerList), 1)
            container = containerList[0]
            self.assertEqual(container.getName(), ccId)
            self.assertEqual(container.getObj("chem_comp").getValue("formula", 0), EXPECTED[ccId]["formula"])
            self.assertEqual(container.getObj("chem_comp_atom").getRowCount(), oem.getMol().NumAtoms())
            self.assertEqual(container.getObj("chem_comp_bond").getRowCount(), oem.getMol().NumBonds())

    def testBuildFromFiles(self) -> None:
        """Test case -  create OE molecules from the input chem comp definition path list."""
        oeU = OeChemCompIoUtils(topCachePath=self.__topCachePath, verbose=self.__verbose, log=self.__lfh)  # pylint: disable=possibly-used-before-assignment
        oemList = oeU.getFromPathList(self.__pathList, use3D=True, coordType="model")
        self.__molListCheck(oemList, ["ATP", "GTP", "ARG"])
        self.__writeAndVerify(oemList, "FF_")
        self.assertEqual(sorted(os.listdir(self.__testoutput)), ["FF_ARG.cif", "FF_ATP.cif", "FF_GTP.cif"])

    def testBuildFromIds(self) -> None:
        """Test case -  create OE molecules from the input chem comp id list."""
        oeU = OeChemCompIoUtils(topCachePath=self.__topCachePath, verbose=self.__verbose, log=self.__lfh)
        oemList = oeU.getFromIdList(self.__idList, use3D=True, coordType="model")
        self.__molListCheck(oemList, self.__idList)
        self.__writeAndVerify(oemList, "ID_")
        self.assertEqual(sorted(os.listdir(self.__testoutput)), ["ID_GTP.cif", "ID_MSE.cif", "ID_TRP.cif"])


def suiteBuildTests() -> unittest.TestSuite:
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(PdbxBuildChemCompTests("testBuildFromIds"))
    suiteSelect.addTest(PdbxBuildChemCompTests("testBuildFromFiles"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = suiteBuildTests()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
