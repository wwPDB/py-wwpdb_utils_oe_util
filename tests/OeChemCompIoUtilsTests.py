##
# File:    OeChemCompIoUtilsTests.py
# Author:  J. Westbrook
# Date:    23-Jan-2012
# Version: 0.001
#
# Updated:
#     6-Jun-2016 jdw general cleanup
#
#
##
"""
Test cases for persistent storage of serialized OE molecule objects.

"""

__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.01"


import io
import os
import unittest
from typing import TYPE_CHECKING, Dict, List, Tuple

try:
    from openeye.oechem import OEFloatArray  # noqa: F401 pylint: disable=unused-import

    skiptests = False
except ImportError:
    skiptests = True

if not skiptests:
    if TYPE_CHECKING:
        from wwpdb.utils.oe_util.build.OeBuildMol import OeBuildMol  # noqa: F401 pylint: disable=unused-import
    from wwpdb.utils.oe_util.build.OeChemCompIoUtils import OeChemCompIoUtils

ATP_ISO_SMILES = "c1nc(c2c(n1)n(cn2)[C@H]3[C@@H]([C@@H]([C@H](O3)CO[P@@](=O)(O)O[P@](=O)(O)OP(=O)(O)O)O)O)N"
GTP_ISO_SMILES = "c1nc2c(=O)[nH]c(nc2n1[C@H]3[C@@H]([C@@H]([C@H](O3)CO[P@](=O)(O)O[P@](=O)(O)OP(=O)(O)O)O)O)N"

# Expected (isomeric SMILES, number of atoms) keyed by component id
EXPECTED: Dict[str, Tuple[str, int]] = {
    "ATP": (ATP_ISO_SMILES, 47),
    "GTP": (GTP_ISO_SMILES, 48),
    "ARG": ("C(C[C@@H](C(=O)O)N)CNC(=[NH2+])N", 27),
    "MSE": ("C[Se]CC[C@@H](C(=O)O)N", 20),
    "TRP": ("c1ccc2c(c1)c(c[nH]2)C[C@@H](C(=O)O)N", 27),
    "XXATP": (ATP_ISO_SMILES, 47),
}


@unittest.skipIf(skiptests, "Could not import openeye")
class OeChemCompIoUtilsTests(unittest.TestCase):
    def setUp(self) -> None:
        self.__lfh = io.StringIO()
        self.__verbose = True
        self.__here = os.path.abspath(os.path.dirname(__file__))
        self.__datadir = os.path.join(self.__here, "data")

        self.__topCachePath = os.path.join(self.__here, "ligand-dict-v3")
        self.__pathList: List[str] = [
            os.path.join(self.__datadir, "ATP.cif"),
            os.path.join(self.__datadir, "GTP.cif"),
            os.path.join(self.__datadir, "ARG.cif"),
        ]
        self.__idList: List[str] = ["MSE", "GTP", "TRP", "XXATP"]

    def tearDown(self) -> None:
        pass

    def __checkList(self, oemList: List["OeBuildMol"], ccIdList: List[str]) -> None:
        self.assertEqual([oem.getCcId() for oem in oemList], ccIdList)
        for oem in oemList:
            ccId = oem.getCcId()
            assert ccId is not None
            isoSmiles, numAtoms = EXPECTED[ccId]
            self.assertEqual(oem.getTitle(), ccId)
            self.assertEqual(oem.getIsoSMILES(), isoSmiles)
            self.assertGreater(len(oem.getCanSMILES()), 0)
            self.assertNotIn("@", oem.getCanSMILES())
            self.assertEqual(oem.getMol().NumAtoms(), numAtoms)

    def testMakeFromFiles(self) -> None:
        """Test case -  create OE molecules from the input chem comp definition path list."""
        oeU = OeChemCompIoUtils(topCachePath=self.__topCachePath, verbose=self.__verbose, log=self.__lfh)  # pylint: disable=possibly-used-before-assignment
        oemList = oeU.getFromPathList(self.__pathList, use3D=True, coordType="model")
        self.__checkList(oemList, ["ATP", "GTP", "ARG"])

    def testMakeFromIds(self) -> None:
        """Test case -  create OE molecules from the input chem comp id list."""
        oeU = OeChemCompIoUtils(topCachePath=self.__topCachePath, verbose=self.__verbose, log=self.__lfh)
        oemList = oeU.getFromIdList(self.__idList, use3D=True, coordType="model")
        self.__checkList(oemList, self.__idList)
        # XXATP is a renamed copy of ATP, located with the two-character hash directory
        self.assertEqual(oemList[3].getInChIKey(), "ZKHQWZAMYRWXGA-KQYNXXCUSA-N")


def suite1() -> unittest.TestSuite:
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeChemCompIoUtilsTests("testMakeFromFiles"))
    suiteSelect.addTest(OeChemCompIoUtilsTests("testMakeFromIds"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = suite1()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
