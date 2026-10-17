##
#
# File:    OeBuildMolTests.py
# Author:  jdw
# Date:    25-Sept-2011
# Version: 0.001
#
# Updates:
#  1-Oct-2011 jdw Add test cases for 2D molecule builder
# 19-Feb-2012 jdw Add test cases serialization/deserialization
# 24-Feb-2012 jdw Revised to support iterators and common source data input.
#  6-Jun-2016 jdw General cleanup
#
##
"""
A collection of tests for the OEBuildMol and related classes.

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
from typing import Dict, List

try:
    from openeye.oechem import OEFloatArray  # noqa: F401 pylint: disable=unused-import

    skiptests = False
except ImportError:
    skiptests = True

if not skiptests:
    from mmcif_utils.persist.PdbxPyIoAdapter import PdbxPyIoAdapter as PdbxIoAdapter

    from wwpdb.utils.oe_util.build.OeBuildMol import OeBuildMol

ATP_ISO_SMILES = "c1nc(c2c(n1)n(cn2)[C@H]3[C@@H]([C@@H]([C@H](O3)CO[P@@](=O)(O)O[P@](=O)(O)OP(=O)(O)O)O)O)N"
ATP_CAN_SMILES = "c1nc(c2c(n1)n(cn2)C3C(C(C(O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N"

# Expected (title, canonical SMILES, isomeric SMILES) for each definition in the path list
EXPECTED: Dict[str, Dict[str, str]] = {
    "ATP": {"can": ATP_CAN_SMILES, "iso": ATP_ISO_SMILES},
    "GTP": {
        "can": "c1nc2c(=O)[nH]c(nc2n1C3C(C(C(O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N",
        "iso": "c1nc2c(=O)[nH]c(nc2n1[C@H]3[C@@H]([C@@H]([C@H](O3)CO[P@](=O)(O)O[P@](=O)(O)OP(=O)(O)O)O)O)N",
    },
    "ARG": {"can": "C(CC(C(=O)O)N)CNC(=[NH2+])N", "iso": "C(C[C@@H](C(=O)O)N)CNC(=[NH2+])N"},
}


@unittest.skipIf(skiptests, "Cannot import openeye.oechem for tests")
class OeBuildMolTests(unittest.TestCase):
    def setUp(self) -> None:
        self.__lfh = io.StringIO()
        self.__verbose = True
        self.__here = os.path.abspath(os.path.dirname(__file__))
        self.__tmpDir = tempfile.TemporaryDirectory()  # pylint: disable=consider-using-with
        self.__testoutput = self.__tmpDir.name
        self.__datadir = os.path.join(self.__here, "data")
        self.__sdfFilePath = os.path.join(self.__datadir, "ATP.sdf")
        self.__pathList: List[str] = [
            os.path.join(self.__datadir, "ATP.cif"),
            os.path.join(self.__datadir, "GTP.cif"),
            os.path.join(self.__datadir, "ARG.cif"),
        ]

    def tearDown(self) -> None:
        self.__tmpDir.cleanup()

    def __newBuilder(self) -> "OeBuildMol":
        return OeBuildMol(verbose=self.__verbose, log=self.__lfh)  # pylint: disable=possibly-used-before-assignment

    def __checkMol(self, oem: "OeBuildMol", ccId: str) -> None:
        self.assertEqual(oem.getTitle(), ccId)
        self.assertEqual(oem.getCanSMILES(), EXPECTED[ccId]["can"])
        self.assertEqual(oem.getIsoSMILES(), EXPECTED[ccId]["iso"])

    def __roundTrip(self, oem: "OeBuildMol", ccId: str) -> None:
        oeS = oem.serialize()
        self.assertGreater(len(oeS), 0)
        oemD = self.__newBuilder()
        self.assertTrue(oemD.deserialize(oeS))
        self.assertEqual(oemD.getCcId(), ccId)
        self.assertEqual(oemD.getCanSMILES(), EXPECTED[ccId]["can"])
        self.assertEqual(oemD.getIsoSMILES(), EXPECTED[ccId]["iso"])
        self.assertEqual(oemD.getMol().NumAtoms(), oem.getMol().NumAtoms())

    def __buildAll(self, use3D: bool, serialize: bool, writeCopy: bool) -> None:
        oem = self.__newBuilder()
        if not use3D:
            oem.setDebug(True)
        ccIdL: List[str] = []
        for pth in self.__pathList:
            myReader = PdbxIoAdapter(self.__verbose, self.__lfh)  # pylint: disable=possibly-used-before-assignment
            self.assertTrue(myReader.read(pdbxFilePath=pth))
            if writeCopy:
                outPath = os.path.join(self.__testoutput, "TMP.cif")
                self.assertTrue(myReader.write(pdbxFilePath=outPath))
                self.assertGreater(os.path.getsize(outPath), 0)
            containerList = myReader.getContainerList()
            self.assertEqual(len(containerList), 1)
            for container in containerList:
                ccId = container.getName()
                ccIdL.append(ccId)
                oem.set(
                    ccId,
                    dcChemCompAtom=container.getObj("chem_comp_atom"),
                    dcChemCompBond=container.getObj("chem_comp_bond"),
                )
                if use3D:
                    self.assertTrue(oem.build3D(coordType="model"))
                else:
                    self.assertTrue(oem.build2D())
                self.__checkMol(oem, ccId)
                if serialize:
                    self.__roundTrip(oem, ccId)
        self.assertEqual(ccIdL, ["ATP", "GTP", "ARG"])

    def testBuildMolFromSDF(self) -> None:
        """Test case -  read a test SDF file and build the corresponding OEGraphMol"""
        oem = self.__newBuilder()
        self.assertTrue(oem.importFile(self.__sdfFilePath, type="3D"))
        self.assertEqual(oem.getTitle(), "ATP")
        self.assertEqual(oem.getCanSMILES(), ATP_CAN_SMILES)
        self.assertEqual(oem.getIsoSMILES(), ATP_ISO_SMILES)

    def testBuildFromFiles(self) -> None:
        """Test case -  build OE molecule from definition file source data."""
        self.__buildAll(use3D=True, serialize=False, writeCopy=False)

    def testSerialize2D(self) -> None:
        """Test case -  build OE molecule using 2D data in the definition source data file
        then serialize and deserialize this molecule.
        """
        self.__buildAll(use3D=False, serialize=True, writeCopy=True)
        self.assertIn("OeBuildMol.build2d - atom", self.__lfh.getvalue())

    def testSerialize3D(self) -> None:
        """Test case -  build OE molecule using 3D data in the definition source data file
        then serialize and deserialize this molecule.

        """
        self.__buildAll(use3D=True, serialize=True, writeCopy=True)


def suiteOeBuild() -> unittest.TestSuite:
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeBuildMolTests("testBuildFromFiles"))
    suiteSelect.addTest(OeBuildMolTests("testSerialize2D"))
    suiteSelect.addTest(OeBuildMolTests("testSerialize3D"))
    return suiteSelect


def suiteOeBuildFromOther() -> unittest.TestSuite:
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeBuildMolTests("testBuildMolFromSDF"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = suiteOeBuild()
    unittest.TextTestRunner(verbosity=2).run(mySuite)

    mySuite = suiteOeBuildFromOther()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
