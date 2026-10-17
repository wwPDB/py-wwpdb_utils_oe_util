##
# File:    OeBuildModelMolTests.py
# Author:  jdw
# Date:    29-Nov-2014
# Version: 0.001
#
# Updates:
#     6-Jun-2016  jdw general cleanup -
#
##
"""
A collection of tests for the OEBuildModelMol and related classes.

"""

__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.01"


import io
import os
import unittest
from typing import List, Optional

try:
    from openeye.oechem import OEFloatArray  # noqa: F401 pylint: disable=unused-import

    skiptests = False
except ImportError:
    skiptests = True

if not skiptests:
    from wwpdb.utils.oe_util.build.OeBuildModelMol import OeBuildModelMol

from mmcif_utils.chemcomp.PdbxChemCompModel import PdbxChemCompModelDescriptor
from mmcif_utils.chemcomp.PdbxChemCompModelIo import PdbxChemCompModelIo

MODEL_ID = "MTGL00001"
# Model has no stereo centers so canonical and isomeric SMILES coincide
MODEL_SMILES = "CCCCCCCCCCCCCCCCCC(=O)OCC(COC(=O)CCCCCCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCCCCCC"


@unittest.skipIf(skiptests, "Requires oe library")
class OeBuildModelMolTests(unittest.TestCase):
    def setUp(self) -> None:
        self.__lfh = io.StringIO()
        self.__verbose = True
        HERE = os.path.abspath(os.path.dirname(__file__))
        datadir = os.path.join(HERE, "data")
        self.__modelFilePath = os.path.join(datadir, "MTGL00001.cif")
        self.__modelFilePathList: List[str] = [self.__modelFilePath]

    def tearDown(self) -> None:
        pass

    def testBuildFromModel(self) -> None:
        """Test case -  build OE molecule from model instance"""
        oem = OeBuildModelMol(verbose=self.__verbose, log=self.__lfh)  # pylint: disable=possibly-used-before-assignment

        modelId = oem.setChemCompModelPath(self.__modelFilePath)
        self.assertEqual(modelId, MODEL_ID)
        self.assertTrue(oem.build3D())
        self.assertEqual(oem.getTitle(), MODEL_ID)
        self.assertEqual(oem.getCanSMILES(), MODEL_SMILES)
        self.assertEqual(oem.getIsoSMILES(), MODEL_SMILES)

        pccm = PdbxChemCompModelIo(verbose=self.__verbose, log=self.__lfh)
        self.assertTrue(pccm.setFilePath(self.__modelFilePath))
        dL = pccm.getDescriptorList()
        self.assertGreater(len(dL), 0)
        typeL: List[Optional[str]] = []
        for d in dL:
            # mmcif_utils annotations disagree: getDescriptorList() rows allow None values
            pd = PdbxChemCompModelDescriptor(d, verbose=self.__verbose, log=self.__lfh)  # type: ignore[arg-type]
            typeL.append(pd.getType())
            if pd.getType() == "SMILES_CANNONICAL":
                # The reconstructed molecule must match the model's stored descriptor
                self.assertEqual(pd.getDescriptor(), oem.getIsoSMILES())
        self.assertIn("SMILES_CANNONICAL", typeL)

    def testSerialize3D(self) -> None:
        """Test case -  build OE molecule using 3D data in the chemical component model instance and
        then serialize and deserialize this molecule.

        """
        oem = OeBuildModelMol(verbose=self.__verbose, log=self.__lfh)
        for pth in self.__modelFilePathList:
            modelId = oem.setChemCompModelPath(pth)
            self.assertEqual(modelId, MODEL_ID)
            self.assertTrue(oem.build3D())
            self.assertEqual(oem.getTitle(), MODEL_ID)
            oeS = oem.serialize()
            self.assertGreater(len(oeS), 0)
            oemD = OeBuildModelMol(verbose=self.__verbose, log=self.__lfh)
            self.assertTrue(oemD.deserialize(oeS))
            self.assertEqual(oemD.getTitle(), MODEL_ID)
            self.assertEqual(oemD.getCanSMILES(), oem.getCanSMILES())
            self.assertEqual(oemD.getIsoSMILES(), oem.getIsoSMILES())
            mol = oem.getMol()
            molD = oemD.getMol()
            assert mol is not None and molD is not None
            self.assertEqual(molD.NumAtoms(), mol.NumAtoms())


def suiteOeModelBuild() -> unittest.TestSuite:
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeBuildModelMolTests("testBuildFromModel"))
    suiteSelect.addTest(OeBuildModelMolTests("testSerialize3D"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = suiteOeModelBuild()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
