##
# File:    OeBuildModelMolExtraTests.py
#
# Additional test cases for OeBuildModelMol covering getters, serialization
# round trips, debug output and failure branches.
##
"""
Additional tests for the OeBuildModelMol class.
"""

import io
import os
import tempfile
import unittest

try:
    from openeye.oechem import OEFloatArray  # noqa: F401 pylint: disable=unused-import

    skiptests = False
except ImportError:
    skiptests = True

if not skiptests:
    from wwpdb.utils.oe_util.build.OeBuildModelMol import OeBuildModelMol

HERE = os.path.abspath(os.path.dirname(__file__))
DATADIR = os.path.join(HERE, "data")

MODEL_ID = "MTGL00001"
MODEL_FORMULA = "C57H110O6"
MODEL_INCHIKEY = "DCXXMTOCNZCJGO-UHFFFAOYSA-N"


@unittest.skipIf(skiptests, "Could not import openeye")
class OeBuildModelMolExtraTests(unittest.TestCase):
    def setUp(self) -> None:
        self.__lfh = io.StringIO()
        self.__modelPath = os.path.join(DATADIR, "MTGL00001.cif")
        self.__tmpDir = tempfile.TemporaryDirectory()  # pylint: disable=consider-using-with

    def tearDown(self) -> None:
        self.__tmpDir.cleanup()

    def __newBuilder(self) -> "OeBuildModelMol":
        return OeBuildModelMol(verbose=False, log=self.__lfh)  # pylint: disable=possibly-used-before-assignment

    def testBuildAndGetters(self) -> None:
        """Build the model molecule and check derived descriptors."""
        oem = self.__newBuilder()
        oem.setDebug(True)
        self.assertEqual(oem.setChemCompModelPath(self.__modelPath), MODEL_ID)
        self.assertTrue(oem.build3D())
        self.assertEqual(oem.getTitle(), MODEL_ID)
        self.assertEqual(oem.getFormula(), MODEL_FORMULA)
        self.assertEqual(oem.getInChIKey(), MODEL_INCHIKEY)
        self.assertTrue(oem.getInChI().startswith("InChI=1S/C57H110O6/"))
        mol = oem.getMol()
        assert mol is not None
        self.assertEqual(mol.NumAtoms(), 173)
        self.assertEqual(oem.getElementCounts(), {1: 110, 6: 57, 8: 6})
        self.assertEqual(oem.getCoords(), [])
        self.assertIn("Atom - CA1", self.__lfh.getvalue())

    def testSuppressH(self) -> None:
        """Hydrogens are removed from the model molecule."""
        oem = self.__newBuilder()
        oem.setChemCompModelPath(self.__modelPath)
        self.assertTrue(oem.build3D())
        self.assertEqual(oem.getGraphMolSuppressH().NumAtoms(), 63)

    def testSerializeRoundTrip(self) -> None:
        """Serialize and deserialize the model molecule (with debug output)."""
        oem = self.__newBuilder()
        oem.setDebug(True)
        oem.setChemCompModelPath(self.__modelPath)
        self.assertTrue(oem.build3D())
        oeS = oem.serialize()
        self.assertGreater(len(oeS), 0)

        oemD = self.__newBuilder()
        oemD.setDebug(True)
        self.assertTrue(oemD.deserialize(oeS))
        self.assertEqual(oemD.getTitle(), MODEL_ID)
        self.assertEqual(oemD.getIsoSMILES(), oem.getIsoSMILES())
        self.assertEqual(oemD.getCanSMILES(), oem.getCanSMILES())
        # Element counts are recomputed from the deserialized molecule
        self.assertEqual(oemD.getElementCounts(), {1: 110, 6: 57, 8: 6})
        self.assertIn("deserialize) mols  1", self.__lfh.getvalue())

    def testDeserializeEmpty(self) -> None:
        """Deserializing empty input returns False."""
        oem = self.__newBuilder()
        self.assertFalse(oem.deserialize(b""))
        self.assertEqual(oem.getElementCounts(), {})

    def testMissingModelFile(self) -> None:
        """A missing model file returns None."""
        oem = self.__newBuilder()
        self.assertIsNone(oem.setChemCompModelPath(os.path.join(self.__tmpDir.name, "missing.cif")))
        self.assertIn("setChemCompModelPath) Fails", self.__lfh.getvalue())

    def testBuildBadBond(self) -> None:
        """A bond referencing an unknown atom causes build3D to return False."""
        with open(self.__modelPath, encoding="utf-8") as ifh:
            text = ifh.read()
        badText = text.replace("MTGL00001 CA1 CA2  SING 1", "MTGL00001 CA1 XX99 SING 1")
        self.assertNotEqual(text, badText)
        badPath = os.path.join(self.__tmpDir.name, "BAD.cif")
        with open(badPath, "w", encoding="utf-8") as ofh:
            ofh.write(badText)
        oem = self.__newBuilder()
        self.assertEqual(oem.setChemCompModelPath(badPath), MODEL_ID)
        self.assertFalse(oem.build3D())


if __name__ == "__main__":
    unittest.main()
