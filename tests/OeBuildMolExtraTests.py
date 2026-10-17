##
# File:    OeBuildMolExtraTests.py
#
# Additional test cases for OeBuildMol covering alternate coordinate types,
# SDF/SMILES import, serialization round trips, getters and failure branches.
##
"""
Additional tests for the OeBuildMol class.
"""

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
    from openeye.oechem import OEGraphMol, OEParseSmiles

    from wwpdb.utils.oe_util.build.OeBuildMol import OeBuildMol

HERE = os.path.abspath(os.path.dirname(__file__))
DATADIR = os.path.join(HERE, "data")

ATP_FORMULA = "C10H16N5O13P3"
ATP_INCHIKEY = "ZKHQWZAMYRWXGA-KQYNXXCUSA-N"
ATP_ELEMENTS: Dict[int, int] = {1: 16, 6: 10, 7: 5, 8: 13, 15: 3}

# A small synthetic definition (hexa-2,4-diene) with E and Z bond stereo assignments
HEXADIENE_CIF = """data_HXD
#
_chem_comp.id HXD
_chem_comp.name hexadiene
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.charge
_chem_comp_atom.pdbx_aromatic_flag
_chem_comp_atom.pdbx_leaving_atom_flag
_chem_comp_atom.pdbx_stereo_config
_chem_comp_atom.pdbx_ordinal
HXD C1 C 0 N N N 1
HXD C2 C 0 N N N 2
HXD C3 C 0 N N N 3
HXD C4 C 0 N N N 4
HXD C5 C 0 N N N 5
HXD C6 C 0 N N N 6
#
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
_chem_comp_bond.pdbx_stereo_config
_chem_comp_bond.pdbx_ordinal
HXD C1 C2 SING N N 1
HXD C2 C3 DOUB N E 2
HXD C3 C4 SING N N 3
HXD C4 C5 DOUB N Z 4
HXD C5 C6 SING N N 5
#
"""


@unittest.skipIf(skiptests, "Could not import openeye")
class OeBuildMolExtraTests(unittest.TestCase):
    def setUp(self) -> None:
        self.__lfh = io.StringIO()
        self.__atpPath = os.path.join(DATADIR, "ATP.cif")
        self.__sdfPath = os.path.join(DATADIR, "ATP.sdf")
        self.__tmpDir = tempfile.TemporaryDirectory()  # pylint: disable=consider-using-with

    def tearDown(self) -> None:
        self.__tmpDir.cleanup()

    def __newBuilder(self) -> "OeBuildMol":
        return OeBuildMol(verbose=False, log=self.__lfh)  # pylint: disable=possibly-used-before-assignment

    def testSetChemCompPath(self) -> None:
        """Reading a definition file returns the component id; a missing file returns None."""
        oem = self.__newBuilder()
        self.assertEqual(oem.setChemCompPath(self.__atpPath), "ATP")
        self.assertEqual(oem.getCcId(), "ATP")
        self.assertIsNone(oem.setChemCompPath(os.path.join(self.__tmpDir.name, "missing.cif")))
        self.assertIn("setChemCompPath) Fails", self.__lfh.getvalue())

    def testBuild3DCoordTypes(self) -> None:
        """Build with model and ideal coordinates yields the same chemistry."""
        oem = self.__newBuilder()
        oem.setChemCompPath(self.__atpPath)
        smiL: List[str] = []
        for coordType in ["model", "ideal"]:
            self.assertTrue(oem.build3D(coordType=coordType))
            self.assertEqual(oem.getTitle(), "ATP")
            self.assertEqual(oem.getFormula(), ATP_FORMULA)
            self.assertEqual(oem.getInChIKey(), ATP_INCHIKEY)
            self.assertTrue(oem.getInChI().startswith("InChI=1S/C10H16N5O13P3/"))
            self.assertEqual(oem.getMol().NumAtoms(), 47)
            self.assertEqual(oem.getElementCounts(), ATP_ELEMENTS)
            smiL.append(oem.getIsoSMILES())
        self.assertEqual(smiL[0], smiL[1])

    def testBuild3DNoTitle(self) -> None:
        """Building without setting the title leaves the molecule title empty."""
        oem = self.__newBuilder()
        oem.setChemCompPath(self.__atpPath)
        self.assertTrue(oem.build3D(coordType="ideal", setTitle=False))
        self.assertEqual(oem.getTitle(), "")
        self.assertEqual(oem.getCcId(), "ATP")

    def testBuild3DDebug(self) -> None:
        """Debug output is written to the log during the build."""
        oem = self.__newBuilder()
        oem.setDebug(True)
        oem.setChemCompPath(self.__atpPath)
        self.assertTrue(oem.build3D(coordType="model"))
        self.assertIn("Atom - PG", self.__lfh.getvalue())

    def testBuildFailures(self) -> None:
        """Building without source data categories fails gracefully."""
        oem = self.__newBuilder()
        oem.set("XXX")
        self.assertFalse(oem.build3D())
        self.assertIn("build3D) Failing", self.__lfh.getvalue())
        self.assertFalse(oem.build2D())

    def testBuild2DBondStereo(self) -> None:
        """2D build assigns E/Z bond stereo from the definition."""
        cifPath = os.path.join(self.__tmpDir.name, "HXD.cif")
        with open(cifPath, "w", encoding="utf-8") as ofh:
            ofh.write(HEXADIENE_CIF)
        oem = self.__newBuilder()
        oem.setDebug(True)
        self.assertEqual(oem.setChemCompPath(cifPath), "HXD")
        self.assertTrue(oem.build2D())
        self.assertEqual(oem.getFormula(), "C6")
        self.assertEqual(oem.getElementCounts(), {6: 6})
        # No implicit hydrogens are added: C2=C3 is E (trans) and C4=C5 is Z (cis)
        self.assertEqual(oem.getIsoSMILES(), "[C]/[C]=[C]/[C]=[C]\\[C]")
        self.assertIn("build2d - atom", self.__lfh.getvalue())

    def testSerializeRoundTrip(self) -> None:
        """Serialize and deserialize a 3D molecule (with debug output)."""
        oem = self.__newBuilder()
        oem.setDebug(True)
        oem.setChemCompPath(self.__atpPath)
        self.assertTrue(oem.build3D(coordType="model"))
        oeS = oem.serialize()
        self.assertGreater(len(oeS), 0)

        oemD = self.__newBuilder()
        oemD.setDebug(True)
        self.assertTrue(oemD.deserialize(oeS))
        self.assertEqual(oemD.getCcId(), "ATP")
        self.assertEqual(oemD.getIsoSMILES(), oem.getIsoSMILES())
        self.assertEqual(oemD.getElementCounts(), ATP_ELEMENTS)
        self.assertIn("deserialize) mols  1", self.__lfh.getvalue())

    def testDeserializeEmpty(self) -> None:
        """Deserializing empty input returns False and no element counts."""
        oem = self.__newBuilder()
        self.assertFalse(oem.deserialize(b""))
        self.assertEqual(oem.getElementCounts(), {})

    def testSetOeMol(self) -> None:
        """Load an existing OE molecule and rename the atoms."""
        mol = OEGraphMol()  # pylint: disable=possibly-used-before-assignment
        self.assertTrue(OEParseSmiles(mol, "CCO"))  # pylint: disable=possibly-used-before-assignment
        mol.SetTitle("ETH")
        oem = self.__newBuilder()
        oem.setOeMol(mol, "EOH")
        self.assertEqual(oem.getCcId(), "EOH")
        self.assertEqual(oem.getTitle(), "ETH")
        self.assertEqual(oem.getElementCounts(), {6: 2, 8: 1})
        oem.simpleAtomNames()
        names = sorted(atom.GetName() for atom in oem.getMol().GetAtoms())
        self.assertEqual(names, ["C1", "C2", "O1"])

    def testImportSmiles(self) -> None:
        """Import valid and invalid SMILES."""
        oem = self.__newBuilder()
        self.assertTrue(oem.importSmiles("OCC"))
        self.assertEqual(oem.getCanSMILES(), "CCO")
        self.assertEqual(oem.getFormula(), "C2H6O")
        self.assertFalse(oem.importSmiles("C1CC(("))

    def testGraphMolSuppressH(self) -> None:
        """Hydrogens are removed from the current molecule."""
        oem = self.__newBuilder()
        oem.setChemCompPath(self.__atpPath)
        self.assertTrue(oem.build3D(coordType="model"))
        mol = oem.getGraphMolSuppressH()
        self.assertEqual(mol.NumAtoms(), 31)

    def testImportFile2D(self) -> None:
        """Import an SDF file with 2D perception and check generated coordinate list."""
        oem = self.__newBuilder()
        self.assertEqual(oem.getCoords(), [])
        self.assertTrue(oem.importFile(self.__sdfPath, type="2D"))
        coords = oem.getCoords()
        self.assertEqual(len(coords), 47)
        self.assertEqual(coords[0][2], 15)
        self.assertEqual(coords[0][3], "P1")
        self.assertEqual(coords[1][3], "O1")
        self.assertEqual(len(coords[0]), 8)
        self.assertEqual(oem.getFormula(), ATP_FORMULA)

    def testImportFile3D(self) -> None:
        """Import an SDF file with 3D perception."""
        oem = self.__newBuilder()
        self.assertTrue(oem.importFile(self.__sdfPath, type="3D"))
        self.assertEqual(len(oem.getCoords()), 47)
        self.assertEqual(oem.getInChIKey(), ATP_INCHIKEY)

    def testUpdatePerceptions3D(self) -> None:
        """Re-running 3D perceptions on an imported molecule preserves stereochemistry."""
        oem = self.__newBuilder()
        self.assertTrue(oem.importFile(self.__sdfPath, type="2D"))
        oem.updatePerceptions3D()
        self.assertEqual(oem.getMol().GetDimension(), 3)
        self.assertEqual(oem.getInChIKey(), ATP_INCHIKEY)

    def testImportFileMissing(self) -> None:
        """A missing input file returns False."""
        oem = self.__newBuilder()
        self.assertFalse(oem.importFile(os.path.join(self.__tmpDir.name, "missing.sdf")))


if __name__ == "__main__":
    unittest.main()
