##
# File:    PdbxBuildChemCompExtraTests.py
#
# Additional test cases for PdbxBuildChemComp covering accessors, stereo
# assignment, writing alternate formats and failure branches.
##
"""
Additional tests for the PdbxBuildChemComp class.
"""

import io
import os
import tempfile
import unittest
from typing import Any, Dict, List

try:
    from openeye.oechem import OEFloatArray  # noqa: F401 pylint: disable=unused-import

    skiptests = False
except ImportError:
    skiptests = True

if not skiptests:
    from openeye.oechem import (
        OEAssignAromaticFlags,
        OEFindRingAtomsAndBonds,
        OEGraphMol,
        OEParseSmiles,
        OEPerceiveChiral,
        OETriposAtomNames,
    )

    from wwpdb.utils.oe_util.build.PdbxBuildChemComp import PdbxBuildChemComp


def makeMol(smiles: str, nameAtoms: bool = True) -> "OEGraphMol":
    mol = OEGraphMol()  # pylint: disable=possibly-used-before-assignment
    if not OEParseSmiles(mol, smiles):  # pylint: disable=possibly-used-before-assignment
        raise ValueError("Bad SMILES %s" % smiles)
    OEFindRingAtomsAndBonds(mol)  # pylint: disable=possibly-used-before-assignment
    OEAssignAromaticFlags(mol)  # pylint: disable=possibly-used-before-assignment
    OEPerceiveChiral(mol)  # pylint: disable=possibly-used-before-assignment
    if nameAtoms:
        OETriposAtomNames(mol)  # pylint: disable=possibly-used-before-assignment
    return mol


@unittest.skipIf(skiptests, "Could not import openeye")
class PdbxBuildChemCompExtraTests(unittest.TestCase):
    def setUp(self) -> None:
        self.__lfh = io.StringIO()
        self.__tmpDir = tempfile.TemporaryDirectory()  # pylint: disable=consider-using-with

    def tearDown(self) -> None:
        self.__tmpDir.cleanup()

    def __newBuilder(self) -> "PdbxBuildChemComp":
        return PdbxBuildChemComp(verbose=False, log=self.__lfh)  # pylint: disable=possibly-used-before-assignment

    def __atomColumn(self, ccB: "PdbxBuildChemComp", attr: str) -> List[Any]:
        return [row["_chem_comp_atom." + attr] for row in ccB.getCategory("chem_comp_atom")]

    def testChemCompCategory(self) -> None:
        """chem_comp category values for a molecule without a name."""
        ccB = self.__newBuilder()
        ccB.setOeMol(makeMol("CCO"), "eoh", missingModelXyz=False)
        ccL = ccB.getChemCompDict()
        self.assertEqual(len(ccL), 1)
        cc: Dict[str, Any] = ccL[0]
        self.assertEqual(cc["_chem_comp.id"], "EOH")
        self.assertEqual(cc["_chem_comp.name"], "?")
        self.assertEqual(cc["_chem_comp.formula"], "C2H6O")
        self.assertEqual(cc["_chem_comp.pdbx_formal_charge"], 0)
        self.assertEqual(cc["_chem_comp.three_letter_code"], "EOH")
        self.assertEqual(cc["_chem_comp.pdbx_model_coordinates_missing_flag"], "N")
        self.assertAlmostEqual(float(cc["_chem_comp.formula_weight"]), 46.07, places=2)
        self.assertEqual(ccB.getChemCompCategory(), ccL)
        self.assertEqual(ccB.getCategory(), ccL)

    def testNameAndMissingFlag(self) -> None:
        """Name and missing model coordinate flag are set; charge summed."""
        ccB = self.__newBuilder()
        ccB.setOeMol(makeMol("C[NH3+]"), "MAM_1", name="methylammonium")
        cc = ccB.getChemCompDict()[0]
        self.assertEqual(cc["_chem_comp.name"], "methylammonium")
        self.assertEqual(cc["_chem_comp.pdbx_formal_charge"], 1)
        self.assertEqual(cc["_chem_comp.three_letter_code"], "MAM")
        self.assertEqual(cc["_chem_comp.pdbx_model_coordinates_missing_flag"], "Y")

    def testDescriptors(self) -> None:
        """Descriptor category contains SMILES, InChI and InChIKey rows."""
        ccB = self.__newBuilder()
        ccB.setOeMol(makeMol("C/C=C\\C"), "BUT")
        rows = ccB.getChemCompCategory("pdbx_chem_comp_descriptor")
        dD = {r["_pdbx_chem_comp_descriptor.type"]: r["_pdbx_chem_comp_descriptor.descriptor"] for r in rows}
        self.assertEqual(dD["SMILES_CANONICAL"], "C/C=C\\C")
        self.assertEqual(dD["SMILES"], "CC=CC")
        self.assertEqual(dD["InChI"], "InChI=1S/C4H8/c1-3-4-2/h3-4H,1-2H3/b4-3-")
        self.assertEqual(dD["InChIKey"], "IAQRGUVFOMOMEM-ARJAWSKDSA-N")
        audit = ccB.getCategory("pdbx_chem_comp_audit")
        self.assertEqual(audit[0]["_pdbx_chem_comp_audit.action_type"], "CREATE")

    def testBondStereo(self) -> None:
        """E/Z bond stereo is reported in the bond category."""
        for smiles, expected in [("C/C=C/C", "E"), ("C/C=C\\C", "Z"), ("CC=CC", "N"), ("C/C=C(/C)C", "N")]:
            ccB = self.__newBuilder()
            ccB.setOeMol(makeMol(smiles), "BUT")
            bondL = ccB.getBondList()
            orders = [row[3] for row in bondL]
            stereo = [row[5] for row in bondL]
            dIdx = orders.index("DOUB")
            self.assertEqual(stereo[dIdx], expected, smiles)
            self.assertEqual([s for i, s in enumerate(stereo) if i != dIdx], ["N"] * (len(stereo) - 1))

    def testAtomStereoAndFlags(self) -> None:
        """Atom R/S stereo and default leaving atom flag."""
        for smiles, expected in [("C[C@H](N)O", "R"), ("C[C@@H](N)O", "S"), ("CC(N)O", "N"), ("C[C@H](C)O", "N")]:
            ccB = self.__newBuilder()
            ccB.setOeMol(makeMol(smiles), "CHI")
            self.assertEqual(self.__atomColumn(ccB, "pdbx_stereo_config")[1], expected, smiles)
            self.assertEqual(set(self.__atomColumn(ccB, "pdbx_leaving_atom_flag")), {"N"})

    def testAromaticAndAlign(self) -> None:
        """Aromatic flags on atoms and bonds; alignment for two-letter elements."""
        ccB = self.__newBuilder()
        ccB.setOeMol(makeMol("c1ccccc1Cl"), "CLB")
        arom = self.__atomColumn(ccB, "pdbx_aromatic_flag")
        symbols = self.__atomColumn(ccB, "type_symbol")
        align = self.__atomColumn(ccB, "pdbx_align")
        self.assertEqual(arom.count("Y"), 6)
        self.assertEqual(align[symbols.index("Cl")], 0)
        self.assertEqual(align[symbols.index("C")], 1)
        bondArom = [row[4] for row in ccB.getBondList()]
        self.assertEqual(bondArom.count("Y"), 6)

    def testUnnamedAtomsSkipBonds(self) -> None:
        """Bonds between unnamed atoms are omitted."""
        ccB = self.__newBuilder()
        ccB.setOeMol(makeMol("CCO", nameAtoms=False), "EOH")
        self.assertEqual(len(ccB.getCategory("chem_comp_atom")), 3)
        self.assertEqual(ccB.getBondList(), [])

    def testBlocksUpdateWrite(self) -> None:
        """Block selection, attribute updates and writing a definition file."""
        ccB = self.__newBuilder()
        ccB.setOeMol(makeMol("CCO"), "EOH")
        self.assertFalse(ccB.setBlock("NOPE"))
        self.assertTrue(ccB.setBlock("EOH"))
        self.assertTrue(ccB.update("chem_comp", "name", "ethanol"))
        self.assertEqual(ccB.getChemCompDict()[0]["_chem_comp.name"], "ethanol")
        self.assertIsInstance(ccB.complyStyle(), bool)
        self.assertTrue(ccB.newBlock("NEW"))

        fp = os.path.join(self.__tmpDir.name, "EOH.cif")
        self.assertTrue(ccB.write(fp))
        with open(fp, encoding="utf-8") as ifh:
            text = ifh.read()
        self.assertIn("data_EOH", text)
        self.assertIn("ethanol", text)

    def testWriteOther(self) -> None:
        """Write the molecule in other formats and handle failure."""
        ccB = self.__newBuilder()
        mol = makeMol("CCO")
        sdfPath = os.path.join(self.__tmpDir.name, "EOH.sdf")
        self.assertTrue(ccB.writeOther(mol, sdfPath, title="EOH"))
        with open(sdfPath, encoding="utf-8") as ifh:
            self.assertEqual(ifh.readline().strip(), "EOH")
        mol2Path = os.path.join(self.__tmpDir.name, "EOH.mol2")
        self.assertTrue(ccB.writeOther(mol, mol2Path, title="EOH", constantMol=True))
        self.assertGreater(os.path.getsize(mol2Path), 0)
        self.assertFalse(ccB.writeOther(None, os.path.join(self.__tmpDir.name, "bad.sdf")))
        self.assertIn("writeOther FAILING", self.__lfh.getvalue())


if __name__ == "__main__":
    unittest.main()
