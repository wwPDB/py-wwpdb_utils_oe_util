##
# File:    OeDescriptorUtilsExtraTests.py
#
# Additional behaviour tests for OeDescriptorUtils.standardizeSmiles()
##
"""
Additional tests for standardizing SMILES descriptors.
"""

import logging
import unittest
from typing import List

try:
    from openeye.oechem import OEFloatArray  # noqa: F401 pylint: disable=unused-import

    skiptests = False
except ImportError:
    skiptests = True

if not skiptests:
    from wwpdb.utils.oe_util.build.OeDescriptorUtils import OeDescriptorUtils

logger = logging.getLogger(__name__)

_LOGNAME = "wwpdb.utils.oe_util.build.OeDescriptorUtils"


@unittest.skipIf(skiptests, "Requires openeye library")
class OeDescriptorUtilsExtraTests(unittest.TestCase):
    def setUp(self) -> None:
        # GTP with stereo in two different input orderings
        self.__gtpStereoList: List[str] = [
            "NC1=Nc2n(cnc2C(=O)N1)[C@@H]3O[C@H](CO[P@@](O)(=O)O[P@@](O)(=O)O[P](O)(O)=O)[C@@H](O)[C@H]3O",
            "c1nc2c(n1[C@H]3[C@@H]([C@@H]([C@H](O3)CO[P@](=O)(O)O[P@](=O)(O)OP(=O)(O)O)O)O)N=C(NC2=O)N",
        ]
        # GTP without stereo, different input orderings
        self.__gtpFlatList: List[str] = [
            "O=P(O)(O)OP(=O)(O)OP(=O)(O)OCC3OC(n2cnc1c2N=C(N)NC1=O)C(O)C3O",
            "c1nc2c(n1C3C(C(C(O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N=C(NC2=O)N",
        ]
        self.__oedu = OeDescriptorUtils()  # pylint: disable=possibly-used-before-assignment

    def testDefaultTypeIsIsomeric(self) -> None:
        """Default type is ISOMERIC"""
        smi = "C[C@H](N)C(=O)O"
        self.assertEqual(self.__oedu.standardizeSmiles(smi), self.__oedu.standardizeSmiles(smi, type="ISOMERIC"))

    def testSimpleSmiles(self) -> None:
        """Achiral molecule gives same SMILES for both types"""
        self.assertEqual(self.__oedu.standardizeSmiles("CC", type="ISOMERIC"), "CC")
        self.assertEqual(self.__oedu.standardizeSmiles("CC", type="CANNONICAL"), "CC")

    def testIsomericKeepsStereoCanonicalDropsIt(self) -> None:
        """ISOMERIC retains stereo marks while CANNONICAL strips them"""
        smi = "C[C@H](N)C(=O)O"
        iso = self.__oedu.standardizeSmiles(smi, type="ISOMERIC")
        can = self.__oedu.standardizeSmiles(smi, type="CANNONICAL")
        self.assertIsNotNone(iso)
        self.assertIsNotNone(can)
        assert iso is not None and can is not None
        self.assertIn("@", iso)
        self.assertNotIn("@", can)
        # Enantiomer is distinguished by ISOMERIC but not CANNONICAL
        ent = "C[C@@H](N)C(=O)O"
        self.assertNotEqual(iso, self.__oedu.standardizeSmiles(ent, type="ISOMERIC"))
        self.assertEqual(can, self.__oedu.standardizeSmiles(ent, type="CANNONICAL"))

    def testGtpNormalization(self) -> None:
        """Different GTP inputs normalize to the same output"""
        canList = [
            self.__oedu.standardizeSmiles(s, type="CANNONICAL") for s in self.__gtpStereoList + self.__gtpFlatList
        ]
        self.assertTrue(all(c is not None for c in canList))
        self.assertEqual(len(set(canList)), 1)
        canonical = canList[0]
        assert canonical is not None
        self.assertNotIn("@", canonical)

        isoStereo = [self.__oedu.standardizeSmiles(s, type="ISOMERIC") for s in self.__gtpStereoList]
        self.assertEqual(len(set(isoStereo)), 1)
        isoS = isoStereo[0]
        assert isoS is not None
        self.assertIn("@", isoS)
        self.assertNotEqual(isoS, canonical)

        isoFlat = [self.__oedu.standardizeSmiles(s, type="ISOMERIC") for s in self.__gtpFlatList]
        # Without stereo input, isomeric equals canonical
        self.assertEqual(set(isoFlat), {canonical})

    def testIdempotent(self) -> None:
        """Standardizing an already standardized SMILES gives the same result"""
        for typ in ("ISOMERIC", "CANNONICAL"):
            for smi in self.__gtpStereoList:
                once = self.__oedu.standardizeSmiles(smi, type=typ)
                assert once is not None
                self.assertEqual(self.__oedu.standardizeSmiles(once, type=typ), once)

    def testUnknownType(self) -> None:
        """Unrecognized type returns None"""
        self.assertIsNone(self.__oedu.standardizeSmiles("CC", type="CANONICAL"))
        self.assertIsNone(self.__oedu.standardizeSmiles("CC", type="bogus"))

    def testInvalidSmiles(self) -> None:
        """Unparsable SMILES returns None and logs an error"""
        with self.assertLogs(_LOGNAME, level="ERROR") as cm:
            self.assertIsNone(self.__oedu.standardizeSmiles("C1CC(((", type="ISOMERIC"))
        self.assertTrue(any("Unable to parse input SMILES" in m for m in cm.output))

    def testNoneSmiles(self) -> None:
        """None input triggers exception handling and returns None"""
        with self.assertLogs(_LOGNAME, level="ERROR") as cm:
            self.assertIsNone(self.__oedu.standardizeSmiles(None))  # type: ignore[arg-type]
        self.assertTrue(any("occured" in m for m in cm.output))


def suiteExtra() -> unittest.TestSuite:
    return unittest.TestLoader().loadTestsFromTestCase(OeDescriptorUtilsExtraTests)


if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suiteExtra())
