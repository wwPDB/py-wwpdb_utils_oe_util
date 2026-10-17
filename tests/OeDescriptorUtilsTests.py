##
# File:    OeDescriptorUtilsTests.py
# Author:  J. Westbrook
# Date:    21-Mar-2017
# Version: 0.001
#
# Update:
#
##
"""
Tests for standardizing descriptors -

"""

__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.01"


import logging
import time
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

GTP_ISO_OE = "c1nc2c(=O)[nH]c(nc2n1[C@H]3[C@@H]([C@@H]([C@H](O3)CO[P@](=O)(O)O[P@](=O)(O)OP(=O)(O)O)O)O)N"
GTP_CAN_OE = "c1nc2c(=O)[nH]c(nc2n1C3C(C(C(O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N"


@unittest.skipIf(skiptests, "Requires openeye library")
class OeDescriptorUtilsTests(unittest.TestCase):
    def setUp(self) -> None:
        self.__gtpSmilesList: List[str] = [
            "O=P(O)(O)OP(=O)(O)OP(=O)(O)OCC3OC(n2cnc1c2N=C(N)NC1=O)C(O)C3O",
            "NC1=Nc2n(cnc2C(=O)N1)[C@@H]3O[C@H](CO[P@@](O)(=O)O[P@@](O)(=O)O[P](O)(O)=O)[C@@H](O)[C@H]3O",
            "NC1=Nc2n(cnc2C(=O)N1)[CH]3O[CH](CO[P](O)(=O)O[P](O)(=O)O[P](O)(O)=O)[CH](O)[CH]3O",
            "c1nc2c(n1[C@H]3[C@@H]([C@@H]([C@H](O3)CO[P@](=O)(O)O[P@](=O)(O)OP(=O)(O)O)O)O)N=C(NC2=O)N",
            "c1nc2c(n1C3C(C(C(O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N=C(NC2=O)N",
        ]
        # Inputs at these indices carry stereo; the others do not
        self.__stereoIndices: List[int] = [1, 3]

    def tearDown(self) -> None:
        pass

    def testSmilesConvert(self) -> None:
        """Test case - SMILES Conversion-"""
        startTime = time.time()
        oedu = OeDescriptorUtils()  # pylint: disable=possibly-used-before-assignment

        for ii, smiles in enumerate(self.__gtpSmilesList):
            smiIso = oedu.standardizeSmiles(smiles, type="ISOMERIC")
            smiCan = oedu.standardizeSmiles(smiles, type="CANNONICAL")
            logger.debug("SMILES (ISO) %s", smiIso)
            logger.debug("SMILES (CAN) %s", smiCan)
            self.assertEqual(smiCan, GTP_CAN_OE)
            if ii in self.__stereoIndices:
                self.assertEqual(smiIso, GTP_ISO_OE)
            else:
                self.assertEqual(smiIso, GTP_CAN_OE)

        logger.debug("Completed (%.2f seconds)", time.time() - startTime)


def testSmilesConversionSuite() -> unittest.TestSuite:
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeDescriptorUtilsTests("testSmilesConvert"))
    return suiteSelect


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    mySuite = testSmilesConversionSuite()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
