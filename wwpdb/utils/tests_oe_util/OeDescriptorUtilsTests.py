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

import time
import unittest

try:
    from openeye.oechem import OEFloatArray  # noqa: F401 pylint: disable=unused-import

    skiptests = False
except ImportError:
    skiptests = True

if not skiptests:
    from wwpdb.utils.oe_util.build.OeDescriptorUtils import OeDescriptorUtils

import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()


@unittest.skipIf(skiptests, "Requires openeye library")
class OeDescriptorUtilsTests(unittest.TestCase):
    def setUp(self):
        #
        self.__gtpSmilesList = [
            "O=P(O)(O)OP(=O)(O)OP(=O)(O)OCC3OC(n2cnc1c2N=C(N)NC1=O)C(O)C3O",
            "NC1=Nc2n(cnc2C(=O)N1)[C@@H]3O[C@H](CO[P@@](O)(=O)O[P@@](O)(=O)O[P](O)(O)=O)[C@@H](O)[C@H]3O",
            "NC1=Nc2n(cnc2C(=O)N1)[CH]3O[CH](CO[P](O)(=O)O[P](O)(=O)O[P](O)(O)=O)[CH](O)[CH]3O",
            "c1nc2c(n1[C@H]3[C@@H]([C@@H]([C@H](O3)CO[P@](=O)(O)O[P@](=O)(O)OP(=O)(O)O)O)O)N=C(NC2=O)N",
            "c1nc2c(n1C3C(C(C(O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N=C(NC2=O)N",
        ]
        self.__smilesGtpIsoOE = "c1nc2c(n1[C@H]3[C@@H]([C@@H]([C@H](O3)CO[P@](=O)(O)O[P@](=O)(O)OP(=O)(O)O)O)O)N=C(NC2=O)N"
        self.__smilesGtpCanOE = "c1nc2c(n1C3C(C(C(O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N=C(NC2=O)N"

    def tearDown(self):
        pass

    def testSmilesConvert(self):
        """Test case - SMILES Conversion-"""
        logger.info("Starting")
        startTime = time.time()
        try:
            oedu = OeDescriptorUtils()
            for smiles in self.__gtpSmilesList:
                smiIso = oedu.standardizeSmiles(smiles, type="ISOMERIC")
                smiCan = oedu.standardizeSmiles(smiles, type="CANNONICAL")
                logger.info("SMILES (ISO) %s", smiIso)
                logger.info("SMILES (CAN) %s", smiCan)

        except Exception as e:
            logger.exception("Error '%s' occured. Arguments %s.", str(e), e.args)
            self.fail()
        #
        endTime = time.time()
        logger.info("Completed at %s (%.2f seconds)", time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime)


def testSmilesConversionSuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeDescriptorUtilsTests("testSmilesConvert"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = testSmilesConversionSuite()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
