##
# File:    simple-depict-example.py
# Author:  jdw
# Date:    7-Sep-2013
# Version: 0.001
#
# Updates:
#      6-Jun-2016 jdw general cleanup -
##
"""
An example of generating 2D chemical diagrams in svg format from chemical component definition.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.01"


import sys
import traceback
import unittest

from wwpdb.utils.oe_util.build.OeChemCompIoUtils import OeChemCompIoUtils
from wwpdb.utils.oe_util.oedepict.OeDepict import OeDepict


class OeDepictTests(unittest.TestCase):
    def setUp(self):
        self.__lfh = sys.stderr
        self.__verbose = True
        # ---------------------------------------------------------------------------------
        # The following needs to be set to the path of a checked out copy of the CVS ligand repository.
        ###
        self.__topCachePath = "/data/components/ligand-dict-v3"
        ##
        # example list of 3-letter-codes -
        ##
        self.__idList = ["atp", "gtp", "A", "C", "G", "DG"]

    def tearDown(self):
        pass

    def testDepictIdList(self):
        """Test case -  get, read, build OE molecule, and depict the molecule."""
        self.__lfh.write("\nStarting OeDepictTests testDepictIDList\n")
        try:
            oeMolTitleList = self.__makeFromIdList(idList=self.__idList)
            if self.__verbose:
                self.__lfh.write("molTitleList length is %d\n" % len(oeMolTitleList))
            #
            for ccId, mol, title in oeMolTitleList:
                imagePath = ccId + ".svg"
                oed = OeDepict(verbose=self.__verbose, log=self.__lfh)
                oed.setMolTitleList([(ccId, mol, title)])
                oed.setDisplayOptions(imageX=250, imageY=250, labelAtomName=True, labelAtomCIPStereo=True, labelAtomIndex=False, labelBondIndex=False, bondDisplayWidth=0.5)
                oed.setGridOptions(rows=1, cols=1)
                oed.prepare()
                oed.write(imagePath)
        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def __makeFromIdList(self, idList):
        """Create OE molecules from the input chemical component definition id list."""
        self.__lfh.write("\nStarting OeDepictTests __makeFromIdList\n")
        try:
            oemList = []
            oeU = OeChemCompIoUtils(topCachePath=self.__topCachePath, verbose=self.__verbose, log=self.__lfh)
            oemList = oeU.getFromIdList(idList, use3D=True, coordType="model")
            for oem in oemList:
                self.__lfh.write("Title              = %s\n" % oem.getTitle())
                self.__lfh.write("SMILES (canonical) = %s\n" % oem.getCanSMILES())
                self.__lfh.write("SMILES (isomeric)  = %s\n" % oem.getIsoSMILES())
            return zip(idList, oemList, idList)
        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()


def suiteDepict():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeDepictTests("testDepictIdList"))
    return suiteSelect


if __name__ == "__main__":
    mySuite2 = suiteDepict()
    unittest.TextTestRunner(verbosity=2).run(mySuite2)
    #
