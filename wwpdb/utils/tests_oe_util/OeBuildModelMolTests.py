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

import os
import sys
import unittest
import traceback


try:
    from wwpdb.utils.oe_util.build.OeBuildModelMol import OeBuildModelMol
    skiptests = False
except ImportError as e:
    skiptests = True
    
from mmcif_utils.chemcomp.PdbxChemCompModelIo import PdbxChemCompModelIo
from mmcif_utils.chemcomp.PdbxChemCompModel import PdbxChemCompModelDescriptor

@unittest.skipIf(skiptests, "Requires oe library")
class OeBuildModelMolTests(unittest.TestCase):

    def setUp(self):
        self.__lfh = sys.stderr
        self.__verbose = True
        HERE = os.path.abspath(os.path.dirname(__file__))
        TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
        mockTopPath = os.path.join(TOPDIR, 'wwpdb', 'mock-data')
        self.__modelFilePath = os.path.join(mockTopPath, 'CCD', 'MTGL00001.cif')
        self.__modelFilePathList = [self.__modelFilePath]

    def tearDown(self):
        pass

    def testBuildFromModel(self):
        """Test case -  build OE molecule from model instance
        """
        self.__lfh.write("\nStarting %s %s\n" % (self.__class__.__name__, sys._getframe().f_code.co_name))
        try:
            oem = OeBuildModelMol(verbose=self.__verbose, log=self.__lfh)
            modelId = oem.setChemCompModelPath(self.__modelFilePath)
            self.__lfh.write("Model              = %s\n" % modelId)
            oem.build3D()
            self.__lfh.write("Title              = %s\n" % oem.getTitle())
            self.__lfh.write("SMILES (canonical) = %s\n" % oem.getCanSMILES())
            self.__lfh.write("SMILES (isomeric)  = %s\n" % oem.getIsoSMILES())
            #
            pccm = PdbxChemCompModelIo(verbose=self.__verbose, log=self.__lfh)
            pccm.setFilePath(self.__modelFilePath)
            dL = pccm.getDescriptorList()
            for d in dL:
                pd = PdbxChemCompModelDescriptor(d, verbose=self.__verbose, log=self.__lfh)
                print(pd.getType())
                if pd.getType() == 'SMILES_CANNONICAL':
                    sm = pd.getDescriptor()
                    if (sm == oem.getIsoSMILES()):
                        self.__lfh.write("+testBuildFromModel. SMILES MATCH for %s\n" % modelId)
                    else:
                        self.__lfh.write("+testBuildFromModel. SMILES MISMATCH for %s\n" % modelId)
                        self.__lfh.write("+testBuildFromModel. SMILES (model)               %s\n" % sm)
                        self.__lfh.write("+testBuildFromModel. SMILES (reconstructed model) %s\n" % oem.getIsoSMILES())
        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testSerialize3D(self):
        """Test case -  build OE molecule using 3D data in the chemical component model instance and
           then serialize and deserialize this molecule.

        """
        self.__lfh.write("\nStarting %s %s\n" % (self.__class__.__name__, sys._getframe().f_code.co_name))
        try:
            oem = OeBuildModelMol(verbose=self.__verbose, log=self.__lfh)
            for pth in self.__modelFilePathList:
                modelId = oem.setChemCompModelPath(pth)
                self.__lfh.write("Model              = %s\n" % modelId)
                oem.build3D()
                self.__lfh.write("Title              = %s\n" % oem.getTitle())
                self.__lfh.write("SMILES (canonical) = %s\n" % oem.getCanSMILES())
                self.__lfh.write("SMILES (isomeric)  = %s\n" % oem.getIsoSMILES())
                oeS = oem.serialize()
                #
                self.__lfh.write("Serialized string length = %d\n" % len(oeS))
                #
                oemD = OeBuildModelMol(verbose=self.__verbose, log=self.__lfh)
                ok = oemD.deserialize(oeS)
                self.__lfh.write("Deserialized status = %d\n" % ok)
                self.__lfh.write("Deserialized SMILES (canonical) = %s\n" % oemD.getCanSMILES())
                self.__lfh.write("Deserialized SMILES (isomeric)  = %s\n" % oemD.getIsoSMILES())
        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()


def suiteOeModelBuild():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeBuildModelMolTests("testBuildFromModel"))
    suiteSelect.addTest(OeBuildModelMolTests("testSerialize3D"))
    return suiteSelect

if __name__ == '__main__':
    if (True):
        mySuite = suiteOeModelBuild()
        unittest.TextTestRunner(verbosity=2).run(mySuite)
