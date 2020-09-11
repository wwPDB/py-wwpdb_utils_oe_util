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

import sys
import unittest
import traceback
import platform
import os

try:
    from openeye.oechem import OEFloatArray  # noqa: F401 pylint: disable=unused-import

    skiptests = False
except ImportError:
    skiptests = True

if not skiptests:
    from wwpdb.utils.oe_util.build.OeBuildMol import OeBuildMol
    from mmcif_utils.persist.PdbxPyIoAdapter import PdbxPyIoAdapter as PdbxIoAdapter


@unittest.skipIf(skiptests, "Cannot import openeye.oechem for tests")
class OeBuildMolTests(unittest.TestCase):
    def setUp(self):
        self.__lfh = sys.stderr
        self.__verbose = True
        self.__here = os.path.abspath(os.path.dirname(__file__))
        self.__testoutput = os.path.join(self.__here, "test-output", platform.python_version())
        if not os.path.exists(self.__testoutput):
            os.makedirs(self.__testoutput)
        self.__datadir = os.path.join(self.__here, "data")
        self.__sdfFilePath = os.path.join(self.__datadir, "ATP.sdf")
        self.__pathList = [os.path.join(self.__datadir, "ATP.cif"), os.path.join(self.__datadir, "GTP.cif"), os.path.join(self.__datadir, "ARG.cif")]

    def tearDown(self):
        pass

    def testBuildMolFromSDF(self):
        """Test case -  read a test SDF file and build the corresponding OEGraphMol"""
        self.__lfh.write("\nStarting OeBuildMolTests testBuildolFromSDF\n")
        try:
            oem = OeBuildMol(verbose=self.__verbose, log=self.__lfh)
            if oem.importFile(self.__sdfFilePath, type="3D"):
                self.__lfh.write("Title              = %s\n" % oem.getTitle())
                self.__lfh.write("SMILES (canonical) = %s\n" % oem.getCanSMILES())
                self.__lfh.write("SMILES (isomeric)  = %s\n" % oem.getIsoSMILES())
            else:
                self.__lfh.write("SDF read failed for %s\n" % self.__sdfFilePath)
        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testBuildFromFiles(self):
        """Test case -  build OE molecule from definition file source data."""
        self.__lfh.write("\nStarting OeBuildMolTests testBuildFromFiles\n")
        try:
            oem = OeBuildMol(verbose=self.__verbose, log=self.__lfh)
            for pth in self.__pathList:
                myReader = PdbxIoAdapter(self.__verbose, self.__lfh)
                ok = myReader.read(pdbxFilePath=pth)
                self.assertTrue(ok)
                # myReader.write(pdbxFilePath="TMP.cif")
                for container in myReader.getContainerList():
                    oem.set(container.getName(), dcChemCompAtom=container.getObj("chem_comp_atom"), dcChemCompBond=container.getObj("chem_comp_bond"))
                    oem.build3D(coordType="model")
                    self.__lfh.write("Title              = %s\n" % oem.getTitle())
                    self.__lfh.write("SMILES (canonical) = %s\n" % oem.getCanSMILES())
                    self.__lfh.write("SMILES (isomeric)  = %s\n" % oem.getIsoSMILES())

        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testSerialize2D(self):
        """Test case -  build OE molecule using 2D data in the definition source data file
        then serialize and deserialize this molecule.
        """
        self.__lfh.write("\nStarting OeBuildMolTests testSerialize2D\n")
        try:
            oem = OeBuildMol(verbose=self.__verbose, log=self.__lfh)
            oem.setDebug(True)
            for pth in self.__pathList:
                myReader = PdbxIoAdapter(self.__verbose, self.__lfh)
                ok = myReader.read(pdbxFilePath=pth)
                myReader.write(pdbxFilePath=os.path.join(self.__testoutput, "TMP.cif"))
                for container in myReader.getContainerList():
                    oem.set(container.getName(), dcChemCompAtom=container.getObj("chem_comp_atom"), dcChemCompBond=container.getObj("chem_comp_bond"))
                    oem.build2D()
                    self.__lfh.write("Title              = %s\n" % oem.getTitle())
                    self.__lfh.write("SMILES (canonical) = %s\n" % oem.getCanSMILES())
                    self.__lfh.write("SMILES (isomeric)  = %s\n" % oem.getIsoSMILES())
                    oeS = oem.serialize()
                    #
                    self.__lfh.write("Serialized string length = %d\n" % len(oeS))
                    #
                    oemD = OeBuildMol(verbose=self.__verbose, log=self.__lfh)
                    ok = oemD.deserialize(oeS)
                    self.__lfh.write("Deserialized status = %d\n" % ok)
                    self.__lfh.write("Deserialized SMILES (canonical) = %s\n" % oemD.getCanSMILES())
                    self.__lfh.write("Deserialized SMILES (isomeric)  = %s\n" % oemD.getIsoSMILES())

        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testSerialize3D(self):
        """Test case -  build OE molecule using 3D data in the definition source data file
        then serialize and deserialize this molecule.

        """
        self.__lfh.write("\nStarting OeBuildMolTests testSerialize3D\n")
        try:
            oem = OeBuildMol(verbose=self.__verbose, log=self.__lfh)
            for pth in self.__pathList:
                myReader = PdbxIoAdapter(self.__verbose, self.__lfh)
                ok = myReader.read(pdbxFilePath=pth)
                myReader.write(pdbxFilePath=os.path.join(self.__testoutput, "TMP.cif"))
                for container in myReader.getContainerList():
                    oem.set(container.getName(), dcChemCompAtom=container.getObj("chem_comp_atom"), dcChemCompBond=container.getObj("chem_comp_bond"))
                    oem.build3D(coordType="model")
                    self.__lfh.write("Title              = %s\n" % oem.getTitle())
                    self.__lfh.write("SMILES (canonical) = %s\n" % oem.getCanSMILES())
                    self.__lfh.write("SMILES (isomeric)  = %s\n" % oem.getIsoSMILES())
                    oeS = oem.serialize()
                    #
                    self.__lfh.write("Serialized string length = %d\n" % len(oeS))
                    #
                    oemD = OeBuildMol(verbose=self.__verbose, log=self.__lfh)
                    ok = oemD.deserialize(oeS)
                    self.__lfh.write("Deserialized status = %d\n" % ok)
                    self.__lfh.write("Deserialized SMILES (canonical) = %s\n" % oemD.getCanSMILES())
                    self.__lfh.write("Deserialized SMILES (isomeric)  = %s\n" % oemD.getIsoSMILES())

        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()


def suiteOeBuild():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeBuildMolTests("testBuildFromFiles"))
    suiteSelect.addTest(OeBuildMolTests("testSerialize2D"))
    suiteSelect.addTest(OeBuildMolTests("testSerialize3D"))
    return suiteSelect


def suiteOeBuildFromOther():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeBuildMolTests("testBuildMolFromSDF"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = suiteOeBuild()
    unittest.TextTestRunner(verbosity=2).run(mySuite)

    mySuite = suiteOeBuildFromOther()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
