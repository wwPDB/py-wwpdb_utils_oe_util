##
#
# File:    OePersistTests.py
# Author:  J. Westbrook
# Date:    23-Jan-2012
# Version: 0.001
#
# Updated:
#      6-Jun-2016  jdw general cleanup
#
##
"""
Test cases for persistent storage of serialized OE molecule objects.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.01"

import sys
import unittest
import platform
import os
import traceback
import time

try:
    from openeye.oechem import OEFloatArray  # noqa: F401 pylint: disable=unused-import

    skiptests = False
except ImportError:
    skiptests = True

if not skiptests:
    from wwpdb.utils.oe_util.build.OePersist import OePersist
    from wwpdb.utils.oe_util.build.OeBuildMol import OeBuildMol
    from mmcif_utils.persist.PdbxPyIoAdapter import PdbxPyIoAdapter as PdbxIoAdapter


@unittest.skipIf(skiptests, "Cannot import openeye.oechem for tests")
class OePersistTests(unittest.TestCase):
    def setUp(self):
        self.__lfh = sys.stdout
        self.__verbose = True
        self.__here = os.path.abspath(os.path.dirname(__file__))
        self.__examples = os.path.join(self.__here, "examples")
        self.__datadir = os.path.join(self.__here, "data")
        self.__testoutput = os.path.join(self.__here, "test-output", platform.python_version())
        if not os.path.exists(self.__testoutput):
            os.makedirs(self.__testoutput)
        self.__pathList = [os.path.join(self.__datadir, "ATP.cif"), os.path.join(self.__datadir, "GTP.cif"), os.path.join(self.__datadir, "ARG.cif")]
        self.__storePath = os.path.join(self.__testoutput, "oe-store.db")

    def tearDown(self):
        pass

    def testCreateStore(self):
        """Test case -  build store from a selection of serialized OE molecules."""
        self.__lfh.write("\nStarting OePersistTests testCreateStore\n")
        try:
            molList = []
            myPersist = OePersist(self.__verbose, self.__lfh)
            oem = OeBuildMol(verbose=self.__verbose, log=self.__lfh)
            for pth in self.__pathList:
                myReader = PdbxIoAdapter(self.__verbose, self.__lfh)
                ok = myReader.read(pdbxFilePath=pth)
                self.assertTrue(ok)
                for container in myReader.getContainerList():
                    name = container.getName()
                    oem.set(name, dcChemCompAtom=container.getObj("chem_comp_atom"), dcChemCompBond=container.getObj("chem_comp_bond"))
                    oem.build2D()
                    self.__lfh.write("Title              = %s\n" % oem.getTitle())
                    self.__lfh.write("SMILES (canonical) = %s\n" % oem.getCanSMILES())
                    self.__lfh.write("SMILES (isomeric)  = %s\n" % oem.getIsoSMILES())
                    molD = {}
                    molD["name"] = name
                    molD["oeb"] = oem.serialize()
                    molList.append(molD)

            myPersist.setMoleculeList(moleculeList=molList)
            myPersist.store(dbFileName=self.__storePath)
            mL = myPersist.getIndex(dbFileName=self.__storePath)
            self.__lfh.write("OePersistTests(testCreateStore) molecule list %r\n" % mL)

        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testFetchAll(self):
        """Test case -  fetch all of the molecules in the persistent store"""
        startTime = time.time()
        self.__lfh.write("\nStarting OePersistTests testFetchAll at %s\n" % time.strftime("%Y %m %d %H:%M:%S", time.localtime()))
        try:
            myPersist = OePersist(self.__verbose, self.__lfh)
            myPersist.open(dbFileName=self.__storePath)
            moleculeNameList = myPersist.getStoreMoleculeIndex()
            #
            oem = OeBuildMol(verbose=self.__verbose, log=self.__lfh)
            for ccId in moleculeNameList:
                molD = myPersist.fetchMolecule(moleculeName=ccId)
                # name = molD['name']
                ok = oem.deserialize(molD["oeb"])
                self.__lfh.write("Deserialized status = %d\n" % ok)
                self.__lfh.write("Deserialized SMILES (canonical) = %s\n" % oem.getCanSMILES())
                self.__lfh.write("Deserialized SMILES (isomeric)  = %s\n" % oem.getIsoSMILES())

            myPersist.close()

        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()

        endTime = time.time()
        self.__lfh.write("\nCompleted OePersistTests testFetchAll at %s (%d seconds)\n" % (time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime))

    def testFetchOne(self):
        """Test case -  fetch all of the molecules in the persistent store one by one"""
        startTime = time.time()
        self.__lfh.write("\nStarting OePersistTests testFetchOne at %s\n" % time.strftime("%Y %m %d %H:%M:%S", time.localtime()))
        try:
            myPersist = OePersist(self.__verbose, self.__lfh)
            moleculeNameList = myPersist.getIndex(dbFileName=self.__storePath)
            #
            oem = OeBuildMol(verbose=self.__verbose, log=self.__lfh)
            for ccId in moleculeNameList:
                molD = myPersist.fetchOneMolecule(dbFileName=self.__storePath, moleculeName=ccId)
                # name = molD['name']
                ok = oem.deserialize(molD["oeb"])
                self.__lfh.write("Deserialized status = %d\n" % ok)
                self.__lfh.write("Deserialized SMILES (canonical) = %s\n" % oem.getCanSMILES())
                self.__lfh.write("Deserialized SMILES (isomeric)  = %s\n" % oem.getIsoSMILES())

            myPersist.close()

        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()

        endTime = time.time()
        self.__lfh.write("\nCompleted OePersistTests testFetchOne at %s (%d seconds)\n" % (time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime))

    def testUpdateStore(self):
        """Test case -  update store from a selection of serialized OE molecules."""
        self.__lfh.write("\nStarting OePersistTests testUpdateStore\n")
        try:
            myPersist = OePersist(self.__verbose, self.__lfh)
            oem = OeBuildMol(verbose=self.__verbose, log=self.__lfh)
            for pth in self.__pathList:
                myReader = PdbxIoAdapter(self.__verbose, self.__lfh)
                ok = myReader.read(pdbxFilePath=pth)
                self.assertTrue(ok)
                for container in myReader.getContainerList():
                    name = container.getName()
                    oem.set(name, dcChemCompAtom=container.getObj("chem_comp_atom"), dcChemCompBond=container.getObj("chem_comp_bond"))
                    oem.build2D()
                    self.__lfh.write("Title              = %s\n" % oem.getTitle())
                    self.__lfh.write("SMILES (canonical) = %s\n" % oem.getCanSMILES())
                    self.__lfh.write("SMILES (isomeric)  = %s\n" % oem.getIsoSMILES())
                    molD = {}
                    molD["name"] = name
                    molD["oeb"] = oem.serialize()
                    myPersist.updateOneMolecule(molD, dbFileName=self.__storePath)

            mL = myPersist.getIndex(dbFileName=self.__storePath)
            self.__lfh.write("OePersistTests(testCreateStore) molecule list %r\n" % mL)

        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()


def suite1():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OePersistTests("testCreateStore"))
    suiteSelect.addTest(OePersistTests("testUpdateStore"))
    suiteSelect.addTest(OePersistTests("testFetchAll"))
    suiteSelect.addTest(OePersistTests("testFetchOne"))

    return suiteSelect


if __name__ == "__main__":
    #
    mySuite = suite1()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
