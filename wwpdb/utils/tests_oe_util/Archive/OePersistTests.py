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
import traceback
import time

from wwpdb.utils.oe_util.build.OePersist import OePersist
from wwpdb.utils.oe_util.build.OeBuildMol import OeBuildMol
from mmcif_utils.persist.PdbxPyIoAdapter import PdbxPyIoAdapter as PdbxIoAdapter


class OePersistTests(unittest.TestCase):

    def setUp(self):
        self.__lfh = sys.stdout
        self.__verbose = True
        self.__pathList = ['../data/ATP.cif', '../data/GTP.cif', '../data/ARG.cif']
        self.__storePath = 'oe-store.db'

    def tearDown(self):
        pass

    def testCreateStore(self):
        """Test case -  build store from a selection of serialized OE molecules.
        """
        self.__lfh.write("\nStarting %s %s\n" % (self.__class__.__name__,
                                                 sys._getframe().f_code.co_name))
        try:
            molList = []
            myPersist = OePersist(self.__verbose, self.__lfh)
            oem = OeBuildMol(verbose=self.__verbose, log=self.__lfh)
            for pth in self.__pathList:
                myReader = PdbxIoAdapter(self.__verbose, self.__lfh)
                ok = myReader.read(pdbxFilePath=pth)
                for container in myReader.getContainerList():
                    name = container.getName()
                    oem.set(name,
                            dcChemCompAtom=container.getObj("chem_comp_atom"),
                            dcChemCompBond=container.getObj("chem_comp_bond"))
                    oem.build2D()
                    self.__lfh.write("Title              = %s\n" % oem.getTitle())
                    self.__lfh.write("SMILES (canonical) = %s\n" % oem.getCanSMILES())
                    self.__lfh.write("SMILES (isomeric)  = %s\n" % oem.getIsoSMILES())
                    molD = {}
                    molD['name'] = name
                    molD['oeb'] = oem.serialize()
                    molList.append(molD)

            myPersist.setMoleculeList(moleculeList=molList)
            myPersist.store(dbFileName=self.__storePath)
            mL = myPersist.getIndex(dbFileName=self.__storePath)
            self.__lfh.write("OePersistTests(testCreateStore) molecule list %r\n" % mL)

        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testFetchAll(self):
        """Test case -  fetch all of the molecules in the persistent store
        """
        startTime = time.time()
        self.__lfh.write("\nStarting %s %s at %s\n" % (self.__class__.__name__,
                                                       sys._getframe().f_code.co_name,
                                                       time.strftime("%Y %m %d %H:%M:%S", time.localtime())))
        try:
            myPersist = OePersist(self.__verbose, self.__lfh)
            myPersist.open(dbFileName=self.__storePath)
            moleculeNameList = myPersist.getStoreMoleculeIndex()
            #
            oem = OeBuildMol(verbose=self.__verbose, log=self.__lfh)
            for ccId in moleculeNameList:
                molD = myPersist.fetchMolecule(moleculeName=ccId)
                name = molD['name']
                ok = oem.deserialize(molD['oeb'])
                self.__lfh.write("Deserialized status = %d\n" % ok)
                self.__lfh.write("Deserialized SMILES (canonical) = %s\n" % oem.getCanSMILES())
                self.__lfh.write("Deserialized SMILES (isomeric)  = %s\n" % oem.getIsoSMILES())

            myPersist.close()

        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()

        endTime = time.time()
        self.__lfh.write("\nCompleted %s %s at %s (%d seconds)\n" % (self.__class__.__name__,
                                                                     sys._getframe().f_code.co_name,
                                                                     time.strftime("%Y %m %d %H:%M:%S", time.localtime()),
                                                                     endTime - startTime))

    def testFetchOne(self):
        """Test case -  fetch all of the molecules in the persistent store one by one
        """
        startTime = time.time()
        self.__lfh.write("\nStarting %s %s at %s\n" % (self.__class__.__name__,
                                                       sys._getframe().f_code.co_name,
                                                       time.strftime("%Y %m %d %H:%M:%S", time.localtime())))
        try:
            myPersist = OePersist(self.__verbose, self.__lfh)
            moleculeNameList = myPersist.getIndex(dbFileName=self.__storePath)
            #
            oem = OeBuildMol(verbose=self.__verbose, log=self.__lfh)
            for ccId in moleculeNameList:
                molD = myPersist.fetchOneMolecule(dbFileName=self.__storePath, moleculeName=ccId)
                name = molD['name']
                ok = oem.deserialize(molD['oeb'])
                self.__lfh.write("Deserialized status = %d\n" % ok)
                self.__lfh.write("Deserialized SMILES (canonical) = %s\n" % oem.getCanSMILES())
                self.__lfh.write("Deserialized SMILES (isomeric)  = %s\n" % oem.getIsoSMILES())

            myPersist.close()

        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()

        endTime = time.time()
        self.__lfh.write("\nCompleted %s %s at %s (%d seconds)\n" % (self.__class__.__name__,
                                                                     sys._getframe().f_code.co_name,
                                                                     time.strftime("%Y %m %d %H:%M:%S", time.localtime()),
                                                                     endTime - startTime))

    def testUpdateStore(self):
        """Test case -  update store from a selection of serialized OE molecules.
        """
        self.__lfh.write("\nStarting %s %s\n" % (self.__class__.__name__,
                                                 sys._getframe().f_code.co_name))
        try:
            molList = []
            myPersist = OePersist(self.__verbose, self.__lfh)
            oem = OeBuildMol(verbose=self.__verbose, log=self.__lfh)
            for pth in self.__pathList:
                myReader = PdbxIoAdapter(self.__verbose, self.__lfh)
                ok = myReader.read(pdbxFilePath=pth)
                for container in myReader.getContainerList():
                    name = container.getName()
                    oem.set(name,
                            dcChemCompAtom=container.getObj("chem_comp_atom"),
                            dcChemCompBond=container.getObj("chem_comp_bond"))
                    oem.build2D()
                    self.__lfh.write("Title              = %s\n" % oem.getTitle())
                    self.__lfh.write("SMILES (canonical) = %s\n" % oem.getCanSMILES())
                    self.__lfh.write("SMILES (isomeric)  = %s\n" % oem.getIsoSMILES())
                    molD = {}
                    molD['name'] = name
                    molD['oeb'] = oem.serialize()
                    myPersist.updateOneMolecule(molD, dbFileName=self.__storePath)

            mL = myPersist.getIndex(dbFileName=self.__storePath)
            self.__lfh.write("OePersistTests(testCreateStore) molecule list %r\n" % mL)

        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()


def suite1():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OePersistTests("testCreateStore"))
    suiteSelect.addTest(OePersistTests("testUpdateStore"))
    suiteSelect.addTest(OePersistTests("testFetchAll"))
    suiteSelect.addTest(OePersistTests("testFetchOne"))

    return suiteSelect

if __name__ == '__main__':
    #
    mySuite = suite1()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
