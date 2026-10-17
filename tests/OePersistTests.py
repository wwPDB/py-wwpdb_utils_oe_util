# pylint: disable=possibly-used-before-assignment
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

import os
import sys
import tempfile
import time
import unittest
from typing import Any, Dict, List

try:
    from openeye.oechem import OEFloatArray  # noqa: F401 pylint: disable=unused-import

    skiptests = False
except ImportError:
    skiptests = True

if not skiptests:
    from mmcif_utils.persist.PdbxPyIoAdapter import PdbxPyIoAdapter as PdbxIoAdapter

    from wwpdb.utils.oe_util.build.OeBuildMol import OeBuildMol
    from wwpdb.utils.oe_util.build.OePersist import OePersist


@unittest.skipIf(skiptests, "Cannot import openeye.oechem for tests")
class OePersistTests(unittest.TestCase):
    def setUp(self) -> None:
        self.__lfh = sys.stdout
        self.__verbose = True
        self.__here = os.path.abspath(os.path.dirname(__file__))
        self.__datadir = os.path.join(self.__here, "data")
        self.__tmpDir = tempfile.TemporaryDirectory()
        self.__testoutput = self.__tmpDir.name
        self.__pathList: List[str] = [
            os.path.join(self.__datadir, "ATP.cif"),
            os.path.join(self.__datadir, "GTP.cif"),
            os.path.join(self.__datadir, "ARG.cif"),
        ]
        self.__expectedNames: List[str] = ["ATP", "GTP", "ARG"]
        self.__storePath = os.path.join(self.__testoutput, "oe-store.db")
        # name -> isomeric SMILES computed before serialization
        self.__smilesD: Dict[str, str] = {}
        self.__molList: List[Dict[str, Any]] = self.__buildMoleculeList()

    def tearDown(self) -> None:
        self.__tmpDir.cleanup()

    def __buildMoleculeList(self) -> List[Dict[str, Any]]:
        """Build serialized OE molecules from the test chemical component definitions."""
        molList: List[Dict[str, Any]] = []
        oem = OeBuildMol(verbose=self.__verbose, log=self.__lfh)
        for pth in self.__pathList:
            myReader = PdbxIoAdapter(self.__verbose, self.__lfh)
            ok = myReader.read(pdbxFilePath=pth)
            self.assertTrue(ok, "read failed for %s" % pth)
            for container in myReader.getContainerList():
                name = container.getName()
                oem.set(
                    name,
                    dcChemCompAtom=container.getObj("chem_comp_atom"),
                    dcChemCompBond=container.getObj("chem_comp_bond"),
                )
                self.assertTrue(oem.build2D(), "build2D failed for %s" % name)
                self.assertEqual(oem.getTitle(), name)
                self.__lfh.write("Title              = %s\n" % oem.getTitle())
                self.__lfh.write("SMILES (canonical) = %s\n" % oem.getCanSMILES())
                self.__lfh.write("SMILES (isomeric)  = %s\n" % oem.getIsoSMILES())
                self.__smilesD[name] = oem.getIsoSMILES()
                oeb = oem.serialize()
                self.assertTrue(len(oeb) > 0)
                molList.append({"name": name, "oeb": oeb})
        self.assertEqual([molD["name"] for molD in molList], self.__expectedNames)
        return molList

    def __createStore(self, molList: List[Dict[str, Any]]) -> "OePersist":
        """Write the input molecule list to the test store and return the persist object used."""
        myPersist = OePersist(self.__verbose, self.__lfh)
        self.assertTrue(myPersist.setMoleculeList(moleculeList=molList))
        self.assertTrue(myPersist.store(dbFileName=self.__storePath))
        return myPersist

    def __assertRoundTrip(self, oem: "OeBuildMol", molD: Any, name: str) -> None:
        """Check that a fetched molecule deserializes to the molecule originally stored."""
        self.assertIsNotNone(molD, "fetch failed for %s" % name)
        self.assertEqual(molD["name"], name)
        ok = oem.deserialize(molD["oeb"])
        self.__lfh.write("Deserialized status = %s\n" % ok)
        self.assertTrue(ok, "deserialize failed for %s" % name)
        self.__lfh.write("Deserialized SMILES (canonical) = %s\n" % oem.getCanSMILES())
        self.__lfh.write("Deserialized SMILES (isomeric)  = %s\n" % oem.getIsoSMILES())
        self.assertEqual(oem.getTitle(), name)
        self.assertEqual(oem.getIsoSMILES(), self.__smilesD[name])

    def testCreateStore(self) -> None:
        """Test case -  build store from a selection of serialized OE molecules."""
        self.__lfh.write("\nStarting OePersistTests testCreateStore\n")
        myPersist = self.__createStore(self.__molList)
        mL = myPersist.getIndex(dbFileName=self.__storePath)
        self.__lfh.write("OePersistTests(testCreateStore) molecule list %r\n" % mL)
        self.assertEqual(mL, self.__expectedNames)
        #
        # Recover the full store into a new instance and compare with the input
        rPersist = OePersist(self.__verbose, self.__lfh)
        self.assertTrue(rPersist.recover(dbFileName=self.__storePath))
        self.assertEqual(rPersist.getMoleculeNameList(), self.__expectedNames)
        self.assertEqual(rPersist.getMoleculeList(), self.__molList)

    def testFetchAll(self) -> None:
        """Test case -  fetch all of the molecules in the persistent store"""
        startTime = time.time()
        self.__lfh.write(
            "\nStarting OePersistTests testFetchAll at %s\n" % time.strftime("%Y %m %d %H:%M:%S", time.localtime())
        )
        self.__createStore(self.__molList)
        myPersist = OePersist(self.__verbose, self.__lfh)
        self.assertTrue(myPersist.open(dbFileName=self.__storePath))
        try:
            moleculeNameList = myPersist.getStoreMoleculeIndex()
            self.assertEqual(moleculeNameList, self.__expectedNames)
            oem = OeBuildMol(verbose=self.__verbose, log=self.__lfh)
            for ccId in moleculeNameList:
                molD = myPersist.fetchMolecule(moleculeName=ccId)
                self.__assertRoundTrip(oem, molD, ccId)
        finally:
            self.assertTrue(myPersist.close())

        endTime = time.time()
        self.__lfh.write(
            "\nCompleted OePersistTests testFetchAll at %s (%d seconds)\n"
            % (time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime)
        )

    def testFetchOne(self) -> None:
        """Test case -  fetch all of the molecules in the persistent store one by one"""
        startTime = time.time()
        self.__lfh.write(
            "\nStarting OePersistTests testFetchOne at %s\n" % time.strftime("%Y %m %d %H:%M:%S", time.localtime())
        )
        self.__createStore(self.__molList)
        myPersist = OePersist(self.__verbose, self.__lfh)
        moleculeNameList = myPersist.getIndex(dbFileName=self.__storePath)
        self.assertEqual(moleculeNameList, self.__expectedNames)
        oem = OeBuildMol(verbose=self.__verbose, log=self.__lfh)
        for ccId in moleculeNameList:
            fetchedD = myPersist.fetchOneMolecule(dbFileName=self.__storePath, moleculeName=ccId)
            self.__assertRoundTrip(oem, fetchedD, ccId)

        endTime = time.time()
        self.__lfh.write(
            "\nCompleted OePersistTests testFetchOne at %s (%d seconds)\n"
            % (time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime)
        )

    def testUpdateStore(self) -> None:
        """Test case -  update store from a selection of serialized OE molecules."""
        self.__lfh.write("\nStarting OePersistTests testUpdateStore\n")
        # Start from a store holding only the first molecule so that updates exercise
        # both the replace and the append paths.
        self.__createStore(self.__molList[:1])
        myPersist = OePersist(self.__verbose, self.__lfh)
        self.assertEqual(myPersist.getIndex(dbFileName=self.__storePath), self.__expectedNames[:1])
        for molD in self.__molList:
            self.assertTrue(myPersist.updateOneMolecule(molD, dbFileName=self.__storePath))

        mL = myPersist.getIndex(dbFileName=self.__storePath)
        self.__lfh.write("OePersistTests(testUpdateStore) molecule list %r\n" % mL)
        self.assertEqual(mL, self.__expectedNames)
        oem = OeBuildMol(verbose=self.__verbose, log=self.__lfh)
        for ccId in self.__expectedNames:
            fetchedD = myPersist.fetchOneMolecule(dbFileName=self.__storePath, moleculeName=ccId)
            self.__assertRoundTrip(oem, fetchedD, ccId)


def suite1() -> unittest.TestSuite:
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OePersistTests("testCreateStore"))
    suiteSelect.addTest(OePersistTests("testUpdateStore"))
    suiteSelect.addTest(OePersistTests("testFetchAll"))
    suiteSelect.addTest(OePersistTests("testFetchOne"))

    return suiteSelect


if __name__ == "__main__":
    mySuite = suite1()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
