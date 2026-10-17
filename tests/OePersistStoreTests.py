##
# File:    OePersistStoreTests.py
#
# Updates:
#
##
"""
Test cases for OePersist store operations using plain molecule dictionaries.

These tests exercise the shelve-based persistence layer directly and do not require OpenEye.
"""

import io
import os
import shutil
import tempfile
import unittest
from typing import Any, Dict, List

from wwpdb.utils.oe_util.build.OePersist import OePersist


def _makeMolecule(name: str, payload: str = "") -> Dict[str, Any]:
    return {"name": name, "oeb": payload or "oeb-%s" % name}


class OePersistStoreTests(unittest.TestCase):
    def setUp(self) -> None:
        self.__tmpDir = tempfile.mkdtemp()
        self.__storePath = os.path.join(self.__tmpDir, "oe-store.db")
        self.__lfh = io.StringIO()
        self.__molList: List[Dict[str, Any]] = [_makeMolecule("ATP"), _makeMolecule("GTP"), _makeMolecule("ARG")]

    def tearDown(self) -> None:
        shutil.rmtree(self.__tmpDir, ignore_errors=True)

    def __newPersist(self) -> OePersist:
        return OePersist(verbose=True, log=self.__lfh)

    def __createStore(self) -> None:
        myPersist = self.__newPersist()
        self.assertTrue(myPersist.setMoleculeList(moleculeList=list(self.__molList)))
        self.assertTrue(myPersist.store(dbFileName=self.__storePath))

    # ---- in-memory list handling ----

    def testInMemoryList(self) -> None:
        """Test case - set, append and look up molecules without a store."""
        myPersist = self.__newPersist()
        self.assertTrue(myPersist.setMoleculeList(moleculeList=self.__molList[:2]))
        self.assertEqual(myPersist.getMoleculeNameList(), ["ATP", "GTP"])
        self.assertTrue(myPersist.appendMoleculeList(moleculeList=self.__molList[2:]))
        self.assertEqual(myPersist.getMoleculeNameList(), ["ATP", "GTP", "ARG"])
        self.assertEqual(len(myPersist.getMoleculeList()), 3)
        self.assertEqual(myPersist.getMolecule("GTP"), self.__molList[1])
        self.assertIsNone(myPersist.getMolecule("XXX"))

    def testSetMoleculeListInvalid(self) -> None:
        """Test case - setting a None molecule list reports failure."""
        myPersist = self.__newPersist()
        self.assertFalse(myPersist.setMoleculeList(moleculeList=None))
        self.assertIn("setMoleculeList() initialization failed", self.__lfh.getvalue())

    def testAppendMoleculeListInvalid(self) -> None:
        """Test case - appending a None molecule list reports failure."""
        myPersist = self.__newPersist()
        self.assertFalse(myPersist.appendMoleculeList(moleculeList=None))
        self.assertIn("appendMoleculeList() initialization failed", self.__lfh.getvalue())

    # ---- whole-store operations ----

    def testStoreAndRecover(self) -> None:
        """Test case - store a molecule list and recover it in a new instance."""
        self.__createStore()
        myPersist = self.__newPersist()
        self.assertEqual(myPersist.getIndex(dbFileName=self.__storePath), ["ATP", "GTP", "ARG"])
        self.assertTrue(myPersist.recover(dbFileName=self.__storePath))
        self.assertEqual(myPersist.getMoleculeList(), self.__molList)
        self.assertEqual(myPersist.getMolecule("ARG"), self.__molList[2])
        # Lock file must be released after each locked operation
        self.assertFalse(os.path.exists(self.__storePath + ".lock"))

    def testRecoverMissingStore(self) -> None:
        """Test case - recovering a nonexistent store fails."""
        myPersist = self.__newPersist()
        self.assertFalse(myPersist.recover(dbFileName=self.__storePath))
        self.assertIn("shelve recover failed", self.__lfh.getvalue())

    def testGetIndexMissingStore(self) -> None:
        """Test case - index of a nonexistent store is empty."""
        myPersist = self.__newPersist()
        self.assertFalse(myPersist.getIndex(dbFileName=self.__storePath))
        self.assertIn("shelve index failed", self.__lfh.getvalue())

    def testStoreInvalidMolecule(self) -> None:
        """Test case - storing a molecule without a name fails."""
        myPersist = self.__newPersist()
        self.assertTrue(myPersist.setMoleculeList(moleculeList=[{"oeb": "no-name"}]))
        self.assertFalse(myPersist.store(dbFileName=self.__storePath))
        self.assertIn("shelve store failed", self.__lfh.getvalue())

    # ---- single molecule fetch ----

    def testFetchOneMolecule(self) -> None:
        """Test case - fetch single molecules atomically from the store."""
        self.__createStore()
        myPersist = self.__newPersist()
        self.assertEqual(myPersist.fetchOneMolecule(dbFileName=self.__storePath, moleculeName="GTP"), self.__molList[1])
        self.assertIsNone(myPersist.fetchOneMolecule(dbFileName=self.__storePath, moleculeName="XXX"))
        self.assertIn("shelve fetch failed", self.__lfh.getvalue())

    def testOpenFetchClose(self) -> None:
        """Test case - open a store and fetch several molecules before closing."""
        self.__createStore()
        myPersist = self.__newPersist()
        self.assertTrue(myPersist.open(dbFileName=self.__storePath))
        self.assertTrue(os.path.exists(self.__storePath + ".lock"))
        self.assertEqual(myPersist.getStoreMoleculeIndex(), ["ATP", "GTP", "ARG"])
        for mol in self.__molList:
            self.assertEqual(myPersist.fetchMolecule(moleculeName=mol["name"]), mol)
        self.assertIsNone(myPersist.fetchMolecule(moleculeName="XXX"))
        self.assertTrue(myPersist.close())
        self.assertFalse(os.path.exists(self.__storePath + ".lock"))

    def testOpenMissingStore(self) -> None:
        """Test case - opening a nonexistent store read-only fails and releases the lock."""
        myPersist = self.__newPersist()
        self.assertFalse(myPersist.open(dbFileName=self.__storePath, flag="r"))
        self.assertIn("open() write failed", self.__lfh.getvalue())
        self.assertFalse(os.path.exists(self.__storePath + ".lock"))

    def testUseWithoutOpen(self) -> None:
        """Test case - store access methods fail cleanly when no store is open."""
        myPersist = self.__newPersist()
        self.assertEqual(myPersist.getStoreMoleculeIndex(), [])
        self.assertIsNone(myPersist.fetchMolecule(moleculeName="ATP"))
        self.assertFalse(myPersist.close())

    # ---- updates ----

    def testUpdateOneMolecule(self) -> None:
        """Test case - update an existing molecule and append a new one."""
        self.__createStore()
        myPersist = self.__newPersist()
        replacement = _makeMolecule("GTP", "new-gtp")
        self.assertTrue(myPersist.updateOneMolecule(replacement, dbFileName=self.__storePath))
        self.assertTrue(myPersist.updateOneMolecule(_makeMolecule("TRP"), dbFileName=self.__storePath))
        self.assertEqual(myPersist.getIndex(dbFileName=self.__storePath), ["ATP", "GTP", "ARG", "TRP"])
        self.assertEqual(myPersist.fetchOneMolecule(dbFileName=self.__storePath, moleculeName="GTP"), replacement)

    def testUpdateOneMoleculeMissingStore(self) -> None:
        """Test case - updating a nonexistent store fails."""
        myPersist = self.__newPersist()
        self.assertFalse(myPersist.updateOneMolecule(_makeMolecule("TRP"), dbFileName=self.__storePath))
        self.assertIn("update failed", self.__lfh.getvalue())

    def testUpdateMoleculeList(self) -> None:
        """Test case - update and append a list of molecules in an existing store."""
        self.__createStore()
        myPersist = self.__newPersist()
        updateList = [_makeMolecule("ATP", "new-atp"), _makeMolecule("TRP"), _makeMolecule("MSE")]
        self.assertTrue(myPersist.updateMoleculeList(dbFileName=self.__storePath, moleculeList=updateList))
        self.assertEqual(myPersist.getIndex(dbFileName=self.__storePath), ["ATP", "GTP", "ARG", "TRP", "MSE"])
        self.assertTrue(myPersist.recover(dbFileName=self.__storePath))
        self.assertEqual(myPersist.getMolecule("ATP"), updateList[0])
        self.assertEqual(myPersist.getMolecule("MSE"), updateList[2])

    def testUpdateMoleculeListInvalid(self) -> None:
        """Test case - updating with a None molecule list fails."""
        self.__createStore()
        myPersist = self.__newPersist()
        self.assertFalse(myPersist.updateMoleculeList(dbFileName=self.__storePath, moleculeList=None))
        self.assertIn("updateMoleculeListShelve() update failed", self.__lfh.getvalue())

    # ---- move ----

    def testMoveStore(self) -> None:
        """Test case - move a file to a new destination under a destination lock."""
        srcPath = os.path.join(self.__tmpDir, "src.txt")
        dstPath = os.path.join(self.__tmpDir, "dst.txt")
        with open(srcPath, "w", encoding="utf-8") as ofh:
            ofh.write("content")
        myPersist = self.__newPersist()
        self.assertTrue(myPersist.moveStore(srcPath, dstPath))
        self.assertFalse(os.path.exists(srcPath))
        with open(dstPath, encoding="utf-8") as ifh:
            self.assertEqual(ifh.read(), "content")
        self.assertFalse(os.path.exists(dstPath + ".lock"))

    def testMoveStoreMissingSource(self) -> None:
        """Test case - moving a nonexistent source fails."""
        myPersist = self.__newPersist()
        self.assertFalse(
            myPersist.moveStore(os.path.join(self.__tmpDir, "none.db"), os.path.join(self.__tmpDir, "dst.db"))
        )
        self.assertIn("move failed", self.__lfh.getvalue())


if __name__ == "__main__":
    unittest.main()
