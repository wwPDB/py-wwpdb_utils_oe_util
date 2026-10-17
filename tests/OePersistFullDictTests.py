# pylint: disable=possibly-used-before-assignment
##
#
# File:    OePersistFullDictTests.py
# Author:  J. Westbrook
# Date:    24-Jan-2012
# Version: 0.001
#
# Updated:
#    25-Feb-2012  jdw add formula index search example
#     6-Jun-2016  jdw general cleanup
#
##
"""
Test cases for persistent storage of serialized OE molecule objects
for the the full chemical component dictionary + PRD chemical components

"""

__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.01"

import fnmatch
import os
import os.path
import sys
import tempfile
import time
import unittest
from typing import Any, ClassVar, Dict, List, Optional, TextIO, Tuple

try:
    from openeye.oechem import OEFloatArray  # noqa: F401 pylint: disable=unused-import

    skiptests = False
except ImportError:
    skiptests = True

if not skiptests:
    from mmcif_utils.persist.PdbxPersist import PdbxPersist
    from wwpdb.utils.cc_dict_util.persist.PdbxChemCompDictIndex import PdbxChemCompDictIndex
    from wwpdb.utils.cc_dict_util.persist.PdbxChemCompDictUtil import PdbxChemCompDictUtil

    from wwpdb.utils.oe_util.build.OeBuildMol import OeBuildMol
    from wwpdb.utils.oe_util.build.OePersist import OePersist

import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()

HERE = os.path.abspath(os.path.dirname(__file__))
##
# Set these as appropriate to checked out versions of the CC and PRDCC CVS repositories
##
PATH_CHEM_COMP_CVS = os.path.join(HERE, "ligand-dict-v3")
PATH_PRD_CHEM_COMP_CVS = os.path.join(HERE, "prdcc-v3")
EXCLUDE_DIRS = ["CVS", "REMOVED", "FULL"]


def getPathList(
    topPath: str, pattern: str = "*", excludeDirs: Optional[List[str]] = None, recurse: bool = True
) -> List[str]:
    """Return a list of file paths in the input topPath which satisfy the input search criteria.

    This version does not follow symbolic links.
    """
    pathList: List[str] = []
    if excludeDirs is None:
        excludeDirs = []
    try:
        names = os.listdir(topPath)
    except OSError:
        return pathList

    # expand pattern
    pattern = pattern or "*"
    patternList = str.split(pattern, ";")

    for name in names:
        fullname = os.path.normpath(os.path.join(topPath, name))
        # check for matching files
        for pat in patternList:
            if fnmatch.fnmatch(name, pat):
                if os.path.isfile(fullname):
                    pathList.append(fullname)
                    continue
        if recurse:
            # recursively scan directories
            if os.path.isdir(fullname) and not os.path.islink(fullname) and (name not in excludeDirs):
                pathList.extend(
                    getPathList(topPath=fullname, pattern=pattern, excludeDirs=excludeDirs, recurse=recurse)
                )

    return pathList


def getContainerNameList(storePath: str, verbose: bool, lfh: TextIO) -> List[str]:
    """Return the container names recorded in the index of a chemical component persistent store."""
    myPersist = PdbxPersist(verbose, lfh)
    indexD = myPersist.getIndex(dbFileName=storePath)
    return [str(tup[0]) for tup in indexD.get("__containers__", [])]


def makeChemCompStore(storePath: str, verbose: bool, lfh: TextIO) -> Tuple[bool, List[str]]:
    """Create a chemical component persistent store from the CC repository.  Return status and the input path list."""
    ccPathList = getPathList(topPath=PATH_CHEM_COMP_CVS, pattern="*.cif", excludeDirs=EXCLUDE_DIRS)
    if verbose:
        lfh.write("Pathlist length is %d\n" % len(ccPathList))
    dUtil = PdbxChemCompDictUtil(verbose=verbose, log=lfh)
    ok = dUtil.makeStoreFromPathList(pathList=ccPathList, storePath=storePath)
    return ok, ccPathList


def makeOeStore(
    ccStorePath: str, oeStorePath: str, verbose: bool, lfh: TextIO, debug: bool = False
) -> Tuple[bool, Dict[str, str]]:
    """Build a persistent store of OE molecules from the contents of a chemical component persistent store.

    Return the store status and a dictionary of isomeric SMILES for each stored molecule.
    """
    myPersist = PdbxPersist(verbose, lfh)
    myPersist.open(dbFileName=ccStorePath)
    containerNameList = myPersist.getStoreContainerIndex()

    oem = OeBuildMol(verbose=verbose, log=lfh)
    molList: List[Dict[str, Any]] = []
    smilesD: Dict[str, str] = {}
    for ccId in containerNameList:
        ccAt = myPersist.fetchObject(containerName=ccId, objectName="chem_comp_atom")
        ccBnd = myPersist.fetchObject(containerName=ccId, objectName="chem_comp_bond")
        if ccAt is None or ccBnd is None:
            continue
        oem.set(ccId, dcChemCompAtom=ccAt, dcChemCompBond=ccBnd)
        oem.build3D()
        if debug:
            lfh.write("Title              = %s\n" % oem.getTitle())
            lfh.write("SMILES (canonical) = %s\n" % oem.getCanSMILES())
            lfh.write("SMILES (isomeric)  = %s\n" % oem.getIsoSMILES())
        smilesD[ccId] = oem.getIsoSMILES()
        molList.append({"name": ccId, "oeb": oem.serialize()})

    myPersist.close()
    oeP = OePersist(verbose, lfh)
    ok = oeP.setMoleculeList(moleculeList=molList) and oeP.store(dbFileName=oeStorePath)
    return ok, smilesD


@unittest.skipIf(skiptests, "Cannot import openeye.oechem for tests")
@unittest.skipUnless(
    os.path.isdir(PATH_CHEM_COMP_CVS), "Chemical component repository %s not available" % PATH_CHEM_COMP_CVS
)
class OePersistFullDictTests(unittest.TestCase):
    # Shared, read-only stores used by the fetch tests - built once in setUpClass
    sharedTmpDir: "ClassVar[Optional[tempfile.TemporaryDirectory[str]]]" = None
    sharedStorePath: ClassVar[str] = ""
    sharedSmilesD: ClassVar[Dict[str, str]] = {}

    @classmethod
    def setUpClass(cls) -> None:
        cls.sharedTmpDir = tempfile.TemporaryDirectory()
        ccStorePath = os.path.join(cls.sharedTmpDir.name, "chemcomp-store.db")
        cls.sharedStorePath = os.path.join(cls.sharedTmpDir.name, "oe-store.db")
        ok, _ = makeChemCompStore(ccStorePath, True, sys.stdout)
        if not ok:
            raise RuntimeError("Failed to build shared chemical component store")
        ok, cls.sharedSmilesD = makeOeStore(ccStorePath, cls.sharedStorePath, True, sys.stdout)
        if not ok or not cls.sharedSmilesD:
            raise RuntimeError("Failed to build shared OE molecule store")

    @classmethod
    def tearDownClass(cls) -> None:
        if cls.sharedTmpDir is not None:
            cls.sharedTmpDir.cleanup()
            cls.sharedTmpDir = None

    def setUp(self) -> None:
        self.__lfh = sys.stdout
        self.__verbose = True
        self.__debug = False
        ##
        self.__tmpDir = tempfile.TemporaryDirectory()
        self.__testoutput = self.__tmpDir.name
        ##
        # file names for persistent stores -
        self.__persistStorePathCC = os.path.join(self.__testoutput, "chemcomp-store.db")
        self.__indexPathCC = os.path.join(self.__testoutput, "chemcomp-index.pic")
        self.__storePath = os.path.join(self.__testoutput, "oe-store.db")

    def tearDown(self) -> None:
        self.__tmpDir.cleanup()

    def __createChemCompStore(self) -> List[str]:
        """Create the per-test chemical component store and return the list of stored container names."""
        ok, ccPathList = makeChemCompStore(self.__persistStorePathCC, self.__verbose, self.__lfh)
        self.assertTrue(ok, "makeStoreFromPathList failed")
        self.assertTrue(len(ccPathList) > 0)
        containerNameList = getContainerNameList(self.__persistStorePathCC, self.__verbose, self.__lfh)
        expectedIds = sorted(os.path.splitext(os.path.basename(pth))[0] for pth in ccPathList)
        self.assertEqual(sorted(containerNameList), expectedIds)
        return containerNameList

    def testCreateChemCompStore(self) -> None:
        """Test case -  create persistent store from a path list of chemical component defintions.

        Extract the path list by searching the file system of the CVS repository
        """
        startTime = time.time()
        self.__lfh.write(
            "\nStarting OePersistFullDictTests testCreateChemCompStore at %s\n"
            % (time.strftime("%Y %m %d %H:%M:%S", time.localtime()))
        )
        containerNameList = self.__createChemCompStore()
        self.assertIn("ATP", containerNameList)

        endTime = time.time()
        self.__lfh.write(
            "\nCompleted OePersistFullDictTests testCreateChemCompStore at %s (%d seconds)\n"
            % (time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime)
        )

    @unittest.skipUnless(os.path.isdir(PATH_PRD_CHEM_COMP_CVS), "PRD chemical component repository not available")
    def testUpdateChemCompStoreWithPrd(self) -> None:
        """Test case -  update persistent store from a path list of PRD chemical component defintions.

        Extract the path list from the file system of the PRD CVS repository.
        """
        startTime = time.time()
        self.__lfh.write(
            "\nStarting OePersistFullDictTests testUpdateChemCompStoreWithPrd at %s\n"
            % (time.strftime("%Y %m %d %H:%M:%S", time.localtime()))
        )
        ccNameList = self.__createChemCompStore()
        ccPathList = getPathList(topPath=PATH_PRD_CHEM_COMP_CVS, pattern="*.cif", excludeDirs=EXCLUDE_DIRS)
        if self.__verbose:
            self.__lfh.write("PRDCC pathlist length is %d\n" % len(ccPathList))
        self.assertTrue(len(ccPathList) > 0)
        dUtil = PdbxChemCompDictUtil(verbose=self.__verbose, log=self.__lfh)
        ok = dUtil.updateStoreByFile(pathList=ccPathList, storePath=self.__persistStorePathCC)
        self.assertTrue(ok, "updateStoreByFile failed")

        updatedNameList = getContainerNameList(self.__persistStorePathCC, self.__verbose, self.__lfh)
        prdIds = sorted(os.path.splitext(os.path.basename(pth))[0] for pth in ccPathList)
        self.assertEqual(sorted(updatedNameList), sorted(ccNameList + prdIds))

        endTime = time.time()
        self.__lfh.write(
            "\nCompleted OePersistFullDictTests testUpdateChemCompStoreWithPrd at %s (%d seconds)\n"
            % (time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime)
        )

    def testCreateChemCompIndex(self) -> None:
        """Test case -  create search index from chemical component persistent store"""
        startTime = time.time()
        self.__lfh.write(
            "\nStarting OePersistFullDictTests testCreateChemCompIndex at %s\n"
            % (time.strftime("%Y %m %d %H:%M:%S", time.localtime()))
        )
        containerNameList = self.__createChemCompStore()
        dIndx = PdbxChemCompDictIndex(verbose=self.__verbose, log=self.__lfh)
        ret = dIndx.makeIndex(storePath=self.__persistStorePathCC, indexPath=self.__indexPathCC)
        self.assertTrue(len(ret) != 0, "Index failed")
        self.assertTrue(os.path.exists(self.__indexPathCC))
        self.assertTrue(set(ret.keys()) <= set(containerNameList))
        self.assertIn("ATP", ret)
        # The index written to disk matches the returned index
        self.assertEqual(sorted(dIndx.readIndex(indexPath=self.__indexPathCC).keys()), sorted(ret.keys()))

        endTime = time.time()
        self.__lfh.write(
            "\nCompleted OePersistFullDictTests testCreateChemCompIndex at %s (%d seconds)\n"
            % (time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime)
        )

    def testCreateStoreOE(self) -> None:
        """Test case -  build persistent store of OE molecules using the contents of the persistent store
        of chemical component defintions.
        """
        startTime = time.time()
        self.__lfh.write(
            "\nStarting OePersistFullDictTests testCreateStoreOE at %s\n"
            % (time.strftime("%Y %m %d %H:%M:%S", time.localtime()))
        )
        containerNameList = self.__createChemCompStore()
        ok, smilesD = makeOeStore(
            self.__persistStorePathCC, self.__storePath, self.__verbose, self.__lfh, debug=self.__debug
        )
        self.assertTrue(ok, "OE molecule store creation failed")
        self.assertTrue(len(smilesD) > 0)
        oeP = OePersist(self.__verbose, self.__lfh)
        mL = oeP.getIndex(dbFileName=self.__storePath)
        if self.__debug:
            self.__lfh.write("OePersistTests(testCreateStore) molecule list %r\n" % mL)
        self.assertEqual(sorted(mL), sorted(smilesD.keys()))
        self.assertTrue(set(mL) <= set(containerNameList))
        self.assertIn("ATP", mL)

        endTime = time.time()
        self.__lfh.write(
            "\nCompleted OePersistFullDictTests testCreateStoreOE at %s (%d seconds)\n"
            % (time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime)
        )

    def __assertRoundTrip(self, oem: "OeBuildMol", molD: Any, ccId: str) -> None:
        """Check that a fetched molecule deserializes to the molecule originally stored."""
        self.assertIsNotNone(molD, "fetch failed for %s" % ccId)
        self.assertEqual(molD["name"], ccId)
        ok = oem.deserialize(molD["oeb"])
        if not ok:
            self.__lfh.write("Deserialized status %s = %s\n" % (ccId, ok))
        self.assertTrue(ok, "deserialize failed for %s" % ccId)
        if self.__debug:
            self.__lfh.write("Deserialized SMILES (canonical) = %s\n" % oem.getCanSMILES())
            self.__lfh.write("Deserialized SMILES (isomeric)  = %s\n" % oem.getIsoSMILES())
        self.assertEqual(oem.getIsoSMILES(), self.sharedSmilesD[ccId], "SMILES mismatch for %s" % ccId)

    def testFetchAll(self) -> None:
        """Test case -  fetch all of the molecules in an 'open' persistent store"""
        startTime = time.time()
        self.__lfh.write(
            "\nStarting OePersistFullDictTests testFetchAll at %s\n"
            % (time.strftime("%Y %m %d %H:%M:%S", time.localtime()))
        )
        myPersist = OePersist(self.__verbose, self.__lfh)
        self.assertTrue(myPersist.open(dbFileName=self.sharedStorePath))
        try:
            moleculeNameList = myPersist.getStoreMoleculeIndex()
            self.assertEqual(sorted(moleculeNameList), sorted(self.sharedSmilesD.keys()))
            oem = OeBuildMol(verbose=self.__verbose, log=self.__lfh)
            for ccId in moleculeNameList:
                molD = myPersist.fetchMolecule(moleculeName=ccId)
                self.__assertRoundTrip(oem, molD, ccId)
        finally:
            self.assertTrue(myPersist.close())

        endTime = time.time()
        self.__lfh.write(
            "\nCompleted OePersistFullDictTests testFetchAll at %s (%d seconds)\n"
            % (time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime)
        )

    def testFetchOne(self) -> None:
        """Test case -  fetch all of the molecules in the persistent store one by one.  Each fetch reopens
        the store.
        """
        startTime = time.time()
        self.__lfh.write(
            "\nStarting OePersistFullDictTests testFetchOne at %s\n"
            % (time.strftime("%Y %m %d %H:%M:%S", time.localtime()))
        )
        myPersist = OePersist(self.__verbose, self.__lfh)
        moleculeNameList = myPersist.getIndex(dbFileName=self.sharedStorePath)
        self.assertEqual(sorted(moleculeNameList), sorted(self.sharedSmilesD.keys()))
        oem = OeBuildMol(verbose=self.__verbose, log=self.__lfh)
        for ccId in moleculeNameList:
            molD = myPersist.fetchOneMolecule(dbFileName=self.sharedStorePath, moleculeName=ccId)
            self.__assertRoundTrip(oem, molD, ccId)

        endTime = time.time()
        self.__lfh.write(
            "\nCompleted OePersistFullDictTests testFetchOne at %s (%d seconds)\n"
            % (time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime)
        )


def suiteCreateStoreCC() -> unittest.TestSuite:
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OePersistFullDictTests("testCreateChemCompStore"))
    #    suiteSelect.addTest(OePersistFullDictTests("testUpdateChemCompStoreWithPrd"))
    return suiteSelect


def suiteCreateIndexCC() -> unittest.TestSuite:
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OePersistFullDictTests("testCreateChemCompIndex"))
    return suiteSelect


def suiteCreateStoreOE() -> unittest.TestSuite:
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OePersistFullDictTests("testCreateStoreOE"))
    #    suiteSelect.addTest(OePersistFullDictTests("testFetchAll"))
    #    suiteSelect.addTest(OePersistFullDictTests("testFetchOne"))

    return suiteSelect


def suiteSearchStoreOE() -> unittest.TestSuite:
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OePersistFullDictTests("testFetchAll"))
    suiteSelect.addTest(OePersistFullDictTests("testFetchOne"))
    return suiteSelect


if __name__ == "__main__":
    # Each test builds the stores it needs in a temporary directory, so the suites are independent
    mySuite = suiteCreateStoreCC()
    unittest.TextTestRunner(verbosity=2).run(mySuite)

    # mySuite = suiteCreateIndexCC()
    # unittest.TextTestRunner(verbosity=2).run(mySuite)

    # mySuite = suiteCreateStoreOE()
    # unittest.TextTestRunner(verbosity=2).run(mySuite)

    # mySuite = suiteSearchStoreOE()
    # unittest.TextTestRunner(verbosity=2).run(mySuite)
