##
#
# File:    OeShapeSearchtests.py
# Author:  J. Westbrook
# Date:    22-Feb-2012
# Version: 0.001
#
# Updated:
#    25-Feb-2012  jdw add formula index search example
#     6-Jun-2016  jdw general cleanup
#
##
"""
Test cases for persistent storage of serialized OE molecule objects
for the the full chemical component dictionary + PRD chemical components.

Examples here illustrate building all of the persistent supporting data stores
and index files required to support the search from scratch.   Each of the
stores can also be updated incrementally.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.01"

import sys
import unittest
import traceback
import fnmatch
import time
import platform
import os
import os.path

from wwpdb.utils.cc_dict_util.persist.PdbxChemCompDictUtil import PdbxChemCompDictUtil
from wwpdb.utils.cc_dict_util.persist.PdbxChemCompDictIndex import PdbxChemCompDictIndex

try:
    from openeye.oechem import OEFloatArray  # noqa: F401 pylint: disable=unused-import

    skiptests = False
except ImportError:
    skiptests = True

if not skiptests:
    from wwpdb.utils.oe_util.build.OePersist import OePersist
    from wwpdb.utils.oe_util.build.OeBuildMol import OeBuildMol
    from wwpdb.utils.oe_util.search.OeShapeSearch import OeShapeSearch

    from mmcif_utils.persist.PdbxPersist import PdbxPersist
    from mmcif_utils.persist.PdbxPyIoAdapter import PdbxPyIoAdapter as PdbxIoAdapter


@unittest.skipIf(skiptests, "Cannot import openeye.oechem for tests")
class OeShapeSearchtests(unittest.TestCase):
    def setUp(self):
        self.__lfh = sys.stdout
        self.__verbose = True
        self.__debug = False
        #
        self.__here = os.path.abspath(os.path.dirname(__file__))
        self.__examples = os.path.join(self.__here, "examples")
        self.__datadir = os.path.join(self.__here, "data")
        self.__testoutput = os.path.join(self.__here, "test-output", platform.python_version())
        if not os.path.exists(self.__testoutput):
            os.makedirs(self.__testoutput)
        ##
        # Set these as appropriate to checked out versions of the CC and PRDCC CVS repositories
        ##
        self.__pathChemCompCVS = os.path.join(self.__here, "ligand-dict-v3")
        self.__pathPrdChemCompCVS = os.path.join(self.__here, "prdcc-v3")
        ##
        # file names for persistent stores -
        self.__persistStorePathCC = os.path.join(self.__testoutput, "chemcomp-store.db")
        self.__indexPathCC = os.path.join(self.__testoutput, "chemcomp-index.pic")
        self.__storePath = os.path.join(self.__testoutput, "oe-store.db")
        #
        # Test list of PRD molecule ids  -- In future check the convention on the identifiers
        #
        #
        self.__prdIdListO = [
            "PRDCC_000009",
            "PRDCC_000109",
            "PRDCC_000159",
            "PRDCC_000199",
            "PRDCC_000209",
            "PRDCC_000219",
            "PRDCC_000239",
            "PRDCC_000259",
            "PRDCC_000289",
            "PRDCC_000299",
            "PRDCC_000309",
            "PRDCC_000319",
            "PRDCC_000339",
            "PRDCC_000379",
            "PRDCC_000409",
            "PRDCC_000419",
        ]

        #
        self.__prdIdList = [
            "PRD_000009",
            "PRD_000109",
            "PRD_000159",
            "PRD_000199",
            "PRD_000209",
            "PRD_000219",
            "PRD_000239",
            "PRD_000259",
            "PRD_000289",
            "PRD_000299",
            "PRD_000309",
            "PRD_000319",
            "PRD_000339",
            "PRD_000379",
            "PRD_000409",
            "PRD_000419",
        ]

    def tearDown(self):
        pass

    def getPathList(self, topPath, pattern="*", excludeDirs=None, recurse=True):
        """Return a list of file paths in the input topPath which satisfy the input search criteria.

        This version does not follow symbolic links.
        """
        pathList = []
        if excludeDirs is None:
            excludeDirs = []
        #
        try:
            names = os.listdir(topPath)
        except os.error:
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
                    pathList.extend(self.getPathList(topPath=fullname, pattern=pattern, excludeDirs=excludeDirs, recurse=recurse))

        return pathList

    def testCreateChemCompStore(self):
        """Test case -  create persistent store from a path list of chemical component defintions.

        Extract the path list by searching the file system of the CVS repository
        """
        startTime = time.time()
        self.__lfh.write("\nStarting OeShapeSearchtests testCreateChemStore at %s\n" % time.strftime("%Y %m %d %H:%M:%S", time.localtime()))
        try:
            ccPathList = self.getPathList(topPath=self.__pathChemCompCVS, pattern="*.cif", excludeDirs=["CVS", "REMOVED", "FULL"])
            if self.__verbose:
                self.__lfh.write("Pathlist length is %d\n" % len(ccPathList))
            dUtil = PdbxChemCompDictUtil(verbose=self.__verbose, log=self.__lfh)
            dUtil.makeStoreFromPathList(pathList=ccPathList, storePath=self.__persistStorePathCC)
        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()

        endTime = time.time()
        self.__lfh.write(
            "\nCompleted OeShapeSearchtests testCreateChemCompStore at %s (%d seconds)\n" % (time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime)
        )

    def testUpdateChemCompStoreWithPrd(self):
        """Test case -  update persistent store from a path list of PRD chemical component defintions.

        Extract the path list from the fie system of the PRD CVS repository.
        """
        startTime = time.time()
        self.__lfh.write("\nStarting OeShapeSearchtests testUpdateChemCompStorePrd at %s\n" % time.strftime("%Y %m %d %H:%M:%S", time.localtime()))
        try:
            ccPathList = self.getPathList(topPath=self.__pathPrdChemCompCVS, pattern="*.cif", excludeDirs=["CVS", "REMOVED", "FULL"])
            if self.__verbose:
                self.__lfh.write("Pathlist length is %d\n" % len(ccPathList))
            #
            dUtil = PdbxChemCompDictUtil(verbose=self.__verbose, log=self.__lfh)
            dUtil.updateStoreByFile(pathList=ccPathList, storePath=self.__persistStorePathCC)
            #
        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()

        endTime = time.time()
        self.__lfh.write(
            "\nCompleted OeShapeSearchtests testUpdateChemCompStoreWithPrd at %s (%d seconds)\n" % (time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime)
        )

    def testCreateChemCompIndex(self):
        """Test case -  create search index from chemical component persistent store"""
        startTime = time.time()
        self.__lfh.write("\nStarting OeShapeSearchtests testCreateChemCompIndex at %s\n" % time.strftime("%Y %m %d %H:%M:%S", time.localtime()))
        try:
            dIndx = PdbxChemCompDictIndex(verbose=self.__verbose, log=self.__lfh)
            dIndx.makeIndex(storePath=self.__persistStorePathCC, indexPath=self.__indexPathCC)
        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()

        endTime = time.time()
        self.__lfh.write("\nCompleted OeShapeSearchteststestCreateChemCompIndex at %s (%d seconds)\n" % (time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime))

    def testCreateStoreOE(self):
        """Test case -  build persistent store of serialized OE molecules using the contents of the chemical
        dictionary persistent store containing all of chemical component defintions.

        ***NOTE - This will display diagnostics from OE toolkit for molecules with
           problematic features.  These issues will not impact the shape search.

        """
        startTime = time.time()
        self.__lfh.write("\nStarting OeShapeSearchtests testCreateCoreOE at %s\n" % time.strftime("%Y %m %d %H:%M:%S", time.localtime()))
        try:
            # Get handle for
            myPersist = PdbxPersist(self.__verbose, self.__lfh)
            indexD = myPersist.getIndex(dbFileName=self.__persistStorePathCC)  # noqa: F841 pylint: disable=unused-variable
            myPersist.open(dbFileName=self.__persistStorePathCC)
            containerNameList = myPersist.getStoreContainerIndex()
            oem = OeBuildMol(verbose=self.__verbose, log=self.__lfh)
            molList = []
            for ccId in containerNameList:
                ccAt = myPersist.fetchObject(containerName=ccId, objectName="chem_comp_atom")
                ccBnd = myPersist.fetchObject(containerName=ccId, objectName="chem_comp_bond")
                if ccAt is None or ccBnd is None:
                    continue
                #
                oem.set(ccId, dcChemCompAtom=ccAt, dcChemCompBond=ccBnd)
                ok = oem.build3D()
                if ok:
                    if self.__debug:
                        self.__lfh.write("Title              = %s\n" % oem.getTitle())
                        self.__lfh.write("SMILES (canonical) = %s\n" % oem.getCanSMILES())
                        self.__lfh.write("SMILES (isomeric)  = %s\n" % oem.getIsoSMILES())
                    molD = {}
                    molD["name"] = ccId
                    molD["oeb"] = oem.serialize()
                    molList.append(molD)
                else:
                    self.__lfh.write("+WARN - failed to build OE molecule for %s\n" % ccId)

            myPersist.close()
            #
            oeP = OePersist(self.__verbose, self.__lfh)
            oeP.setMoleculeList(moleculeList=molList)
            oeP.store(dbFileName=self.__storePath)
            mL = oeP.getIndex(dbFileName=self.__storePath)
            #
            if self.__debug:
                self.__lfh.write("OePersistTests(testCreateStore) molecule list %r\n" % mL)

        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()

        endTime = time.time()
        self.__lfh.write("\nCompleted OeShapeSearchtests testCreateOEStore at %s (%d seconds)\n" % (time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime))

    def testSimpleShapeSearch(self):
        """Test case -  build OE molecule from chemical component definitions and perform shape search"""
        startTime = time.time()
        self.__lfh.write("\nStarting OeShapeSearchtests testSimpleShapeSearch at %s\n" % time.strftime("%Y %m %d %H:%M:%S", time.localtime()))
        try:
            #
            # Get the path list of PRD CC  -
            ccPathList = self.getPathList(topPath=self.__pathPrdChemCompCVS, pattern="*.cif", excludeDirs=["CVS", "REMOVED", "FULL"])
            for ii, ccPath in enumerate(ccPathList):
                self.__lfh.write(" %d %s\n" % (ii, ccPath))

            #
            # Initialize the shape search and set the reference molecule
            #
            oeMolRef = self.__getOEMol(ccPath=ccPathList[5])
            oeShape = OeShapeSearch(verbose=self.__verbose, log=self.__lfh)
            oeShape.setRefMol(oeMolRef)
            retList = []
            for ccPath in ccPathList[279:290]:
                oeMolFit = self.__getOEMol(ccPath=ccPath)
                #
                # Test each of the molecule in the path list -
                # Invoke the shape search here and save the score result list
                #
                rD = oeShape.setFitMol(oeMolFit)
                retList.append(rD)
                #
            # Output the list of shape scores --
            for r in retList:
                self.__lfh.write("Reference %s match %r\n" % (oeMolRef.GetTitle(), r.items()))
        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()

        endTime = time.time()
        self.__lfh.write("\nCompleted OeShapeSearchtests testSimpleShapeSearch at %s (%d seconds)\n" % (time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime))

    def __getOEMol(self, ccPath):
        """Test case -  build OE molecule using 3D data in the definition source data"""
        self.__lfh.write("\nStarting OeShapeSearchtests __getOEMol\n")
        try:
            oem = OeBuildMol(verbose=self.__verbose, log=self.__lfh)
            myReader = PdbxIoAdapter(self.__verbose, self.__lfh)
            ok = myReader.read(pdbxFilePath=ccPath)
            self.assertTrue(ok)
            #
            # myReader.write(pdbxFilePath="TMP.cif")
            #
            for container in myReader.getContainerList():
                oem.set(container.getName(), dcChemCompAtom=container.getObj("chem_comp_atom"), dcChemCompBond=container.getObj("chem_comp_bond"))
                oem.build3D(coordType="model")
                if self.__debug:
                    self.__lfh.write("Title              = %s\n" % oem.getTitle())
                    self.__lfh.write("SMILES (canonical) = %s\n" % oem.getCanSMILES())
                    self.__lfh.write("SMILES (isomeric)  = %s\n" % oem.getIsoSMILES())

            return oem.getMol()

        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()


def suiteCreateStoreCC():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeShapeSearchtests("testCreateChemCompStore"))
    suiteSelect.addTest(OeShapeSearchtests("testUpdateChemCompStoreWithPrd"))
    return suiteSelect


def suiteCreateIndexCC():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeShapeSearchtests("testCreateChemCompIndex"))
    return suiteSelect


def suiteCreateStoreOE():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeShapeSearchtests("testCreateStoreOE"))
    return suiteSelect


def suiteSearchStoreOE():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeShapeSearchtests("testBoundedFormulaShapeSearch"))
    return suiteSelect


def suiteSearchSimple():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeShapeSearchtests("testSimpleShapeSearch"))
    return suiteSelect


if __name__ == "__main__":
    #

    if not os.access("chemcomp-store.db", os.F_OK):
        mySuite = suiteCreateStoreCC()
        unittest.TextTestRunner(verbosity=2).run(mySuite)

    if not os.access("chemcomp-index.pic", os.F_OK):
        mySuite = suiteCreateIndexCC()
        unittest.TextTestRunner(verbosity=2).run(mySuite)

    if not os.access("oe-store.db", os.F_OK):
        mySuite = suiteCreateStoreOE()
        unittest.TextTestRunner(verbosity=2).run(mySuite)

    mySuite = suiteSearchStoreOE()
    unittest.TextTestRunner(verbosity=2).run(mySuite)

    mySuite = suiteSearchSimple()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
