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

import sys
import unittest
import traceback
import fnmatch
import time
import os
import os.path
import string

try:
    from wwpdb.utils.cc_dict_util.persist.PdbxChemCompDictUtil import PdbxChemCompDictUtil
    from wwpdb.utils.cc_dict_util.persist.PdbxChemCompDictIndex import PdbxChemCompDictIndex

    from wwpdb.utils.oe_util.build.OePersist import OePersist
    from wwpdb.utils.oe_util.build.OeBuildMol import OeBuildMol

    from mmcif_utils.persist.PdbxPersist import PdbxPersist
    skiptests = False
except ImportError:
    skiptests = True


@unittest.skipIf(skiptests, "Cannot import openeye.oechem for tests")
class OePersistFullDictTests(unittest.TestCase):

    def setUp(self):
        self.__lfh = sys.stdout
        self.__verbose = True
        self.__debug = False
        ##
        # Set these as appropriate to checked out versions of the CC and PRDCC CVS repositories
        ##
        self.__pathChemCompCVS = "../../../../../reference/components/ligand-dict-v3"
        self.__pathPrdChemCompCVS = "../../../../../reference/components/prdcc-v3"
        ##
        # file names for persistent stores -
        self.__persistStorePathCC = "chemcomp-store.db"
        self.__indexPathCC = "chemcomp-index.pic"
        self.__storePath = "oe-store.db"
        #
        # Test list of PRD molecule ids
        #
        self.__prdIdListO = ['PRDCC_000009', 'PRDCC_000109', 'PRDCC_000159', 'PRDCC_000199', 'PRDCC_000209', 'PRDCC_000219',
                             'PRDCC_000239', 'PRDCC_000259', 'PRDCC_000289', 'PRDCC_000299', 'PRDCC_000309', 'PRDCC_000319', 'PRDCC_000339',
                             'PRDCC_000379', 'PRDCC_000409', 'PRDCC_000419']

        #
        self.__prdIdList = ['PRD_000009', 'PRD_000109', 'PRD_000159', 'PRD_000199', 'PRD_000209', 'PRD_000219',
                            'PRD_000239', 'PRD_000259', 'PRD_000289', 'PRD_000299', 'PRD_000309', 'PRD_000319', 'PRD_000339',
                            'PRD_000379', 'PRD_000409', 'PRD_000419']

    def tearDown(self):
        pass

    def getPathList(self, topPath, pattern='*', excludeDirs=[], recurse=True):
        """ Return a list of file paths in the input topPath which satisfy the input search criteria.

            This version does not follow symbolic links.
        """
        pathList = []
        #
        try:
            names = os.listdir(topPath)
        except os.error:
            return pathList

        # expand pattern
        pattern = pattern or '*'
        patternList = string.splitfields(pattern, ';')

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
        self.__lfh.write("\nStarting %s %s at %s\n" % (self.__class__.__name__,
                                                       sys._getframe().f_code.co_name,
                                                       time.strftime("%Y %m %d %H:%M:%S", time.localtime())))
        try:
            ccPathList = self.getPathList(topPath=self.__pathChemCompCVS, pattern="*.cif", excludeDirs=['CVS', 'REMOVED', 'FULL'])
            if (self.__verbose):
                self.__lfh.write("Pathlist length is %d\n" % len(ccPathList))
            dUtil = PdbxChemCompDictUtil(verbose=self.__verbose, log=self.__lfh)
            dUtil.makeStoreFromPathList(pathList=ccPathList, storePath=self.__persistStorePathCC)
        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()

        endTime = time.time()
        self.__lfh.write("\nCompleted %s %s at %s (%d seconds)\n" % (self.__class__.__name__,
                                                                     sys._getframe().f_code.co_name,
                                                                     time.strftime("%Y %m %d %H:%M:%S", time.localtime()),
                                                                     endTime - startTime))

    def testUpdateChemCompStoreWithPrd(self):
        """Test case -  update persistent store from a path list of PRD chemical component defintions.

                        Extract the path list from the file system of the PRD CVS repository.
        """
        startTime = time.time()
        self.__lfh.write("\nStarting %s %s at %s\n" % (self.__class__.__name__,
                                                       sys._getframe().f_code.co_name,
                                                       time.strftime("%Y %m %d %H:%M:%S", time.localtime())))
        try:
            ccPathList = self.getPathList(topPath=self.__pathPrdChemCompCVS, pattern="*.cif", excludeDirs=['CVS', 'REMOVED', 'FULL'])
            if (self.__verbose):
                self.__lfh.write("PRDCC pathlist length is %d\n" % len(ccPathList))
            #
            dUtil = PdbxChemCompDictUtil(verbose=self.__verbose, log=self.__lfh)
            dUtil.updateStoreByFile(pathList=ccPathList, storePath=self.__persistStorePathCC)
            #
        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()

        endTime = time.time()
        self.__lfh.write("\nCompleted %s %s at %s (%d seconds)\n" % (self.__class__.__name__,
                                                                     sys._getframe().f_code.co_name,
                                                                     time.strftime("%Y %m %d %H:%M:%S", time.localtime()),
                                                                     endTime - startTime))

    def testCreateChemCompIndex(self):
        """Test case -  create search index from chemical component persistent store
        """
        startTime = time.time()
        self.__lfh.write("\nStarting %s %s at %s\n" % (self.__class__.__name__,
                                                       sys._getframe().f_code.co_name,
                                                       time.strftime("%Y %m %d %H:%M:%S", time.localtime())))
        try:
            dIndx = PdbxChemCompDictIndex(verbose=self.__verbose, log=self.__lfh)
            dIndx.makeIndex(storePath=self.__persistStorePathCC, indexPath=self.__indexPathCC)
        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()

        endTime = time.time()
        self.__lfh.write("\nCompleted %s %s at %s (%d seconds)\n" % (self.__class__.__name__,
                                                                     sys._getframe().f_code.co_name,
                                                                     time.strftime("%Y %m %d %H:%M:%S", time.localtime()),
                                                                     endTime - startTime))

    def testCreateStoreOE(self):
        """Test case -  build persistent store of OE molecules using the contents of the persistent store
                        of chemical component defintions.
        """
        startTime = time.time()
        self.__lfh.write("\nStarting %s %s at %s\n" % (self.__class__.__name__,
                                                       sys._getframe().f_code.co_name,
                                                       time.strftime("%Y %m %d %H:%M:%S", time.localtime())))
        try:
            myPersist = PdbxPersist(self.__verbose, self.__lfh)
            indexD = myPersist.getIndex(dbFileName=self.__persistStorePathCC)
            myPersist.open(dbFileName=self.__persistStorePathCC)
            containerNameList = myPersist.getStoreContainerIndex()

            oem = OeBuildMol(verbose=self.__verbose, log=self.__lfh)
            molList = []
            for ccId in containerNameList:
                ccAt = myPersist.fetchObject(containerName=ccId, objectName='chem_comp_atom')
                ccBnd = myPersist.fetchObject(containerName=ccId, objectName='chem_comp_bond')
                if ccAt is None or ccBnd is None:
                    continue
                #
                oem.set(ccId, dcChemCompAtom=ccAt, dcChemCompBond=ccBnd)
                oem.build3D()
                if (self.__debug):
                    self.__lfh.write("Title              = %s\n" % oem.getTitle())
                    self.__lfh.write("SMILES (canonical) = %s\n" % oem.getCanSMILES())
                    self.__lfh.write("SMILES (isomeric)  = %s\n" % oem.getIsoSMILES())
                molD = {}
                molD['name'] = ccId
                molD['oeb'] = oem.serialize()
                molList.append(molD)

            myPersist.close()
            #
            oeP = OePersist(self.__verbose, self.__lfh)
            oeP.setMoleculeList(moleculeList=molList)
            oeP.store(dbFileName=self.__storePath)
            mL = oeP.getIndex(dbFileName=self.__storePath)
            #
            if (self.__debug):
                self.__lfh.write("OePersistTests(testCreateStore) molecule list %r\n" % mL)

        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()

        endTime = time.time()
        self.__lfh.write("\nCompleted %s %s at %s (%d seconds)\n" % (self.__class__.__name__,
                                                                     sys._getframe().f_code.co_name,
                                                                     time.strftime("%Y %m %d %H:%M:%S", time.localtime()),
                                                                     endTime - startTime))

    def testFetchAll(self):
        """Test case -  fetch all of the molecules in an 'open' persistent store
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
                if not ok:
                    self.__lfh.write("Deserialized status %s = %d\n" % (ccId, ok))
                    continue
                if (self.__debug):
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
        """Test case -  fetch all of the molecules in the persistent store one by one.  Each fetch reopens
           the store.
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
                if not ok:
                    self.__lfh.write("Deserialized status %s = %d\n" % (ccId, ok))
                    continue

                if (self.__debug):
                    self.__lfh.write("Deserialized SMILES (canonical) = %s\n" % oem.getCanSMILES())
                    self.__lfh.write("Deserialized SMILES (isomeric)  = %s\n" % oem.getIsoSMILES())

        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()

        endTime = time.time()
        self.__lfh.write("\nCompleted %s %s at %s (%d seconds)\n" % (self.__class__.__name__,
                                                                     sys._getframe().f_code.co_name,
                                                                     time.strftime("%Y %m %d %H:%M:%S", time.localtime()),
                                                                     endTime - startTime))

    def testBoundedFormulaSearch(self):
        """Test case -  for the set of PRD tests cases, perform a bounded formula search on the index file
           and then extract the candidate matching molecules from the OE persistent store.

           This approximates a typical search scenario.

        """
        startTime = time.time()
        self.__lfh.write("\nStarting %s %s at %s\n" % (self.__class__.__name__,
                                                       sys._getframe().f_code.co_name,
                                                       time.strftime("%Y %m %d %H:%M:%S", time.localtime())))
        try:
            myPersist = OePersist(self.__verbose, self.__lfh)
            oem = OeBuildMol(verbose=self.__verbose, log=self.__lfh)
            #
            dIndx = PdbxChemCompDictIndex(verbose=self.__verbose, log=self.__lfh)
            ccIdx = dIndx.readIndex(indexPath=self.__indexPathCC)
            self.__lfh.write("Index %s length = %d\n" % (self.__indexPathCC, len(ccIdx)))
            #
            for prdId in self.__prdIdList:
                if prdId in ccIdx:
                    d = ccIdx[prdId]
                    eCounts = d['typeCounts']
                    mL = dIndx.searchBoundedFormula(elementCounts=eCounts, indexPath=self.__indexPathCC)
                    self.__lfh.write("Formula filter count for %s = %d\n" % (prdId, len(mL)))
                    for ccId in mL:
                        molD = myPersist.fetchOneMolecule(dbFileName=self.__storePath, moleculeName=ccId)
                        name = molD['name']
                        ok = oem.deserialize(molD['oeb'])
                        if not ok:
                            self.__lfh.write("Deserialized status %s = %d\n" % (ccId, ok))
                            continue
                        if (self.__debug):
                            self.__lfh.write("Deserialized SMILES (canonical) = %s\n" % oem.getCanSMILES())
                            self.__lfh.write("Deserialized SMILES (isomeric)  = %s\n" % oem.getIsoSMILES())
                else:
                    self.__lfh.write("+WARNING - Index %s missing id %s \n" % (self.__indexPathCC, prdId))

        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()

        endTime = time.time()
        self.__lfh.write("\nCompleted %s %s at %s (%d seconds)\n" % (self.__class__.__name__,
                                                                     sys._getframe().f_code.co_name,
                                                                     time.strftime("%Y %m %d %H:%M:%S", time.localtime()),
                                                                     endTime - startTime))


def suiteCreateStoreCC():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OePersistFullDictTests("testCreateChemCompStore"))
    suiteSelect.addTest(OePersistFullDictTests("testUpdateChemCompStoreWithPrd"))
    return suiteSelect


def suiteCreateIndexCC():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OePersistFullDictTests("testCreateChemCompIndex"))
    return suiteSelect


def suiteCreateStoreOE():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OePersistFullDictTests("testCreateStoreOE"))
    #    suiteSelect.addTest(OePersistFullDictTests("testFetchAll"))
    #    suiteSelect.addTest(OePersistFullDictTests("testFetchOne"))

    return suiteSelect


def suiteSearchStoreOE():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OePersistFullDictTests("testFetchAll"))
    suiteSelect.addTest(OePersistFullDictTests("testFetchOne"))
    suiteSelect.addTest(OePersistFullDictTests("testBoundedFormulaSearch"))
    return suiteSelect

if __name__ == '__main__':
    #
    if (not os.access("chemcomp-store.db", os.F_OK)):
        mySuite = suiteCreateStoreCC()
        unittest.TextTestRunner(verbosity=2).run(mySuite)

    if (not os.access("chemcomp-index.pic", os.F_OK)):
        mySuite = suiteCreateIndexCC()
        unittest.TextTestRunner(verbosity=2).run(mySuite)

    if (not os.access("oe-store.db", os.F_OK)):
        mySuite = suiteCreateStoreOE()
        unittest.TextTestRunner(verbosity=2).run(mySuite)

    mySuite = suiteSearchStoreOE()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
