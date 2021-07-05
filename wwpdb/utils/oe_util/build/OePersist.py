##
# File: OePersist.py
# Date: 22-Feb-2012  John Westbrook
#
# Updates:
#
#      6-Jun-2016  jdw  - general cleanup
##
"""
Class supporting persistant storage of serialized OE molecules.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.01"


import shelve
import shutil
import sys
import traceback

try:
    import cPickle as pickle
except ImportError:
    import pickle

from mmcif_utils.persist.LockFile import LockFile


class OePersist(object):
    """ Persistent storage for instances of OE molecules.

    """

    def __init__(self, verbose=True, log=sys.stderr):
        self.__verbose = verbose
        self.__debug = True
        self.__lfh = log
        #
        self.__moleculeList = []

        self.__moleculeNameList = []
        #
        # Parameters to tune lock file management --
        self.__timeoutSeconds = 10
        self.__retrySeconds = 0.2
        # placeholder for LockFile object for open/close methods
        self.__lockObj = None
        self.__isExpired = True
        #
        self.__db = None

    def __reIndex(self):
        """ Rebuild name and type lists from molecule object list.
        """
        #
        self.__moleculeNameList = []
        for molecule in self.__moleculeList:
            self.__moleculeNameList.append(molecule['name'])
        #
        self.__isExpired = False

    def __getMoleculeIndex(self, name):
        try:
            if (self.__isExpired):
                self.__reIndex()
            return self.__moleculeNameList.index(name)
        except:  # noqa: E722 pylint: disable=bare-except
            return None

    def setMoleculeList(self, moleculeList=None):
        """ Initialize molecule data in the internal molecule to be persisted to the data store.
        """
        try:
            self.__moleculeList = moleculeList
            self.__isExpired = True
            if (self.__debug):
                self.__lfh.write("+OePersist.setMoleculeList() - Container list length %d\n" % len(self.__moleculeList))
            return True
        except Exception as e:
            if (self.__verbose):
                self.__lfh.write("+ERROR- OePersist.setMoleculeList() initialization failed %s\n" % str(e))
            if (self.__debug):
                traceback.print_exc(file=self.__lfh)
            return False

    def appendMoleculeList(self, moleculeList=None):
        """ Append container data to the internal moleculeList to be persisted to the data store.
        """
        try:
            self.__moleculeList.extend(moleculeList)
            self.__isExpired = True
            if (self.__debug):
                self.__lfh.write("+OePersist.appendMoleculeList() - Container list length %d\n" % len(self.__moleculeList))
            return True
        except Exception as e:
            if (self.__verbose):
                self.__lfh.write("+ERROR- OePersist.appendMoleculeList() initialization failed %s\n" % str(e))
            if (self.__debug):
                traceback.print_exc(file=self.__lfh)
            return False

    def getMolecule(self, moleculeName):
        """ Return the molecule from internal MoleculeList corresponding to the input molecule name.

            **Note that this method is NOT directly fetching data from the internal store.  First recover the
            the store to internal storage before calling this method.
        """
        idx = self.__getMoleculeIndex(moleculeName)
        if idx is not None:
            return self.__moleculeList[idx]
        else:
            return None

    def getMoleculeList(self):
        """ Return the current molecule list.
        """
        return self.__moleculeList

    def getMoleculeNameList(self):
        """ Return the current list of molecule names.
        """
        if (self.__isExpired):
            self.__reIndex()
        return self.__moleculeNameList

    #
    def open(self, dbFileName="my.db", flag='r'):
        """ Open the persistent store with the input database file name and access mode.

            flag = 'r read-only'
                   'c/n new/create'
                   'w read/write'
        """
        try:
            self.__lockObj = LockFile(dbFileName, timeoutSeconds=self.__timeoutSeconds, retrySeconds=self.__retrySeconds,
                                      verbose=self.__verbose, log=self.__lfh)
            self.__lockObj.acquire()
            self.__db = shelve.open(dbFileName, flag=flag, protocol=pickle.HIGHEST_PROTOCOL)
            return True
        except Exception as e:
            if self.__lockObj is not None:
                self.__lockObj.release()
            if (self.__verbose):
                self.__lfh.write("+ERROR- OePersist.open() write failed for file %s %s\n" % (dbFileName, str(e)))
            if (self.__debug):
                traceback.print_exc(file=self.__lfh)
            return False

    def getStoreMoleculeIndex(self):
        try:
            return self.__db['__index__']
        except:  # noqa: E722 pylint: disable=bare-except
            return []

    def close(self):
        """ Close the persistent store
        """
        try:
            self.__db.close()
            if self.__lockObj is not None:
                self.__lockObj.release()
            self.__db = None
            return True
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__lockObj is not None:
                self.__lockObj.release()
            return False

    def moveStore(self, srcDbFilePath, dstDbFilePath):
        """ Move source store to destination store with destination locking.

            No locking is performed on the source file.
        """
        with LockFile(dstDbFilePath, timeoutSeconds=self.__timeoutSeconds, retrySeconds=self.__retrySeconds,
                      verbose=self.__verbose, log=self.__lfh) as lf:  # noqa: F841 pylint: disable=unused-variable
            retVal = self.__moveStore(srcDbFilePath, dstDbFilePath)
        return retVal

    def store(self, dbFileName="my.db"):
        """ Store the current molecule list in persistent database.
        """
        with LockFile(dbFileName, timeoutSeconds=self.__timeoutSeconds, retrySeconds=self.__retrySeconds,
                      verbose=self.__verbose, log=self.__lfh) as lf:  # noqa: F841 pylint: disable=unused-variable
            retVal = self.__storeShelve(dbFileName)
        return retVal

    def recover(self, dbFileName="my.db"):
        """ Recover the stored state to the current in-memory molecule representation.
        """
        with LockFile(dbFileName, timeoutSeconds=self.__timeoutSeconds, retrySeconds=self.__retrySeconds,
                      verbose=self.__verbose, log=self.__lfh) as lf:  # noqa: F841 pylint: disable=unused-variable
            retVal = self.__recoverShelve(dbFileName)
        return retVal

    def getIndex(self, dbFileName="my.db"):
        """ Recover the index of the persistent store.
        """
        with LockFile(dbFileName, timeoutSeconds=self.__timeoutSeconds, retrySeconds=self.__retrySeconds,
                      verbose=self.__verbose, log=self.__lfh) as lf:  # noqa: F841  pylint: disable=unused-variable
            retVal = self.__indexShelve(dbFileName)
        return retVal

    def updateOneMolecule(self, inputMolecule, dbFileName="my.db"):
        """ Update or append a molecule to an existing store.
        """
        with LockFile(dbFileName, timeoutSeconds=self.__timeoutSeconds, retrySeconds=self.__retrySeconds,
                      verbose=self.__verbose, log=self.__lfh) as lf:  # noqa: F841 pylint: disable=unused-variable
            retVal = self.__updateMoleculeShelve(dbFileName=dbFileName, inputMolecule=inputMolecule)
        return retVal

    def updateMoleculeList(self, dbFileName="my.db", moleculeList=None):
        """ Update or append the contents of the input molecule list to an existing store.
        """
        with LockFile(dbFileName, timeoutSeconds=self.__timeoutSeconds, retrySeconds=self.__retrySeconds,
                      verbose=self.__verbose, log=self.__lfh) as lf:  # noqa: F841 pylint: disable=unused-variable
            retVal = self.__updateMoleculeListShelve(dbFileName=dbFileName, updateMoleculeList=moleculeList)
        return retVal

    #
    def fetchOneMolecule(self, dbFileName="my.db", moleculeName=None):
        """  Fetch a single object from a named molecule.  This is atomic operation with respect to the store.
        """
        with LockFile(dbFileName, timeoutSeconds=self.__timeoutSeconds, retrySeconds=self.__retrySeconds,
                      verbose=self.__verbose, log=self.__lfh) as lf:  # noqa: F841  pylint: disable=unused-variable
            retVal = self.__fetchOneMoleculeShelve(dbFileName, moleculeName)
        return retVal

    def fetchMolecule(self, moleculeName=None):
        """ Fetch molecule from an open store.

            Use this method to extract multiple molecules from an open store.
        """
        return self.__fetchMoleculeShelve(moleculeName)

    def __moveStore(self, src, dst):
        """  Internal method to perform file move.
        """
        try:
            shutil.move(src, dst)
            return True
        except Exception as e:
            if (self.__verbose):
                self.__lfh.write("+ERROR- OePersist.__moveStore() move failed for file %s %s\n" % (src, str(e)))
            if (self.__debug):
                traceback.print_exc(file=self.__lfh)
            return False

    def __storeShelve(self, dbFileName="my.db"):
        """  Create a new persistent store using the internal molecule list.
        """
        try:
            if (self.__isExpired):
                self.__reIndex()
            db = shelve.open(dbFileName, flag='c', protocol=pickle.HIGHEST_PROTOCOL)
            db['__index__'] = self.__moleculeNameList
            if (self.__debug):
                self.__lfh.write("+OePersist.__storeShelve() - Molecule list length  %d\n" % len(self.__moleculeList))
            for molecule in self.__moleculeList:
                ky = molecule['name']
                db[ky] = molecule
            db.close()
            return True
        except Exception as e:
            if (self.__verbose):
                self.__lfh.write("+ERROR- OePersist.__storeShelve() shelve store failed for file %s %s\n" % (dbFileName, str(e)))
            if (self.__debug):
                traceback.print_exc(file=self.__lfh)
            return False

    def __indexShelve(self, dbFileName="my.db"):
        """  Recover the index of molecules in the persistent store.
        """
        try:
            db = shelve.open(dbFileName, flag='r', protocol=pickle.HIGHEST_PROTOCOL)
            self.__moleculeNameList = db['__index__']
            db.close()
            self.__isExpired = False
            return self.__moleculeNameList

        except Exception as e:
            if (self.__verbose):
                self.__lfh.write("+ERROR- OePersist.__indexShelve() shelve index failed for file %s %s\n" % (dbFileName, str(e)))
            if (self.__debug):
                traceback.print_exc(file=self.__lfh)
            return {}

    def __recoverShelve(self, dbFileName="my.db"):
        """  Recover the list of molecules from the persistent store.
        """
        try:
            self.__moleculeList = []
            db = shelve.open(dbFileName, flag='r', protocol=pickle.HIGHEST_PROTOCOL)
            self.__moleculeNameList = db['__index__']
            if (self.__debug):
                self.__lfh.write("+OePersist.__recoverShelve() - Molecule name list %r\n" % self.__moleculeNameList)

            for moleculeName in self.__moleculeNameList:
                self.__moleculeList.append(db[moleculeName])

            db.close()
            return True
        except Exception as e:
            if (self.__verbose):
                self.__lfh.write("+ERROR- OePersist.__recoverShelve() shelve recover failed for file %s %s\n" % (dbFileName, str(e)))
            if (self.__debug):
                traceback.print_exc(file=self.__lfh)
            return False

    def __fetchOneMoleculeShelve(self, dbFileName="my.db", moleculeName=None,):
        """ Recover the molecule from from persistent store corresponding to the
            input molecule name.
        """
        try:
            #
            db = shelve.open(dbFileName, flag='r', protocol=pickle.HIGHEST_PROTOCOL)
            ky = moleculeName
            molecule = db[ky]
            db.close()
            return molecule
        except Exception as e:
            if (self.__verbose):
                self.__lfh.write("+ERROR- OePersist.__fetchOneMoleculeShelve() shelve fetch failed for file %s %s %s\n"
                                 % (dbFileName, moleculeName, str(e)))
            if (self.__debug):
                traceback.print_exc(file=self.__lfh)
            return None

    def __fetchMoleculeShelve(self, moleculeName=None):
        """ Recover the molecule from from persistent store corresponding to the
            input molecule name.

            shelve store must be opened from prior call.
        """
        try:
            #
            ky = moleculeName
            molecule = self.__db[ky]
            return molecule
        except Exception as e:
            if (self.__debug):
                self.__lfh.write("+ERROR- OePersist.__fetchObjectShelve() shelve fetch failed %s %s\n"
                                 % (moleculeName, str(e)))
            if (self.__debug):
                traceback.print_exc(file=self.__lfh)
            #
            return None

    def __updateMoleculeShelve(self, dbFileName="my.db", inputMolecule=None):
        """ Update/append the single input molecule object into the persistent store.
        """
        try:
            #
            moleculeName = inputMolecule['name']
            #
            db = shelve.open(dbFileName, flag='w', protocol=pickle.HIGHEST_PROTOCOL)
            #
            # get molecule index -
            #
            moleculeNameList = db['__index__']
            #
            # Modify or append input molecule --
            #
            if moleculeName in moleculeNameList:
                # idx = moleculeNameList.index(moleculeName)
                pass
            else:
                moleculeNameList.append(moleculeName)
            #
            # store updated molecule index lists -
            #
            db['__index__'] = moleculeNameList
            #
            # ------------------------------------------
            #    Update the object contents --
            #
            ky = moleculeName
            db[ky] = inputMolecule
            db.close()
            return True
        except Exception as e:
            if (self.__verbose):
                self.__lfh.write("+ERROR- OePersist.__updateObjectShelve() update failed for file %s %s %s\n" % (dbFileName, moleculeName, str(e)))
            if (self.__debug):
                traceback.print_exc(file=self.__lfh)
            return False

    def __updateMoleculeListShelve(self, dbFileName="my.db", updateMoleculeList=None):
        """ Update/append the contents of the input molecule list into an existing persistent store.
        """
        try:
            db = shelve.open(dbFileName, flag='w', protocol=pickle.HIGHEST_PROTOCOL)
            #
            # get molecule index -
            #
            moleculeNameList = db['__index__']

            for updateMolecule in updateMoleculeList:
                #
                updateMoleculeName = updateMolecule['name']
                #
                # Modify or append input molecule --
                #
                if updateMoleculeName in moleculeNameList:
                    # idx = moleculeNameList.index(updateMoleculeName)
                    pass
                else:
                    moleculeNameList.append(updateMoleculeName)
                #
                # Store updated molecule index lists -
                #
                db['__index__'] = moleculeNameList

                ky = updateMoleculeName
                db[ky] = updateMolecule
            #
            db.close()
            return True
        except Exception as e:
            if (self.__verbose):
                self.__lfh.write("+ERROR- OePersist.__updateMoleculeListShelve() update failed for store %s %s\n" % (dbFileName, str(e)))
            if (self.__debug):
                traceback.print_exc(file=self.__lfh)
            return False

    # def __attributePart(self, name):
    #     i = name.find(".")
    #     if i == -1:
    #         return None
    #     else:
    #         return name[i + 1:]
