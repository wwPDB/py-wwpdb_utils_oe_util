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


import pickle  # noqa: S403
import shelve  # noqa: S403
import shutil
import sys
import traceback
from typing import Any, Dict, List, Literal, Optional, TextIO, Union

from mmcif_utils.persist.LockFile import LockFile


class OePersist:
    """Persistent storage for instances of OE molecules."""

    def __init__(self, verbose: bool = True, log: TextIO = sys.stderr) -> None:
        self.__verbose = verbose
        self.__debug = True
        self.__lfh = log
        self.__moleculeList: List[Dict[str, Any]] = []

        self.__moleculeNameList: List[str] = []
        #
        # Parameters to tune lock file management --
        self.__timeoutSeconds = 10
        self.__retrySeconds = 0.2
        # placeholder for LockFile object for open/close methods
        self.__lockObj: Optional[LockFile] = None
        self.__isExpired = True
        self.__db: Optional[shelve.Shelf[Any]] = None

    def __reIndex(self) -> None:
        """Rebuild name and type lists from molecule object list."""
        self.__moleculeNameList = []
        for molecule in self.__moleculeList:
            self.__moleculeNameList.append(molecule["name"])
        self.__isExpired = False

    def __getMoleculeIndex(self, name: str) -> Optional[int]:
        try:
            if self.__isExpired:
                self.__reIndex()
            return self.__moleculeNameList.index(name)
        except:  # noqa: E722 pylint: disable=bare-except
            return None

    def setMoleculeList(self, moleculeList: Optional[List[Dict[str, Any]]] = None) -> bool:
        """Initialize molecule data in the internal molecule to be persisted to the data store."""
        try:
            # A None argument is stored as-is; len() below then raises and False is returned.
            self.__moleculeList = moleculeList  # type: ignore[assignment]  # None kept as-is
            self.__isExpired = True
            if self.__debug:
                self.__lfh.write("+OePersist.setMoleculeList() - Container list length %d\n" % len(self.__moleculeList))
            return True
        except Exception as e:  # noqa: BLE001
            if self.__verbose:
                self.__lfh.write("+ERROR- OePersist.setMoleculeList() initialization failed %s\n" % str(e))
            if self.__debug:
                traceback.print_exc(file=self.__lfh)
            return False

    def appendMoleculeList(self, moleculeList: Optional[List[Dict[str, Any]]] = None) -> bool:
        """Append container data to the internal moleculeList to be persisted to the data store."""
        try:
            if moleculeList is None:
                # Same failure extend(None) would raise
                raise TypeError("'NoneType' object is not iterable")
            self.__moleculeList.extend(moleculeList)
            self.__isExpired = True
            if self.__debug:
                self.__lfh.write(
                    "+OePersist.appendMoleculeList() - Container list length %d\n" % len(self.__moleculeList)
                )
            return True
        except Exception as e:  # noqa: BLE001
            if self.__verbose:
                self.__lfh.write("+ERROR- OePersist.appendMoleculeList() initialization failed %s\n" % str(e))
            if self.__debug:
                traceback.print_exc(file=self.__lfh)
            return False

    def getMolecule(self, moleculeName: str) -> Optional[Dict[str, Any]]:
        """Return the molecule from internal MoleculeList corresponding to the input molecule name.

        **Note that this method is NOT directly fetching data from the internal store.  First recover the
        the store to internal storage before calling this method.
        """
        idx = self.__getMoleculeIndex(moleculeName)
        if idx is not None:
            return self.__moleculeList[idx]
        return None

    def getMoleculeList(self) -> List[Dict[str, Any]]:
        """Return the current molecule list."""
        return self.__moleculeList

    def getMoleculeNameList(self) -> List[str]:
        """Return the current list of molecule names."""
        if self.__isExpired:
            self.__reIndex()
        return self.__moleculeNameList

    def open(self, dbFileName: str = "my.db", flag: Literal["r", "w", "c", "n"] = "r") -> bool:
        """Open the persistent store with the input database file name and access mode.

        flag = 'r read-only'
               'c/n new/create'
               'w read/write'
        """
        try:
            self.__lockObj = LockFile(
                dbFileName,
                timeoutSeconds=self.__timeoutSeconds,
                retrySeconds=self.__retrySeconds,
                verbose=self.__verbose,
                log=self.__lfh,
            )
            self.__lockObj.acquire()
            self.__db = shelve.open(dbFileName, flag=flag, protocol=pickle.HIGHEST_PROTOCOL)  # noqa: S301
            return True
        except Exception as e:  # noqa: BLE001
            if self.__lockObj is not None:
                self.__lockObj.release()
            if self.__verbose:
                self.__lfh.write("+ERROR- OePersist.open() write failed for file %s %s\n" % (dbFileName, str(e)))
            if self.__debug:
                traceback.print_exc(file=self.__lfh)
            return False

    def getStoreMoleculeIndex(self) -> List[str]:
        try:
            if self.__db is None:
                # Same failure subscripting None would raise
                raise TypeError("'NoneType' object is not subscriptable")
            storeIndex: List[str] = self.__db["__index__"]
            return storeIndex
        except:  # noqa: E722 pylint: disable=bare-except
            return []

    def close(self) -> bool:
        """Close the persistent store"""
        try:
            if self.__db is None:
                # Same failure None.close() would raise
                raise AttributeError("'NoneType' object has no attribute 'close'")
            self.__db.close()
            if self.__lockObj is not None:
                self.__lockObj.release()
            self.__db = None
            return True
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__lockObj is not None:
                self.__lockObj.release()
            return False

    def moveStore(self, srcDbFilePath: str, dstDbFilePath: str) -> bool:
        """Move source store to destination store with destination locking.

        No locking is performed on the source file.
        """
        with LockFile(
            dstDbFilePath,
            timeoutSeconds=self.__timeoutSeconds,
            retrySeconds=self.__retrySeconds,
            verbose=self.__verbose,
            log=self.__lfh,
        ) as lf:  # noqa: F841 pylint: disable=unused-variable
            retVal = self.__moveStore(srcDbFilePath, dstDbFilePath)
        return retVal

    def store(self, dbFileName: str = "my.db") -> bool:
        """Store the current molecule list in persistent database."""
        with LockFile(
            dbFileName,
            timeoutSeconds=self.__timeoutSeconds,
            retrySeconds=self.__retrySeconds,
            verbose=self.__verbose,
            log=self.__lfh,
        ) as lf:  # noqa: F841 pylint: disable=unused-variable
            retVal = self.__storeShelve(dbFileName)
        return retVal

    def recover(self, dbFileName: str = "my.db") -> bool:
        """Recover the stored state to the current in-memory molecule representation."""
        with LockFile(
            dbFileName,
            timeoutSeconds=self.__timeoutSeconds,
            retrySeconds=self.__retrySeconds,
            verbose=self.__verbose,
            log=self.__lfh,
        ) as lf:  # noqa: F841 pylint: disable=unused-variable
            retVal = self.__recoverShelve(dbFileName)
        return retVal

    def getIndex(self, dbFileName: str = "my.db") -> Union[List[str], Dict[Any, Any]]:
        """Recover the index of the persistent store."""
        with LockFile(
            dbFileName,
            timeoutSeconds=self.__timeoutSeconds,
            retrySeconds=self.__retrySeconds,
            verbose=self.__verbose,
            log=self.__lfh,
        ) as lf:  # noqa: F841  pylint: disable=unused-variable
            retVal = self.__indexShelve(dbFileName)
        return retVal

    def updateOneMolecule(self, inputMolecule: Dict[str, Any], dbFileName: str = "my.db") -> bool:
        """Update or append a molecule to an existing store."""
        with LockFile(
            dbFileName,
            timeoutSeconds=self.__timeoutSeconds,
            retrySeconds=self.__retrySeconds,
            verbose=self.__verbose,
            log=self.__lfh,
        ) as lf:  # noqa: F841 pylint: disable=unused-variable
            retVal = self.__updateMoleculeShelve(dbFileName=dbFileName, inputMolecule=inputMolecule)
        return retVal

    def updateMoleculeList(
        self, dbFileName: str = "my.db", moleculeList: Optional[List[Dict[str, Any]]] = None
    ) -> bool:
        """Update or append the contents of the input molecule list to an existing store."""
        with LockFile(
            dbFileName,
            timeoutSeconds=self.__timeoutSeconds,
            retrySeconds=self.__retrySeconds,
            verbose=self.__verbose,
            log=self.__lfh,
        ) as lf:  # noqa: F841 pylint: disable=unused-variable
            retVal = self.__updateMoleculeListShelve(dbFileName=dbFileName, updateMoleculeList=moleculeList)
        return retVal

    def fetchOneMolecule(
        self, dbFileName: str = "my.db", moleculeName: Optional[str] = None
    ) -> Optional[Dict[str, Any]]:
        """Fetch a single object from a named molecule.  This is atomic operation with respect to the store."""
        with LockFile(
            dbFileName,
            timeoutSeconds=self.__timeoutSeconds,
            retrySeconds=self.__retrySeconds,
            verbose=self.__verbose,
            log=self.__lfh,
        ) as lf:  # noqa: F841  pylint: disable=unused-variable
            retVal = self.__fetchOneMoleculeShelve(dbFileName, moleculeName)
        return retVal

    def fetchMolecule(self, moleculeName: Optional[str] = None) -> Optional[Dict[str, Any]]:
        """Fetch molecule from an open store.

        Use this method to extract multiple molecules from an open store.
        """
        return self.__fetchMoleculeShelve(moleculeName)

    def __moveStore(self, src: str, dst: str) -> bool:
        """Internal method to perform file move."""
        try:
            shutil.move(src, dst)
            return True
        except Exception as e:  # noqa: BLE001
            if self.__verbose:
                self.__lfh.write("+ERROR- OePersist.__moveStore() move failed for file %s %s\n" % (src, str(e)))
            if self.__debug:
                traceback.print_exc(file=self.__lfh)
            return False

    def __storeShelve(self, dbFileName: str = "my.db") -> bool:
        """Create a new persistent store using the internal molecule list."""
        try:
            if self.__isExpired:
                self.__reIndex()
            db = shelve.open(dbFileName, flag="c", protocol=pickle.HIGHEST_PROTOCOL)  # noqa: S301
            db["__index__"] = self.__moleculeNameList
            if self.__debug:
                self.__lfh.write("+OePersist.__storeShelve() - Molecule list length  %d\n" % len(self.__moleculeList))
            for molecule in self.__moleculeList:
                ky = molecule["name"]
                db[ky] = molecule
            db.close()
            return True
        except Exception as e:  # noqa: BLE001
            if self.__verbose:
                self.__lfh.write(
                    "+ERROR- OePersist.__storeShelve() shelve store failed for file %s %s\n" % (dbFileName, str(e))
                )
            if self.__debug:
                traceback.print_exc(file=self.__lfh)
            return False

    def __indexShelve(self, dbFileName: str = "my.db") -> Union[List[str], Dict[Any, Any]]:
        """Recover the index of molecules in the persistent store."""
        try:
            db = shelve.open(dbFileName, flag="r", protocol=pickle.HIGHEST_PROTOCOL)  # noqa: S301
            self.__moleculeNameList = db["__index__"]
            db.close()
            self.__isExpired = False
            return self.__moleculeNameList

        except Exception as e:  # noqa: BLE001
            if self.__verbose:
                self.__lfh.write(
                    "+ERROR- OePersist.__indexShelve() shelve index failed for file %s %s\n" % (dbFileName, str(e))
                )
            if self.__debug:
                traceback.print_exc(file=self.__lfh)
            return {}

    def __recoverShelve(self, dbFileName: str = "my.db") -> bool:
        """Recover the list of molecules from the persistent store."""
        try:
            self.__moleculeList = []
            db = shelve.open(dbFileName, flag="r", protocol=pickle.HIGHEST_PROTOCOL)  # noqa: S301
            self.__moleculeNameList = db["__index__"]
            if self.__debug:
                self.__lfh.write("+OePersist.__recoverShelve() - Molecule name list %r\n" % self.__moleculeNameList)

            for moleculeName in self.__moleculeNameList:
                self.__moleculeList.append(db[moleculeName])

            db.close()
            return True
        except Exception as e:  # noqa: BLE001
            if self.__verbose:
                self.__lfh.write(
                    "+ERROR- OePersist.__recoverShelve() shelve recover failed for file %s %s\n" % (dbFileName, str(e))
                )
            if self.__debug:
                traceback.print_exc(file=self.__lfh)
            return False

    def __fetchOneMoleculeShelve(
        self,
        dbFileName: str = "my.db",
        moleculeName: Optional[str] = None,
    ) -> Optional[Dict[str, Any]]:
        """Recover the molecule from from persistent store corresponding to the
        input molecule name.
        """
        try:
            db = shelve.open(dbFileName, flag="r", protocol=pickle.HIGHEST_PROTOCOL)  # noqa: S301
            ky = moleculeName
            if ky is None:
                # Same failure shelve raises for a None key
                raise AttributeError("'NoneType' object has no attribute 'encode'")
            molecule: Dict[str, Any] = db[ky]
            db.close()
            return molecule
        except Exception as e:  # noqa: BLE001
            if self.__verbose:
                self.__lfh.write(
                    "+ERROR- OePersist.__fetchOneMoleculeShelve() shelve fetch failed for file %s %s %s\n"
                    % (dbFileName, moleculeName, str(e))
                )
            if self.__debug:
                traceback.print_exc(file=self.__lfh)
            return None

    def __fetchMoleculeShelve(self, moleculeName: Optional[str] = None) -> Optional[Dict[str, Any]]:
        """Recover the molecule from from persistent store corresponding to the
        input molecule name.

        shelve store must be opened from prior call.
        """
        try:
            ky = moleculeName
            if self.__db is None:
                # Same failure subscripting None would raise
                raise TypeError("'NoneType' object is not subscriptable")
            if ky is None:
                # Same failure shelve raises for a None key
                raise AttributeError("'NoneType' object has no attribute 'encode'")
            molecule: Dict[str, Any] = self.__db[ky]
            return molecule
        except Exception as e:  # noqa: BLE001
            if self.__debug:
                self.__lfh.write(
                    "+ERROR- OePersist.__fetchObjectShelve() shelve fetch failed %s %s\n" % (moleculeName, str(e))
                )
            if self.__debug:
                traceback.print_exc(file=self.__lfh)
            return None

    def __updateMoleculeShelve(self, dbFileName: str = "my.db", inputMolecule: Optional[Dict[str, Any]] = None) -> bool:
        """Update/append the single input molecule object into the persistent store."""
        try:
            if inputMolecule is None:
                # Same failure subscripting None would raise
                raise TypeError("'NoneType' object is not subscriptable")
            moleculeName = inputMolecule["name"]
            db = shelve.open(dbFileName, flag="w", protocol=pickle.HIGHEST_PROTOCOL)  # noqa: S301
            #
            # get molecule index -
            #
            moleculeNameList = db["__index__"]
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
            db["__index__"] = moleculeNameList
            #
            # ------------------------------------------
            #    Update the object contents --
            #
            ky = moleculeName
            db[ky] = inputMolecule
            db.close()
            return True
        except Exception as e:  # noqa: BLE001
            if self.__verbose:
                self.__lfh.write(
                    "+ERROR- OePersist.__updateObjectShelve() update failed for file %s %s %s\n"
                    % (dbFileName, moleculeName, str(e))
                )
            if self.__debug:
                traceback.print_exc(file=self.__lfh)
            return False

    def __updateMoleculeListShelve(
        self, dbFileName: str = "my.db", updateMoleculeList: Optional[List[Dict[str, Any]]] = None
    ) -> bool:
        """Update/append the contents of the input molecule list into an existing persistent store."""
        try:
            db = shelve.open(dbFileName, flag="w", protocol=pickle.HIGHEST_PROTOCOL)  # noqa: S301
            #
            # get molecule index -
            #
            moleculeNameList = db["__index__"]

            if updateMoleculeList is None:
                # Same failure iterating None would raise
                raise TypeError("'NoneType' object is not iterable")
            for updateMolecule in updateMoleculeList:
                updateMoleculeName = updateMolecule["name"]
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
                db["__index__"] = moleculeNameList

                ky = updateMoleculeName
                db[ky] = updateMolecule
            db.close()
            return True
        except Exception as e:  # noqa: BLE001
            if self.__verbose:
                self.__lfh.write(
                    "+ERROR- OePersist.__updateMoleculeListShelve() update failed for store %s %s\n"
                    % (dbFileName, str(e))
                )
            if self.__debug:
                traceback.print_exc(file=self.__lfh)
            return False

    # def __attributePart(self, name):
    #     i = name.find(".")
    #     if i == -1:
    #         return None
    #     else:
    #         return name[i + 1:]
