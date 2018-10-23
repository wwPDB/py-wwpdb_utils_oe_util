##
# File:    OeChemCompIoUtils.py
# Author:  jdw
# Date:    04-May-2013
# Version: 0.001
#
# Updates:
#      6-Jun-2016  jdw  general cleanup - use native Python parser on darwin
##
"""
Construct OE mol data structures from chemical component definitions.  Methods are provided
provided to manage the chemical component definition IO

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.01"

import sys
import traceback
import time
import os
import os.path


from oe_util.build.OeBuildMol import OeBuildMol
#
#  Storage adapater can be changed here --
if sys.platform in ['darwin']:
    from mmcif_utils.persist.PdbxPyIoAdapter import PdbxPyIoAdapter as PdbxIoAdapter
else:
    from mmcif_utils.persist.PdbxPyIoAdapter import PdbxPyIoAdapter as PdbxIoAdapter


class OeChemCompIoUtils(object):
    ''' Construct OE mol data structures from chemical component definitions.  Methods are provided
        provided to manage the chemical component definition IO
    '''

    def __init__(self, topCachePath='/data/components/ligand-dict-v3', verbose=True, log=sys.stderr):
        self.__topCachePath = topCachePath
        self.__verbose = verbose
        self.__debug = False
        self.__lfh = log
        #

    def getFromPathList(self, pathList, use3D=True, coordType="model", setTitle=True):
        """ Return a list of OE mols constructed from the input pathList of chemical definitions.
        """
        oemList = []
        try:
            for pth in pathList:
                myReader = PdbxIoAdapter(self.__verbose, self.__lfh)
                ok = myReader.read(pdbxFilePath=pth)
                for container in myReader.getContainerList():
                    oem = OeBuildMol(verbose=self.__verbose, log=self.__lfh)
                    oem.setDebug(self.__debug)
                    oem.set(container.getName(),
                            dcChemCompAtom=container.getObj("chem_comp_atom"),
                            dcChemCompBond=container.getObj("chem_comp_bond"))
                    if use3D:
                        oem.build3D(coordType=coordType, setTitle=setTitle)
                    else:
                        oem.build2D(setTitle=setTitle)
                    #
                    oemList.append(oem)
                    #
                    if self.__debug:
                        self.__lfh.write("+OeChemCompIoUtils.getOeMols() Title              = %s\n" % oem.getTitle())
                        self.__lfh.write("+OeChemCompIoUtils.getOeMols() SMILES (canonical) = %s\n" % oem.getCanSMILES())
                        self.__lfh.write("+OeChemCompIoUtils.getOeMols() SMILES (isomeric)  = %s\n" % oem.getIsoSMILES())
        except:
            if self.__verbose:
                self.__lfh.write("+OeChemCompIoUtils.getOeMols() Failed \n")
                traceback.print_exc(file=self.__lfh)
        return oemList

    def getFromIdList(self, idList, use3D=True, coordType="model", setTitle=True):
        """ Return a list of OE mols constructed from the input Id of chemical definitions.
        """
        pathList = []
        for id in idList:
            if len(id) < 1:
                continue
            idU = str(id).upper()
            if idU.startswith("PRDCC_"):
                hash = idU[-1]
                pth = os.path.join(self.__topCachePath, hash, idU + '.cif')
            else:
                hash = idU[0]
                pth = os.path.join(self.__topCachePath, hash, idU, idU + '.cif')
            pathList.append(pth)
        return self.getFromPathList(pathList, use3D=use3D, coordType=coordType, setTitle=setTitle)

    def __getPathList(self, topPath, pattern='*', excludeDirs=[], recurse=True):
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
