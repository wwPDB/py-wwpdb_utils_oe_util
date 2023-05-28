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

import os
import os.path
import sys
import traceback

from wwpdb.utils.oe_util.build.OeBuildMol import OeBuildMol

#
#  Storage adapater can be changed here --
if sys.platform in ['darwin']:
    from mmcif_utils.persist.PdbxPyIoAdapter import \
        PdbxPyIoAdapter as PdbxIoAdapter
else:
    from mmcif_utils.persist.PdbxPyIoAdapter import \
        PdbxPyIoAdapter as PdbxIoAdapter


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
                ok = myReader.read(pdbxFilePath=pth)  # noqa: F841 pylint: disable=unused-variable
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
        except Exception as e:
            if self.__verbose:
                self.__lfh.write("+OeChemCompIoUtils.getOeMols() Failed %s\n" % str(e))
                traceback.print_exc(file=self.__lfh)
        return oemList

    def getFromIdList(self, idList, use3D=True, coordType="model", setTitle=True):
        """ Return a list of OE mols constructed from the input Id of chemical definitions.
        """
        pathList = []
        for ide in idList:
            if len(ide) < 1:
                continue
            idU = str(ide).upper()
            if idU.startswith("PRDCC_"):
                hashd = idU[-1]
                pth = os.path.join(self.__topCachePath, hashd, idU + '.cif')
            else:
                hashd = self.__getCcdHash(idU)
                pth = os.path.join(self.__topCachePath, hashd, idU, idU + '.cif')
            pathList.append(pth)
        return self.getFromPathList(pathList, use3D=use3D, coordType=coordType, setTitle=setTitle)

    def __getCcdHash(self, idCode):
        """Returns the hash code for a CCD id.  Currently first letter"""
        if not idCode:
            return None

        if len(idCode) > 3:
            hash_key = idCode.upper()[-2:]
        else:
            hash_key = idCode.upper()[0]

        return hash_key
