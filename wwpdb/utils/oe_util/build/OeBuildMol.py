##
# File:    OeBuildMol.py
# Author:  jdw
# Date:    25-Sep-2011
# Version: 0.001
#
# Updates:
#  1-Oct-2011 jdw Add 2D builder, CIP Stereo and SMILES
# 18-Feb-2012 jdw add 3D builder with aromatic and CIP perception.
# 20-Feb-2012 jdw Add serialization and deserialization
# 20-Feb-2012 jdw Add CIP perception for 3D builder
# 23-Feb-2012 jdw Moved to oe_util/build path.
#                 Revised to use chemical component iterator classes
#                 and to support flat file and persist storage models.
#  5-May-2014 jdw add method build molecule from input SMILES importSmiles(smiles)
#  6-Jun-2016 jdw general cleanup and platform selection parsing -  use native python on darwin
#  8-Jun-2016 jdw add method to instantiate using an existing OEMOL and method to return ccId.
# 31-Aug-2020 ZF  in importSmiles() function, added checking OEParseSmiles() return status
#
##
"""
Classes to build OE molecule objects from chemical component definition data.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.01"


import sys
import traceback

from mmcif.io.IoAdapterCore import IoAdapterCore as IoAdapter
# from mmcif.api.PdbxContainers import *
from mmcif_utils.chemcomp.PdbxChemComp import PdbxChemCompConstants
from openeye.oechem import (OE3DToInternalStereo, OEAddExplicitHydrogens,
                            OEAroModelOpenEye, OEAssignAromaticFlags,
                            OECIPAtomStereo_R, OECIPAtomStereo_S,
                            OECIPBondStereo_E, OECIPBondStereo_Z,
                            OECreateCanSmiString, OECreateInChI,
                            OECreateInChIKey, OECreateIsoSmiString,
                            OEFindRingAtomsAndBonds, OEFloatArray,
                            OEFormat_OEB, OEGetAtomicSymbol, OEGraphMol, OEMol,
                            OEMolecularFormula, OEParseSmiles,
                            OEPerceiveChiral, OEPerceiveCIPStereo,
                            OEReadMolecule, OESetCIPStereo,
                            OESuppressHydrogens, OETriposAtomNames,
                            OEWriteMolecule, oemolistream, oemolostream)
from wwpdb.utils.cc_dict_util.persist.PdbxChemCompPersist import (
    PdbxChemCompAtomIt, PdbxChemCompBondIt)


class OeBuildMol(object):
    ''' Utility methods for constructing OEGraphMols from chemical component definition objects.
    '''

    def __init__(self, verbose=True, log=sys.stderr):
        self.__verbose = verbose
        self.__debug = False
        self.__lfh = log
        #
        # File system path to the chemical component dictionary definitions in (CVS checkout organization)
        #
        # Internal storage for current OE molecule
        self.__oeMol = None
        #
        # Component identifier
        #
        self.__ccId = None
        #
        # dictionary of element counts eD[atno]=count
        self.__eD = {}
        #
        # Source data categories objects from chemical component definitions.
        self.__dcChemCompAtom = None
        self.__dcChemCompBond = None
        #
        self.__molXyzL = []

    def setDebug(self, flag):
        self.__debug = flag

    def setChemCompPath(self, ccPath):
        try:
            myReader = IoAdapter(self.__verbose, self.__lfh)
            cL = myReader.readFile(ccPath)
            self.__ccId = cL[0].getName()
            self.__dcChemCompAtom = cL[0].getObj("chem_comp_atom")
            self.__dcChemCompBond = cL[0].getObj("chem_comp_bond")
            return self.__ccId
        except Exception as e:
            self.__lfh.write("OeBuildMol(setChemCompPath) Fails for %s %s\n" % (ccPath, str(e)))
            traceback.print_exc(file=self.__lfh)
        return None

    def setOeMol(self, inpOeMol, ccId):
        """  Load this object with an existing oeMOL()
        """
        self.__clear()
        self.__oeMol = OEMol(inpOeMol)
        self.__ccId = ccId
        self.getElementCounts()

    def set(self, ccId, dcChemCompAtom=None, dcChemCompBond=None):
        """  Assign source data categories -
        """
        self.__ccId = ccId
        self.__dcChemCompAtom = dcChemCompAtom
        self.__dcChemCompBond = dcChemCompBond

    def __clear(self):
        self.__oeMol = None
        self.__eD = {}

    def serialize(self):
        """ Create a string representing the content of the current OE molecule.   This
            serialization uses the OE internal binary format.
        """
        oms = oemolostream()
        oms.SetFormat(OEFormat_OEB)
        oms.openstring()
        OEWriteMolecule(oms, self.__oeMol)
        if (self.__debug):
            self.__lfh.write("OeBuildMol(Serialize) SMILES %s\n" % OECreateCanSmiString(self.__oeMol))
            self.__lfh.write("OeBuildMol(Serialize) atoms = %d\n" % self.__oeMol.NumAtoms())
        return oms.GetString()

    def deserialize(self, oeS):
        """ Reconstruct an OE molecule from the input string serialization (OE binary).

            The deserialized molecule is used to initialize the internal OE molecule
            within this object.

            Returns True for success or False otherwise.
        """
        self.__clear()
        ims = oemolistream()
        ims.SetFormat(OEFormat_OEB)
        ims.openstring(oeS)

        nmol = 0
        mList = []
        # for mol in ims.GetOEGraphMols():
        for mol in ims.GetOEMols():
            if (self.__debug):
                self.__lfh.write("OeBuildMol(deserialize) SMILES %s\n" % OECreateCanSmiString(mol))
                self.__lfh.write("OeBuildMol(deserialize) title  %s\n" % mol.GetTitle())
                self.__lfh.write("OeBuildMol(deserialize) atoms  %d\n" % mol.NumAtoms())
            # mList.append(OEGraphMol(mol))
            mList.append(OEMol(mol))
            nmol += 1
        #
        if nmol >= 1:
            self.__oeMol = mList[0]
            self.__ccId = self.__oeMol.GetTitle()
            #
            if (self.__debug):
                self.__lfh.write("OeBuildMol(deserialize) mols  %d\n" % nmol)
                self.__lfh.write("OeBuildMol(deserialize) id %s\n" % self.__ccId)
                self.__lfh.write("OeBuildMol(deserialize) atoms  %d\n" % self.__oeMol.NumAtoms())
            return True
        else:
            return False

    def simpleAtomNames(self):
        """
        """
        for atom in self.__oeMol.GetAtoms():
            atom.SetIntType(atom.GetAtomicNum())
            atom.SetType(OEGetAtomicSymbol(atom.GetAtomicNum()))
        OETriposAtomNames(self.__oeMol)

    def getElementCounts(self):
        """ Get the dictionary of element counts (eg. eD[iAtNo]=iCount).
        """
        if len(self.__eD) == 0:
            # calculate from current oeMol
            try:
                self.__eD = {}
                for atom in self.__oeMol.GetAtoms():
                    atNo = atom.GetAtomicNum()
                    if atNo not in self.__eD:
                        self.__eD[atNo] = 1
                    else:
                        self.__eD[atNo] += 1
            except:  # noqa: E722 pylint: disable=bare-except
                pass

        return self.__eD

    def build3D(self, coordType="model", setTitle=True):
        try:
            self.__build3D(coordType=coordType, setTitle=setTitle)
            return True
        except Exception as e:
            self.__lfh.write("OeBuildMol(build3D) Failing %s\n" % str(e))
            traceback.print_exc(file=self.__lfh)
        return False

    def __build3D(self, coordType="model", setTitle=True):
        """ Build OE molecule using 3D coordinates and OE stereo perception.
        """
        self.__clear()
        # self.__oeMol=OEGraphMol()
        self.__oeMol = OEMol()
        #
        if setTitle:
            self.__oeMol.SetTitle(self.__ccId)
        aL = []

        # Atom index dictionary
        aD = {}
        i = 1

        atomIt = PdbxChemCompAtomIt(self.__dcChemCompAtom, self.__verbose, self.__lfh)
        for ccAt in atomIt:

            atName = ccAt.getName()
            aD[atName] = i
            i += 1
            atNo = ccAt.getAtNo()
            if atNo not in self.__eD:
                self.__eD[atNo] = 1
            else:
                self.__eD[atNo] += 1

            atType = ccAt.getType()
            fc = ccAt.getFormalCharge()
            chFlag = ccAt.isChiral()
            # arFlag = ccAt.isAromatic()
            isotope = ccAt.getIsotope()
            leavingAtom = ccAt.getLeavingAtomFlag()

            oeAt = self.__oeMol.NewAtom(atNo)
            oeAt.SetName(atName)
            oeAt.SetFormalCharge(fc)
            oeAt.SetStringData("pdbx_leaving_atom_flag", leavingAtom)
            oeAt.SetChiral(chFlag)
            oeAt.SetIsotope(isotope)

            # oeAt.SetAromatic(arFlag)
            # if chFlag:
            #    st=ccAt.getCIPStereo()
            #    if st == 'S' or st == 'R':
            #        oeAt.SetStringData("StereoInfo",st)
            # if (self.__debug):
            #    self.__lfh.write("Atom - %s type %s atno %d isotope %d fc %d chFlag %r\n" % (atName,atType,atNo,isotope,fc,chFlag))

            if ((coordType == 'model') and ccAt.hasModelCoordinates()):
                cTup = ccAt.getModelCoordinates()
                # if (self.__verbose):
                #    self.__lfh.write("CC %s Atom - %s cTup %r\n" % (self.__ccId,atName,cTup))
                self.__oeMol.SetCoords(oeAt, cTup)
            elif ((coordType == 'ideal') and ccAt.hasIdealCoordinates()):
                cTup = ccAt.getIdealCoordinates()
                self.__oeMol.SetCoords(oeAt, cTup)
            else:
                pass

            if (self.__debug):
                self.__lfh.write("Atom - %s type %s atno %d isotope %d fc %d (xyz) %r\n" % (atName, atType, atNo, isotope, fc, cTup))
            aL.append(oeAt)

        bondIt = PdbxChemCompBondIt(self.__dcChemCompBond, self.__verbose, self.__lfh)
        for ccBnd in bondIt:
            (at1, at2) = ccBnd.getBond()
            iat1 = aD[at1] - 1
            iat2 = aD[at2] - 1
            iType = ccBnd.getIntegerType()
            #
            arFlag = ccBnd.isAromatic()  # noqa: F841 pylint: disable=unused-variable
            if (self.__debug):
                self.__lfh.write(" %s %d -- %s %d (%d)\n" % (at1, iat1, at2, iat2, iType))

            oeBnd = self.__oeMol.NewBond(aL[iat1], aL[iat2], iType)  # noqa: F841 pylint: disable=unused-variable

            # oeBnd.SetAromatic(arFlag)
            # if arFlag:
            #    oeBnd.SetIntType(5)
            # st=ccBnd.getStereo()
            # if st == 'E' or st =='Z':
            #    oeBnd.SetStringData("StereoInfo",st)

        #
        # run standard perceptions --
        #
        self.__oeMol.SetDimension(3)
        OE3DToInternalStereo(self.__oeMol)
        OEFindRingAtomsAndBonds(self.__oeMol)
        # Other aromatic models: OEAroModelMDL or OEAroModelDaylight
        OEAssignAromaticFlags(self.__oeMol, OEAroModelOpenEye)
        self.updateCIPStereoOE()

    def updatePerceptions3D(self):
        self.__oeMol.SetDimension(3)
        OE3DToInternalStereo(self.__oeMol)
        OEFindRingAtomsAndBonds(self.__oeMol)
        # Other aromatic models: OEAroModelMDL or OEAroModelDaylight
        OEAssignAromaticFlags(self.__oeMol, OEAroModelOpenEye)
        self.updateCIPStereoOE()

    def updateCIPStereoOE(self):
        """ OE perception of CIP stereo -
        """
        for atom in self.__oeMol.GetAtoms():
            OEPerceiveCIPStereo(self.__oeMol, atom)

        for bond in self.__oeMol.GetBonds():
            if (bond.GetOrder() == 2):
                OEPerceiveCIPStereo(self.__oeMol, bond)

    def build2D(self, setTitle=True):  # pylint: disable=unused-argument
        try:
            self.__build2D(setTitle=True)
            return True
        except:  # noqa: E722 pylint: disable=bare-except
            return False

    def __build2D(self, setTitle=True):
        """  Build molecule using existing assignments of chemical information in the CC definition.
        """
        self.__clear()
        self.__oeMol = OEGraphMol()
        if setTitle:
            self.__oeMol.SetTitle(self.__ccId)
        aL = []
        i = 1
        # Atom index dictionary
        aD = {}

        atomIt = PdbxChemCompAtomIt(self.__dcChemCompAtom, self.__verbose, self.__lfh)
        for ccAt in atomIt:
            atName = ccAt.getName()
            aD[atName] = i
            i += 1
            atNo = ccAt.getAtNo()
            if atNo not in self.__eD:
                self.__eD[atNo] = 1
            else:
                self.__eD[atNo] += 1
            atType = ccAt.getType()
            fc = ccAt.getFormalCharge()
            chFlag = ccAt.isChiral()
            arFlag = ccAt.isAromatic()
            isotope = ccAt.getIsotope()
            leavingAtom = ccAt.getLeavingAtomFlag()

            oeAt = self.__oeMol.NewAtom(atNo)
            oeAt.SetName(atName)
            oeAt.SetFormalCharge(fc)
            oeAt.SetStringData("pdbx_leaving_atom_flag", leavingAtom)
            oeAt.SetChiral(chFlag)
            oeAt.SetIsotope(isotope)
            oeAt.SetAromatic(arFlag)
            if chFlag:
                st = ccAt.getCIPStereo()
                if st == 'S' or st == 'R':
                    oeAt.SetStringData("StereoInfo", st)
            if (self.__debug):
                self.__lfh.write("Atom - %s type %s atno %d isotope %d fc %d chFlag %r\n" % (atName, atType, atNo, isotope, fc, chFlag))
            aL.append(oeAt)

        bondIt = PdbxChemCompBondIt(self.__dcChemCompBond, self.__verbose, self.__lfh)
        for ccBnd in bondIt:
            (at1, at2) = ccBnd.getBond()
            iat1 = aD[at1] - 1
            iat2 = aD[at2] - 1
            iType = ccBnd.getIntegerType()
            arFlag = ccBnd.isAromatic()
            if (self.__debug):
                self.__lfh.write(" %s %d -- %s %d (%d)\n" % (at1, iat1, at2, iat2, iType))

            oeBnd = self.__oeMol.NewBond(aL[iat1], aL[iat2], iType)
            oeBnd.SetAromatic(arFlag)
            if arFlag:
                oeBnd.SetIntType(5)
            st = ccBnd.getStereo()
            if st == 'E' or st == 'Z':
                oeBnd.SetStringData("StereoInfo", st)

        #
        # run standard perceptions --
        OEFindRingAtomsAndBonds(self.__oeMol)
        OEPerceiveChiral(self.__oeMol)

        for oeAt in self.__oeMol.GetAtoms():
            st = oeAt.GetStringData("StereoInfo")
            if st == 'R':
                OESetCIPStereo(self.__oeMol, oeAt, OECIPAtomStereo_R)
            elif st == 'S':
                OESetCIPStereo(self.__oeMol, oeAt, OECIPAtomStereo_S)

        for oeBnd in self.__oeMol.GetBonds():
            st = oeBnd.GetStringData("StereoInfo")
            if st == 'E':
                OESetCIPStereo(self.__oeMol, oeBnd, OECIPBondStereo_E)
            elif st == 'Z':
                OESetCIPStereo(self.__oeMol, oeBnd, OECIPBondStereo_Z)
        if (self.__debug):
            for ii, atm in enumerate(self.__oeMol.GetAtoms()):
                self.__lfh.write("OeBuildMol.build2d - atom  %d %s\n" % (ii, atm.GetName()))

    def importFile(self, filePath, type='2D'):  # pylint: disable=redefined-builtin
        """  Contruct a OEGraphMol using the content of the input file.  The input
             file must have a file extension recognized by the OE toolkit (e.g. .sdf)
        """
        ifs = oemolistream()
        if not ifs.open(filePath):
            return False
        #
        # self.__oeMol = OEGraphMol()
        self.__oeMol = OEMol()
        OEReadMolecule(ifs, self.__oeMol)
        #        OETriposAtomNames(self.__oeMol)
        if type == '2D':
            # run standard perceptions --
            OEFindRingAtomsAndBonds(self.__oeMol)
            OEPerceiveChiral(self.__oeMol)

            for oeAt in self.__oeMol.GetAtoms():
                st = oeAt.GetStringData("StereoInfo")
                if st == 'R':
                    OESetCIPStereo(self.__oeMol, oeAt, OECIPAtomStereo_R)
                elif st == 'S':
                    OESetCIPStereo(self.__oeMol, oeAt, OECIPAtomStereo_S)

            for oeBnd in self.__oeMol.GetBonds():
                st = oeBnd.GetStringData("StereoInfo")
                if st == 'E':
                    OESetCIPStereo(self.__oeMol, oeBnd, OECIPBondStereo_E)
                elif st == 'Z':
                    OESetCIPStereo(self.__oeMol, oeBnd, OECIPBondStereo_Z)
        elif type == '3D':
            # run standard perceptions --
            #
            self.__oeMol.SetDimension(3)
            OE3DToInternalStereo(self.__oeMol)
            OEFindRingAtomsAndBonds(self.__oeMol)
            # Other aromatic models: OEAroModelMDL or OEAroModelDaylight
            OEAssignAromaticFlags(self.__oeMol, OEAroModelOpenEye)
            self.updateCIPStereoOE()
            OEAddExplicitHydrogens(self.__oeMol)

        self.__molXyzL = []
        aC = {}
        for ii, atm in enumerate(self.__oeMol.GetAtoms()):
            iAtNum = atm.GetAtomicNum()
            if iAtNum in aC:
                aC[iAtNum] += 1
            else:
                aC[iAtNum] = 1
            # Less than idea - should have an API
            atName = PdbxChemCompConstants._periodicTable[iAtNum - 1] + str(aC[iAtNum])  # pylint: disable=protected-access
            atm.SetName(atName)
            #
            xyzL = OEFloatArray(3)
            self.__oeMol.GetCoords(atm, xyzL)
            self.__molXyzL.append((ii, atm.GetIdx(), atm.GetAtomicNum(), atm.GetName(), atm.GetType(), xyzL[0], xyzL[1], xyzL[2]))

        return True

    def importSmiles(self, smiles):
        """  Contruct a OEGraphMol using the input descriptor.
        """
        self.__oeMol = OEGraphMol()
        if OEParseSmiles(self.__oeMol, smiles):
            OEFindRingAtomsAndBonds(self.__oeMol)
            OEPerceiveChiral(self.__oeMol)
            return True
        #
        return False

    def getGraphMolSuppressH(self):
        """ Return the current constructed OE molecule with hydrogens suppressed.
        """
        # OESuppressHydrogens(self.__oeMol, retainPolar=False,retainStereo=True,retainIsotope=True)
        OESuppressHydrogens(self.__oeMol)
        return self.__oeMol

    def getMol(self):
        """ Return the current constructed OE molecule.
        """
        return OEMol(self.__oeMol)

    def getCanSMILES(self):
        """ Return the cannonical SMILES string derived from the current OD molecule.
        """
        return OECreateCanSmiString(self.__oeMol)

    def getIsoSMILES(self):
        """ Return the cannonical stereo SMILES string derived from the current OE molecule.
        """
        return OECreateIsoSmiString(self.__oeMol)

    def getFormula(self):
        """ Return the Hill order formulat  derived from the current OE molecule.
        """
        return OEMolecularFormula(self.__oeMol)

    def getInChIKey(self):
        """ Return the InChI key derived from the current OE molecule.
        """
        return OECreateInChIKey(self.__oeMol)

    def getInChI(self):
        """ Return the InChI string derived from the current OE molecule.
        """
        return OECreateInChI(self.__oeMol)

    def getTitle(self):
        """ Return the title assigned to the current OE molecule
        """
        return self.__oeMol.GetTitle()

    def getCcId(self):
        """ Return the CC id of this object -
        """
        return self.__ccId

    def getCoords(self):
        """  Return coordinate list if a 3D molecule is built -- otherwise an empty list --

        """
        return self.__molXyzL
