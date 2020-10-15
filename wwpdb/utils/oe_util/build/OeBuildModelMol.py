##
# File:    OeBuildModelMol.py
# Author:  jdw
# Date:    28-Nov-2014
# Version: 0.001
#
# Updates:
#   6-Jun-2016  jdw  general cleanup
#
##
"""
Classes to build OE molecule objects from chemical component model instances.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.01"


import sys
import traceback

from mmcif_utils.chemcomp.PdbxChemCompModel import (PdbxChemCompModel,
                                                    PdbxChemCompModelAtom,
                                                    PdbxChemCompModelBond)
from mmcif_utils.chemcomp.PdbxChemCompModelIo import PdbxChemCompModelIo
from openeye.oechem import (OE3DToInternalStereo, OEAroModelOpenEye,
                            OEAssignAromaticFlags, OECreateCanSmiString,
                            OECreateInChI, OECreateInChIKey,
                            OECreateIsoSmiString, OEFindRingAtomsAndBonds,
                            OEFormat_OEB, OEMol, OEMolecularFormula,
                            OEPerceiveCIPStereo, OESuppressHydrogens,
                            OEWriteMolecule, oemolistream, oemolostream)


class OeBuildModelMol(object):
    ''' Utility methods for constructing OEMols from chemical component model instances.
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
        self.__modelId = None
        #
        # dictionary of element counts eD[atno]=count
        self.__eD = {}
        #
        # Source data categories objects from chemical component definitions.
        #
        self.__ccAtomDL = None
        self.__ccBondDL = None
        self.__molXyzL = []

    def setDebug(self, flag):
        self.__debug = flag

    def setChemCompModelPath(self, modelPath):
        try:
            ccm = PdbxChemCompModelIo(verbose=self.__verbose, log=self.__lfh)
            ccm.setFilePath(modelPath)
            ccDL = ccm.getAttribDictList(catName='pdbx_chem_comp_model')
            cmp = PdbxChemCompModel(ccDL[0], self.__verbose, self.__lfh)
            self.__modelId = cmp.getId()
            self.__ccId = cmp.getCompId()

            self.__ccAtomDL = ccm.getAttribDictList(catName='pdbx_chem_comp_model_atom')
            self.__ccBondDL = ccm.getAttribDictList(catName='pdbx_chem_comp_model_bond')
            return self.__modelId
        except Exception as e:
            self.__lfh.write("OeBuildModelMol(setChemCompModelPath) Fails for %s %s\n" % (modelPath, str(e)))
            traceback.print_exc(file=self.__lfh)
        return None

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
            self.__lfh.write("OeBuildModelMol(Serialize) SMILES %s\n" % OECreateCanSmiString(self.__oeMol))
            self.__lfh.write("OeBuildModelMol(Serialize) atoms = %d\n" % self.__oeMol.NumAtoms())
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
                self.__lfh.write("OeBuildModelMol(deserialize) SMILES %s\n" % OECreateCanSmiString(mol))
                self.__lfh.write("OeBuildModelMol(deserialize) title  %s\n" % mol.GetTitle())
                self.__lfh.write("OeBuildModelMol(deserialize) atoms  %d\n" % mol.NumAtoms())
            # mList.append(OEGraphMol(mol))
            mList.append(OEMol(mol))
            nmol += 1
        #
        if nmol >= 1:
            self.__oeMol = mList[0]
            self.__ccId = self.__oeMol.GetTitle()
            #
            if (self.__debug):
                self.__lfh.write("OeBuildModelMol(deserialize) mols  %d\n" % nmol)
                self.__lfh.write("OeBuildModelMol(deserialize) id %s\n" % self.__ccId)
                self.__lfh.write("OeBuildModelMol(deserialize) atoms  %d\n" % self.__oeMol.NumAtoms())
            return True
        else:
            return False

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

    def build3D(self):
        try:
            self.__build3D()
            return True
        except:  # noqa: E722 pylint: disable=bare-except
            return False

    def __build3D(self):
        """ Build OE molecule using model instance 3D coordinates and OE stereo perception.
        """
        self.__clear()
        self.__oeMol = OEMol()
        #
        self.__oeMol.SetTitle(self.__modelId)
        aL = []

        # Atom index dictionary
        aD = {}
        i = 1

        for d in self.__ccAtomDL:
            ccAt = PdbxChemCompModelAtom(d, self.__verbose, self.__lfh)
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
            # chFlag=ccAt.isChiral()
            # arFlag=ccAt.isAromatic()
            isotope = ccAt.getIsotope()
            oeAt = self.__oeMol.NewAtom(atNo)
            oeAt.SetName(atName)
            oeAt.SetFormalCharge(fc)
            # oeAt.SetChiral(chFlag)
            oeAt.SetIsotope(isotope)

            # oeAt.SetAromatic(arFlag)
            # if chFlag:
            #    st=ccAt.getCIPStereo()
            #    if st == 'S' or st == 'R':
            #        oeAt.SetStringData("StereoInfo",st)
            # if (self.__debug):
            #    self.__lfh.write("Atom - %s type %s atno %d isotope %d fc %d chFlag %r\n" % (atName,atType,atNo,isotope,fc,chFlag))

            cTup = ccAt.getModelCoordinates()
            # if (self.__verbose):
            #    self.__lfh.write("CC %s Atom - %s cTup %r\n" % (self.__ccId,atName,cTup))
            self.__oeMol.SetCoords(oeAt, cTup)
            if (self.__debug):
                self.__lfh.write("Atom - %s type %s atno %d isotope %d fc %d (xyz) %r\n" % (atName, atType, atNo, isotope, fc, cTup))
            aL.append(oeAt)

        for d in self.__ccBondDL:
            ccBnd = PdbxChemCompModelBond(d, self.__verbose, self.__lfh)
            (at1, at2) = ccBnd.getBond()
            iat1 = aD[at1] - 1
            iat2 = aD[at2] - 1
            iType = ccBnd.getIntegerType()
            #
            # arFlag=ccBnd.isAromatic()
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

    def updateCIPStereoOE(self):
        """ OE perception of CIP stereo -
        """
        for atom in self.__oeMol.GetAtoms():
            OEPerceiveCIPStereo(self.__oeMol, atom)

        for bond in self.__oeMol.GetBonds():
            if (bond.GetOrder() == 2):
                OEPerceiveCIPStereo(self.__oeMol, bond)

    def getGraphMolSuppressH(self):
        """ Return the current constructed OE molecule with hydrogens suppressed.
        """
        # OESuppressHydrogens(self.__oeMol, retainPolar=False,retainStereo=True,retainIsotope=True)
        OESuppressHydrogens(self.__oeMol)
        return self.__oeMol

    def getMol(self):
        """ Return the current constructed OE molecule.
        """
        return self.__oeMol

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

    def getCoords(self):
        """  Return coordinate list if a 3D molecule is built -- otherwise an empty list --

        """
        return self.__molXyzL
