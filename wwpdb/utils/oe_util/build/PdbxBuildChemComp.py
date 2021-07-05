##
# File:    PdbxBuildChemComp.py
# Author:  jdw
# Date:    9-Jun-2016
# Version: 0.001
#
# Updates:
#
#
##
"""
Utilities to build chemical component definitions from OE molecule objects.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.01"

import sys
import time
import traceback

from mmcif_utils.style.ChemCompCategoryStyle import ChemCompCategoryStyle
from mmcif_utils.style.PdbxStyleIoUtil import PdbxStyleIoUtil
from openeye.oechem import (OECalculateMolecularWeight,
                            OECIPAtomStereo_NotStereo, OECIPAtomStereo_R,
                            OECIPAtomStereo_S, OECIPAtomStereo_UnspecStereo,
                            OECIPBondStereo_E, OECIPBondStereo_NotStereo,
                            OECIPBondStereo_UnspecStereo, OECIPBondStereo_Z,
                            OECreateCanSmiString, OECreateInChI,
                            OECreateInChIKey, OECreateIsoSmiString,
                            OEGetAtomicSymbol, OEMol, OEMolecularFormula,
                            OEPerceiveCIPStereo, OEWriteConstMolecule,
                            OEWriteMolecule, oemolostream)


class PdbxBuildChemComp(PdbxStyleIoUtil):
    ''' Methods to build chemical component definitions from OE molecule objects.

    '''

    def __init__(self, verbose=True, log=sys.stderr):
        super(PdbxBuildChemComp, self).__init__(styleObject=ChemCompCategoryStyle(), verbose=verbose, log=log)

        # self.__verbose = verbose
        # self.__debug = False
        self.__lfh = log
        #

    def getCategory(self, catName='chem_comp'):
        return self.getItemDictList(catName)

    def getChemCompCategory(self, catName='chem_comp'):
        return self.getItemDictList(catName)

    def getChemCompDict(self):
        return self.getItemDictList(catName='chem_comp')

    def getBondList(self):
        return self.getRowDataList(catName='chem_comp_bond')

    def complyStyle(self):
        return self.testStyleComplete(self.__lfh)

    def setBlock(self, blockId):
        return self.setContainer(containerName=blockId)

    def newBlock(self, blockId):
        return self.newContainer(containerName=blockId)

    def update(self, catName, attributeName, value, iRow=0):
        return self.updateAttribute(catName, attributeName, value, iRow=iRow)

    def write(self, filePath):
        return self.writeFile(filePath)

    #
    def setOeMol(self, oeMol, ccId, name=None, missingModelXyz=True):
        ccIdU = str(ccId).strip().upper()
        self.newContainer(containerName=ccIdU, overWrite=True)
        self.setContainer(containerName=ccIdU)
        #
        self.newCategory('chem_comp', container=None, overWrite=True)
        #
        formula = OEMolecularFormula(oeMol)
        fCharge = self.__getFormalCharge(oeMol)
        fWeight = OECalculateMolecularWeight(oeMol)
        # JDW
        rowD = self.__makeChemCompCategory(ccIdU, name=name, formula=formula, charge=fCharge, fW=fWeight, site='RCSB', missingModelXyz=missingModelXyz)
        self.updateRowByAttribute(rowD, 'chem_comp', iRow=0)
        #
        self.newCategory('chem_comp_atom', container=None, overWrite=True)
        rowDL = self.__makeChemCompAtomCategory(ccIdU, oeMol)
        for ii, rowD in enumerate(rowDL):
            self.updateRowByAttribute(rowD, 'chem_comp_atom', iRow=ii)
        #
        self.newCategory('chem_comp_bond', container=None, overWrite=True)
        rowDL = self.__makeChemCompBondCategory(ccIdU, oeMol)
        for ii, rowD in enumerate(rowDL):
            self.updateRowByAttribute(rowD, 'chem_comp_bond', iRow=ii)
        #
        self.newCategory('pdbx_chem_comp_descriptor', container=None, overWrite=True)
        rowDL = self.__makeChemCompDescriptorCategory(ccIdU, oeMol)
        for ii, rowD in enumerate(rowDL):
            self.updateRowByAttribute(rowD, 'pdbx_chem_comp_descriptor', iRow=ii)
        #
        self.newCategory('pdbx_chem_comp_audit', container=None, overWrite=True)
        rowD = self.__makeChemCompAuditRow(ccIdU)
        self.updateRowByAttribute(rowD, 'pdbx_chem_comp_audit', iRow=0)
        #
        #

    def writeOther(self, oeMol, filePath, title='None', constantMol=False):
        try:
            ofs = oemolostream()
            ofs.open(filePath)
            myMol = OEMol(oeMol)
            myMol.SetTitle(title)
            self.__lfh.write("+PdbxBuildChemComp.writeOther writing %s title %s\n" % (filePath, myMol.GetTitle()))
            if constantMol:
                OEWriteConstMolecule(ofs, myMol)
            else:
                OEWriteMolecule(ofs, myMol)
            return True
        except Exception as e:
            self.__lfh.write("+PdbxBuildChemComp.writeOther FAILING %s\n" % str(e))
            traceback.print_exc(file=self.__lfh)
        return False

    def __getFormalCharge(self, oeMol):
        fCharge = 0
        for atom in oeMol.GetAtoms():
            fCharge += atom.GetFormalCharge()
        return fCharge

    def __makeChemCompAuditRow(self, ccId, action='CREATE', date=None, processingSite='RCSB', annotator='?', details='?'):
        '''
            loop_
            _pdbx_chem_comp_audit.comp_id
            _pdbx_chem_comp_audit.action_type
            _pdbx_chem_comp_audit.date
            _pdbx_chem_comp_audit.processing_site
            _pdbx_chem_comp_audit.annotator
            _pdbx_chem_comp_audit.details
            ARG "Create component"  1999-07-08 RCSB ? ?

        '''
        aRow = {}
        aRow['comp_id'] = ccId
        aRow['action_type'] = action
        if date is None:
            date = time.strftime("%Y-%m-%d", time.localtime())
        aRow['date'] = date
        aRow['processing_site'] = processingSite
        aRow['annotator'] = annotator
        aRow['details'] = details
        #
        return aRow

    def __makeChemCompDescriptorCategory(self, ccId, oeMol):
        '''
            loop_
            _pdbx_chem_comp_descriptor.comp_id
            _pdbx_chem_comp_descriptor.type
            _pdbx_chem_comp_descriptor.program
            _pdbx_chem_comp_descriptor.program_version
            _pdbx_chem_comp_descriptor.descriptor
            ARG SMILES           ACDLabs              10.04 "O=C(O)C(N)CCCN\\C(=[NH2+])N"
            ARG SMILES_CANONICAL CACTVS               3.341 "N[C@@H](CCCNC(N)=[NH2+])C(O)=O"
            ARG SMILES           CACTVS               3.341 "N[CH](CCCNC(N)=[NH2+])C(O)=O"
            ARG SMILES_CANONICAL "OpenEye OEToolkits" 1.5.0 "C(C[C@@H](C(=O)O)N)CNC(=[NH2+])N"
            ARG SMILES           "OpenEye OEToolkits" 1.5.0 "C(CC(C(=O)O)N)CNC(=[NH2+])N"
            ARG InChI            InChI                1.03  "InChI=1S/C6H14N4O2/c7-4(5(11)12)2-1-3-1..... "
            ARG InChIKey         InChI                1.03  ODKSFYDXXFIFQN-BYPYZUCNSA-O
            #
        '''
        rowL = []
        #
        aRow = {}
        aRow['comp_id'] = ccId
        aRow['type'] = 'SMILES_CANONICAL'
        aRow['program'] = "OpenEye OEToolkits"
        aRow['program_version'] = '2016.2'
        aRow['descriptor'] = OECreateIsoSmiString(oeMol)
        rowL.append(aRow)
        #
        aRow = {}
        aRow['comp_id'] = ccId
        aRow['type'] = 'SMILES'
        aRow['program'] = "OpenEye OEToolkits"
        aRow['program_version'] = '2016.2'
        aRow['descriptor'] = OECreateCanSmiString(oeMol)
        rowL.append(aRow)
        #
        aRow = {}
        aRow['comp_id'] = ccId
        aRow['type'] = 'InChI'
        aRow['program'] = "OpenEye OEToolkits"
        aRow['program_version'] = '2016.2'
        aRow['descriptor'] = OECreateInChI(oeMol)
        rowL.append(aRow)
        #
        aRow = {}
        aRow['comp_id'] = ccId
        aRow['type'] = 'InChIKey'
        aRow['program'] = "OpenEye OEToolkits"
        aRow['program_version'] = '2016.2'
        aRow['descriptor'] = OECreateInChIKey(oeMol)
        rowL.append(aRow)
        #
        return rowL

    def __makeChemCompCategory(self, ccId, name=None, formula=None, charge=None, fW=None, site='RCSB', missingModelXyz=False):
        #
        lt = time.strftime("%Y-%m-%d", time.localtime())
        ccRow = {}
        ccRow['id'] = ccId
        if name is not None:
            ccRow['name'] = name
        else:
            ccRow['name'] = '?'
        ccRow['type'] = 'NON-POLYMER'
        ccRow['pdbx_type'] = '?'
        if formula is not None:
            ccRow['formula'] = formula
        else:
            ccRow['formula'] = '?'
        ccRow['mon_nstd_parent_comp_id'] = '?'
        ccRow['pdbx_synonyms'] = '?'
        if charge is not None:
            ccRow['pdbx_formal_charge'] = charge
        else:
            ccRow['pdbx_formal_charge'] = '?'
        ccRow['pdbx_ambiguous_flag'] = 'N'
        ccRow['pdbx_initial_date'] = lt
        ccRow['pdbx_modified_date'] = lt
        ccRow['pdbx_release_status'] = 'HOLD'
        ccRow['pdbx_replaced_by'] = '?'
        ccRow['pdbx_replaces'] = '?'
        if fW is not None:
            ccRow['formula_weight'] = fW
        else:
            ccRow['formula_weight'] = '?'
        ccRow['one_letter_code'] = '?'
        tlc = ccId.split("_")[0]
        ccRow['three_letter_code'] = tlc
        ccRow['pdbx_model_coordinates_details'] = '?'
        ccRow['pdbx_ideal_coordinates_details'] = '?'
        if missingModelXyz:
            ccRow['pdbx_model_coordinates_missing_flag'] = 'Y'
        else:
            ccRow['pdbx_model_coordinates_missing_flag'] = 'N'
        ccRow['pdbx_model_coordinates_db_code'] = '?'
        ccRow['pdbx_processing_site'] = site
        ccRow['pdbx_subcomponent_list'] = '?'
        return ccRow

    def __makeChemCompAtomCategory(self, ccId, oeMol):
        """ Populate elements of chemical component definition for atoms and bonds.
        """
        #
        idCode = ccId
        #
        rowL = []
        for ii, atom in enumerate(oeMol.GetAtoms()):
            #
            atRow = {}
            atRow['comp_id'] = idCode
            (x, y, z) = oeMol.GetCoords(atom)
            # pdbx_model_Cartn_x_ideal
            atRow['pdbx_model_Cartn_x_ideal'] = '%0.3f' % x
            atRow['pdbx_model_Cartn_y_ideal'] = '%0.3f' % y
            atRow['pdbx_model_Cartn_z_ideal'] = '%0.3f' % z
            atRow['atom_id'] = atom.GetName().strip()
            atRow['alt_atom_id'] = atom.GetName().strip()
            atRow['type_symbol'] = OEGetAtomicSymbol(atom.GetAtomicNum())
            atRow['charge'] = atom.GetFormalCharge()
            if atom.GetStringData('pdbx_leaving_atom_flag') in ['Y', 'N']:
                atRow['pdbx_leaving_atom_flag'] = atom.GetStringData('pdbx_leaving_atom_flag')
            else:
                atRow['pdbx_leaving_atom_flag'] = 'N'
            if len(atRow['atom_id']) > 3 or len(atRow['type_symbol']) == 2:
                atRow['pdbx_align'] = 0
            else:
                atRow['pdbx_align'] = 1

            if atom.IsAromatic():
                atRow['pdbx_aromatic_flag'] = 'Y'
            else:
                atRow['pdbx_aromatic_flag'] = 'N'
            # oeSt = OEGetCIPStereo(mol, atom)
            oeSt = None
            # if atom.IsChiral():
            cip = OEPerceiveCIPStereo(oeMol, atom)
            if atom.HasStereoSpecified():
                if cip == OECIPAtomStereo_S:
                    oeSt = "S"
                if cip == OECIPAtomStereo_R:
                    oeSt = "R"
                if cip == OECIPAtomStereo_NotStereo:
                    oeSt = None
                if cip == OECIPAtomStereo_UnspecStereo:
                    oeSt = "U"
                # oeSt = OEGetCIPStereo(oeMol, atom)
                # oeSt = atom.GetCIPStereo()
            if (oeSt is not None and (len(oeSt) > 0) and (oeSt == "R" or oeSt == "S")):
                atRow['pdbx_stereo_config'] = oeSt
            else:
                atRow['pdbx_stereo_config'] = 'N'
            atRow['pdbx_component_atom_id'] = atom.GetName().strip()
            atRow['pdbx_component_comp_id'] = idCode
            atRow['pdbx_ordinal'] = str(ii + 1)
            rowL.append(atRow)
        #
        return rowL

    def __makeChemCompBondCategory(self, ccId, oeMol):
        idCode = ccId
        rowL = []
        iBondD = {1: "SING", 2: "DOUB", 3: "TRIP", 4: "QUAD", 5: "AROM", 0: "DELO"}

        for ii, bond in enumerate(oeMol.GetBonds()):
            oeOrder = str(bond.GetOrder()).upper()
            b1 = bond.GetBgn().GetName().strip()
            b2 = bond.GetEnd().GetName().strip()
            if ((b1 is None) or (len(b1) < 1)):
                continue
            if ((b2 is None) or (len(b1) < 1)):
                continue
            if (int(oeOrder) in iBondD):
                bndRow = {}
                bndRow['comp_id'] = idCode
                bndRow['atom_id_1'] = b1
                bndRow['atom_id_2'] = b2
                bndRow['value_order'] = iBondD[int(oeOrder)]

                if bond.IsAromatic():
                    bndRow['pdbx_aromatic_flag'] = 'Y'
                else:
                    bndRow['pdbx_aromatic_flag'] = 'N'

                oeSt = None
                if (bond.GetOrder() == 2):
                    cip = OEPerceiveCIPStereo(oeMol, bond)
                    if bond.HasStereoSpecified():
                        if cip == OECIPBondStereo_E:
                            oeSt = "E"
                        if cip == OECIPBondStereo_Z:
                            oeSt = "Z"
                        if cip == OECIPBondStereo_NotStereo:
                            oeSt = None
                        if cip == OECIPBondStereo_UnspecStereo:
                            oeSt = None
                if (oeSt is not None and (len(oeSt) > 0) and (oeSt == "E" or oeSt == "Z")):
                    bndRow['pdbx_stereo_config'] = oeSt
                else:
                    bndRow['pdbx_stereo_config'] = 'N'

                bndRow['pdbx_ordinal'] = str(ii + 1)
                rowL.append(bndRow)
        return rowL
