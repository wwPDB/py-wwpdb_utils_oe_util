##
# File:  OeMCSS.py
# Date:  29-Nov-2014  J. Westbrook
#
# Updates:
#    6-Jun-2016 jdw general cleanup
#
##
"""
Classes implementing MCSS comparisons --

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.01"


import os.path
import sys
import traceback

from openeye.oechem import (
    OEAddExplicitHydrogens,
    OEExprOpts_AtomicNumber,
    OEExprOpts_DefaultAtoms,
    OEExprOpts_DefaultBonds,
    OEExprOpts_ExactAtoms,
    OEExprOpts_ExactBonds,
    OEFloatArray,
    OEMCSMaxAtoms,
    OEMCSSearch,
    OEMCSType_Approximate,
)
from wwpdb.utils.cc_dict_util.timeout.TimeoutMultiProc import timeout
from wwpdb.utils.oe_util.build.OeBuildMol import OeBuildMol


class OeMCSS(object):
    """Perform MCSS comparisons.  Targets can be chemical component identifiers
    or paths to chemical component definition files.  Inputs can be in the the form of pairs,
    lists, and pair lists of chemical component definitions.

    """

    def __init__(self, verbose=True, log=sys.stderr):
        #
        self.__verbose = verbose
        # self.__debug = False
        self.__lfh = log
        #
        self.__refId = None
        self.__refmol = None
        # self.__refTitle = None
        #
        self.__fitId = None
        self.__fitmol = None
        # self.__fitTitle = None
        #
        self.__pairTupleList = []
        #
        # self.__minAtomMatchFraction = 0.50
        #
        self.__searchType = "default"
        #
        self.__refFD = {}
        self.__fitFD = {}
        #
        self.__refPath = None
        self.__fitPath = None
        self.__mcss = None

    def setSearchType(self, sType="default"):
        self.__searchType = sType
        return self.__searchType

    def setRefId(self, ccId, title=None, suppressHydrogens=False, cachePath="/data/components/ligand-dict-v3"):
        """Set the query molecule for MCSS comparison using the input chemical component ID.
        It is assumed that the definition for this ID can be obtained from the chemical component
        repository.

        Once the reference molecule is built, the MCSS calculation is initialized.

        A title is optionally provided otherwise the component Id will be used.

        The hydrogen flag can be used to perform the MCSS using only heavy atoms.
        """
        self.__refId = ccId
        ccIdU = ccId.upper()
        self.__refPath = os.path.join(cachePath, ccIdU[0], ccIdU, ccIdU + ".cif")

        cid, self.__refmol, self.__refFD = self.__getCCDefFile(self.__refPath, suppressHydrogens=suppressHydrogens)  # pylint: disable=unused-variable
        #
        # Insert title here -
        if title is not None:
            self.__refmol.SetTitle(title)
            # self.__refTitle = title
        else:
            self.__refmol.SetTitle(self.__refId)
            # self.__refTitle = None
        #
        #
        self.__setupMCSS(self.__refmol)

    def setRefPath(self, ccPath, title=None, suppressHydrogens=False, type="CC", importType="2D"):  # pylint: disable=redefined-builtin
        """Set the query molecule for MCSS comparison using the input file path.

        The file type is either ['CC'] for a chemical component definition or another file type
        supported by OE toolkit assumed to have a conventional file extension for this type.

        Once the reference molecule is built, the MCSS calculation is initialized.

        A title is optionally provided otherwise the component Id will be used.

        The hydrogen flag can be used to perform the MCSS using only heavy atoms.
        """
        self.__refPath = ccPath
        if type in ["CC"]:
            (self.__refId, self.__refmol, self.__refFD) = self.__getCCDefFile(ccPath, suppressHydrogens=suppressHydrogens)
        else:
            (self.__refId, self.__refmol, self.__refFD) = self.__getMiscFile(ccPath, suppressHydrogens=suppressHydrogens, importType=importType)

        if self.__verbose:
            self.__lfh.write("Derived ref ID     = %s\n" % self.__refId)
            self.__lfh.write("SMILES (stereo)  = %s\n" % self.__refFD["SMILES_STEREO"])
        #
        # Insert title here -
        if title is not None:
            self.__refmol.SetTitle(title)
            # self.__refTitle = title
        else:
            self.__refmol.SetTitle(self.__refId)
            # self.__refTitle = None

        self.__setupMCSS(self.__refmol)

    def setFitId(self, ccId, title=None, suppressHydrogens=False, cachePath="/data/components/ligand-dict-v3"):
        """Set the ID of the target/library molecule for MCSS comparison."""
        self.__fitId = ccId
        ccIdU = ccId.upper()
        self.__fitPath = os.path.join(cachePath, ccIdU[0], ccIdU, ccIdU + ".cif")
        self.__fitId, self.__fitmol, self.__fitFD = self.__getCCDefFile(self.__fitPath, suppressHydrogens=suppressHydrogens)
        if self.__verbose:
            self.__lfh.write("Fit ID             = %s\n" % self.__fitId)
            self.__lfh.write("SMILES (isomeric)  = %s\n" % self.__fitFD["SMILES_STEREO"])
        if title is not None:
            self.__fitmol.SetTitle(title)
            # self.__fitTitle = title
        else:
            self.__fitmol.SetTitle(self.__fitId)
            # self.__fitTitle = None

    def setFitPath(self, ccPath, title=None, suppressHydrogens=False):
        """Set the path of the target/library molecule for MCSS comparison."""
        (self.__fitId, self.__fitmol, self.__fitFD) = self.__getCCDefFile(ccPath, suppressHydrogens=suppressHydrogens)
        if title is not None:
            self.__fitmol.SetTitle(title)
            # self.__fitTitle = title
        else:
            self.__fitmol.SetTitle(self.__fitId)
            # self.__fitTitle = None

    def setFitIdList(self, ccIdList, cachePath="/data/components/ligand-dict-v3", suppressHydrogens=False):  # pylint: disable=unused-argument
        """Set the list of IDs to be compared with reference molecule by MCSS.

        From the input ID list build the internal pair list of
        tuples  [(refId,refPath,refTitle,fitId,fitPath,fitTitle),(),...]
        """
        self.__pairTupleList = []
        for ccId in ccIdList:
            refId = self.__refId
            refPath = self.__refPath
            #
            fitId = ccId
            ccIdU = ccId.upper()
            fitPath = os.path.join(cachePath, ccIdU[0], ccIdU, ccIdU + ".cif")
            #
            refTitle = refId + "/" + fitId
            fitTitle = fitId + "/" + refId
            self.__pairTupleList.append((refId, refPath, refTitle, fitId, fitPath, fitTitle))

    def setPairIdList(self, pairIdList, cachePath="/data/components/ligand-dict-v3", suppressHydrogens=False):  # pylint: disable=unused-argument
        """Set the list of ID pais to be aligned by MCSS.

        From the input ID list build the internal pair list of
        tuples  [(refId,refPath,refTitle,fitId,fitPath,fitTitle),(),...]
        """

        self.__pairTupleList = []
        for refId, fitId in pairIdList:
            ccIdU = refId.upper()
            refPath = os.path.join(cachePath, ccIdU[0], ccIdU, ccIdU + ".cif")
            #
            ccIdU = fitId.upper()
            fitPath = os.path.join(cachePath, ccIdU[0], ccIdU, ccIdU + ".cif")
            #
            refTitle = refId + "/" + fitId
            fitTitle = fitId + "/" + refId
            self.__pairTupleList.append((refId, refPath, refTitle, fitId, fitPath, fitTitle))

    def __getCCDefFile(self, ccPath, suppressHydrogens=False):
        """Fetch the molecule definition (ccPath) and build OE molecules
        for comparison and depiction.

        """
        #
        oem = OeBuildMol(verbose=self.__verbose, log=self.__lfh)
        ccId = oem.setChemCompPath(ccPath)
        #
        oem.build2D()

        if self.__verbose:
            self.__lfh.write("+OeMCSS.__getCCDefFile() for %s\n" % ccId)
            self.__lfh.write("  Title              = %s\n" % oem.getTitle())
            self.__lfh.write("  SMILES             = %s\n" % oem.getCanSMILES())
            self.__lfh.write("  SMILES (stereo)    = %s\n" % oem.getIsoSMILES())
            self.__lfh.write("  Formula (Hill)     = %s\n" % oem.getFormula())
            self.__lfh.write("  InChI key          = %s\n" % oem.getInChIKey())
            self.__lfh.write("  InChI              = %s\n" % oem.getInChI())

        fD = {}
        fD = {"Formula": oem.getFormula(), "SMILES": oem.getCanSMILES(), "SMILES_STEREO": oem.getIsoSMILES(), "InChI": oem.getInChI(), "InChIKey": oem.getInChIKey()}

        if suppressHydrogens:
            tMol = oem.getGraphMolSuppressH()
        else:
            tMol = oem.getMol()

        fD["OEMOL"] = tMol

        return (ccId, tMol, fD)

    def __getMiscFile(self, ccPath, suppressHydrogens=False, importType="2D"):
        """Fetch a miscellaneous chemical file (ccPath) and build OE molecules
        for comparison and depiction.

        """
        try:
            oem = OeBuildMol(verbose=self.__verbose, log=self.__lfh)
            if oem.importFile(ccPath, type=importType):
                if self.__verbose:
                    self.__lfh.write("+OeMCSS.__getMiscFile()\n")
                    self.__lfh.write("  Title              = %s\n" % oem.getTitle())
                    self.__lfh.write("  SMILES             = %s\n" % oem.getCanSMILES())
                    self.__lfh.write("  SMILES (stereo)    = %s\n" % oem.getIsoSMILES())
                    self.__lfh.write("  Formula (Hill)     = %s\n" % oem.getFormula())
                    self.__lfh.write("  InChI key          = %s\n" % oem.getInChIKey())
                    self.__lfh.write("  InChI              = %s\n" % oem.getInChI())
            else:
                self.__lfh.write("+OeMCSS.__getMiscFile() Read failed for %s\n" % ccPath)
                return None, None, None
            #
            # oem.build2D()
            ccId = oem.getTitle()
            if suppressHydrogens:
                tMol = oem.getGraphMolSuppressH()
            else:
                tMol = oem.getMol()

            molXyzL = []
            if importType == "3D":
                for ii, atm in enumerate(tMol.GetAtoms()):
                    xyzL = OEFloatArray(3)
                    tMol.GetCoords(atm, xyzL)
                    molXyzL.append((ii, atm.GetIdx(), atm.GetAtomicNum(), atm.GetName(), atm.GetType(), "%.3f" % xyzL[0], "%.3f" % xyzL[1], "%.3f" % xyzL[2]))
            fD = {}
            fD = {
                "Formula": oem.getFormula(),
                "SMILES": oem.getCanSMILES(),
                "SMILES_STEREO": oem.getIsoSMILES(),
                "InChI": oem.getInChI(),
                "InChIKey": oem.getInChIKey(),
                "xyz": molXyzL,
            }

            for ii, atm in enumerate(tMol.GetAtoms()):
                xyzL = OEFloatArray(3)
                tMol.GetCoords(atm, xyzL)
                if self.__verbose:
                    self.__lfh.write("OeMCSS.__getMiscFile - atom  %d %s %s %s %s %r\n" % (ii, atm.GetIdx(), atm.GetAtomicNum(), atm.GetName(), atm.GetType(), xyzL))

            fD["OEMOL"] = tMol
            return (ccId, tMol, fD)
        except Exception as e:  # noqa: F841 pylint: disable=unused-variable
            traceback.print_exc(file=self.__lfh)
            # self.fail()

        return None, None, None

    def __setupMCSS(self, refmol):
        """Internal initialization for the MCSS comparison."""
        #
        if self.__searchType == "default":
            self.__mcss = OEMCSSearch(OEMCSType_Approximate)
            atomexpr = OEExprOpts_DefaultAtoms
            bondexpr = OEExprOpts_DefaultBonds
        elif self.__searchType == "relaxed":
            self.__mcss = OEMCSSearch(OEMCSType_Approximate)
            atomexpr = OEExprOpts_AtomicNumber
            bondexpr = 0
        elif self.__searchType == "exact":
            self.__mcss = OEMCSSearch(OEMCSType_Approximate)
            # self.__mcss = OEMCSSearch(OEMCSType_Exhaustive)
            atomexpr = OEExprOpts_ExactAtoms
            bondexpr = OEExprOpts_ExactBonds
            OEAddExplicitHydrogens(refmol)
        else:
            self.__mcss = OEMCSSearch(OEMCSType_Approximate)
            atomexpr = OEExprOpts_DefaultAtoms
            bondexpr = OEExprOpts_DefaultBonds
        #
        # atomexpr = OEExprOpts_AtomicNumber|OEExprOpts_EqAromatic
        # bondexpr = 0
        #
        # atomexpr = OEExprOpts_AtomicNumber|OEExprOpts_Aromaticity
        # bondexpr = OEExprOpts_BondOrder|OEExprOpts_EqNotAromatic
        #

        self.__mcss.Init(refmol, atomexpr, bondexpr)

        #
        # self.__mcss.SetMCSFunc(OEMCSMaxBondsCompleteCycles())
        # self.__mcss.SetMCSFunc(OEMCSMaxAtoms())
        #
        # Half of the reference molecule --
        #
        # nAtomsRef=refmol.NumAtoms()
        # self.__mcss.SetMinAtoms(nAtomsRef/2)

    @timeout(30)
    def doAlign(self, suppressHydrogens=False, unique=True, minFrac=1.0):  # pylint: disable=unused-argument
        """Test the MCSS comparison between current reference and fit molecules -
        Return list of corresponding atoms on success or an empty list otherwise.
        """
        atomMap = []
        #
        nAtomsRef = self.__refmol.NumAtoms()
        nAtomsFit = self.__fitmol.NumAtoms()
        minAtoms = int(min(nAtomsRef, nAtomsFit) * minFrac)
        # self.__mcss.SetMinAtoms( int(minAtoms*self.__minAtomMatchFraction) )
        #
        # -------
        self.__mcss.SetMCSFunc(OEMCSMaxAtoms())
        self.__mcss.SetMinAtoms(minAtoms)
        OEAddExplicitHydrogens(self.__refmol)
        OEAddExplicitHydrogens(self.__fitmol)
        #
        # --------
        miter = self.__mcss.Match(self.__fitmol, unique)
        if miter.IsValid():
            match = miter.Target()
            for mAt in match.GetAtoms():
                atomMap.append(
                    (
                        self.__refId,
                        mAt.pattern.GetIdx(),
                        mAt.pattern.GetAtomicNum(),
                        mAt.pattern.GetName(),
                        self.__fitId,
                        mAt.target.GetIdx(),
                        mAt.target.GetAtomicNum(),
                        mAt.target.GetName(),
                    )
                )
        return (nAtomsRef, self.__refFD, nAtomsFit, self.__fitFD, atomMap)
