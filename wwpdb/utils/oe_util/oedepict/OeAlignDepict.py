##
# File:  OeAlignDepict.py
# Date:  2-Oct-2011  J. Westbrook
#
# Updates:
#  5-Jan-2012 jdw Revise for new OE depiction API
#  5-Mar-2012 jdw refactor
#  7-Nar-2012 jdw adapt to support latest OeBuildMol
#  9-Mar-2012 jdw simplify calling interface
# 10-Jun-2016 jdw add method setFitPathList
#
##
"""
Classes to depict MCS alignments.

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

from openeye.oechem import (
    OEAddExplicitHydrogens,
    OEBlueTint,
    OEExprOpts_AtomicNumber,
    OEExprOpts_DefaultAtoms,
    OEExprOpts_DefaultBonds,
    OEExprOpts_ExactAtoms,
    OEExprOpts_ExactBonds,
    OEFloatArray,
    OEGreenTint,
    OEIsAtomMember,
    OEIsBondMember,
    OEMCSMaxAtoms,
    OEMCSSearch,
    OEMCSType_Approximate,
    OENotAtom,
    OENotBond,
    OEPinkTint,
)
from openeye.oedepict import (
    OE2DMolDisplay,
    OE2DMolDisplayOptions,
    OEAddHighlighting,
    OEAtomStereoStyle_Display_All,
    OEGetMoleculeScale,
    OEHighlightStyle_BallAndStick,
    OEHighlightStyle_Stick,
    OEImage,
    OEImageGrid,
    OEMultiPageImageFile,
    OEPageOrientation_Landscape,
    OEPageOrientation_Portrait,
    OEPageSize_US_Letter,
    OEPrepareAlignedDepiction,
    OEPrepareDepiction,
    OERenderMolecule,
    OEScale_AutoScale,
    OEWriteImage,
    OEWriteMultiPageImage,
)
from wwpdb.utils.cc_dict_util.timeout.TimeoutMultiProc import timeout
from wwpdb.utils.oe_util.build.OeBuildMol import OeBuildMol


class OeDepictMCSAlign(object):

    """Create 2D depictions of MCSS alignments.  Targets can be chemical component identifiers
    or paths to chemical component definition files.  Inputs can be in the the form of pairs,
    lists, and pair lists of chemical component definitions.

    """

    def __init__(self, verbose=True, log=sys.stderr):
        #
        self.__verbose = verbose
        self.__debug = False
        self.__lfh = log
        #
        self.__refId = None
        self.__refmol = None
        self.__refTitle = None
        #
        self.__fitId = None
        self.__fitmol = None
        self.__fitTitle = None
        #
        self.__pairTupleList = []
        #
        self.__minAtomMatchFraction = 0.50
        self.__pageOrientation = "portrait"
        #
        self.__searchType = "default"
        #
        self.__refFD = {}
        self.__fitFD = {}
        #
        self.__refPath = None
        self.__fitPath = None
        self.__mcss = None
        #
        self.__image = None
        self.__grid = None
        self.__opts = None
        self.__gridRows = None
        self.__gridCols = None
        self.__multi = None

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

        _id, self.__refmol, self.__refFD = self.__getCCDefFile(self.__refPath, suppressHydrogens=suppressHydrogens)
        #
        # Insert title here -
        if title is not None:
            self.__refmol.SetTitle(title)
            self.__refTitle = title
        else:
            self.__refmol.SetTitle(self.__refId)
            self.__refTitle = None
        #
        #
        OEPrepareDepiction(self.__refmol)
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
            self.__refTitle = title
        else:
            self.__refmol.SetTitle(self.__refId)
            self.__refTitle = None

        OEPrepareDepiction(self.__refmol)
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
            self.__fitTitle = title
        else:
            self.__fitmol.SetTitle(self.__fitId)
            self.__fitTitle = None

    def setFitPath(self, ccPath, title=None, suppressHydrogens=False):
        """Set the path of the target/library molecule for MCSS comparison."""
        (self.__fitId, self.__fitmol, self.__fitFD) = self.__getCCDefFile(ccPath, suppressHydrogens=suppressHydrogens)
        if title is not None:
            self.__fitmol.SetTitle(title)
            self.__fitTitle = title
        else:
            self.__fitmol.SetTitle(self.__fitId)
            self.__fitTitle = None

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

    def setFitPathList(self, fitPathTupList, suppressHydrogens=False):  # pylint: disable=unused-argument
        """Set the list of paths for the target/library molecules to be compared with reference molecule by MCSS.

        From the input path tuple list build the internal pair list of
        tuples  [(fitId,fitPath,fitTitle),(),...]
        """
        self.__pairTupleList = []
        for fitId, fitPath, fitTitle in fitPathTupList:
            refId = self.__refId
            refPath = self.__refPath
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
        oem.build2D()

        if self.__verbose:
            self.__lfh.write("+OEAlignDepilsct.__getCCDefFile() for %s\n" % ccId)
            self.__lfh.write("  Title              = %s\n" % oem.getTitle())
            self.__lfh.write("  SMILES             = %s\n" % oem.getCanSMILES())
            self.__lfh.write("  SMILES (stereo)    = %s\n" % oem.getIsoSMILES())
            self.__lfh.write("  Formula (Hill)     = %s\n" % oem.getFormula())
            self.__lfh.write("  InChI key          = %s\n" % oem.getInChIKey())
            self.__lfh.write("  InChI              = %s\n" % oem.getInChI())

        fD = {}
        fD = {"Formula": oem.getFormula(), "SMILES": oem.getCanSMILES(), "SMILES_STEREO": oem.getIsoSMILES(), "InChI": oem.getInChI(), "InChIKey": oem.getInChIKey()}

        if suppressHydrogens:
            return (ccId, oem.getGraphMolSuppressH(), fD)
        else:
            return (ccId, oem.getMol(), fD)

    def __getMiscFile(self, ccPath, suppressHydrogens=False, importType="2D"):
        """Fetch a miscellaneous chemical file (ccPath) and build OE molecules
        for comparison and depiction.

        """
        try:
            oem = OeBuildMol(verbose=self.__verbose, log=self.__lfh)
            if oem.importFile(ccPath, type=importType):
                if self.__verbose:
                    self.__lfh.write("+OEAlignDepilsct.__getMiscFile()\n")
                    self.__lfh.write("  Title              = %s\n" % oem.getTitle())
                    self.__lfh.write("  SMILES             = %s\n" % oem.getCanSMILES())
                    self.__lfh.write("  SMILES (stereo)    = %s\n" % oem.getIsoSMILES())
                    self.__lfh.write("  Formula (Hill)     = %s\n" % oem.getFormula())
                    self.__lfh.write("  InChI key          = %s\n" % oem.getInChIKey())
                    self.__lfh.write("  InChI              = %s\n" % oem.getInChI())
            else:
                self.__lfh.write("+OEAlignDepict.__getMiscFile() Read failed for %s\n" % ccPath)
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
                    self.__lfh.write("OeAlignDepict.__getMiscFile - atom  %d %s %s %s %s %r\n" % (ii, atm.GetIdx(), atm.GetAtomicNum(), atm.GetName(), atm.GetType(), xyzL))

            return (ccId, tMol, fD)
        except:  # noqa: E722 pylint: disable=bare-except
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

    @timeout(15)
    def testAlign(self, suppressHydrogens=False, unique=True, minFrac=1.0):  # pylint: disable=unused-argument
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
        # --------
        miter = self.__mcss.Match(self.__fitmol, unique)
        if miter.IsValid():
            match = miter.Target()
            for mAt in match.GetAtoms():
                # atomMap.append( (self.__refId,mAt.pattern.GetName(),self.__fitId,mAt.target.GetName() ))
                atomMap.append((self.__refId, mAt.pattern.GetIdx(), mAt.pattern.GetType(), mAt.pattern.GetName(), self.__fitId, mAt.target.GetIdx(), mAt.target.GetName()))

        return atomMap

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
        # --------
        miter = self.__mcss.Match(self.__fitmol, unique)
        if miter.IsValid():
            match = miter.Target()
            for mAt in match.GetAtoms():
                # atomMap.append( (self.__refId,mAt.pattern.GetName(),self.__fitId,mAt.target.GetName() ))
                atomMap.append((self.__refId, mAt.pattern.GetIdx(), mAt.pattern.GetType(), mAt.pattern.GetName(), self.__fitId, mAt.target.GetIdx(), mAt.target.GetName()))

        return (nAtomsRef, self.__refFD["SMILES_STEREO"], nAtomsFit, self.__fitFD["SMILES_STEREO"], atomMap)

    @timeout(30)
    def doAlignAlt(self, suppressHydrogens=False, unique=True, minFrac=1.0):  # pylint: disable=unused-argument
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

    def __setupImage(self, imageX=400, imageY=400):
        """Internal method to configure a single pair alignment image."""
        #
        self.__image = OEImage(imageX, imageY)
        rows = 1
        cols = 2
        self.__grid = OEImageGrid(self.__image, rows, cols)
        if self.__debug:
            self.__lfh.write("Num columns %d\n" % self.__grid.NumCols())
            self.__lfh.write("Num rows    %d\n" % self.__grid.NumRows())
        self.__opts = OE2DMolDisplayOptions(self.__grid.GetCellWidth(), self.__grid.GetCellHeight(), OEScale_AutoScale)

    def __setupImageMulti(self, gridRows=2, gridCols=2):
        """Internal method to configure a multipage image."""
        #
        self.__gridRows = gridRows
        self.__gridCols = gridCols
        #
        if self.__pageOrientation == "landscape":
            self.__multi = OEMultiPageImageFile(OEPageOrientation_Landscape, OEPageSize_US_Letter)
        else:
            self.__multi = OEMultiPageImageFile(OEPageOrientation_Portrait, OEPageSize_US_Letter)

        # self.__image = self.__multi.NewPage()
        # self.__opts = OE2DMolDisplayOptions()
        self.__newPage()

    def __newPage(self):
        """Internal method to advance to a new page in a multipage configuration."""
        rows = self.__gridRows
        cols = self.__gridCols
        self.__image = self.__multi.NewPage()
        self.__grid = OEImageGrid(self.__image, rows, cols)
        self.__grid.SetCellGap(20)
        self.__grid.SetMargins(20)
        if self.__debug:
            self.__lfh.write("Num columns %d\n" % self.__grid.NumCols())
            self.__lfh.write("Num rows    %d\n" % self.__grid.NumRows())
        self.__opts = OE2DMolDisplayOptions(self.__grid.GetCellWidth(), self.__grid.GetCellHeight(), OEScale_AutoScale)

    def alignPair(self, imagePath="mcs-match.svg", imageX=400, imageY=400):
        """Compare current reference and fit molecules by MCSS.

        Map of corresponding atoms is returned.

        imagePath is the path of the output image.
        """
        atomMap = []
        self.__setupImage(imageX=imageX, imageY=imageY)

        if self.__refTitle is None:
            self.__refTitle = self.__refId + "/" + self.__fitId
        self.__refmol.SetTitle(self.__refTitle)
        OEPrepareDepiction(self.__refmol)

        #
        if self.__fitTitle is None:
            self.__fitTitle = self.__fitId + "/" + self.__refId
        self.__fitmol.SetTitle(self.__fitTitle)
        OEPrepareDepiction(self.__fitmol)

        #
        # opts = OE2DMolDisplayOptions(self.__grid.GetCellWidth(), self.__grid.GetCellHeight(), OEScale_AutoScale)
        #
        refcell = self.__grid.GetCell(1, 1)
        fitcell = self.__grid.GetCell(1, 2)

        self.__opts.SetAtomStereoStyle(OEAtomStereoStyle_Display_All)
        # opts.SetTitleLocation(OETitleLocation_Hidden)

        refscale = OEGetMoleculeScale(self.__refmol, self.__opts)
        fitscale = OEGetMoleculeScale(self.__refmol, self.__opts)
        self.__opts.SetScale(min(refscale, fitscale))

        hstyle = OEHighlightStyle_BallAndStick

        unique = True
        miter = self.__mcss.Match(self.__fitmol, unique)
        if miter.IsValid():
            match = miter.Target()
            OEPrepareAlignedDepiction(self.__fitmol, self.__mcss.GetPattern(), match)

            # Depict reference molecule with MCS highlighting

            refdisp = OE2DMolDisplay(self.__mcss.GetPattern(), self.__opts)
            patoms = OEIsAtomMember(match.GetPatternAtoms())
            OEAddHighlighting(refdisp, OEBlueTint, hstyle, patoms)
            pbonds = OEIsBondMember(match.GetPatternBonds())
            OEAddHighlighting(refdisp, OEBlueTint, hstyle, pbonds)
            OERenderMolecule(refcell, refdisp)

            # Depict fit molecule with MCS highlighting

            fitdisp = OE2DMolDisplay(self.__fitmol, self.__opts)
            tatoms = OEIsAtomMember(match.GetTargetAtoms())
            OEAddHighlighting(fitdisp, OEBlueTint, hstyle, tatoms)
            tbonds = OEIsBondMember(match.GetTargetBonds())
            OEAddHighlighting(fitdisp, OEBlueTint, hstyle, tbonds)
            OERenderMolecule(fitcell, fitdisp)

            for mAt in match.GetAtoms():
                atomMap.append((self.__refId, mAt.pattern.GetName(), self.__fitId, mAt.target.GetName()))

        OEWriteImage(imagePath, self.__image)
        return atomMap

    def alignPairList(
        self,
        imagePath="multi.pdf",
        suppressHydrogens=False,
        gridRows=2,
        gridCols=2,
        optInverseFit=True,
        highLightStyle="stick",
        highLightMatchColor="green",
        highLightNotMatchColor="pink",
    ):
        """Compare molecule pairs in the current __pairTupleList.

        pairTupleList = (refId,refPath,refTitle,fitId,fitPath,fitTitle)

        Map of corresponding atoms is returned.

        Image Output is in multipage layout.
        """
        #
        self.__setupImageMulti(gridRows=gridRows, gridCols=gridCols)
        #
        atomMap = []
        #
        numRows = self.__grid.NumRows()
        rowIdx = 0
        for (refId, refPath, refTitle, fitId, fitPath, fitTitle) in self.__pairTupleList:
            self.setRefPath(refPath, title=refTitle, suppressHydrogens=suppressHydrogens)
            OEPrepareDepiction(self.__refmol)

            #
            _tId, fitmol, _fitFD = self.__getCCDefFile(fitPath, suppressHydrogens=suppressHydrogens)
            fitmol.SetTitle(fitTitle)
            #
            # self.__minAtomMatchFraction %  of the smaller molecule --
            #
            nAtomsRef = self.__refmol.NumAtoms()
            nAtomsFit = fitmol.NumAtoms()
            minAtoms = min(nAtomsRef, nAtomsFit)
            mcssMinAtoms = int(minAtoms * self.__minAtomMatchFraction)
            self.__mcss.SetMinAtoms(mcssMinAtoms)

            OEPrepareDepiction(fitmol)
            rowIdx += 1
            if rowIdx > numRows:
                self.__newPage()
                rowIdx = 1
            refcell = self.__grid.GetCell(rowIdx, 1)
            fitcell = self.__grid.GetCell(rowIdx, 2)

            self.__opts.SetAtomStereoStyle(OEAtomStereoStyle_Display_All)
            # opts.SetTitleLocation(OETitleLocation_Hidden)

            refscale = OEGetMoleculeScale(self.__refmol, self.__opts)
            fitscale = OEGetMoleculeScale(self.__refmol, self.__opts)
            self.__opts.SetScale(min(refscale, fitscale))

            if highLightStyle == "ballAndStick":
                hstyle = OEHighlightStyle_BallAndStick
            elif highLightStyle == "stick":
                hstyle = OEHighlightStyle_Stick
            else:
                hstyle = OEHighlightStyle_BallAndStick

            if highLightMatchColor == "blue":
                myHighLightMatchColor = OEBlueTint
            elif highLightMatchColor == "green":
                myHighLightMatchColor = OEGreenTint
            elif highLightMatchColor == "pink":
                myHighLightMatchColor = OEPinkTint
            else:
                myHighLightMatchColor = OEBlueTint

            if highLightNotMatchColor == "blue":
                myHighLightNotMatchColor = OEBlueTint
            elif highLightNotMatchColor == "green":
                myHighLightNotMatchColor = OEGreenTint
            elif highLightNotMatchColor == "pink":
                myHighLightNotMatchColor = OEPinkTint
            else:
                myHighLightNotMatchColor = OEBlueTint

            optInverseFit = True

            self.__lfh.write("+OeAlignDepict.alignPairList refId %s fitId %s nAtomsRef %d nAtomsFit %s mcssMinAtoms %d\n" % (refId, fitId, nAtomsRef, nAtomsFit, mcssMinAtoms))
            unique = True

            miter = self.__mcss.Match(fitmol, unique)
            if miter.IsValid():
                match = miter.Target()
                OEPrepareAlignedDepiction(fitmol, self.__mcss.GetPattern(), match)

                # Depict reference molecule with MCS highlighting

                refdisp = OE2DMolDisplay(self.__mcss.GetPattern(), self.__opts)
                patoms = OEIsAtomMember(match.GetPatternAtoms())
                OEAddHighlighting(refdisp, myHighLightMatchColor, hstyle, patoms)
                pbonds = OEIsBondMember(match.GetPatternBonds())
                OEAddHighlighting(refdisp, myHighLightMatchColor, hstyle, pbonds)
                OERenderMolecule(refcell, refdisp)

                # Depict fit molecule with MCS highlighting
                if optInverseFit:
                    fitdisp = OE2DMolDisplay(fitmol, self.__opts)
                    tatoms = OEIsAtomMember(match.GetTargetAtoms())
                    OEAddHighlighting(fitdisp, myHighLightNotMatchColor, hstyle, OENotAtom(tatoms))
                    tbonds = OEIsBondMember(match.GetTargetBonds())
                    OEAddHighlighting(fitdisp, myHighLightNotMatchColor, hstyle, OENotBond(tbonds))
                    OERenderMolecule(fitcell, fitdisp)
                else:
                    fitdisp = OE2DMolDisplay(fitmol, self.__opts)
                    tatoms = OEIsAtomMember(match.GetTargetAtoms())
                    OEAddHighlighting(fitdisp, myHighLightNotMatchColor, hstyle, tatoms)
                    tbonds = OEIsBondMember(match.GetTargetBonds())
                    OEAddHighlighting(fitdisp, myHighLightNotMatchColor, hstyle, tbonds)
                    OERenderMolecule(fitcell, fitdisp)

                for mAt in match.GetAtoms():
                    atomMap.append((self.__refId, mAt.pattern.GetName(), fitId, mAt.target.GetName()))

        OEWriteMultiPageImage(imagePath, self.__multi)
        return atomMap
