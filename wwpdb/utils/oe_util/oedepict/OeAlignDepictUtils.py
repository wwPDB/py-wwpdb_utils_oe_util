##
# File:  OeAlignDepictUtils.py
# Date:  29-Jun-2013  J. Westbrook
#
# Updates:
#         1-July-2013 jdw  Add OeDepictMCSAlignSingle(OeDepictAlignBase,OeDepictBase)
#         2-July-2013 jdw  Revise api for OeDepictAlignBase() and add method to extend
#                          comparison list incrementally.
#        12-June-2016 jdw  add  method setFitPathList()  extend the _params[] to support more
#                          display styles in  _setHighlightStyleRef _setHighlightStyleFit
##
"""
Classes to depict aligned 2D chemical diagrams.

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
    OEBlueTint,
    OEExprOpts_AtomicNumber,
    OEExprOpts_DefaultAtoms,
    OEExprOpts_DefaultBonds,
    OEExprOpts_ExactAtoms,
    OEExprOpts_ExactBonds,
    OEGetMCSExhaustiveSearchTruncationLimit,
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
    OEBlackPen,
    OEDrawBorder,
    OEGetMoleculeScale,
    OEHighlightStyle_BallAndStick,
    OEHighlightStyle_Stick,
    OEImage,
    OEImageGrid,
    OEMultiPageImageFile,
    OEPageOrientation_Landscape,
    OEPageOrientation_Portrait,
    OEPageSize_US_Letter,
    OEPen,
    OEPrepareAlignedDepiction,
    OEPrepareDepiction,
    OERenderMolecule,
    OEScale_AutoScale,
    OEWriteImage,
    OEWriteMultiPageImage,
)
from wwpdb.utils.cc_dict_util.timeout.TimeoutMultiProc import TimeoutException, timeout
from wwpdb.utils.oe_util.build.OeChemCompIoUtils import OeChemCompIoUtils
from wwpdb.utils.oe_util.oedepict.OeDepict import OeDepictBase


# Note for pylint testing.  This base class references self._opt and self._param, but it is defined in a child class of a parallel class. We do not init it
class OeDepictAlignBase(object):
    """"""

    # pylint: disable=no-member

    def __init__(self, verbose=True, log=sys.stderr):
        super(OeDepictAlignBase, self).__init__()
        #
        self.__verbose = verbose
        # self.__debug = False
        self.__lfh = log
        #
        self._refId = None
        self._refPath = None
        self._refMol = None
        self._refTitle = None
        self._refImagePath = None
        #
        self._fitId = None
        self._fitPath = None
        self._fitMol = None
        self._fitTitle = None
        self._fitImagePath = None
        #
        self._pairMolList = []
        #
        self._searchType = "default"
        self._minAtomMatchFraction = 0.50
        self._mcss = None
        self._miter = None
        #
        # We cannot define these here - class heirarchy might be init at higher level.
        # self._opts = None
        # self._params = None

    def setRefMol(self, oeMol, ccId, title=None, imagePath=None):
        try:
            self._refId = ccId
            self._refMol = oeMol
            #
            if title is not None:
                self._refMol.SetTitle(title)
                self._refTitle = title
            else:
                self._refMol.SetTitle(self._refId)
                self._refTitle = None
            #
            self._refImagePath = imagePath if imagePath is not None else self._refId + ".svg"
            return True
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                traceback.print_exc(file=self.__lfh)
        return False

    def setRefId(self, ccId, title=None, imagePath=None, suppressHydrogens=False, cachePath="/data/components/ligand-dict-v3"):
        """Set the query reference molecule for MCSS comparison using the input chemical component ID.
        It is assumed that the definition for this ID can be obtained from the chemical component
        repository (in cachePath).

        Once the reference molecule is built, the MCSS calculation is initialized.

        A title and imagePath are optionally provided otherwise the component Id will be used for title and
        the basis of the image file name.

        The hydrogen flag can be used to perform the MCSS using only heavy atoms.
        """
        try:
            self._refId = ccId
            ccIdU = ccId.upper()
            self._refPath = os.path.join(cachePath, ccIdU[0], ccIdU, ccIdU + ".cif")

            _id, self._refMol, _refIsoSMILES = self.__makeMolfromCCDef(path=self._refPath, suppressHydrogens=suppressHydrogens)
            #
            # Insert title here -
            if title is not None:
                self._refMol.SetTitle(title)
                self._refTitle = title
            else:
                self._refMol.SetTitle(self._refId)
                self._refTitle = None
            #
            self._refImagePath = imagePath if imagePath is not None else self._refId + ".svg"
            return True
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                traceback.print_exc(file=self.__lfh)
        return False

    def setRefPath(self, refId, ccPath, title=None, imagePath=None, suppressHydrogens=False):
        """Set the query reference molecule for MCSS comparison using the input file path.

        Once the reference molecule is built, the MCSS calculation is initialized.

        A title is optionally provided otherwise the component Id will be used.

        The hydrogen flag can be used to perform the MCSS using only heavy atoms.
        """
        try:
            self._refPath = ccPath
            self._refId = refId
            (_rId, self._refMol, _refIsoSMILES) = self.__makeMolfromCCDef(path=self._refPath, suppressHydrogens=suppressHydrogens)

            # Insert title here -
            if title is not None:
                self._refMol.SetTitle(title)
                self._refTitle = title
            else:
                self._refMol.SetTitle(self._refId)
                self._refTitle = self._refId

            self._refImagePath = imagePath if imagePath is not None else self._refId + ".svg"
            return True
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                traceback.print_exc(file=self.__lfh)
        return False

    def setFitMol(self, oeMol, ccId, title=None, imagePath=None):
        try:
            self._fitId = ccId
            self._fitMol = oeMol
            #
            if title is not None:
                self._fitMol.SetTitle(title)
                self._fitTitle = title
            else:
                self._fitMol.SetTitle(self._fitId)
                self._fitTitle = None
            #
            self._fitImagePath = imagePath if imagePath is not None else self._fitId + ".svg"
            return True
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                traceback.print_exc(file=self.__lfh)
        return False

    def setFitId(self, ccId, title=None, imagePath=None, suppressHydrogens=False, cachePath="/data/components/ligand-dict-v3"):
        """Set the test/fit molecule for MCSS comparison using the input chemical component ID.
        It is assumed that the definition for this ID can be obtained from the chemical component
        repository (in cachePath).

        Once the test/fit molecule is built, the MCSS calculation is initialized.

        A title and imagePath are optionally provided otherwise the component Id will be used for title and
        the basis of the image file name.

        The hydrogen flag can be used to perform the MCSS using only heavy atoms.
        """
        try:
            self._fitId = ccId
            ccIdU = ccId.upper()
            self._fitPath = os.path.join(cachePath, ccIdU[0], ccIdU, ccIdU + ".cif")
            self._fitId, self._fitMol, _fitIsoSMILES = self.__makeMolfromCCDef(path=self._fitPath, suppressHydrogens=suppressHydrogens)
            if title is not None:
                self._fitMol.SetTitle(title)
                self._fitTitle = title
            else:
                self._fitMol.SetTitle(self._fitId)
                self._fitTitle = ""
            self._fitImagePath = imagePath if imagePath is not None else self._fitId + ".svg"
            return True
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                traceback.print_exc(file=self.__lfh)
        return False

    def setFitPath(self, fitId, ccPath, title=None, imagePath=None, suppressHydrogens=False):  # pylint: disable=unused-argument
        """Set the test/fit molecule for MCSS comparison using the input file path.

        A title and imagePath are optionally provided otherwise the component Id will be used for title and
        the basis of the image file name.

        The hydrogen flag can be used to perform the MCSS using only heavy atoms.
        """
        try:
            (self._fitId, self._fitMol, _fitIsoSMILES) = self.__makeMolfromCCDef(path=ccPath, suppressHydrogens=suppressHydrogens)
            if title is not None:
                self._fitMol.SetTitle(title)
                self._fitTitle = title
            else:
                self._fitMol.SetTitle(self._fitId)
                self._fitTitle = None
            self._fitImagePath = imagePath if imagePath is not None else self._fitId + ".svg"
            return True
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                traceback.print_exc(file=self.__lfh)
        return False

    def setFitIdList(self, ccIdList, cachePath="/data/components/ligand-dict-v3", suppressHydrogens=False, imageDirPath="."):
        """Set the list of IDs to be compared with reference molecule by MCSS.

        It is assumed that each definition on the list can be obtained from the chemical component
        repository (in cachePath).

        From the input ID list build the internal pair list (self._pairMolList) of
        tuples  [(refId,refPath,refTitle,refImagePath,fitId,fitPath,fitTitle,fitImagePath),(),...]
        """
        self._pairMolList = []
        try:
            for ccId in ccIdList:
                (refId, refMol, _refIsoSMILES) = self.__makeMolfromCCDef(path=self._refPath, suppressHydrogens=suppressHydrogens)
                refTitle = refId
                refImagePath = os.path.join(imageDirPath, refId.upper() + ".svg")
                #
                ccIdU = ccId.upper()
                fitPath = os.path.join(cachePath, ccIdU[0], ccIdU, ccIdU + ".cif")
                (fitId, fitMol, _fitIsoSMILES) = self.__makeMolfromCCDef(path=fitPath, suppressHydrogens=suppressHydrogens)
                fitTitle = fitId + "/" + refId
                fitImagePath = os.path.join(imageDirPath, ccIdU + ".svg")
                self._pairMolList.append((refId, refMol, refTitle, refImagePath, fitId, fitMol, fitTitle, fitImagePath))
            return True
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                traceback.print_exc(file=self.__lfh)
        return False

    def addFitPath(self, fitId, ccPath, title=None, suppressHydrogens=False, imagePath=None):
        """Add ID/path to the list of molecules to be compared with reference molecule by MCSS.

        The input ID/path is appended to the internal pair list (self._pairMolList).

        """
        if self._refMol is None:
            return False

        try:
            fitId = fitId.upper()
            fitPath = ccPath
            (_fId, fitMol, _fitIsoSMILES) = self.__makeMolfromCCDef(path=fitPath, suppressHydrogens=suppressHydrogens)
            if title is not None:
                fitMol.SetTitle(title)
                fitTitle = title
            else:
                fitMol.SetTitle(fitId)
                fitTitle = fitId + "/" + self._refId

            fitImagePath = imagePath if imagePath is not None else fitId + ".svg"

            self._pairMolList.append((self._refId, self._refMol, self._refTitle, self._refImagePath, fitId, fitMol, fitTitle, fitImagePath))

            return True
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                traceback.print_exc(file=self.__lfh)
        return False

    def setPairIdList(self, pairIdList, cachePath="/data/components/ligand-dict-v3", suppressHydrogens=False, imageDirPath="."):
        """Set the list of ID pairs to be aligned by MCSS.

        It is assumed that each definition on the list can be obtained from the chemical component
        repository (in cachePath).

        From the input pair ID list build the internal pair list (self._pairMolList).
        """
        self._pairMolList = []
        #
        try:
            for refId, fitId in pairIdList:
                ccIdU = refId.upper()
                refPath = os.path.join(cachePath, ccIdU[0], ccIdU, ccIdU + ".cif")
                (_rId, refMol, _refIsoSMILES) = self.__makeMolfromCCDef(path=refPath, suppressHydrogens=suppressHydrogens)
                refImagePath = os.path.join(imageDirPath, ccIdU + ".svg")
                #
                ccIdU = fitId.upper()
                fitPath = os.path.join(cachePath, ccIdU[0], ccIdU, ccIdU + ".cif")
                (_fId, fitMol, _fitIsoSMILES) = self.__makeMolfromCCDef(path=fitPath, suppressHydrogens=suppressHydrogens)
                #
                refTitle = refId + "/" + fitId
                fitTitle = fitId + "/" + refId
                fitImagePath = os.path.join(imageDirPath, ccIdU + ".svg")

                self._pairMolList.append((refId, refMol, refTitle, refImagePath, fitId, fitMol, fitTitle, fitImagePath))
            return True
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                traceback.print_exc(file=self.__lfh)
        return False

    def setFitPathList(self, fitPathTupList, suppressHydrogens=False, imagePath=None):
        """Add ID/path to the list of molecules to be compared with reference molecule by MCSS.

        The input ID/path is appended to the internal pair list (self._pairMolList).

        """
        if self._refMol is None:
            return False

        try:
            for fitId, fitPath, title in fitPathTupList:
                fitId = fitId.upper()
                (_fId, fitMol, _fitIsoSMILES) = self.__makeMolfromCCDef(path=fitPath, suppressHydrogens=suppressHydrogens)
                if title is not None:
                    fitMol.SetTitle(title)
                    fitTitle = title
                else:
                    fitMol.SetTitle(fitId)
                    fitTitle = fitId + "/" + self._refId

                fitImagePath = imagePath if imagePath is not None else fitId + ".svg"

                self._pairMolList.append((self._refId, self._refMol, self._refTitle, self._refImagePath, fitId, fitMol, fitTitle, fitImagePath))

            return True
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                traceback.print_exc(file=self.__lfh)
        return False

    def __makeMolfromCCDef(self, path=None, suppressHydrogens=False, use3D=False, coordType="model"):
        """Create OE molecule from the chem comp definition in the input path.

        Construction is done using the connectivity information only (i.e. no coordinate perception).
        """
        try:
            oemList = []
            oeU = OeChemCompIoUtils(verbose=self.__verbose, log=self.__lfh)
            oemList = oeU.getFromPathList([path], use3D=use3D, coordType=coordType)
            oem = oemList[0]
            ccId = oem.getTitle()
            isoSmiles = oem.getIsoSMILES()
            if self.__verbose:
                self.__lfh.write("Creating OEMOL for %s  SMILES: %s  PATH: %s\n" % (ccId, isoSmiles, path))
            if suppressHydrogens:
                return (ccId, oem.getGraphMolSuppressH(), isoSmiles)
            else:
                return (ccId, oem.getMol(), isoSmiles)
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                traceback.print_exc(file=self.__lfh)
        return (None, None, None)
        #

    def setSearchType(self, sType="default"):
        self._searchType = sType
        return self._searchType

    def _setupMCSS(self, refmol):
        """Internal initialization for the MCSS comparison."""
        self._mcss = OEMCSSearch(OEMCSType_Approximate)
        # self._mcss = OEMCSSearch(OEMCSType_Exhaustive)
        #
        if self._searchType == "default":
            atomexpr = OEExprOpts_DefaultAtoms
            bondexpr = OEExprOpts_DefaultBonds
        elif self._searchType == "relaxed":
            atomexpr = OEExprOpts_AtomicNumber
            bondexpr = 0
        elif self._searchType == "exact":
            atomexpr = OEExprOpts_ExactAtoms
            bondexpr = OEExprOpts_ExactBonds
            # OEAddExplicitHydrogens(refmol)
        else:
            atomexpr = OEExprOpts_DefaultAtoms
            bondexpr = OEExprOpts_DefaultBonds
        #
        # atomexpr = OEExprOpts_AtomicNumber|OEExprOpts_EqAromatic
        # bondexpr = 0
        #
        # atomexpr = OEExprOpts_AtomicNumber|OEExprOpts_Aromaticity
        # bondexpr = OEExprOpts_BondOrder|OEExprOpts_EqNotAromatic
        #
        self._mcss.Init(refmol, atomexpr, bondexpr)
        #
        # self._mcss.SetMCSFunc(OEMCSMaxBondsCompleteCycles())
        # self._mcss.SetMCSFunc(OEMCSMaxAtoms())
        #
        # Half of the reference molecule --
        #
        # nAtomsRef=refmol.NumAtoms()
        # self._mcss.SetMinAtoms(nAtomsRef/2)

    def _setHighlightStyleRef(self, matchObj, refMol):
        refdisp = OE2DMolDisplay(refMol, self._opts)
        #
        optInverseFit = False
        if self._params["highlightStyleRef"] == "ballAndStick":
            hstyle = OEHighlightStyle_BallAndStick
        elif self._params["highlightStyleRef"] == "stick":
            hstyle = OEHighlightStyle_Stick
        elif self._params["highlightStyleRef"] == "stickInverse":
            hstyle = OEHighlightStyle_Stick
            optInverseFit = True
        elif self._params["highlightStyleRef"] == "ballAndStickInverse":
            hstyle = OEHighlightStyle_BallAndStick
            optInverseFit = True
        else:
            hstyle = OEHighlightStyle_BallAndStick

        if self._params["highLightMatchColorRef"] == "blue":
            myHighLightMatchColor = OEBlueTint
        elif self._params["highLightMatchColorRef"] == "green":
            myHighLightMatchColor = OEGreenTint
        elif self._params["highLightMatchColorRef"] == "pink":
            myHighLightMatchColor = OEPinkTint
        else:
            myHighLightMatchColor = OEBlueTint

        if self._params["highLightNotMatchColorRef"] == "blue":
            myHighLightNotMatchColor = OEBlueTint
        elif self._params["highLightNotMatchColorRef"] == "green":
            myHighLightNotMatchColor = OEGreenTint
        elif self._params["highLightNotMatchColorRef"] == "pink":
            myHighLightNotMatchColor = OEPinkTint
        else:
            myHighLightNotMatchColor = OEBlueTint

        #
        matchedatoms = OEIsAtomMember(matchObj.GetPatternAtoms())
        matchedbonds = OEIsBondMember(matchObj.GetPatternBonds())
        if optInverseFit:
            OEAddHighlighting(refdisp, myHighLightNotMatchColor, hstyle, OENotAtom(matchedatoms), OENotBond(matchedbonds))
        else:
            OEAddHighlighting(refdisp, myHighLightMatchColor, hstyle, matchedatoms, matchedbonds)

        return refdisp

    def _setHighlightStyleFit(self, matchObj, fitMol):
        fitdisp = OE2DMolDisplay(fitMol, self._opts)

        optInverseFit = False
        if self._params["highlightStyleFit"] == "ballAndStick":
            hstyle = OEHighlightStyle_BallAndStick
        elif self._params["highlightStyleFit"] == "stick":
            hstyle = OEHighlightStyle_Stick
        elif self._params["highlightStyleFit"] == "stickInverse":
            hstyle = OEHighlightStyle_Stick
            optInverseFit = True
        elif self._params["highlightStyleFit"] == "ballAndStickInverse":
            hstyle = OEHighlightStyle_BallAndStick
            optInverseFit = True
        else:
            hstyle = OEHighlightStyle_BallAndStick

        if self._params["highLightMatchColorFit"] == "blue":
            myHighLightMatchColor = OEBlueTint
        elif self._params["highLightMatchColorFit"] == "green":
            myHighLightMatchColor = OEGreenTint
        elif self._params["highLightMatchColorFit"] == "pink":
            myHighLightMatchColor = OEPinkTint
        else:
            myHighLightMatchColor = OEBlueTint

        if self._params["highLightNotMatchColorFit"] == "blue":
            myHighLightNotMatchColor = OEBlueTint
        elif self._params["highLightNotMatchColorFit"] == "green":
            myHighLightNotMatchColor = OEGreenTint
        elif self._params["highLightNotMatchColorFit"] == "pink":
            myHighLightNotMatchColor = OEPinkTint
        else:
            myHighLightNotMatchColor = OEBlueTint

        #
        # matchedatoms = OEIsAtomMember(matchObj.GetPatternAtoms())
        # matchedbonds = OEIsBondMember(matchObj.GetPatternBonds())
        matchedatoms = OEIsAtomMember(matchObj.GetTargetAtoms())
        matchedbonds = OEIsBondMember(matchObj.GetTargetBonds())

        if optInverseFit:
            OEAddHighlighting(fitdisp, myHighLightNotMatchColor, hstyle, OENotAtom(matchedatoms), OENotBond(matchedbonds))
        else:
            OEAddHighlighting(fitdisp, myHighLightMatchColor, hstyle, matchedatoms, matchedbonds)

        return fitdisp


class OeDepictMCSAlignMulti(OeDepictAlignBase, OeDepictBase):
    """Create 2D depictions of MCSS alignments.  Targets can be chemical component identifiers
    or paths to chemical component definition files.  Inputs can be in the the form of pairs,
    lists, and pair lists of chemical component definitions.

    Output can span multiple pages each with grid layout in image formats supporting pagination.

    """

    def __init__(self, verbose=False, log=sys.stderr):
        super(OeDepictMCSAlignMulti, self).__init__(verbose=verbose, log=log)
        #
        self.__debug = False
        self.__verbose = verbose
        self.__lfh = log
        #
        self.__grid = None
        self.__gridRows = None
        self.__gridCols = None
        self.__multi = None
        self.__image = None
        self.__citer = None

    def alignPairListMulti(self, imagePath="multi.pdf"):
        aM = []
        self._params["gridCols"] = 2
        try:
            aM = self.__alignListMultiWorker(imagePath=imagePath, layout="pairs")
        except TimeoutException:
            self.__lfh.write("+OeDepictMCSAlignMulti.alignPairListMulti timeout exception\n")
            self.__lfh.write("+OeDepictMCSAlignMulti.alignPairListMulti no output written to %r\n" % imagePath)
        except Exception as e:
            self.__lfh.write("+OeDepictMCSAlignMulti.alignPairListMulti general exception %s\n" % str(e))
            if self.__verbose:
                traceback.print_exc(file=self.__lfh)
        return aM

    def alignOneWithListMulti(self, imagePath="multi.pdf"):
        aM = []
        try:
            aM = self.__alignListMultiWorker(imagePath=imagePath, layout="list")
        except TimeoutException:
            self.__lfh.write("+OeDepictMCSAlignMulti.alignOneWithListMulti timeout exception\n")
            self.__lfh.write("+OeDepictMCSAlignMulti.alignOneWithListMulti no output written to %r\n" % imagePath)
        except Exception as e:
            self.__lfh.write("+OeDepictMCSAlignMulti.alignOneWithListMulti general exception here %s\n" % str(e))
            if self.__verbose:
                traceback.print_exc(file=self.__lfh)
        return aM

    def __setupImageMulti(self):
        """Internal method to configure a multipage image."""
        #
        self.__gridRows = self._params["gridRows"]
        self.__gridCols = self._params["gridCols"]
        #
        if self._params["pageOrientation"] == "landscape":
            self.__multi = OEMultiPageImageFile(OEPageOrientation_Landscape, OEPageSize_US_Letter)
        else:
            self.__multi = OEMultiPageImageFile(OEPageOrientation_Portrait, OEPageSize_US_Letter)

        self.__newPage()

    def __newPage(self):
        """Internal method to advance to a new page in a multipage configuration."""
        rows = self.__gridRows
        cols = self.__gridCols
        self.__image = self.__multi.NewPage()
        self.__grid = OEImageGrid(self.__image, rows, cols)
        self.__grid.SetCellGap(self._params["cellGap"])
        self.__grid.SetMargins(self._params["cellMargin"])
        if self.__debug:
            self.__lfh.write("Num columns %d\n" % self.__grid.NumCols())
            self.__lfh.write("Num rows    %d\n" % self.__grid.NumRows())
        self._opts = OE2DMolDisplayOptions(self.__grid.GetCellWidth(), self.__grid.GetCellHeight(), OEScale_AutoScale)
        self._assignDisplayOptions()
        self.__citer = self.__grid.GetCells()

    # @timeout(500)
    def __alignListMultiWorker(self, imagePath="multi.pdf", layout="pairs"):
        """Working method comparing a reference molecule with a list of fit molecules.

        pairMolList = (refId,refMol,refTitle,fitId,fitMol,fitTitle)

        Map of corresponding atoms is returned.

        Image Output is in multipage layout.
        """
        #
        self.__setupImageMulti()
        #
        atomMap = []
        firstOne = True
        iCount = 0
        for (refId, refMol, _refTitle, _refImagePath, fitId, fitMol, fitTitle, _fitImagePath) in self._pairMolList:
            iCount += 1
            self.__lfh.write(
                "+OeDepictMCSAlignMulti.__alignListMultiWorker Starting match pair refId %s fitId %s count %d of %d\n" % (refId, fitId, iCount, len(self._pairMolList))
            )

            #
            OEPrepareDepiction(refMol)
            self._setupMCSS(refMol)
            #
            fitMol.SetTitle(fitTitle)
            OEPrepareDepiction(fitMol)
            #
            nAtomsRef = refMol.NumAtoms()
            nAtomsFit = fitMol.NumAtoms()
            minAtoms = min(nAtomsRef, nAtomsFit)
            mcssMinAtoms = int(minAtoms * self._minAtomMatchFraction)
            self._mcss.SetMinAtoms(mcssMinAtoms)

            # scaling
            refscale = OEGetMoleculeScale(refMol, self._opts)
            fitscale = OEGetMoleculeScale(refMol, self._opts)
            self._opts.SetScale(min(refscale, fitscale))

            unique = True
            self.__lfh.write(
                "+OeDepictMCSAlignMulti.__alignListMultiWorker refId %s fitId %s nAtomsRef %d nAtomsFit %d mcssMinAtoms %d\n" % (refId, fitId, nAtomsRef, nAtomsFit, mcssMinAtoms)
            )
            self._miter = self._mcss.Match(fitMol, unique)

            self.__lfh.write("+OeDepictMCSAlignMulti.__alignListMultiWorker mcss match completed for refId %s fitId %s\n" % (refId, fitId))
            if self._miter.IsValid():
                match = self._miter.Target()
                OEPrepareAlignedDepiction(fitMol, self._mcss.GetPattern(), match)

                # Depict reference molecule with MCS highlighting
                if (firstOne and (layout in ["list"])) or (layout in ["pairs"]):
                    firstOne = False
                    refdisp = self._setHighlightStyleRef(matchObj=match, refMol=self._mcss.GetPattern())
                    if not self.__citer.IsValid():
                        self.__newPage()
                    cell = self.__citer.Target()
                    self.__citer.Next()
                    OERenderMolecule(cell, refdisp)
                    if self._params["cellBorders"]:
                        OEDrawBorder(cell, OEPen(OEBlackPen))

                # Depict fit molecule with MCS highlighting
                fitdisp = self._setHighlightStyleFit(matchObj=match, fitMol=fitMol)
                if not self.__citer.IsValid():
                    self.__newPage()
                cell = self.__citer.Target()
                self.__citer.Next()
                OERenderMolecule(cell, fitdisp)

                if self._params["cellBorders"]:
                    OEDrawBorder(cell, OEPen(OEBlackPen))

                for mAt in match.GetAtoms():
                    atomMap.append((refId, mAt.pattern.GetName(), fitId, mAt.target.GetName()))

                self.__lfh.write("+OeDepictMCSAlignMulti.__alignListMultiWorker mcss match completed for refId %s fitId %s total map length %d \n" % (refId, fitId, len(atomMap)))

        self.__lfh.write("+OeDepictMCSAlignMulti.__alignListMultiWorker writing image %s\n" % imagePath)
        OEWriteMultiPageImage(imagePath, self.__multi)
        self.__lfh.write("+OeDepictMCSAlignMulti.__alignListMultiWorker completed with map lenth %d\n" % len(atomMap))
        #
        return atomMap


class OeDepictMCSAlign(OeDepictAlignBase, OeDepictBase):
    """Create 2D depictions of MCSS alignments.  Targets can be chemical component identifiers
    or paths to chemical component definition files.  Inputs can be in the the form of pairs,
    lists, and pair lists of chemical component definitions.

    Output is to a single page image with grid layout.
    """

    def __init__(self, verbose=False, log=sys.stderr):
        super(OeDepictMCSAlign, self).__init__(verbose=verbose, log=log)
        #
        self.__debug = True
        # self.__verbose = verbose
        self.__lfh = log
        #
        self.__image = None
        self.__grid = None
        self.__citer = None

    def __setupImage(self):
        """Internal method to configure a single pair alignment image."""
        #
        self.__image = OEImage(self._params["imageSizeX"], self._params["imageSizeY"])
        rows = self._params["gridRows"]
        cols = self._params["gridCols"]
        self.__grid = OEImageGrid(self.__image, rows, cols)
        self.__grid.SetCellGap(self._params["cellGap"])
        self.__grid.SetMargins(self._params["cellMargin"])
        if self.__debug:
            self.__lfh.write("Num columns %d\n" % self.__grid.NumCols())
            self.__lfh.write("Num rows    %d\n" % self.__grid.NumRows())
        self._opts = OE2DMolDisplayOptions(self.__grid.GetCellWidth(), self.__grid.GetCellHeight(), OEScale_AutoScale)
        self._assignDisplayOptions()
        self.__citer = self.__grid.GetCells()

    #

    def alignPair(self, imagePath="single-pair.png"):
        """Depict a single aligned ref/fit molecule pair or the first ref/fit molecule pair on the
        current _pairMolList.  Display options set for a single grid row with two columns.
        """
        self._params["gridCols"] = 2
        self._params["gridRows"] = 1
        if len(self._pairMolList) < 1:
            self._pairMolList.append((self._refId, self._refMol, self._refTitle, None, self._fitId, self._fitMol, self._fitTitle, None))
        return self.__alignListWorker(imagePath=imagePath, layout="pairs")

    def alignPairList(self, imagePath="single-pair-list.png"):
        self._params["gridCols"] = 2
        return self.__alignListWorker(imagePath=imagePath, layout="pairs")

    def alignOneWithList(self, imagePath="single-list.png"):
        return self.__alignListWorker(imagePath=imagePath, layout="list")

    # @timeout(15)
    def __alignListWorker(self, imagePath="single.pdf", layout="pairs"):
        """Working method comparing a reference molecule with a list of fit molecules.

        pairMolList = (refId,refMol,refTitle,refImgPath,fitId,fitMol,fitTitle,fitImgPath)

        Map of corresponding atoms is returned.

        Output image is a single-page with grid layout.
        """
        #
        self.__setupImage()
        #
        atomMap = []

        firstOne = True

        for (refId, refMol, _refTitle, _refImagePath, fitId, fitMol, fitTitle, _fitImagePath) in self._pairMolList:
            #
            OEPrepareDepiction(refMol)
            self._setupMCSS(refMol)
            #
            fitMol.SetTitle(fitTitle)
            OEPrepareDepiction(fitMol)
            #
            #
            nAtomsRef = refMol.NumAtoms()
            nAtomsFit = fitMol.NumAtoms()
            minAtoms = min(nAtomsRef, nAtomsFit)
            self._mcss.SetMinAtoms(int(minAtoms * self._minAtomMatchFraction))

            # scaling
            refscale = OEGetMoleculeScale(refMol, self._opts)
            fitscale = OEGetMoleculeScale(refMol, self._opts)
            self._opts.SetScale(min(refscale, fitscale))

            unique = True
            self._miter = self._mcss.Match(fitMol, unique)

            if self._miter.IsValid():
                match = self._miter.Target()
                OEPrepareAlignedDepiction(fitMol, self._mcss.GetPattern(), match)

                # Depict reference molecule with MCS highlighting
                if (firstOne and (layout in ["list"])) or (layout in ["pairs"]):
                    firstOne = False
                    refdisp = self._setHighlightStyleRef(matchObj=match, refMol=self._mcss.GetPattern())
                    if not self.__citer.IsValid():
                        break
                    cell = self.__citer.Target()
                    self.__citer.Next()
                    OERenderMolecule(cell, refdisp)
                    if self._params["cellBorders"]:
                        OEDrawBorder(cell, OEPen(OEBlackPen))

                # Depict fit molecule with MCS highlighting
                fitdisp = self._setHighlightStyleFit(matchObj=match, fitMol=fitMol)
                if not self.__citer.IsValid():
                    break
                cell = self.__citer.Target()
                self.__citer.Next()
                OERenderMolecule(cell, fitdisp)

                if self._params["cellBorders"]:
                    OEDrawBorder(cell, OEPen(OEBlackPen))

                for mAt in match.GetAtoms():
                    atomMap.append((refId, mAt.pattern.GetName(), fitId, mAt.target.GetName()))

        OEWriteImage(imagePath, self.__image)
        return atomMap


class OeDepictMCSAlignSingle(OeDepictAlignBase, OeDepictBase):
    """Create 2D depictions of MCSS alignments.  Targets can be chemical component identifiers
    or paths to chemical component definition files.  Inputs can be in the the form of pairs,
    lists, and pair lists of chemical component definitions.

    Output is to single image files with a single diagram per file.
    """

    def __init__(self, verbose=False, log=sys.stderr):
        super(OeDepictMCSAlignSingle, self).__init__(verbose=verbose, log=log)
        #
        # self.__debug = True
        # self.__verbose = verbose
        # self.__lfh = log
        self.__imageRef = None
        self.__imageFit = None
        #

    def __setupImage(self):
        """Internal method to configure a single page image."""
        #
        self.__imageRef = OEImage(self._params["imageSizeX"], self._params["imageSizeY"])
        self.__imageFit = OEImage(self._params["imageSizeX"], self._params["imageSizeY"])
        self._opts = OE2DMolDisplayOptions(self.__imageRef.GetWidth(), self.__imageRef.GetHeight(), OEScale_AutoScale)
        self._assignDisplayOptions()

    def alignPair(self):
        """Depict a single aligned ref/fit molecule pair or the first ref/fit molecule pair on the
        current _pairMolList.  Display options set for a single grid row with two columns.
        """
        if len(self._pairMolList) < 1:
            self._pairMolList.append((self._refId, self._refMol, self._refTitle, self._refImagePath, self._fitId, self._fitMol, self._fitTitle, self._fitImagePath))
        return self.__alignListWorker()

    def alignPairList(self):
        return self.__alignListWorker()

    def alignOneWithList(self):
        return self.__alignListWorker()

    @timeout(15)
    def __alignListWorker(self):
        """Working method comparing a reference molecule with a list of fit molecules.

        pairMolList = (refId,refMol,refTitle,fitId,fitMol,fitTitle)

        Map of corresponding atoms is returned.

        Output image is a single-page with grid layout.
        """
        #
        atomMap = []
        firstOne = True

        for (refId, refMol, _refTitle, refImagePath, fitId, fitMol, _fitTitle, fitImagePath) in self._pairMolList:
            #
            self.__setupImage()
            #
            OEPrepareDepiction(refMol)
            self._setupMCSS(refMol)

            OEPrepareDepiction(fitMol)
            #
            nAtomsRef = refMol.NumAtoms()
            nAtomsFit = fitMol.NumAtoms()
            minAtoms = min(nAtomsRef, nAtomsFit)
            self._mcss.SetMinAtoms(int(minAtoms * self._minAtomMatchFraction))

            # scaling
            refscale = OEGetMoleculeScale(refMol, self._opts)
            fitscale = OEGetMoleculeScale(refMol, self._opts)
            self._opts.SetScale(min(refscale, fitscale))

            unique = True
            self._miter = self._mcss.Match(fitMol, unique)
            if self._miter.IsValid():
                match = self._miter.Target()
                OEPrepareAlignedDepiction(fitMol, self._mcss.GetPattern(), match)

                # Depict reference molecule with MCS highlighting
                if firstOne:
                    firstOne = False
                    refdisp = self._setHighlightStyleRef(matchObj=match, refMol=self._mcss.GetPattern())
                    OERenderMolecule(self.__imageRef, refdisp)
                    OEWriteImage(refImagePath, self.__imageRef)

                # Depict fit molecule with MCS highlighting
                fitdisp = self._setHighlightStyleFit(matchObj=match, fitMol=fitMol)
                OERenderMolecule(self.__imageFit, fitdisp)
                OEWriteImage(fitImagePath, self.__imageFit)

                for mAt in match.GetAtoms():
                    atomMap.append((refId, mAt.pattern.GetName(), fitId, mAt.target.GetName()))

        return atomMap


class OeTestMCSAlign(OeDepictAlignBase):
    """Create 2D depictions of MCSS alignments.  Targets can be chemical component identifiers
    or paths to chemical component definition files.  Inputs can be in the the form of pairs,
    lists, and pair lists of chemical component definitions.
    """

    def __init__(self, verbose=False, log=sys.stderr, maxMatches=2048):
        super(OeTestMCSAlign, self).__init__(verbose=verbose, log=log)
        #
        # self.__debug = True
        # self.__verbose = verbose
        self.__lfh = log
        self.__maxMatches = maxMatches

    def doAlign(self):
        """Test the MCSS comparison between the current reference and fit molecules -

        Return list of corresponding atoms on success or an empty list otherwise.
        """
        atomMap = []
        self._setupMCSS(self._refMol)
        #
        nAtomsRef = self._refMol.NumAtoms()
        nAtomsFit = self._fitMol.NumAtoms()
        minAtoms = min(nAtomsRef, nAtomsFit)
        self._mcss.SetMinAtoms(int(minAtoms * self._minAtomMatchFraction))
        unique = True
        self._miter = self._mcss.Match(self._fitMol, unique)
        if self._miter.IsValid():
            match = self._miter.Target()
            for mAt in match.GetAtoms():
                atomMap.append((self._refId, mAt.pattern.GetName(), self._fitId, mAt.target.GetName()))

        return atomMap

    def __getAtomIndex(self, oeMol):
        atD = {}
        try:
            for ii, atom in enumerate(oeMol.GetAtoms()):
                atD[atom.GetName().strip()] = ii
        except:  # noqa: E722 pylint: disable=bare-except
            pass
        return atD

    def __getChargeIndex(self, oeMol):
        atD = {}
        try:
            for _ii, atom in enumerate(oeMol.GetAtoms()):
                atD[atom.GetName().strip()] = atom.GetFormalCharge()
        except:  # noqa: E722 pylint: disable=bare-except
            pass
        return atD

    def __getNeighbors(self, oeMol, atNameList):
        nL = []
        try:
            for _ii, atom in enumerate(oeMol.GetAtoms()):
                aN = atom.GetName().strip()
                if aN in atNameList:
                    tL = []
                    for nbr in atom.GetAtoms():
                        tL.append(nbr.GetName())
                        nL.append((aN, tL))
        except:  # noqa: E722 pylint: disable=bare-except
            pass
        return nL

    def doTestAlign(self):
        """Test the MCSS comparison between the current reference and fit molecules -

        Return list of corresponding atoms on success or an empty list otherwise.
        """
        self._minAtomMatchFraction = 0.9
        unique = True
        atomMap = []
        self._setupMCSS(self._refMol)
        #
        nAtomsRef = self._refMol.NumAtoms()
        nAtomsFit = self._fitMol.NumAtoms()
        #
        #
        minAtoms = min(nAtomsRef, nAtomsFit)
        self._mcss.SetMCSFunc(OEMCSMaxAtoms())
        self._mcss.SetMinAtoms(int(minAtoms * self._minAtomMatchFraction))
        unique = False

        self._miter = self._mcss.Match(self._fitMol, unique)
        if self._miter.IsValid():
            match = self._miter.Target()
            for mAt in match.GetAtoms():
                atomMap.append((self._refId, mAt.pattern.GetName(), self._fitId, mAt.target.GetName()))

        self.__lfh.write("+OeTestMCSAlign.doTestAlign nAtomsRef %d nAtomsFit %d match length %d \n" % (nAtomsRef, nAtomsFit, len(atomMap)))

        return atomMap

    def doAlignWithAnal(self):
        """Test the MCSS comparison between the current reference and fit molecules -

        Return list of corresponding atoms on success or an empty list otherwise.
        """
        self._minAtomMatchFraction = 0.9
        unique = True
        atomMap = []
        self._setupMCSS(self._refMol)
        #
        nAtomsRef = self._refMol.NumAtoms()
        nAtomsFit = self._fitMol.NumAtoms()
        #
        atomRefD = self.__getAtomIndex(self._refMol)
        atomFitD = self.__getAtomIndex(self._fitMol)
        #
        chgRefD = self.__getChargeIndex(self._refMol)
        chgFitD = self.__getChargeIndex(self._fitMol)
        #
        minAtoms = min(nAtomsRef, nAtomsFit)
        #
        self._mcss.SetMaxMatches(self.__maxMatches)
        nlim = OEGetMCSExhaustiveSearchTruncationLimit()
        self.__lfh.write("+OeTestMCSAlign.doAlignWithAnal search limit %d max matches %d \n" % (nlim, self.__maxMatches))
        self._mcss.SetMCSFunc(OEMCSMaxAtoms())
        self._mcss.SetMinAtoms(int(minAtoms * self._minAtomMatchFraction))
        unique = False

        atomRefMatchD = {}
        atomFitMatchD = {}
        self._miter = self._mcss.Match(self._fitMol, unique)
        if self._miter.IsValid():
            match = self._miter.Target()
            for mAt in match.GetAtoms():
                atomMap.append((self._refId, mAt.pattern.GetName(), self._fitId, mAt.target.GetName()))
                atomRefMatchD[mAt.pattern.GetName()] = mAt.target.GetName().strip()
                atomFitMatchD[mAt.target.GetName()] = mAt.pattern.GetName().strip()
        self.__lfh.write("+OeTestMCSAlign.doAlignWithAnal nAtomsRef %d nAtomsFit %d match length %d \n" % (nAtomsRef, nAtomsFit, len(atomMap)))

        #
        # Get unmapped reference and fit atoms -
        #
        uRefAtomList = []
        uRefNList = []
        for atN in atomRefD:
            if atN not in atomRefMatchD:
                uRefAtomList.append(atN)
        #
        if len(uRefAtomList) > 0:
            uRefNList = self.__getNeighbors(self._refMol, uRefAtomList)
        #
        uFitAtomList = []
        uFitNList = []
        for atN in atomFitD:
            if atN not in atomFitMatchD:
                uFitAtomList.append(atN)
        #
        #
        if len(uFitAtomList) > 0:
            uFitNList = self.__getNeighbors(self._fitMol, uFitAtomList)
        #
        chgDifRefD = {}
        for atN0 in chgRefD:
            if atN0 in atomRefMatchD:
                atN = atomRefMatchD[atN0]
                # if ((atN not in chgFitD) or (chgRefD[atN0] != chgFitD[atN])):
                if chgRefD[atN0] != chgFitD[atN]:
                    chgDifRefD[atN] = chgRefD[atN0]
        #
        chgDifFitD = {}
        for atN0 in chgFitD:
            if atN0 in atomFitMatchD:
                atN = atomFitMatchD[atN0]
                # if ((atN not in chgRefD) or (chgFitD[atN0] != chgRefD[atN])):
                if chgFitD[atN0] != chgRefD[atN]:
                    chgDifFitD[atN] = chgFitD[atN0]

        return atomMap, uRefAtomList, uRefNList, chgDifRefD, uFitAtomList, uFitNList, chgDifFitD

    def doAlignList(self):
        """Test MCSS comparison between the current reference molecule with a list of fit molecules.

        pairMolList = (refId,refMol,refTitle,fitId,fitMol,fitTitle)

        Map of corresponding atoms is returned.

        """
        atomMap = []
        for (refId, refMol, _refTitle, _refImagePath, fitId, fitMol, _fitTitle, _fitImagePath) in self._pairMolList:
            #
            self._setupMCSS(refMol)
            #
            #
            nAtomsRef = refMol.NumAtoms()
            nAtomsFit = fitMol.NumAtoms()
            minAtoms = min(nAtomsRef, nAtomsFit)
            self._mcss.SetMinAtoms(int(minAtoms * self._minAtomMatchFraction))

            unique = True
            self._miter = self._mcss.Match(fitMol, unique)

            if self._miter.IsValid():
                match = self._miter.Target()
                for mAt in match.GetAtoms():
                    atomMap.append((refId, mAt.pattern.GetName(), fitId, mAt.target.GetName()))

        return atomMap
