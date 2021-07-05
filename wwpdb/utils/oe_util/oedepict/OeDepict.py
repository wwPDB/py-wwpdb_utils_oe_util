##
# File:  OeDepict.py
# Date:  3-May-2013  J. Westbrook
#
# Updates:
# 05-May-2013  jdw strip out IO methods -
# 21-Jan-2014  rps setGridOptions updated to support option for explicitly setting cellBorders
#  6-Jun-2016  jdw general cleanup
# 12-Jun-2016  jdw updated display parameter defaults
#
##
"""
Classes to build 2D OE molecule depictions from ChemComp definitions.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.01"


import sys

from openeye.oechem import OEBlack, OEDarkBlue, OEDarkGreen
from openeye.oedepict import (
    OE2DMolDisplay,
    OE2DMolDisplayOptions,
    OEAtomStereoStyle_Display_All,
    OEBlackPen,
    OEDisplayAtomIdx,
    OEDisplayAtomPropBase,
    OEDisplayBondIdx,
    OEDrawBorder,
    OEFill_On,
    OEFont,
    OEImage,
    OEImageGrid,
    OEMultiPageImageFile,
    OEPageOrientation_Landscape,
    OEPageOrientation_Portrait,
    OEPageSize_US_Letter,
    OEPen,
    OEPrepareDepiction,
    OERenderMolecule,
    OEScale_AutoScale,
    OEWriteImage,
    OEWriteMultiPageImage,
)


class LabelAtoms(OEDisplayAtomPropBase):
    def __init__(self):
        OEDisplayAtomPropBase.__init__(self)

    def __call__(self, atom):
        return atom.GetName()

    def CreateCopy(self):
        # __disown__ is required to allow C++ to take
        # ownsership of this object and its memory
        copy = LabelAtoms()
        return copy.__disown__()


class OeDepictBase(object):

    """Base class for 2D depictions in single and multi-page format."""

    def __init__(self, verbose=True, log=sys.stderr):  # pylint: disable=unused-argument
        super(OeDepictBase, self).__init__()
        # self.__verbose = verbose
        # self.__debug = False
        # self.__lfh = log
        #
        self._molTitleList = []
        self._opts = None
        #
        # internal dictionary of display parameters -
        #
        self._params = {
            "imageSizeX": 500,
            "imageSizeY": 500,
            "cellBorders": True,
            "suppressHydrogens": False,
            "labelAtomName": False,
            "labelAtomCIPStereo": False,
            "labelAtomIndex": False,
            "labelBondIndex": False,
            "bondDisplayWidth": None,
            "gridRows": 2,
            "gridCols": 2,
            "cellGap": 5,
            "cellMargin": 10,
            "pageOrientation": "landscape",
            "highlightStyleRef": "none",
            "highlightStyleFit": "ballAndStick",
            "highLightMatchColorRef": "blue",
            "highLightMatchColorFit": "blue",
            "highLightNotMatchColorRef": "pink",
            "highLightNotMatchColorFit": "pink",
        }

    def setDisplayOptions(self, **kwargs):
        self._params.update(kwargs)

    def setGridOptions(self, rows=1, cols=1, cellBorders=True):
        self._params["gridRows"] = rows
        self._params["gridCols"] = cols
        self._params["cellBorders"] = cellBorders

    def setMolTitleList(self, oeMolTitleList):
        """Set the list of OE Mols to be depicted as a list of tuples containing
        [(ccId,oeMol,titleString),(ccId,oeMol,titleString),...]
        """
        self._molTitleList = oeMolTitleList

    def _assignDisplayOptions(self):
        if self._params["labelAtomCIPStereo"]:
            self._opts.SetAtomStereoStyle(OEAtomStereoStyle_Display_All)

        if self._params["labelAtomIndex"]:
            self._opts.SetAtomPropertyFunctor(OEDisplayAtomIdx())
            self._opts.SetAtomPropLabelFont(OEFont(OEDarkGreen))
            self._opts.SetAtomPropLabelFontScale(0.75)

        if self._params["labelBondIndex"]:
            self._opts.SetBondPropertyFunctor(OEDisplayBondIdx())
            self._opts.SetBondPropLabelFont(OEFont(OEDarkBlue))

        if self._params["labelAtomName"]:
            atomlabel = LabelAtoms()
            self._opts.SetAtomPropertyFunctor(atomlabel)
            self._opts.SetAtomPropLabelFont(OEFont(OEDarkGreen))
            self._opts.SetAtomPropLabelFontScale(0.75)

        if self._params["bondDisplayWidth"] is not None:
            pen = OEPen(OEBlack, OEBlack, OEFill_On, self._params["bondDisplayWidth"])
            self._opts.SetDefaultBondPen(pen)
            # remove for the moment.  not supported on all platforms
            # self._opts.SetBondWidthScaling(False)
        #
        # 5.0 is minimum size -
        self._opts.SetTitleHeight(5.0)
        #


class OeDepictMultiPage(OeDepictBase):

    """Create 2D depictions in multipage format from a list of OE molecules and title strings"""

    def __init__(self, verbose=True, log=sys.stderr, useTitle=True):
        super(OeDepictMultiPage, self).__init__(verbose=verbose, log=log)
        # self.__verbose = verbose
        # self.__debug = False
        # self.__lfh = log
        self.__useTitle = useTitle
        #
        self.__multi = None
        self.__image = None

    def __setupImage(self):
        if self._params["pageOrientation"] == "landscape":
            self.__multi = OEMultiPageImageFile(OEPageOrientation_Landscape, OEPageSize_US_Letter)
        else:
            self.__multi = OEMultiPageImageFile(OEPageOrientation_Portrait, OEPageSize_US_Letter)
        self.__image = self.__multi.NewPage()
        self._opts = OE2DMolDisplayOptions()

    def prepare(self):
        self.__setupImage()
        rows = self._params["gridRows"]
        cols = self._params["gridCols"]
        grid = OEImageGrid(self.__image, rows, cols)

        citer = grid.GetCells()

        for _ccId, oeMol, title in self._molTitleList:
            if not citer.IsValid():
                # go to next page
                self.__image = self.__multi.NewPage()
                grid = OEImageGrid(self.__image, rows, cols)
                grid.SetCellGap(self._params["cellGap"])
                grid.SetMargins(self._params["cellMargin"])
                citer = grid.GetCells()

            cell = citer.Target()
            #
            if self._params["suppressHydrogens"]:
                mol = oeMol.getGraphMolSuppressH()
            else:
                mol = oeMol.getMol()

            if self.__useTitle:
                mol.SetTitle(title)
                self._opts.SetTitleHeight(5.0)
            else:
                mol.SetTitle("")
            #
            #
            OEPrepareDepiction(mol)
            self._opts.SetDimensions(cell.GetWidth(), cell.GetHeight(), OEScale_AutoScale)
            self._assignDisplayOptions()

            disp = OE2DMolDisplay(mol, self._opts)
            OERenderMolecule(cell, disp)
            OEDrawBorder(cell, OEPen(OEBlackPen))

            citer.Next()

    def write(self, imagePath):
        OEWriteMultiPageImage(imagePath, self.__multi)


class OeDepict(OeDepictBase):

    """Create 2D depictions in single-page format from a list of OE molecules & title strings"""

    def __init__(self, verbose=True, log=sys.stderr, useTitle=True):
        super(OeDepict, self).__init__(verbose=verbose, log=log)
        # self.__verbose = verbose
        self.__debug = False
        self.__lfh = log
        self.__useTitle = useTitle
        #
        self.__image = None
        self.__grid = None

    def __setupImage(self):
        """Internal method to configure a single page image."""
        #
        self.__image = OEImage(self._params["imageSizeX"], self._params["imageSizeY"])
        self.__grid = OEImageGrid(self.__image, self._params["gridRows"], self._params["gridCols"])
        self.__grid.SetCellGap(self._params["cellGap"])
        self.__grid.SetMargins(self._params["cellMargin"])
        self._opts = OE2DMolDisplayOptions(self.__grid.GetCellWidth(), self.__grid.GetCellHeight(), OEScale_AutoScale)
        #
        if self.__debug:
            self.__lfh.write("Num columns %d\n" % self.__grid.NumCols())
            self.__lfh.write("Num rows    %d\n" % self.__grid.NumRows())

    def prepare(self):
        self.__setupImage()
        # We convert to list as may come in as a python 3 zipobject with is not indexable
        mollist = list(self._molTitleList)
        for idx, cell in enumerate(self.__grid.GetCells()):
            _ccId, oeMol, title = mollist[idx]
            #
            if self._params["suppressHydrogens"]:
                mol = oeMol.getGraphMolSuppressH()
            else:
                mol = oeMol.getMol()
            if self.__useTitle:
                mol.SetTitle(title)
                self._opts.SetTitleHeight(5.0)
            else:
                mol.SetTitle("")
            #
            #
            OEPrepareDepiction(mol)
            self._opts.SetDimensions(cell.GetWidth(), cell.GetHeight(), OEScale_AutoScale)

            self._assignDisplayOptions()

            disp = OE2DMolDisplay(mol, self._opts)
            OERenderMolecule(cell, disp)
            if self._params["cellBorders"]:
                OEDrawBorder(cell, OEPen(OEBlackPen))

    def write(self, imagePath):
        OEWriteImage(imagePath, self.__image)
