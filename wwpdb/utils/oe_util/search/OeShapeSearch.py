##
# File:    OeShapeSearch.py
# Author:  J. Westbrook
# Date:    24-Jan-2012
# Version: 0.001
#
# Updated:
#    6-Jun-2016 jdw general cleanup
#
#
##
"""
Grid based shape search including chemical property coloring.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.01"

import sys

from openeye.oeshape import OEBestOverlay, OEBestOverlayScoreIter, OECalcVolume, OEColorFFType_ImplicitMillsDean, OEHighestTanimotoCombo, OESortOverlayScores


class OeShapeSearch(object):
    def __init__(self, verbose=False, log=sys.stderr):  # pylint: disable=unused-argument
        # self.__verbose = verbose
        # self.__lfh = log
        self.__refMol = None
        self.__fitMol = None
        self.__bestSearch = None
        self.__scoreiter = None

    def setRefMol(self, refMol):
        self.__refMol = refMol
        self.__bestSearch = OEBestOverlay()
        self.__bestSearch.SetRefMol(self.__refMol)
        self.__bestSearch.SetColorForceField(OEColorFFType_ImplicitMillsDean)
        self.__bestSearch.SetColorOptimize(True)
        self.__scoreiter = OEBestOverlayScoreIter()

    def setFitMol(self, fitMol):
        self.__fitMol = fitMol
        OESortOverlayScores(self.__scoreiter, self.__bestSearch.Overlay(self.__fitMol), OEHighestTanimotoCombo())
        retD = {}
        rL = []
        for score in self.__scoreiter:
            sD = {}
            sD["tanimotoCombo"] = score.GetTanimotoCombo()
            sD["tanimoto"] = score.GetShapeTanimoto()
            sD["tanimotoColor"] = score.GetColorTanimoto()
            sD["fitVolume"] = OECalcVolume(self.__fitMol, False)
            sD["refVolume"] = OECalcVolume(self.__refMol, False)
            rL.append(sD)
        retD[fitMol.GetTitle()] = rL
        #
        return retD
