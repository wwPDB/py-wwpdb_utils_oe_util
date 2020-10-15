##
#
# File:    OeDepictTests.py
# Author:  jdw
# Date:    5-May-2013
# Version: 0.001
#
# Updates:
#  4-May-2014 jdw add example for depiction from SMILES input
#  6-Jun-2016 jdw general cleanup
##
"""
A collection of tests for the OEDepict and related classes.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.01"


import sys
import unittest
import platform
import traceback
import os
import os.path
import fnmatch

try:
    from openeye.oechem import OEFloatArray  # noqa: F401 pylint: disable=unused-import

    skiptests = False
except ImportError:
    skiptests = True

if not skiptests:
    from wwpdb.utils.oe_util.oedepict.OeDepict import OeDepict, OeDepictMultiPage
    from wwpdb.utils.oe_util.build.OeChemCompIoUtils import OeChemCompIoUtils
    from wwpdb.utils.oe_util.build.OeBuildMol import OeBuildMol


@unittest.skipIf(skiptests, "Cannot import openeye.oechem for tests")
class OeDepictTests(unittest.TestCase):
    def setUp(self):
        self.__lfh = sys.stderr
        self.__verbose = True

        self.__here = os.path.abspath(os.path.dirname(__file__))
        self.__examples = os.path.join(self.__here, "examples")
        self.__datadir = os.path.join(self.__here, "data")
        self.__testoutput = os.path.join(self.__here, "test-output", platform.python_version())
        if not os.path.exists(self.__testoutput):
            os.makedirs(self.__testoutput)

        self.__sdfFilePath = os.path.join(self.__datadir, "ATP.sdf")
        self.__topCachePath = os.path.join(self.__here, "ligand-dict-v3")
        self.__idList = ["atp", "gtp", "A", "C", "G", "DG", "UAR"]
        self.__pathList = [
            os.path.join(self.__topCachePath, "H/H2U/H2U.cif"),
            os.path.join(self.__topCachePath, "A/A23/A23.cif"),
            os.path.join(self.__topCachePath, "A/A/A.cif"),
            os.path.join(self.__topCachePath, "A/A2L/A2L.cif"),
            os.path.join(self.__topCachePath, "A/A5O/A5O.cif"),
            os.path.join(self.__topCachePath, "A/AET/AET.cif"),
            os.path.join(self.__topCachePath, "A/ATP/ATP.cif"),
            os.path.join(self.__topCachePath, "D/DG/DG.cif"),
            os.path.join(self.__topCachePath, "A/A/A.cif"),
        ]

        self.__pathList2 = ["../data/PRD_000027.cif"]
        self.__pathPrdChemCompCVS = os.path.join(self.__here, "prdcc-v3")
        self.__pathChemCompCVS = self.__topCachePath

    def tearDown(self):
        pass

    def __getIdsFromFile(self, fPath):
        ifh = open(fPath, "r")
        idList = []
        for line in ifh:
            idList.append(line[:-1])
        ifh.close()
        return idList

    def __testMakeFromFiles(self, pathList=None):
        """Test case -  create OE molecules from the input chem comp definition path list."""
        self.__lfh.write("\nStarting OeDepictTests __testMakeFromFiles\n")
        if pathList is None:
            pathList = []
        try:
            oemList = []
            idList = []
            oeU = OeChemCompIoUtils(topCachePath=self.__pathChemCompCVS, verbose=self.__verbose, log=self.__lfh)
            oemList = oeU.getFromPathList(pathList, use3D=True, coordType="model")
            for oem in oemList:
                title = oem.getTitle()
                idList.append(title)
                self.__lfh.write("Title              = %s\n" % title)
                self.__lfh.write("SMILES (canonical) = %s\n" % oem.getCanSMILES())
                self.__lfh.write("SMILES (isomeric)  = %s\n" % oem.getIsoSMILES())
            return zip(idList, oemList, pathList)
        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def __testMakeFromIds(self, idList):
        """Test case -  create OE molecules from the input chem comp definition id list."""
        self.__lfh.write("\nStarting OeDepictTests __testMakeFromIds\n")
        try:
            oemList = []
            oeU = OeChemCompIoUtils(topCachePath=self.__pathChemCompCVS, verbose=self.__verbose, log=self.__lfh)
            oemList = oeU.getFromIdList(idList, use3D=True, coordType="ideal")
            for oem in oemList:
                self.__lfh.write("Title              = %s\n" % oem.getTitle())
                self.__lfh.write("SMILES (canonical) = %s\n" % oem.getCanSMILES())
                self.__lfh.write("SMILES (isomeric)  = %s\n" % oem.getIsoSMILES())
            return zip(idList, oemList, idList)
        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def __getPathList(self, topPath, pattern="*", excludeDirs=None, recurse=True):
        """Return a list of file paths in the input topPath which satisfy the input search criteria.

        This version does not follow symbolic links.
        """
        if excludeDirs is None:
            excludeDirs = []
        pathList = []
        #
        try:
            names = os.listdir(topPath)
        except os.error:
            return pathList

        # expand pattern
        pattern = pattern or "*"
        patternList = str.split(pattern, ";")

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
                    pathList.extend(self.__getPathList(topPath=fullname, pattern=pattern, excludeDirs=excludeDirs, recurse=recurse))

        return pathList

    def testDepictCCIdList(self):
        """Test case -  get, read, build OE molecule, and depict the molecule."""
        self.__lfh.write("\nStarting OeDepictTests testDepicCCIdList\n")
        try:
            # idList = self.__getIdsFromFile('IDLIST.list')
            idList = self.__idList
            oeMolTitleZip = self.__testMakeFromIds(idList)
            oeMolTitleList = list(oeMolTitleZip)
            if self.__verbose:
                self.__lfh.write("molTitleList length is %d\n" % len(oeMolTitleList))
            #
            for ccId, mol, title in oeMolTitleList:
                # dirPath, fName = os.path.split(title)
                imagePath = os.path.join(self.__testoutput, ccId + ".svg")
                oed = OeDepict(verbose=self.__verbose, log=self.__lfh)
                title = ""
                oed.setMolTitleList([(ccId, mol, title)])
                oed.setDisplayOptions(labelAtomName=False, labelAtomCIPStereo=True, labelAtomIndex=False, labelBondIndex=False, cellBorders=False, bondDisplayWidth=0.5)
                oed.setGridOptions(rows=1, cols=1)
                oed.prepare()
                oed.write(imagePath)
        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testDepictPrdCCPathList(self):
        """Test case -  get, read, build OE molecule, and depict the molecule."""
        self.__lfh.write("\nStarting OeDepictTests testDepictPrdCCPathList\n")
        try:
            ccPathList = self.__getPathList(topPath=self.__pathPrdChemCompCVS, pattern="*.cif", excludeDirs=["CVS", "REMOVED", "FULL"])
            oeMolTitleZip = self.__testMakeFromFiles(pathList=ccPathList)
            oeMolTitleList = list(oeMolTitleZip)
            if self.__verbose:
                self.__lfh.write("molTitleList length is %d\n" % len(oeMolTitleList))
            #
            for ccId, mol, title in oeMolTitleList:
                _dirPath, fName = os.path.split(title)
                imagePath = os.path.join(self.__testoutput, fName[:-4] + ".svg")
                oed = OeDepict(verbose=self.__verbose, log=self.__lfh)
                oed.setMolTitleList([(ccId, mol, title)])
                oed.setDisplayOptions(labelAtomName=False, labelAtomCIPStereo=True, labelAtomIndex=False, labelBondIndex=False, bondDisplayWidth=0.5)
                oed.setGridOptions(rows=1, cols=1)
                oed.prepare()
                oed.write(imagePath)
        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testDepictIdList(self):
        """Test case -  get, read, build OE molecule, and depict the molecule in a single image."""
        self.__lfh.write("\nStarting OeDepictTests testDepictIdList\n")
        try:
            oeMolTitleList = self.__testMakeFromIds(self.__idList)
            oed = OeDepict(verbose=self.__verbose, log=self.__lfh)
            oed.setMolTitleList(oeMolTitleList)
            oed.setDisplayOptions(imageX=1000, imageY=1000, labelAtomName=True, labelAtomCIPStereo=True, labelAtomIndex=False, labelBondIndex=False, bondDisplayWidth=0.5)
            oed.setGridOptions(rows=2, cols=2)
            oed.prepare()
            oed.write(os.path.join(self.__testoutput, "myIdListtest.png"))
        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testDepictPathList(self):
        """Test case -  get, read, build OE molecule, and depict the molecule."""
        self.__lfh.write("\nStarting OeDepictTests testDepictPathList\n")
        try:
            oeMolTitleList = self.__testMakeFromFiles(self.__pathList)
            oed = OeDepict(verbose=self.__verbose, log=self.__lfh)
            oed.setMolTitleList(oeMolTitleList)
            oed.setDisplayOptions(imageX=1500, imageY=1500, labelAtomName=True, labelAtomCIPStereo=True, labelAtomIndex=False, labelBondIndex=False, bondDisplayWidth=0.5)
            oed.setGridOptions(rows=3, cols=3)
            oed.prepare()
            oed.write(os.path.join(self.__testoutput, "pathListtest.png"))
        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testDepictIdListMulti(self):
        """Test case -  get, read, build OE molecule, and depict the molecule."""
        self.__lfh.write("\nStarting OeDepictTests testDepictIdListMulti\n")
        try:
            oeMolTitleList = self.__testMakeFromIds(self.__idList)
            oed = OeDepictMultiPage(verbose=self.__verbose, log=self.__lfh)
            oed.setMolTitleList(oeMolTitleList)
            oed.prepare()
            oed.write(os.path.join(self.__testoutput, "mulitIdListtest.pdf"))
        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testDepictPathListMulti(self):
        """Test case -  get, read, build OE molecule, and depict the molecule."""
        self.__lfh.write("\nStarting OeDepictTests testDepictPathListMulti\n")
        try:
            oeMolTitleList = self.__testMakeFromFiles(self.__pathList)
            oed = OeDepictMultiPage(verbose=self.__verbose, log=self.__lfh)
            oed.setMolTitleList(oeMolTitleList)
            oed.setDisplayOptions(pageOrientation="Portrait", labelAtomName=True, labelAtomCIPStereo=True, labelAtomIndex=False, labelBondIndex=False, bondDisplayWidth=0.5)
            oed.setGridOptions(rows=2, cols=1)
            oed.prepare()
            oed.write(os.path.join(self.__testoutput, "multiPathListtest.pdf"))
        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testDepictWithErrorsMulti(self):
        """Test case -  depicting missing bits --"""
        self.__lfh.write("\nStarting OeDepictTests testDepictWithErrorsMulti\n")
        try:
            oeMolTitleList = self.__testMakeFromFiles(self.__pathList2)
            oed = OeDepictMultiPage(verbose=self.__verbose, log=self.__lfh)
            oed.setMolTitleList(oeMolTitleList)
            oed.setDisplayOptions(labelAtomName=True, labelAtomCIPStereo=True, labelAtomIndex=False, labelBondIndex=False, bondDisplayWidth=0.5)
            oed.prepare()
            oed.write(os.path.join(self.__testoutput, "myErrortest.pdf"))
        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testDepictOneSDF(self):
        """Test case -  get, read, build OE molecule from SDF file, and depict the molecule."""
        self.__lfh.write("\nStarting OeDepictTests testDepictOneSDF\n")
        try:
            oem = OeBuildMol(verbose=self.__verbose, log=self.__lfh)
            if oem.importFile(self.__sdfFilePath, type="3D"):
                self.__lfh.write("Title              = %s\n" % oem.getTitle())
            #
            imagePath = os.path.join(self.__testoutput, "ATP.svg")
            oed = OeDepict(verbose=self.__verbose, log=self.__lfh)
            oed.setMolTitleList([("ATP", oem, "Title for ATP")])
            oed.setDisplayOptions(labelAtomName=True, labelAtomCIPStereo=True, labelAtomIndex=False, labelBondIndex=False, bondDisplayWidth=0.5)
            oed.setGridOptions(rows=1, cols=1)
            oed.prepare()
            oed.write(imagePath)
        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testDepictSMILES(self):
        """Test case -  create depiction from SMILES descriptor."""
        self.__lfh.write("\nStarting OeDepictTests testDepictSMILES\n")
        try:
            imagePath = os.path.join(self.__testoutput, "benzene.svg")
            oem = OeBuildMol(verbose=self.__verbose, log=self.__lfh)
            ok = oem.importSmiles("c1ccccc1")
            self.assertTrue(ok)

            oed = OeDepict(verbose=self.__verbose, log=self.__lfh)
            oed.setMolTitleList([("benzene", oem, "Title for benzene")])
            oed.setDisplayOptions(labelAtomName=False, labelAtomCIPStereo=True, labelAtomIndex=False, labelBondIndex=False, bondDisplayWidth=1.0)
            oed.setGridOptions(rows=1, cols=1)
            oed.prepare()
            oed.write(imagePath)
        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()


def suiteSmiles():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeDepictTests("testDepictSMILES"))
    return suiteSelect


def suiteDepict():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeDepictTests("testDepictIdList"))
    suiteSelect.addTest(OeDepictTests("testDepictPathList"))
    return suiteSelect


def suiteDepictMulti():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeDepictTests("testDepictIdListMulti"))
    suiteSelect.addTest(OeDepictTests("testDepictPathListMulti"))
    # suiteSelect.addTest(OeDepictTests("testDepictWithErrorsMulti"))
    return suiteSelect


def suiteDepictPrd():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeDepictTests("testDepictPrdCCPathList"))
    return suiteSelect


def suiteDepictCC():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeDepictTests("testDepictCCIdList"))
    return suiteSelect


def suiteDepictSDF():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeDepictTests("testDepictOneSDF"))
    return suiteSelect


if __name__ == "__main__":
    # unittest.main()
    mySuite1 = suiteDepictMulti()
    unittest.TextTestRunner(verbosity=2).run(mySuite1)
    #
    mySuite2 = suiteDepict()
    unittest.TextTestRunner(verbosity=2).run(mySuite2)
    #
    mySuite3 = suiteDepictPrd()
    unittest.TextTestRunner(verbosity=2).run(mySuite3)

    mySuite3 = suiteDepictCC()
    unittest.TextTestRunner(verbosity=2).run(mySuite3)

    mySuite5 = suiteDepictSDF()
    unittest.TextTestRunner(verbosity=2).run(mySuite5)

    mySuite5 = suiteSmiles()
    unittest.TextTestRunner(verbosity=2).run(mySuite5)

    mySuite3 = suiteDepictCC()
    unittest.TextTestRunner(verbosity=2).run(mySuite3)
