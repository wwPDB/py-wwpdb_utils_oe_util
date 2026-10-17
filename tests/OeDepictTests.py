# pylint: disable=possibly-used-before-assignment
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


import fnmatch
import os
import os.path
import sys
import tempfile
import unittest
from typing import Any, List, Optional, Tuple

try:
    from openeye.oechem import OEFloatArray  # noqa: F401 pylint: disable=unused-import

    skiptests = False
except ImportError:
    skiptests = True

if not skiptests:
    from wwpdb.utils.oe_util.build.OeBuildMol import OeBuildMol
    from wwpdb.utils.oe_util.build.OeChemCompIoUtils import OeChemCompIoUtils
    from wwpdb.utils.oe_util.oedepict.OeDepict import OeDepict, OeDepictMultiPage

MolTitleList = List[Tuple[str, Any, str]]


@unittest.skipIf(skiptests, "Cannot import openeye.oechem for tests")
class OeDepictTests(unittest.TestCase):
    def setUp(self) -> None:
        self.__lfh = sys.stderr
        self.__verbose = True

        self.__here = os.path.abspath(os.path.dirname(__file__))
        self.__datadir = os.path.join(self.__here, "data")
        self.__tmpDir = tempfile.TemporaryDirectory()  # pylint: disable=consider-using-with
        self.__testoutput = self.__tmpDir.name

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

        # A PRD definition which is not present in the test data
        self.__pathList2 = [os.path.join(self.__datadir, "PRD_000027.cif")]
        self.__pathPrdChemCompCVS = os.path.join(self.__here, "prdcc-v3")
        self.__pathChemCompCVS = self.__topCachePath

    def tearDown(self) -> None:
        self.__tmpDir.cleanup()

    def __assertFile(self, path: str, magic: Optional[bytes] = None) -> None:
        """Assert that the output file exists, is non-empty and optionally starts with magic"""
        self.assertTrue(os.path.isfile(path), path)
        self.assertGreater(os.path.getsize(path), 0, path)
        if magic is not None:
            with open(path, "rb") as ifh:
                self.assertEqual(ifh.read(len(magic)), magic, path)

    def __makeFromFiles(self, pathList: List[str]) -> MolTitleList:
        """Create OE molecules from the input chem comp definition path list.

        Returns a list of (title, molecule, path) tuples.
        """
        self.__lfh.write("\nStarting OeDepictTests __makeFromFiles\n")
        oeU = OeChemCompIoUtils(topCachePath=self.__pathChemCompCVS, verbose=self.__verbose, log=self.__lfh)
        oemList = oeU.getFromPathList(pathList, use3D=True, coordType="model")
        idList: List[str] = []
        for oem in oemList:
            title = oem.getTitle()
            idList.append(title)
            self.__lfh.write("Title              = %s\n" % title)
            self.__lfh.write("SMILES (canonical) = %s\n" % oem.getCanSMILES())
            self.__lfh.write("SMILES (isomeric)  = %s\n" % oem.getIsoSMILES())
        return list(zip(idList, oemList, pathList))

    def __makeFromIds(self, idList: List[str]) -> MolTitleList:
        """Create OE molecules from the input chem comp definition id list.

        Returns a list of (id, molecule, id) tuples.
        """
        self.__lfh.write("\nStarting OeDepictTests __makeFromIds\n")
        oeU = OeChemCompIoUtils(topCachePath=self.__pathChemCompCVS, verbose=self.__verbose, log=self.__lfh)
        oemList = oeU.getFromIdList(idList, use3D=True, coordType="ideal")
        self.assertEqual(len(oemList), len(idList))
        for ccId, oem in zip(idList, oemList):
            self.assertEqual(oem.getCcId(), ccId.upper())
            self.assertTrue(oem.getCanSMILES())
            self.__lfh.write("Title              = %s\n" % oem.getTitle())
            self.__lfh.write("SMILES (canonical) = %s\n" % oem.getCanSMILES())
            self.__lfh.write("SMILES (isomeric)  = %s\n" % oem.getIsoSMILES())
        return list(zip(idList, oemList, idList))

    def __getPathList(
        self, topPath: str, pattern: str = "*", excludeDirs: Optional[List[str]] = None, recurse: bool = True
    ) -> List[str]:
        """Return a list of file paths in the input topPath which satisfy the input search criteria.

        This version does not follow symbolic links.
        """
        if excludeDirs is None:
            excludeDirs = []
        pathList: List[str] = []
        try:
            names = os.listdir(topPath)
        except OSError:
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
                    pathList.extend(
                        self.__getPathList(topPath=fullname, pattern=pattern, excludeDirs=excludeDirs, recurse=recurse)
                    )

        return pathList

    def testDepictCCIdList(self) -> None:
        """Test case -  get, read, build OE molecule, and depict the molecule."""
        self.__lfh.write("\nStarting OeDepictTests testDepicCCIdList\n")
        oeMolTitleList = self.__makeFromIds(self.__idList)
        self.assertEqual(len(oeMolTitleList), 7)
        for ccId, mol, _title in oeMolTitleList:
            imagePath = os.path.join(self.__testoutput, ccId + ".svg")
            oed = OeDepict(verbose=self.__verbose, log=self.__lfh)
            oed.setMolTitleList([(ccId, mol, "")])
            oed.setDisplayOptions(
                labelAtomName=False,
                labelAtomCIPStereo=True,
                labelAtomIndex=False,
                labelBondIndex=False,
                cellBorders=False,
                bondDisplayWidth=0.5,
            )
            oed.setGridOptions(rows=1, cols=1)
            oed.prepare()
            oed.write(imagePath)
            self.__assertFile(imagePath, b"<?xml")
        self.assertEqual(sorted(os.listdir(self.__testoutput)), sorted(ccId + ".svg" for ccId in self.__idList))

    def testDepictPrdCCPathList(self) -> None:
        """Test case -  get, read, build OE molecule, and depict the molecule."""
        self.__lfh.write("\nStarting OeDepictTests testDepictPrdCCPathList\n")
        ccPathList = self.__getPathList(
            topPath=self.__pathPrdChemCompCVS, pattern="*.cif", excludeDirs=["CVS", "REMOVED", "FULL"]
        )
        self.assertEqual(len(ccPathList), 9)
        oeMolTitleList = self.__makeFromFiles(pathList=ccPathList)
        self.assertEqual(len(oeMolTitleList), len(ccPathList))
        expected: List[str] = []
        for ccId, mol, title in oeMolTitleList:
            _dirPath, fName = os.path.split(title)
            # The molecule title is the PRDCC data block name which matches the file name
            self.assertEqual(ccId, fName[:-4])
            self.assertTrue(ccId.startswith("PRDCC_"))
            imagePath = os.path.join(self.__testoutput, fName[:-4] + ".svg")
            expected.append(fName[:-4] + ".svg")
            oed = OeDepict(verbose=self.__verbose, log=self.__lfh)
            oed.setMolTitleList([(ccId, mol, title)])
            oed.setDisplayOptions(
                labelAtomName=False,
                labelAtomCIPStereo=True,
                labelAtomIndex=False,
                labelBondIndex=False,
                bondDisplayWidth=0.5,
            )
            oed.setGridOptions(rows=1, cols=1)
            oed.prepare()
            oed.write(imagePath)
            self.__assertFile(imagePath, b"<?xml")
        self.assertEqual(sorted(os.listdir(self.__testoutput)), sorted(expected))

    def testDepictIdList(self) -> None:
        """Test case -  get, read, build OE molecule, and depict the molecule in a single image."""
        self.__lfh.write("\nStarting OeDepictTests testDepictIdList\n")
        oeMolTitleList = self.__makeFromIds(self.__idList)
        oed = OeDepict(verbose=self.__verbose, log=self.__lfh)
        oed.setMolTitleList(oeMolTitleList)
        oed.setDisplayOptions(
            imageX=1000,
            imageY=1000,
            labelAtomName=True,
            labelAtomCIPStereo=True,
            labelAtomIndex=False,
            labelBondIndex=False,
            bondDisplayWidth=0.5,
        )
        oed.setGridOptions(rows=2, cols=2)
        oed.prepare()
        imagePath = os.path.join(self.__testoutput, "myIdListtest.png")
        oed.write(imagePath)
        self.__assertFile(imagePath, b"\x89PNG\r\n\x1a\n")

    def testDepictPathList(self) -> None:
        """Test case -  get, read, build OE molecule, and depict the molecule."""
        self.__lfh.write("\nStarting OeDepictTests testDepictPathList\n")
        oeMolTitleList = self.__makeFromFiles(self.__pathList)
        self.assertEqual(len(oeMolTitleList), len(self.__pathList))
        self.assertEqual(
            [title for title, _mol, _path in oeMolTitleList],
            ["H2U", "A23", "A", "A2L", "A5O", "AET", "ATP", "DG", "A"],
        )
        oed = OeDepict(verbose=self.__verbose, log=self.__lfh)
        oed.setMolTitleList(oeMolTitleList)
        oed.setDisplayOptions(
            imageX=1500,
            imageY=1500,
            labelAtomName=True,
            labelAtomCIPStereo=True,
            labelAtomIndex=False,
            labelBondIndex=False,
            bondDisplayWidth=0.5,
        )
        oed.setGridOptions(rows=3, cols=3)
        oed.prepare()
        imagePath = os.path.join(self.__testoutput, "pathListtest.png")
        oed.write(imagePath)
        self.__assertFile(imagePath, b"\x89PNG\r\n\x1a\n")

    def testDepictIdListMulti(self) -> None:
        """Test case -  get, read, build OE molecule, and depict the molecule."""
        self.__lfh.write("\nStarting OeDepictTests testDepictIdListMulti\n")
        oeMolTitleList = self.__makeFromIds(self.__idList)
        oed = OeDepictMultiPage(verbose=self.__verbose, log=self.__lfh)
        oed.setMolTitleList(oeMolTitleList)
        oed.prepare()
        imagePath = os.path.join(self.__testoutput, "mulitIdListtest.pdf")
        oed.write(imagePath)
        self.__assertFile(imagePath, b"%PDF")

    def testDepictPathListMulti(self) -> None:
        """Test case -  get, read, build OE molecule, and depict the molecule."""
        self.__lfh.write("\nStarting OeDepictTests testDepictPathListMulti\n")
        oeMolTitleList = self.__makeFromFiles(self.__pathList)
        self.assertEqual(len(oeMolTitleList), len(self.__pathList))
        oed = OeDepictMultiPage(verbose=self.__verbose, log=self.__lfh)
        oed.setMolTitleList(oeMolTitleList)
        oed.setDisplayOptions(
            pageOrientation="Portrait",
            labelAtomName=True,
            labelAtomCIPStereo=True,
            labelAtomIndex=False,
            labelBondIndex=False,
            bondDisplayWidth=0.5,
        )
        oed.setGridOptions(rows=2, cols=1)
        oed.prepare()
        imagePath = os.path.join(self.__testoutput, "multiPathListtest.pdf")
        oed.write(imagePath)
        self.__assertFile(imagePath, b"%PDF")

    def testDepictWithErrorsMulti(self) -> None:
        """Test case -  depicting missing bits --

        The input definition does not exist, so no molecules are built and an
        empty multi-page document is written.
        """
        self.__lfh.write("\nStarting OeDepictTests testDepictWithErrorsMulti\n")
        self.assertFalse(os.path.exists(self.__pathList2[0]))
        oeMolTitleList = self.__makeFromFiles(self.__pathList2)
        self.assertEqual(oeMolTitleList, [])
        oed = OeDepictMultiPage(verbose=self.__verbose, log=self.__lfh)
        oed.setMolTitleList(oeMolTitleList)
        oed.setDisplayOptions(
            labelAtomName=True,
            labelAtomCIPStereo=True,
            labelAtomIndex=False,
            labelBondIndex=False,
            bondDisplayWidth=0.5,
        )
        oed.prepare()
        imagePath = os.path.join(self.__testoutput, "myErrortest.pdf")
        oed.write(imagePath)
        self.__assertFile(imagePath, b"%PDF")

    def testDepictOneSDF(self) -> None:
        """Test case -  get, read, build OE molecule from SDF file, and depict the molecule."""
        self.__lfh.write("\nStarting OeDepictTests testDepictOneSDF\n")
        oem = OeBuildMol(verbose=self.__verbose, log=self.__lfh)
        self.assertTrue(oem.importFile(self.__sdfFilePath, type="3D"))
        self.assertEqual(oem.getTitle(), "ATP")
        self.__lfh.write("Title              = %s\n" % oem.getTitle())
        imagePath = os.path.join(self.__testoutput, "ATP.svg")
        oed = OeDepict(verbose=self.__verbose, log=self.__lfh)
        oed.setMolTitleList([("ATP", oem, "Title for ATP")])
        oed.setDisplayOptions(
            labelAtomName=True,
            labelAtomCIPStereo=True,
            labelAtomIndex=False,
            labelBondIndex=False,
            bondDisplayWidth=0.5,
        )
        oed.setGridOptions(rows=1, cols=1)
        oed.prepare()
        oed.write(imagePath)
        self.__assertFile(imagePath, b"<?xml")
        with open(imagePath, encoding="utf-8") as ifh:
            self.assertIn("Title for ATP", ifh.read())

    def testDepictSMILES(self) -> None:
        """Test case -  create depiction from SMILES descriptor."""
        self.__lfh.write("\nStarting OeDepictTests testDepictSMILES\n")
        imagePath = os.path.join(self.__testoutput, "benzene.svg")
        oem = OeBuildMol(verbose=self.__verbose, log=self.__lfh)
        ok = oem.importSmiles("c1ccccc1")
        self.assertTrue(ok)

        oed = OeDepict(verbose=self.__verbose, log=self.__lfh)
        oed.setMolTitleList([("benzene", oem, "Title for benzene")])
        oed.setDisplayOptions(
            labelAtomName=False,
            labelAtomCIPStereo=True,
            labelAtomIndex=False,
            labelBondIndex=False,
            bondDisplayWidth=1.0,
        )
        oed.setGridOptions(rows=1, cols=1)
        oed.prepare()
        oed.write(imagePath)
        self.__assertFile(imagePath, b"<?xml")
        with open(imagePath, encoding="utf-8") as ifh:
            self.assertIn("Title for benzene", ifh.read())


def suiteSmiles() -> unittest.TestSuite:
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeDepictTests("testDepictSMILES"))
    return suiteSelect


def suiteDepict() -> unittest.TestSuite:
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeDepictTests("testDepictIdList"))
    suiteSelect.addTest(OeDepictTests("testDepictPathList"))
    return suiteSelect


def suiteDepictMulti() -> unittest.TestSuite:
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeDepictTests("testDepictIdListMulti"))
    suiteSelect.addTest(OeDepictTests("testDepictPathListMulti"))
    # suiteSelect.addTest(OeDepictTests("testDepictWithErrorsMulti"))
    return suiteSelect


def suiteDepictPrd() -> unittest.TestSuite:
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeDepictTests("testDepictPrdCCPathList"))
    return suiteSelect


def suiteDepictCC() -> unittest.TestSuite:
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeDepictTests("testDepictCCIdList"))
    return suiteSelect


def suiteDepictSDF() -> unittest.TestSuite:
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeDepictTests("testDepictOneSDF"))
    return suiteSelect


if __name__ == "__main__":
    # unittest.main()
    mySuite1 = suiteDepictMulti()
    unittest.TextTestRunner(verbosity=2).run(mySuite1)
    mySuite2 = suiteDepict()
    unittest.TextTestRunner(verbosity=2).run(mySuite2)
    mySuite3 = suiteDepictPrd()
    unittest.TextTestRunner(verbosity=2).run(mySuite3)

    mySuite3 = suiteDepictCC()
    unittest.TextTestRunner(verbosity=2).run(mySuite3)

    mySuite5 = suiteDepictSDF()
    unittest.TextTestRunner(verbosity=2).run(mySuite5)

    mySuite5 = suiteSmiles()
    unittest.TextTestRunner(verbosity=2).run(mySuite5)
