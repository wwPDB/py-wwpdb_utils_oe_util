# pylint: disable=possibly-used-before-assignment
##
# File:    OeAlignDepictExtraTests.py
# Date:    13-Oct-2026
#
##
"""
Additional tests for OeAlignDepict covering setters, alternate search and
highlighting options, foreign file import and failure branches.
"""

import io
import os
import os.path
import tempfile
import unittest
from unittest import mock

try:
    from openeye.oechem import OEFloatArray  # noqa: F401 pylint: disable=unused-import

    skiptests = False
except ImportError:
    skiptests = True

if not skiptests:
    from wwpdb.utils.oe_util.oedepict.OeAlignDepict import OeDepictMCSAlign

HERE = os.path.abspath(os.path.dirname(__file__))
CACHE_PATH = os.path.join(HERE, "ligand-dict-v3")
DATA_DIR = os.path.join(HERE, "data")
MODULE = "wwpdb.utils.oe_util.oedepict.OeAlignDepict"


@unittest.skipIf(skiptests, "Could not import openeye")
class OeAlignDepictExtraTests(unittest.TestCase):
    def setUp(self) -> None:
        self.__tmpDir = tempfile.TemporaryDirectory()  # pylint: disable=consider-using-with
        self.__outDir = self.__tmpDir.name
        self.__log = io.StringIO()

    def tearDown(self) -> None:
        self.__tmpDir.cleanup()

    def __out(self, fileName: str) -> str:
        return os.path.join(self.__outDir, fileName)

    def __assertImage(self, path: str) -> None:
        self.assertTrue(os.path.exists(path), path)
        self.assertGreater(os.path.getsize(path), 0)

    def testSetIdsWithTitles(self) -> None:
        """setRefId/setFitId with explicit titles and verbose logging"""
        oed = OeDepictMCSAlign(verbose=True, log=self.__log)
        oed.setRefId("A", title="Adenosine", cachePath=CACHE_PATH, suppressHydrogens=True)
        oed.setFitId("ATP", title="ATP title", cachePath=CACHE_PATH, suppressHydrogens=True)
        fName = self.__out("titled-pair.svg")
        aML = oed.alignPair(imagePath=fName, imageX=300, imageY=300)
        self.assertGreater(len(aML), 0)
        self.__assertImage(fName)
        logText = self.__log.getvalue()
        self.assertIn("Fit ID             = ATP", logText)
        self.assertIn("SMILES (isomeric)", logText)

    def testSetFitPathNoTitle(self) -> None:
        """setFitPath without a title uses the component id"""
        oed = OeDepictMCSAlign(verbose=False, log=self.__log)
        oed.setRefPath(os.path.join(DATA_DIR, "A.cif"))
        oed.setFitPath(os.path.join(DATA_DIR, "T.cif"))
        aML = oed.alignPair(imagePath=self.__out("a-t.svg"))
        self.assertGreater(len(aML), 0)
        self.assertEqual({(r, f) for r, _ra, f, _fa in aML}, {("A", "T")})

    def testSetFitPathList(self) -> None:
        """setFitPathList builds pairs used by alignPairList"""
        oed = OeDepictMCSAlign(verbose=False, log=self.__log)
        oed.setRefPath(os.path.join(DATA_DIR, "ATP.cif"))
        fitTupList = [
            ("A", os.path.join(DATA_DIR, "A.cif"), None),
            ("GTP", os.path.join(DATA_DIR, "GTP.cif"), "GTP title"),
        ]
        oed.setFitPathList(fitTupList)
        fName = self.__out("fit-path-list.pdf")
        aML = oed.alignPairList(imagePath=fName, gridRows=1)
        self.__assertImage(fName)
        self.assertEqual({f for _r, _ra, f, _fa in aML}, {"A", "GTP"})

    def testAlignPairListHighlightOptions(self) -> None:
        """alignPairList with alternate highlight styles and colors"""
        optList = [
            ("ballAndStick", "blue", "green"),
            ("other", "pink", "blue"),
            ("stick", "other", "other"),
        ]
        for style, matchColor, notMatchColor in optList:
            oed = OeDepictMCSAlign(verbose=False, log=self.__log)
            oed.setPairIdList([("c", "dg"), ("a", "atp")], cachePath=CACHE_PATH)
            fName = self.__out("pair-list-%s-%s.pdf" % (style, matchColor))
            aML = oed.alignPairList(
                imagePath=fName,
                highLightStyle=style,
                highLightMatchColor=matchColor,
                highLightNotMatchColor=notMatchColor,
                optInverseFit=False,
            )
            self.__assertImage(fName)
            self.assertEqual({r for r, _ra, _f, _fa in aML}, {"C", "A"})

    def testUnknownSearchType(self) -> None:
        """An unrecognized search type falls back to default matching"""
        oed = OeDepictMCSAlign(verbose=False, log=self.__log)
        self.assertEqual(oed.setSearchType(sType="bogus"), "bogus")
        oed.setRefId("ATP", cachePath=CACHE_PATH)
        oed.setFitId("GTP", cachePath=CACHE_PATH)
        aML = oed.testAlign(minFrac=0.5)
        self.assertGreater(len(aML), 0)
        self.assertEqual(len(aML[0]), 7)

    def testDoAlign(self) -> None:
        """doAlign returns atom counts, SMILES and the atom map"""
        oed = OeDepictMCSAlign(verbose=False, log=self.__log)
        oed.setRefId("ATP", cachePath=CACHE_PATH)
        oed.setFitId("ATP", cachePath=CACHE_PATH)
        nRef, refSmiles, nFit, fitSmiles, aML = oed.doAlign(unique=True, minFrac=0.5)
        self.assertGreater(nRef, 0)
        self.assertGreater(nFit, 0)
        self.assertEqual(refSmiles, fitSmiles)
        self.assertGreater(len(aML), 0)
        refId, _rIdx, _rType, _rName, fitId, _fIdx, _fName = aML[0]
        self.assertEqual(refId, "ATP")
        self.assertEqual(fitId, "ATP")

    def testDoAlignAlt(self) -> None:
        """doAlignAlt returns atom counts, feature dictionaries and the atom map"""
        oed = OeDepictMCSAlign(verbose=False, log=self.__log)
        oed.setRefId("GTP", cachePath=CACHE_PATH)
        oed.setFitId("GTP", cachePath=CACHE_PATH)
        nRef, refFD, nFit, fitFD, aML = oed.doAlignAlt(minFrac=0.5)
        assert refFD is not None
        self.assertGreater(nRef, 0)
        self.assertGreater(nFit, 0)
        for ky in ["Formula", "SMILES", "SMILES_STEREO", "InChI", "InChIKey"]:
            self.assertIn(ky, refFD)
            self.assertEqual(refFD[ky], fitFD[ky])
        self.assertGreater(len(aML), 0)
        self.assertEqual(len(aML[0]), 8)
        # Atomic numbers of mapped atoms agree
        for tup in aML:
            self.assertEqual(tup[2], tup[6])

    def testDoAlignNoMatch(self) -> None:
        """doAlign returns an empty map when the minimum match size cannot be reached"""
        oed = OeDepictMCSAlign(verbose=False, log=self.__log)
        oed.setRefId("ATP", cachePath=CACHE_PATH)
        oed.setFitId("A", cachePath=CACHE_PATH)
        _nRef, _rS, _nFit, _fS, aML = oed.doAlign(minFrac=10.0)
        self.assertEqual(aML, [])

    def testRefPathSdf3DVerbose(self) -> None:
        """setRefPath with a foreign SDF file imported in 3D with verbose logging"""
        oed = OeDepictMCSAlign(verbose=True, log=self.__log)
        oed.setRefPath(os.path.join(DATA_DIR, "ATP.sdf"), title="ATP-sdf", type="SDF", importType="3D")
        oed.setFitPath(os.path.join(DATA_DIR, "ATP.cif"), title="ATP")
        aML = oed.testAlign(unique=False, minFrac=0.5)
        self.assertGreater(len(aML), 0)
        logText = self.__log.getvalue()
        self.assertIn("__getMiscFile()", logText)
        self.assertIn("OeAlignDepict.__getMiscFile - atom", logText)

    def testRefPathMiscReadFailure(self) -> None:
        """An unreadable foreign file is logged; setRefPath then fails on the missing molecule"""
        oed = OeDepictMCSAlign(verbose=False, log=self.__log)
        with self.assertRaises((AttributeError, TypeError)):
            oed.setRefPath(os.path.join(self.__outDir, "missing.sdf"), type="SDF")
        self.assertIn("Read failed for", self.__log.getvalue())

    def testRefPathMiscException(self) -> None:
        """An exception while importing a foreign file is caught and logged"""
        oed = OeDepictMCSAlign(verbose=False, log=self.__log)
        with mock.patch(MODULE + ".OeBuildMol") as mockBuild:
            mockBuild.return_value.importFile.side_effect = RuntimeError("import boom")
            with self.assertRaises((AttributeError, TypeError)):
                oed.setRefPath(os.path.join(DATA_DIR, "ATP.sdf"), type="SDF")
        self.assertIn("import boom", self.__log.getvalue())


if __name__ == "__main__":
    unittest.main()
