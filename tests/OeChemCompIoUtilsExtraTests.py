##
# File:    OeChemCompIoUtilsExtraTests.py
#
# Additional behaviour tests for OeChemCompIoUtils
##
"""
Additional tests for building OE molecules from chemical component definitions.
"""

import importlib
import io
import os
import sys
import unittest
from typing import List, TextIO
from unittest import mock

try:
    from openeye.oechem import OEFloatArray  # noqa: F401 pylint: disable=unused-import

    skiptests = False
except ImportError:
    skiptests = True

if not skiptests:
    from mmcif_utils.persist.PdbxPyIoAdapter import PdbxPyIoAdapter

    import wwpdb.utils.oe_util.build.OeChemCompIoUtils as ioModule
    from wwpdb.utils.oe_util.build.OeChemCompIoUtils import OeChemCompIoUtils

_MODNAME = "wwpdb.utils.oe_util.build.OeChemCompIoUtils"


@unittest.skipIf(skiptests, "Could not import openeye")
class OeChemCompIoUtilsExtraTests(unittest.TestCase):
    def setUp(self) -> None:
        self.__here = os.path.abspath(os.path.dirname(__file__))
        self.__datadir = os.path.join(self.__here, "data")
        self.__topCachePath = os.path.join(self.__here, "ligand-dict-v3")
        self.__prdCachePath = os.path.join(self.__here, "prdcc-v3")
        self.__lfh = io.StringIO()

    def __getUtil(self, topCachePath: str, verbose: bool = True) -> "OeChemCompIoUtils":
        return OeChemCompIoUtils(topCachePath=topCachePath, verbose=verbose, log=self.__lfh)  # pylint: disable=possibly-used-before-assignment

    def testPathListContent(self) -> None:
        """Molecules built from paths have expected ids, titles and SMILES"""
        oeU = self.__getUtil(self.__topCachePath)
        pathList = [os.path.join(self.__datadir, "ATP.cif"), os.path.join(self.__datadir, "GTP.cif")]
        oemList = oeU.getFromPathList(pathList, use3D=True, coordType="model")
        self.assertEqual([o.getCcId() for o in oemList], ["ATP", "GTP"])
        self.assertEqual([o.getTitle() for o in oemList], ["ATP", "GTP"])
        for oem in oemList:
            self.assertTrue(len(oem.getCanSMILES()) > 0)
            self.assertIn("@", oem.getIsoSMILES())
            self.assertNotIn("@", oem.getCanSMILES())

    def testBuild2DvsBuild3D(self) -> None:
        """2D and 3D builds give the same molecule graph"""
        oeU = self.__getUtil(self.__topCachePath)
        pathList = [os.path.join(self.__datadir, "ATP.cif")]
        m3 = oeU.getFromPathList(pathList, use3D=True, coordType="ideal")
        m2 = oeU.getFromPathList(pathList, use3D=False)
        self.assertEqual(len(m3), 1)
        self.assertEqual(len(m2), 1)
        self.assertEqual(m3[0].getCanSMILES(), m2[0].getCanSMILES())
        self.assertEqual(m3[0].getMol().NumAtoms(), m2[0].getMol().NumAtoms())

    def testSetTitleFalse3D(self) -> None:
        """setTitle=False leaves the molecule title empty (3D build)"""
        oeU = self.__getUtil(self.__topCachePath)
        oemList = oeU.getFromPathList([os.path.join(self.__datadir, "ATP.cif")], use3D=True, setTitle=False)
        self.assertEqual(len(oemList), 1)
        self.assertEqual(oemList[0].getTitle(), "")
        self.assertEqual(oemList[0].getCcId(), "ATP")

    @unittest.expectedFailure
    def testSetTitleFalse2D(self) -> None:
        """setTitle=False leaves the molecule title empty (2D build)

        Known issue: OeBuildMol.build2D() ignores its setTitle argument and always sets the title.
        """
        oeU = self.__getUtil(self.__topCachePath)
        oemList = oeU.getFromPathList([os.path.join(self.__datadir, "ATP.cif")], use3D=False, setTitle=False)
        self.assertEqual(len(oemList), 1)
        self.assertEqual(oemList[0].getTitle(), "")

    def testEmptyLists(self) -> None:
        """Empty inputs and empty ids yield empty results"""
        oeU = self.__getUtil(self.__topCachePath)
        self.assertEqual(oeU.getFromPathList([]), [])
        self.assertEqual(oeU.getFromIdList([]), [])
        self.assertEqual(oeU.getFromIdList(["", ""]), [])

    def testIdNotInCache(self) -> None:
        """Ids missing from the cache are skipped without failure"""
        oeU = self.__getUtil(self.__topCachePath)
        oemList = oeU.getFromIdList(["ZZZ", "atp", "QQQQ"], use3D=False)
        self.assertEqual([o.getCcId() for o in oemList], ["ATP"])

    def testMissingPath(self) -> None:
        """A non-existent path yields no molecules"""
        oeU = self.__getUtil(self.__topCachePath)
        self.assertEqual(oeU.getFromPathList([os.path.join(self.__datadir, "NOSUCHFILE.cif")]), [])

    def testPrdccId(self) -> None:
        """PRDCC ids are located via last-character hash directory"""
        oeU = self.__getUtil(self.__prdCachePath)
        oemList = oeU.getFromIdList(["prdcc_000010"], use3D=False)
        self.assertEqual(len(oemList), 1)
        self.assertEqual(oemList[0].getCcId(), "PRDCC_000010")

    def testIdPathConstruction(self) -> None:
        """Paths constructed from ids follow the cache layout conventions"""
        top = "/top"
        oeU = self.__getUtil(top)
        with mock.patch.object(oeU, "getFromPathList", return_value=[]) as mGet:
            ret = oeU.getFromIdList(["atp", "", "PRDCC_000123", "XXATP", "A"], use3D=False, coordType="ideal")
        self.assertEqual(ret, [])
        mGet.assert_called_once()
        pathList: List[str] = mGet.call_args[0][0]
        self.assertEqual(
            pathList,
            [
                os.path.join(top, "A", "ATP", "ATP.cif"),
                os.path.join(top, "3", "PRDCC_000123.cif"),
                os.path.join(top, "TP", "XXATP", "XXATP.cif"),
                os.path.join(top, "A", "A", "A.cif"),
            ],
        )
        self.assertEqual(mGet.call_args[1], {"use3D": False, "coordType": "ideal", "setTitle": True})

    def testReaderExceptionVerbose(self) -> None:
        """Exceptions from the reader are caught and logged when verbose"""
        oeU = self.__getUtil(self.__topCachePath, verbose=True)
        with mock.patch(_MODNAME + ".PdbxIoAdapter", side_effect=RuntimeError("boom")):
            oemList = oeU.getFromPathList([os.path.join(self.__datadir, "ATP.cif")])
        self.assertEqual(oemList, [])
        out = self.__lfh.getvalue()
        self.assertIn("Failed boom", out)
        self.assertIn("Traceback", out)

    def testReaderExceptionQuiet(self) -> None:
        """Exceptions from the reader are caught silently when not verbose"""
        oeU = self.__getUtil(self.__topCachePath, verbose=False)
        with mock.patch(_MODNAME + ".PdbxIoAdapter", side_effect=RuntimeError("boom")):
            oemList = oeU.getFromPathList([os.path.join(self.__datadir, "ATP.cif")])
        self.assertEqual(oemList, [])
        self.assertEqual(self.__lfh.getvalue(), "")

    def testPartialResultsOnException(self) -> None:
        """Molecules built before a failure are still returned"""
        oeU = self.__getUtil(self.__topCachePath)
        realAdapter = PdbxPyIoAdapter  # pylint: disable=possibly-used-before-assignment
        calls: List[int] = []

        def factory(verbose: bool, log: TextIO) -> "PdbxPyIoAdapter":
            calls.append(1)
            if len(calls) > 1:
                raise RuntimeError("second fails")
            return realAdapter(verbose, log)

        pathList = [os.path.join(self.__datadir, "ATP.cif"), os.path.join(self.__datadir, "GTP.cif")]
        with mock.patch(_MODNAME + ".PdbxIoAdapter", side_effect=factory):
            oemList = oeU.getFromPathList(pathList, use3D=False)
        self.assertEqual([o.getCcId() for o in oemList], ["ATP"])

    def testNonDarwinImportBranch(self) -> None:
        """Module import selects the Python IO adapter on non-darwin platforms"""
        try:
            with mock.patch.object(sys, "platform", "linux"):
                mod = importlib.reload(ioModule)  # pylint: disable=possibly-used-before-assignment
                self.assertIs(vars(mod)["PdbxIoAdapter"], PdbxPyIoAdapter)
        finally:
            importlib.reload(ioModule)
        # PdbxIoAdapter is an imported (not re-exported) name, so read it from the module namespace
        self.assertIs(vars(ioModule)["PdbxIoAdapter"], PdbxPyIoAdapter)


def suiteExtra() -> unittest.TestSuite:
    return unittest.TestLoader().loadTestsFromTestCase(OeChemCompIoUtilsExtraTests)


if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suiteExtra())
