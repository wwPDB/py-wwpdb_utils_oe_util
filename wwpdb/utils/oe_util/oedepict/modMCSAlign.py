#!/usr/bin/env python
# @ <SNIPPET-MCS-ALIGN>
from openeye.oechem import (
    OEExprOpts_DefaultAtoms,
    OEExprOpts_DefaultBonds,
    OEGraphMol,
    OEIsAtomMember,
    OEIsBondMember,
    OEMCSMaxBondsCompleteCycles,
    OEMCSSearch,
    OEMCSType_Approximate,
    OENotAtom,
    OENotBond,
    OEParseSmiles,
    OEPinkTint,
)
from openeye.oedepict import (
    OE2DMolDisplay,
    OE2DMolDisplayOptions,
    OEAddHighlighting,
    OEGetMoleculeScale,
    OEHighlightStyle_BallAndStick,
    OEImage,
    OEImageGrid,
    OEPrepareAlignedDepiction,
    OEPrepareDepiction,
    OERenderMolecule,
    OEScale_AutoScale,
    OETitleLocation_Hidden,
    OEWriteImage,
)

refmol = OEGraphMol()
OEParseSmiles(refmol, "c1cc(c2cc(cnc2c1)CCCO)C(=O)CCO")
OEPrepareDepiction(refmol)

fitmol = OEGraphMol()
OEParseSmiles(fitmol, "c1cc2ccc(cc2c(c1)C(=O)O)CCO")
OEPrepareDepiction(fitmol)

mcss = OEMCSSearch(OEMCSType_Approximate)
atomexpr = OEExprOpts_DefaultAtoms
bondexpr = OEExprOpts_DefaultBonds
mcss.Init(refmol, atomexpr, bondexpr)
mcss.SetMCSFunc(OEMCSMaxBondsCompleteCycles())

image = OEImage(400, 200)

rows = 1
cols = 2
grid = OEImageGrid(image, rows, cols)
refcell = grid.GetCell(1, 1)
fitcell = grid.GetCell(1, 2)

opts = OE2DMolDisplayOptions(grid.GetCellWidth(), grid.GetCellHeight(), OEScale_AutoScale)
opts.SetTitleLocation(OETitleLocation_Hidden)

refscale = OEGetMoleculeScale(refmol, opts)
fitscale = OEGetMoleculeScale(fitmol, opts)
opts.SetScale(min(refscale, fitscale))

hstyle = OEHighlightStyle_BallAndStick

unique = True
miter = mcss.Match(fitmol, unique)
if miter.IsValid():
    match = miter.Target()
    OEPrepareAlignedDepiction(fitmol, mcss.GetPattern(), match)

    # depict reference molecule with MCS highlighting
    refdisp = OE2DMolDisplay(mcss.GetPattern(), opts)
    matchedatoms = OEIsAtomMember(match.GetPatternAtoms())
    matchedbonds = OEIsBondMember(match.GetPatternBonds())
    OEAddHighlighting(refdisp, OEPinkTint, hstyle, OENotAtom(matchedatoms), OENotBond(matchedbonds))

    # abset = OEAtomBondSet(match.GetPatternAtoms(), match.GetPatternBonds())
    # OEAddHighlighting(refdisp, OEBlueTint, hstyle, abset)
    OERenderMolecule(refcell, refdisp)

    # depict fit molecule with MCS highlighting
    fitdisp = OE2DMolDisplay(fitmol, opts)

    matchedatoms = OEIsAtomMember(match.GetTargetAtoms())
    matchedbonds = OEIsBondMember(match.GetTargetBonds())
    OEAddHighlighting(fitdisp, OEPinkTint, hstyle, OENotAtom(matchedatoms), OENotBond(matchedbonds))

    # abset = OEAtomBondSet(match.GetTargetAtoms(), match.GetTargetBonds())
    # OEAddHighlighting(fitdisp, OEBlueTint, hstyle, abset)
    OERenderMolecule(fitcell, fitdisp)

OEWriteImage("modMCSAlign.png", image)
