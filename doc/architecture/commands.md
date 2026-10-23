# gtBurst commands and GUI workflows

Part of the [pyBurstAnalysisGUI architecture](README.md).

Each module in `python/GtBurst/commands/` defines a module-level
`thisCommand` (a `commandDefiner.Command`) and a `run(**kwargs)` function.
The same object is used in two ways:

- **From the GUI.** `CommandPipeline` builds a form from the command's
  parameters, fills `DATASETFILE` parameters from the current `Dataset`,
  and calls `command.run(**params)` in-process. The returned
  `('key', value, ...)` pairs are written back into the `Dataset`.
- **From a shell.** Run `python -m GtBurst.commands.<name> key=value ...`.
  With no arguments it prints the parameter table. The conda console
  scripts of the same names are broken; see the README.

## GUI workflows (`gtburst.py`)

| Menu item | Method | Dataset filter | Steps |
|-----------|--------|----------------|-------|
| Tasks → Make spectra for XSPEC | `commandInterface` | GBM and LLE (`detector != "LAT"`) | `gtllesrcbindef` → `gtllebkgGUI` → `gtllesrc`, then `writeXspecScript` writes `.xcm` files |
| Tasks → Make likelihood analysis | `likelihoodAnalysis` | `detector == "LAT"` | `gtdocountsmap` → `gtbuildxmlmodel` → `gteditxmlmodel` → `gtdolike` |
| Tasks → Find source with TS map | `tsmap` | LAT | `gtdocountsmap` → `gtbuildxmlmodel` → `gteditxmlmodel` → `gtdotsmap` |
| *(commented out)* Interactively recenter ROI | `recenterROI` | LAT | `gtdocountsmap` (the `gtinteractiveRaDec` step is commented out) |
| *(commented out)* Simulate GRB observation | `simulateObservation` | LAT | `gteditxmlmodelsim` → `gtconvertxmlmodel` → `gtdosimulation` |
| Tools → Make navigation plots | `makeNavigationPlots` | LAT or LLE | `dataHandling.makeNavigationPlots(ft2, ra, dec, trigtime)` |
| File → Download trigger datasets… | `downloadDataSet` / `downloadDataSetFromFTP` | — | `GBMdataCollector`, `LLEdataCollector`, `downloadTransientData` |
| File → Download point source dataset… | `downloadDataFGLsource` | — | `FGLSelector`, then `downloadTransientData` |

## Commands

The *Main code path* column names functions in `dataHandling.py` unless
it says otherwise.

| Command | Purpose (from `shortDescription`, corrected where wrong) | Main code path |
|---------|---------------------------------------------------------|----------------|
| `gtdownloadGBMdata` | Download GBM trigger data (TTE, CSPEC, CTIME, RSP). The description wrongly says "LAT data". | `getGBMfiles.GBMdataCollector` (FTPS) |
| `gtdownloadLLEdata` | Download LAT LLE trigger data. The description wrongly says "LAT data". | `getLLEfiles.LLEdataCollector` (FTPS) |
| `gtdownloadLATdata` | Download LAT FT1/FT2 around a trigger. | `downloadTransientData` (FSSC query) |
| `gtllebin` | Rebin LLE/GBM TTE into a CSPEC, either linearly or by cloning the bins of another CSPEC. | `LLEData.binByEnergyAndTime` (`gtselect`, `gtbin`, `gtbindef`) |
| `gtllebindef`, `gtllesrcbindef`, `gtllebkgbindef` | Define time intervals (wrappers around `gtbindef`). The `src`/`bkg` variants are interactive in the GUI. | `GtApp('gtbindef')` |
| `gtllebkg` / `gtllebkgGUI` | Fit a polynomial background to each channel and write background spectra. | `CspecBackground` |
| `gtllesrc` | Extract source spectra for the chosen intervals and compute or attach the response. | `Spectra`, `RSPweight` (uses `addrmf` if present) |
| `gtdocountsmap` | Standard LAT cuts and a counts map. | `LATData.performStandardCut` (`gtselect`, `gtmktime`), `doSkyMap` (`gtbin`) |
| `gtbuildxmlmodel` | Build the likelihood XML: GRB point source, 4FGL sources and diffuse templates (or the BKGE estimator). | `LikelihoodComponent` |
| `gteditxmlmodel` | Tk editor for the XML model. The `.py` suffix is part of the console-script name. | `xmlModelGUI` |
| `gteditxmlmodelsim` | Same editor, for the simulation model. | `xmlModelGUI` |
| `gtdolike` | Unbinned or binned likelihood fit, TS, upper limits, optional position optimisation, spectral files and model image. | `LATData.doUnbinnedLikelihoodAnalysis`/`doBinnedLikelihoodAnalysis` (`gtltcube`, `gtexpmap`/`gtexpcube2`, `gtdiffrsp`, `gtsrcmaps`, pyLikelihood), `optimizeSourcePosition` (`gtfindsrc`), `doSpectralFiles` (`gtbin`, `gtrspgen`, `gtbkg`) |
| `gtdotsmap` | TS map around the source. | `LATData.makeTSmap` (`gttsmap`, through `gttsmap_mp` with several CPUs) |
| `gtdosed` | SED from energy-resolved likelihood fits. | `LATData` per energy bin |
| `gtconvertxmlmodel` | Convert a gtlike XML model to gtobssim format. Needs `uw` (pointlike), so it is unavailable in the Fermitools. | `dataHandling.convertXML` |
| `gtdosimulation` | Simulate LAT data. | `LATData.makeSimulation`, `Simulation` (`gtobssim`) |
| `gtinteractiveRaDec` | Choose RA/Dec by clicking on a counts map. Its description is a copy of `gtdocountsmap`'s. | `InteractiveFt1Display` (needs `aplpy`) |
| `fits2png` | Render a FITS image to PNG. It has no `commandDefiner` interface; it reads `sys.argv` directly. | `aplpy.FITSFigure` |
| `fromLogToScript` | Turn a `gtburst.log` into a Python script that replays the commands. | regular expressions over the log |
| `gtburstInstTest` | Environment smoke test. It runs **at import**, because it has no `__main__` guard. | `GtApp` construction for 16 tools, template variables, GUI import |

## Scripts without launchers (`python/GtBurst/scripts/`)

| Script | Purpose |
|--------|---------|
| `doTimeResolvedLike.py` | Time-resolved LAT likelihood over a list of intervals. It calls `gtdocountsmap`/`gtbuildxmlmodel`/`gtdolike` for each interval and can do an `emcee` posterior. |
| `applyBayesianBlocks.py` | Bayesian Blocks on an event list using `GtBurst.BayesianBlocks`. |
| `like_profile_fit.py` | Profile-likelihood fitting. Needs `threeML`, which is not part of the Fermitools. |
| `rawdata2package.py` | Package raw FT1/FT2 into the gtBurst dataset layout. |

Run these as `python -m GtBurst.scripts.<name>` or with
`python $CONDA_PREFIX/lib/python3.X/site-packages/fermitools/GtBurst/scripts/<name>.py`.
