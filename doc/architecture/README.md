---
package: pyBurstAnalysisGUI
summary: gtBurst, a Tkinter/matplotlib GUI (`gtburst`) and a set of Python "commands" for downloading and analysing Fermi GBM/LLE/LAT data on GRBs and solar flares by driving the installed gt tools and pyLikelihood.
layer: python
depends_on: [sane, pyLikelihood, tkinter, matplotlib, numpy, scipy, astropy, aplpy]
used_by: []
last_reviewed_commit: c065852
---

# pyBurstAnalysisGUI Architecture

**pyBurstAnalysisGUI** is *gtBurst*, written by G. Vianello and N. Omodei
around 2012–2014 as a stand-alone package and later folded into the
Fermitools. It is pure Python (82 modules, about 20,500 lines) plus about
29 MB of data. It provides:

1. The **`gtburst` GUI** (`python/gtburst.py`). It downloads trigger
   data, builds GBM/LLE background and source spectra for XSPEC, and runs
   LAT likelihood, TS-map and SED workflows step by step.
2. **Commands** (`python/GtBurst/commands/*.py`). Each step of those
   workflows is also a command-line program, with `key=value` arguments
   and its own small parameter system (`commandDefiner`). The GUI runs
   the same command objects in-process.
3. A **task layer** (`GtBurst/dataHandling.py` and friends). It wraps the
   compiled gt tools (through `sane`'s `GtApp`) and pyLikelihood
   (`UnbinnedAnalysis`, `BinnedAnalysis`).

It is not an `st_app` tool. It has no `.par` files of its own, and the
parameter, exit-status and help conventions in `src/st_app/doc/architecture/`
do not apply to it. They do apply to the gt tools it runs.

`README`, `README.md` and `doc/README` are the original stand-alone install
notes, and they are stale. They mention Python 2.7, `pywcs`, `PYTHONPATH`
setup, and downloading Pass 7 background files. `doc/release.notes` stops
in 2014. The version string is `python/GtBurst/version.py`
(`03-01-00p2`).

## Purpose and scope

**Does:**
- **Downloads data** for a GBM trigger. GBM TTE/CSPEC/CTIME/RSP and LAT LLE
  files come from the HEASARC FTPS archive. LAT FT1/FT2 come from the FSSC
  LAT data server. LAT data around a 4FGL source can also be requested.
- **GBM/LLE spectral preparation:** custom binning (`gtllebin`), a
  polynomial background fit per channel (`gtllebkg`, interactive in
  `gtllebkgGUI`), source spectra and responses (`gtllesrc`), and
  `.xcm` scripts for XSPEC. It writes the scripts but never runs XSPEC.
- **LAT analysis:** event selection (`gtselect`/`gtmktime`), counts maps,
  XML model building from 4FGL plus diffuse templates, unbinned or binned
  likelihood with upper limits, position optimisation (`gtfindsrc`), TS maps,
  SEDs, time-resolved likelihood (`GtBurst/scripts/doTimeResolvedLike.py`),
  and navigation plots from FT2.
- **Simulation:** `gtobssim` via `gtdosimulation`. Its GUI menu entry is
  commented out (see below).

**Does not:**
- Implement any instrument response or likelihood itself. Everything
  quantitative is done by the gt tools and pyLikelihood.
- Validate results. There are no tests in this repository.

## Position in the system

| Direction | Package | Why |
|-----------|---------|-----|
| Depends on (Python import) | `sane` | `GtApp` (wraps each gt tool's par file and runs it), `gt_apps` (in `gtapps_mp`), `readXml` (in `FuncFactory`). |
| Depends on (Python import) | `pyLikelihood` | `UnbinnedAnalysis`, `BinnedAnalysis` and `UpperLimits` for the fits. The `pyLikelihood` SWIG module provides `SkyDir` (`cutout.py`) and `JulianDate`/`SolarSystem` (`sunpos.py`). `dataHandling` imports `BinnedAnalysis`/`UnbinnedAnalysis` at module load, so nothing in the task layer imports without pyLikelihood. |
| Runs (executables) | `dataSubselector` | `gtselect`, `gtmktime`. |
| Runs (executables) | `evtbin` | `gtbin`, `gtbindef` (LLE/GBM binning, counts maps). |
| Runs (executables) | `Likelihood` | `gtltcube`, `gtltsum`, `gtexpmap`, `gtexpcube2`, `gtdiffrsp`, `gtsrcmaps`, `gtlike`, `gtmodel`, `gtfindsrc`, `gttsmap`, `gtbkg`, `gtsrcprob`. |
| Runs (executables) | `rspgen` | `gtrspgen` (LAT RSP for XSPEC). |
| Runs (executables) | `observationSim` | `gtobssim` (simulation only). |
| Runtime environment | `facilities`, conda activation | `INST_DIR`, `FERMI_DIR`, `CALDB`, `PFILES` (see *Environment variables*). |
| Help text | `fermitools-fhelp` | `fhelp_files/gtburst.txt` describes the GUI only. |
| External (required) | `tkinter`, `matplotlib` (TkAgg), `numpy`, `scipy`, `astropy` (`io.fits`, `wcs`) | GUI, plotting, fitting and FITS I/O. `scipy` is not listed in the conda run requirements; it arrives transitively (installed in the FT-2.5.1 environment). |
| External (feature) | `aplpy` | Sky-map and TS-map figures (`gtdotsmap`, `gtdolike` with `showmodelimage`, `fits2png`, `InteractiveFt1Display`). A conda run requirement only `[osx or x86_64]`. |
| External (optional, guarded) | `uw` (pointlike), `BKGE_interface`, `ROOT`, `emcee`, `threeML`, `numexpr`, `feedparser`/`chardet` | `uw` is needed by `gtconvertxmlmodel` and `dataHandling.convertXML`, which disable themselves without it. `BKGE_interface` is the LAT background estimator (`bkge.active`). `ROOT` is used only on the `test=True` path of the background polynomial fit. None of these ship with the Fermitools. |
| External (programs) | `curl` | Fetches LAT query results. |
| External (programs) | `addrmf` (HEASoft) | Optional. Without it, GBM uses `.rsp` instead of time-dependent `.rsp2`, and the GUI shows a warning. |
| Used by | — | No other package imports `GtBurst` or `gtburst`. |

The Level 1 dependency table says "— (pure Python; uses the installed
tools at runtime)". That is incomplete. The package imports `sane` and
`pyLikelihood` Python modules at run time, so it needs both to be
installed. CMake does not express this.

## Key abstractions

```
 gtburst.py  GUI(metaclass=MetaForExceptions)  ── menus: File / Tasks / Tools / Update
    │   datasets: [Dataset]   object: TransientSource   configuration: Configuration
    │
    ├─ CommandPipeline(frame, commands=[Command,...], datasets, ...)   (wizard)
    │     builds an EntryPoint form per Command; on "Run":
    │       params = dataset values + user entries → command.run(**params)
    │       returned ('key', value, ...) pairs are written back into the Dataset
    │
    ├─ commands/<name>.py   module-level  thisCommand = commandDefiner.Command(...)
    │       thisCommand.addParameter(..., partype=DATASETFILE|INPUTFILE|OUTPUTFILE|
    │                                      PYTHONONLY|HIDDEN, possiblevalues=...)
    │       def run(**kwargs) → dataHandling ...
    │       thisCommand.run = run
    │       if __name__ == '__main__':  key=value argv → run()
    │
    ├─ dataHandling.py   (task layer, 4,000 lines)
    │       LLEData (GBM/LLE): gtselect/gtbin/gtbindef wrappers
    │       LATData(LLEData): performStandardCut, makeLivetimeCube, makeExposureMap,
    │                         doUnbinned/BinnedLikelihoodAnalysis, makeTSmap,
    │                         optimizeSourcePosition, doSpectralFiles, makeSimulation …
    │       CspecBackground (polynomial background per channel), Spectrum/Spectra (PHA I/O)
    │       multiprocessScienceTools → my_gtltcube, my_gtexpmap, my_gtdiffrsp, my_gttsmap
    │
    ├─ LikelihoodComponent.py   XML model assembly: 4FGL point/extended sources,
    │                            Galactic/isotropic templates, GRB point source
    ├─ IRFS.py                  PROC_VER → event classes/IRFs/diffuse templates
    ├─ dataCollector → GBMdataCollector / LLEdataCollector   (FTPS downloads)
    ├─ downloadTransientData    (LAT data server query + poll + curl)
    └─ updater                  (self-update from GitHub, see gotchas)
```

| Module | Role |
|--------|------|
| `commandDefiner.Command`, `.Parameter` | A small, self-contained parameter system. Each parameter has a name, description, mandatory flag, default, `possiblevalues` and a type. The types tell the GUI what to do: `DATASETFILE` is filled from the current `Dataset`, `INPUTFILE`/`OUTPUTFILE` get a file browser, `PYTHONONLY` is not settable from the CLI, `HIDDEN` is set by the GUI. `getHelp()` prints the parameter table. |
| `Dataset` | A `dict` for one detector's files (`eventfile`, `ft2file`, `rspfile`, `cspecfile`, …) plus trigger name, time and detector. Command outputs are added to it. |
| `TransientSource`, `Configuration` | Object name/RA/Dec/date. Settings: data directory, FTP root and CPU count (from environment variables or defaults). |
| `CommandPipeline`, `EntryPoint`, `SubWindow`, `ConsoleText` | GUI plumbing. `ConsoleText` mirrors stdout into the GUI (with `-gui`) and writes `gtburst.log` in the current directory. |
| `IRFS` | Tables of event classes, IRF names and diffuse template file names by processing version: `PROCS['202','203']` (P7REP), `'302'` (P8R2), `'305'` (P8R3). The `PROC_VER` keyword of the FT1 file selects the row. |
| `gtapps_mp/` | J. Perkins' multi-process wrappers (`gtltcube_mp`, `gtexpmap_mp`, `gtdiffrsp_mp`, `gttsmap_mp`), run as separate executables. |
| `BayesianBlocks.py` | A private NumPy implementation for `scripts/applyBayesianBlocks.py`. It does not use the `BayesianBlocks` package. |
| `html2text.py` | Vendored third-party module, used to show HTML help. |

## Data flow and control flow

### GUI session

```
gtburst [-gui]
  → fermitools.gtburst:main()
      (module import already parsed argv with argparse, imported every command,
       selected the TkAgg backend)
      GUI() → Configuration() (creates ~/FermiData if missing) → getDataPath() → main(): Tk window
  File → Download trigger datasets…
       TriggerSelector (HEASARC W3Browse fermigtrig table, optional)
       → GBMdataCollector / LLEdataCollector: FTPS heasarc.gsfc.nasa.gov
           /fermi/data/{gbm|lat}/triggers/<year>/bn<trigger>/current
       → downloadTransientData: POST to FSSC LATDataQuery.cgi, poll result page,
           curl -LO https://fermi.gsfc.nasa.gov/FTP/fermi/data/lat/queries/<file>
       → files land in $GTBURST_DATA/bn<trigger>/ → registerDatasets → [Dataset]
  Tasks → Make likelihood analysis
       CommandPipeline([gtdocountsmap, gtbuildxmlmodel, gteditxmlmodel, gtdolike])
         gtdocountsmap.run → LATData.performStandardCut (gtselect, gtmktime) → doSkyMap (gtbin)
         gtbuildxmlmodel.run → LikelihoodComponent (4FGL from data/gll_psc_v17.xml,
                               templates from $FERMI_DIR/refdata/fermi/galdiffuse)
         gteditxmlmodel.run → xmlModelGUI (Tk editor)
         gtdolike.run → LATData.doUnbinnedLikelihoodAnalysis
                        (gtltcube, gtexpmap, gtdiffrsp via GtApp; fit via UnbinnedAnalysis)
       → afterLikelihood: light curves redrawn
```

Every `GUI` method is wrapped by `MetaForExceptions`. A `GtBurstException`
becomes an error dialog. Any other exception becomes a "Snap!" dialog, and
the traceback goes to the console and `gtburst.log`. The GUI then keeps
running.

### Command line

```
python -m GtBurst.commands.gtdolike filteredeventfile=… rspfile=… ft2file=… xmlmodel=…
  thisCommand.greetings()
  dict(arg.split('=') for arg in argv[1:]) → gtdolike(**args) → run(**args)
     no arguments          → print help, return (exit 0)
     missing mandatory     → RuntimeError (exit 1)
     unknown key           → KeyError caught, help printed, return (exit 0)
```

### How gt tools are run

`GtApp(name)` (from `sane`) loads the tool's par file through `PFILES`,
stores parameter values, and runs the tool in a shell. With the default
`catchError="at the top level:"`, `GtApp.run` decides that a tool failed
by **scanning its output for the `st_app` error banner**, not by checking
the exit status. It then raises `RuntimeError`. (See
`src/st_app/doc/architecture/` for the banner and exit-status contract.)
`dataHandling.runShellCommand` (used for `curl`) ignores non-zero exit
codes altogether.

## Public interfaces

### Applications

The launchers are generated by conda from `recipe/meta.yaml`
`build.entry_points` in the umbrella repo. CMake installs no executables.

| Script | Entry point | Status |
|--------|-------------|--------|
| `gtburst` | `fermitools.gtburst:main` | Works. Accepts `-gui`, which mirrors console output into the GUI; per its own banner this can hang on Linux. Without it, output stays in the terminal. |
| `fits2png`, `fromLogToScript`, `gtbuildxmlmodel`, `gtburstInstTest`, `gtconvertxmlmodel`, `gtdocountsmap`, `gtdolike`, `gtdosed`, `gtdosimulation`, `gtdotsmap`, `gtdownloadGBMdata`, `gtdownloadLATdata`, `gtdownloadLLEdata`, `gteditxmlmodel.py`, `gteditxmlmodelsim`, `gtinteractiveRaDec`, `gtllebin`, `gtllebindef`, `gtllebkg`, `gtllebkgbindef`, `gtllebkgGUI`, `gtllesrc`, `gtllesrcbindef` (23) | `fermitools.GtBurst.commands.<name>:__main__` | **All broken.** No module defines a callable named `__main__`, so every script exits 1 with `ImportError: cannot import name '__main__'` (confirmed at runtime for `gtdolike`, `fits2png`, `gtllebin` in the FT-2.5.1 environment). Run them as `python -m GtBurst.commands.<name> key=value …`, which works (confirmed at runtime). |

`gteditxmlmodel.py` keeps a `.py` suffix in its script name, unlike the
others. `gtburstInstTest` has no `__main__` guard: importing it runs the
installation test. The four `GtBurst/scripts/*.py`
(`doTimeResolvedLike.py`, `applyBayesianBlocks.py`, `like_profile_fit.py`,
`rawdata2package.py`) are installed, but they have no launcher.

For each command's purpose, see [commands.md](commands.md).

### Python

- `python/` is installed verbatim into `site-packages/fermitools/`, giving
  `fermitools/gtburst.py` and the `fermitools/GtBurst/` package (with
  `commands`, `gtapps_mp`, `scripts` and an empty `tasks` subpackage).
- `fermitools.pth` also puts `site-packages/fermitools` itself on
  `sys.path`. That is why the code's own `from GtBurst import …` works, and
  why `GtApp`, `UnbinnedAnalysis` and the rest import as top-level names.
- The command `run(**kwargs)` functions are the scriptable API.
  `fromLogToScript` turns a `gtburst.log` into a script that calls them.

### Environment variables

| Variable | Read in | Effect |
|----------|---------|--------|
| `GTBURST_DATA` | `Configuration` | Download/data directory. Default `~/FermiData`, which is created on first use. |
| `GTBURST_FTP` | `Configuration` | Archive root. Default `ftps://heasarc.gsfc.nasa.gov/fermi/data`. |
| `GTBURST_NCPUS` | `Configuration` | Maximum processes for `gtapps_mp`. Default: all CPUs. |
| `GTBURST_TEMPLATE_PATH` | `LikelihoodComponent` | Overrides the directory searched for Galactic and isotropic templates. |
| `FERMI_DIR`, `INST_DIR`, `GLAST_EXT` | `LikelihoodComponent`, `getDataPath`, `updater` | Locate `refdata/fermi/galdiffuse` and the package data (`$INST_DIR/data/pyBurstAnalysisGUI`). Set by conda activation. See "Runtime environment" in ScienceTools' `docs/architecture/README.md`. |
| `DIFFUSE_VER` | `LikelihoodComponent` | Old SLAC-internal template layout only. |
| `CONDA_PREFIX` | `updater` | Chooses the data path. |
| `DISPLAY` | `dataHandling` | If set, forces the TkAgg backend at import. |
| `SIMDIR` | `dataHandling` (written) | Set to the working directory before `gtobssim`. |
| `THISBKGE` | `LikelihoodComponent` (written and read) | Passes the BKGE template path internally. |
| `ISOTROPIC_TEMPLATE`, `GALACTIC_DIFFUSE_TEMPLATE` | `gtburstInstTest` only | That test exits with -1 unless both are set. No other code uses them, and conda activation does not set them. |

### Network endpoints (hard-coded)

Checked on 2026-10-22 with one request each; all were reachable.

| Purpose | Endpoint | Code |
|---------|----------|------|
| GBM and LLE trigger data | `ftps://heasarc.gsfc.nasa.gov/fermi/data/{gbm,lat}/triggers/<yyyy>/bn<trig>/current` (anonymous FTP over TLS, `ftplib.FTP_TLS`) | `Configuration`, `dataCollector` |
| LAT photon data | `https://fermi.gsfc.nasa.gov/cgi-bin/ssc/LAT/LATDataQuery.cgi`. The code scrapes the "event database currently holds … MET … to … seconds" sentence, then polls the results page and downloads with `curl -LO` from `https://fermi.gsfc.nasa.gov/FTP/fermi/data/lat/queries/`. | `downloadTransientData` |
| Trigger list | HEASARC W3Browse `fermigtrig` batch query (`https://heasarc.gsfc.nasa.gov/db-perl/W3Browse/w3table.pl?...`) | `TriggerSelector` |
| Self-update | `https://raw.githubusercontent.com/fermi-lat/pyBurstAnalysisGUI/master/` | `updater` |
| Browser links | `http://heasarc.gsfc.nasa.gov/...`, `http://fermi.gsfc.nasa.gov/...` (plain HTTP, help only) | `gtburst.py` |

All of these depend on the format of HTML pages, so a change on the
server side breaks downloads without any change in this repo.

### File formats

- **Reads:** GBM TTE/CSPEC/CTIME (PHA2), RSP/RSP2, LLE FT1-like files, LAT
  FT1/FT2, `gtburst.log`.
- **Writes:** PHA1/PHA2 source and background spectra, `.xcm` XSPEC
  scripts, likelihood XML models and results, counts/TS maps (FITS and PNG),
  `gtburst.log`.
- **Package data** (`data/`, installed to `share/fermitools/data/pyBurstAnalysisGUI`):
  - `gll_psc_v17.xml`: an initial 4FGL release, 5099 sources. This is the
    catalog actually used.
  - `gll_psc_v16.xml` (3FGL) and `3fgl_extendedData.txt`: not referenced
    by the code.
  - `4fgl_extendedData.txt`.
  - `templates/`: 38 extended-source FITS templates, 12 MB.
  - `tcl_extensions/`: `fsdialog` and `msgcat`, used for the file dialog.
  - Logos and icons.

## Build and test

- **CMake:** two `install(DIRECTORY …)` rules: `python/` to
  `${FERMI_INSTALL_PYTHON}` and `data/` to
  `${FERMI_INSTALL_DATADIR}/pyBurstAnalysisGUI`
  (`share/fermitools/data/…`, *not* `refdata/fermi/…`). Nothing is built.
- **Data lookup:** `GtBurst/getDataPath.getDataPath()` does not use
  `facilities`. It tries four layouts (stand-alone, SCons, public ST
  `refdata`, and finally `$INST_DIR/data/pyBurstAnalysisGUI`), checking
  each for `glast_logo.png`. In a conda install only the last one matches,
  so `INST_DIR` must be set.
- **Executable bits:** `install(DIRECTORY)` does not keep source
  permissions, but `dataHandling` runs `gtapps_mp/gttsmap_mp.py` and
  `gtexpmap_mp.py` directly as programs. The umbrella's
  `recipe/activate.sh`/`activate.csh` work around this with a `chmod u+x`
  on **every** activation. They hard-code
  `lib/python3.12/site-packages/fermitools/GtBurst/gtapps_mp`.
- **Not built or stale:**
  - `src/gtburst.cxx`: a 2013 C++ launcher that ran
    `python $INST_DIR/python/gtburst.py`. Nothing compiles it.
  - `SConscript`: the old SCons build; it lists a `GtBurst/aplpy/`
    directory that no longer exists.
  - `azure-pipelines.yml`: the old CI.
  - The `pyBurstAnalysisGUI/spaceHolder` and `doc/placeHolder` files.
- **Developer helpers:**
  - `install.sh` copies `python/` over the `GtBurst` installed in the
    active conda environment. It is a quick way to try changes without a
    rebuild.
  - `uploadToStable.csh` bumps `version.py`'s `patch` and regenerates
    `__file_list` (the MD5 list the self-updater reads). It uses macOS
    `/sbin/md5`.
- **Tests:** none in this repository. The `recipe/tests/ST-unit-test` CI
  harness defines `check_import_gtburst` (`import gtburst`). It is listed
  only in `unitTestsWithRoot`. In the default `unitTestsWithoutRoot` list
  it is commented out ("Throws AttributeError 10/25 JSA"), so CI does not
  exercise this package at all. `recipe/meta.yaml` `test.imports` does
  not include it either. `gtburstInstTest` is a manual smoke test of the
  environment (`GtApp` construction for 16 tools, template variables, GUI
  import).
- **Branches:** `remove_gtobssim` (1 commit, 2022) is not merged. It
  makes `GtApp('gtobssim')` optional in `LATData` and `gtburstInstTest`,
  so that gtBurst still works in a build without `gtobssim`; on `master`,
  constructing `LATData` fails if `gtobssim.par` cannot be found.
  `ftc-143` (2 commits, 2023: `np.float` → `float`) is superseded, because
  `master` no longer uses `np.float`. The other remote branches have no
  extra commits.

## Invariants and gotchas

**Installation and launch**
- The 23 command console scripts do not work (see *Applications*). This
  goes beyond the `fermitools-fhelp` note that they lack help text. The
  fix is in `recipe/meta.yaml` (point at a real function) and/or in the
  command modules (add a `main()`).
- `import gtburst` has side effects. It runs `argparse` on `sys.argv`,
  which exits with status 2 on unknown arguments; this is a likely reason
  the CI import test was disabled. It also imports every command module,
  selects TkAgg, and writes progress dots to stderr. `GtBurst/updater.py`,
  imported by `gtburst`, calls `socket.setdefaulttimeout(60)` for the
  whole process.
- The package is importable under two names: `GtBurst.X`, used internally
  via `fermitools.pth`, and `fermitools.GtBurst.X`, used by the conda
  entry points. Anything imported both ways is loaded twice, with separate
  module state (for example the `thisCommand` objects).
- `Configuration()` creates `~/FermiData` (or `$GTBURST_DATA`) as a side
  effect, even for command-line use.

**Self-updater (*Update → Update to the latest version*)**
- This menu item changes directory into the **installed** `site-packages`
  directory. It then downloads `__file_list` from GitHub `master` and
  overwrites every installed file whose MD5 differs from the list, which
  includes the updater itself; it can also create new files. Finally it
  restarts the process with `os.execl`. The installed package can end up
  a mix of the conda build and GitHub `master`, outside the control of
  conda and the umbrella's submodule pin. If the environment is
  read-only, the update fails.
- `__file_list` is out of date: 13 of its 138 MD5s do not match the files
  on `master` (checked locally). An "update" would therefore re-download
  files that are already current and report them as updated (from
  reading the code; not run).

**Commands and parameters**
- CLI parsing is `dict(arg.split('=') for arg in sys.argv[1:])`. A value
  that contains `=`, or an argument without one, raises `ValueError`.
- Commands only partly signal failure through exit status. A missing
  mandatory parameter raises an exception (exit 1), but an unrecognised
  key only prints help, and so does running with no arguments (exit 0).
- Each command module creates a **module-level `thisCommand`**, and the
  GUI changes it in place. For example, `likelihoodAnalysis` sets
  `gtdocountsmap.definedParameters['irf'].possibleValues` and marks
  `ra`/`dec` as `HIDDEN`. That state persists for the rest of the process.

**Analysis defaults and expiring data**
- `IRFS.PROCS` knows only processing versions 202/203 (P7REP), 302 (P8R2)
  and 305 (P8R3). An FT1 file with any other `PROC_VER` gives a
  `KeyError` in `likelihoodAnalysis`. For FSSC data (`ORIGIN` contains
  `FSSC`), the `TRANSIENT100*` classes are removed from the choices.
- P8R3 template choices are the `galactic` list (`template_4years_P8_V2_scaled.fits`,
  `gll_iem_v07.fits`) and `iso_P8R3_*_V3_v1.txt`, all looked up in
  `$FERMI_DIR/refdata/fermi/galdiffuse`.
- Point sources come from the bundled initial 4FGL XML (`gll_psc_v17`),
  not from a current 4FGL-DR4 catalog. They are not read from the
  `fermitools-data` package.
- The `simirfs` mapping in `dataHandling` covers only P7 class names.

**Multiprocessing**
- `my_gtltcube.multiproc_run` and `my_gtdiffrsp.multiproc_run` return
  `singleproc_run()` straight away, so in practice only `gtexpmap` and
  `gttsmap` use `gtapps_mp`. The rest of those two methods is dead code.
- The `gtapps_mp` scripts launch `python …` and copy files with `cp`
  through `subprocess`.

**GUI**
- `-gui` (capture console output in the GUI) "can hang on Linux" during
  likelihood, according to the program's own banner. The default (no
  flag) leaves output in the terminal.
- *Tasks → Simulate GRB observation* and *Interactively recenter ROI* are
  commented out of the menu. Their code (`simulateObservation`,
  `recenterROI`) remains.
- XSPEC is never run. The GUI only writes `.xcm` scripts. The `fhelp`
  text's "internal calls to XSPEC" is inaccurate.

**Optional features that silently disappear**
- Without `uw` (pointlike), which is not part of the Fermitools,
  `gtconvertxmlmodel` and the gtobssim model conversion are unavailable.
  Without `BKGE_interface`, the background estimator is off. Without
  `aplpy`, which is not installed on linux-aarch64, sky-map plotting
  fails.

## Extension points

- **New command/step:**
  1. Create `GtBurst/commands/<name>.py` with a module-level
     `thisCommand = commandDefiner.Command(...)` and `addParameter` calls.
  2. Write `run(**kwargs)` that returns a flat `('key', value, …)` tuple.
     Keys that match `Dataset` keys feed later steps.
  3. Set `thisCommand.run = run` and add the `__main__` block.
  4. To use it in the GUI, import `thisCommand` in `gtburst.py` and append
     it to a `CommandPipeline` command list.
  5. For a console script, add an entry in `recipe/meta.yaml`. It must
     point at a real callable (not `:__main__`).
- **New processing version / IRFs:** add `IRF(...)` entries to
  `IRFS.IRFS` with the new `PROC_VER` string, the event class bit and the
  template names, so that `PROCS` gains the key.
- **New catalog:** replace `data/gll_psc_v17.xml` and update the file names in
  `LikelihoodComponent.py` and `FGLSelector.py`, plus the matching
  `*_extendedData.txt` and `templates/`.
- **Wrapping another gt tool:** use `GtApp('<tool>')` in `dataHandling`.
  To run it in parallel, subclass `multiprocessScienceTools`.

Check changes manually: run `gtburst` with a known trigger, and run
`python -m GtBurst.commands.<name>` for each changed command. Nothing
automated covers this package.

## Further reading

- [commands.md](commands.md): the command modules and the GUI workflows
  that use them.
- `src/fermitools-fhelp/fhelp_files/gtburst.txt`: user help for the GUI.
- `src/sane/doc/architecture/` (`GtApp`) and
  `src/pyLikelihood/doc/architecture/`: the layers this package calls.
- `src/st_app/doc/architecture/`: the error banner that `GtApp` scans for.
- `python/GtBurst/gtapps_mp/README`: the multi-process scripts.
- FSSC "gtburst" analysis thread (fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/gtburst.html).
