# Status and limitations

BlueMolecule is [vibe-coded](vibe-coding.md) and has been tested almost entirely on synthetic files. This page lists,
as honestly as possible, what has been checked, what has not, and what the program does not do yet.

## What has been checked

- **About 32,000 automated checks in 32 test suites** (22,000 of them cut example outputs as a running job would
  leave them), all passing for this release
  ([how to run them](developer-guide.md#tests)). They cover:
  - every reader and writer, including format detection and the units and conventions of each program;
  - the chemistry routines;
  - the selection language;
  - the launcher, including its local-server security checks;
  - the packaging;
  - the interface itself, in a headless browser.
- **Large files:** synthetic 600 MB Gaussian and ORCA logs, a 1.2 GB log, a 1.1 GB DCD trajectory, and synthetic
  5 GB and 10 GB Gaussian ONIOM optimisation and frequency logs, with `tests/run_all.py --big`.
- **Amber NetCDF:** besides files written by SciPy, real trajectories written by pmemd 16, sander 9 and cpptraj 6.4
  (public test files of the MDAnalysis project, not included here) were read and compared with SciPy's reading:
  the same coordinates, boxes and times.
- **Script injection:** crafted files with HTML in titles, atom and residue names, labels and file names. None of it
  may become markup in the page (`tests/browser_xss_test.js`).
- **Screenshots** of the interface, produced by the browser tests, were inspected during development.

## Not yet verified

- **Real output files.** The example files were generated to imitate each program's output. Real outputs from other
  versions, print levels, keywords and unusual job types may contain layouts the readers do not handle. In
  particular, these have not been checked against real files:
  - ORCA QM/MM output, whose wording differs between versions;
  - Gaussian ONIOM microiterations;
  - CP2K QM/MM;
  - xtb ONIOM.
- **Real huge outputs.** Condensed reading was tested on synthetic logs of up to 10 GB and on every example output.
  A real multi-gigabyte QM/MM output has not been opened with it yet.
- **The macOS app** has never been compiled or run on a Mac. It was written and checked on Linux.
- **The Windows version** has never been run on Windows. Its launcher (C# 5) and installer (Windows PowerShell)
  were written and checked on Linux. The checks cover syntax limits, registry symmetry and the page they produce,
  which is opened in Chromium the way Edge would open it. The launcher itself has not been compiled yet.
- **The Linux Qt window** has been tested against a simulated Qt in the test suite, and the app-window fallbacks in
  a headless browser. It still needs confirmation on real desktops.
- **The macromolecule views** were tested on generated protein, membrane, glycan and DNA structures, not yet on
  structures from the PDB or from simulation set-up tools.

## Known limitations

| Area | Limitation |
|---|---|
| Orbitals | Orbitals, densities and potentials are **not computed** from fchk or Molden basis sets. Isosurfaces need a cube or XSF grid from your program (`cubegen`, `orca_plot`) |
| Surfaces | No molecular surfaces (solvent-accessible or solvent-excluded) and no electrostatic potential mapped on a surface |
| Secondary structure | A simplified DSSP: α-helices and β-ladders only (no 3₁₀ or π helices, no turns), computed once per structure rather than per trajectory frame. Records in PDB or mmCIF files are used when present |
| Cartoons | Rebuilt during trajectory playback only for systems up to 20,000 atoms |
| Glycans | Identification by shape works for six-membered (pyranose) rings; unknown five-membered (furanose) rings get white symbols |
| Trajectories | Amber NetCDF 4 files (HDF5 inside) and text (ASCII) `mdcrd` trajectories are not read; convert them to NetCDF 3. The loading dialog appears from 1 GB; smaller long trajectories are thinned automatically when they would need more than about 512 MB |
| Huge outputs | Condensed reading (only what the readers use) covers Gaussian and ORCA outputs and XYZ trajectories. Huge outputs of other programs are read whole, within the JavaScript engine's limit of about 4 GB |
| Readers | No readers for Molpro or Dalton outputs (Molpro inputs can be written) |
| Gaussian inputs | ONIOM layer columns in `.gjf` files are skipped; layers are read from outputs |
| Clean-up | A generic force field for tidying sketches; not a published force field and not a physical method |
| Building | No conformer search, crystal builder or download from PubChem or the PDB |
| Windows | "Open with" hands over files up to 64 MB (larger ones: Open… or drag and drop in the window); no Open Babel yet |
| Browser-only use | The single HTML file has no Open Babel and is limited by the browser's memory |
| Language | The interface and documentation are in English only |

## Ideas for future versions

Not promised, and in no particular order:
- orbitals and densities computed from fchk and Molden files;
- molecular surfaces and electrostatic potentials;
- ASCII Amber trajectories (`mdcrd`);
- Molpro and Dalton outputs;
- conformer search and crystal tools.

Validating the readers against real outputs comes before all of these. [Reports with example files](../CONTRIBUTING.md)
help most.
