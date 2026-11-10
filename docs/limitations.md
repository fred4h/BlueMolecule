# Status and limitations

BlueMolecule is [vibe-coded](vibe-coding.md) and has been tested almost entirely on synthetic files. This page lists,
as honestly as possible, what has been checked, what has not, and what the program does not do yet.

## What has been checked

- **About 8,000 automated checks in 25 test suites**, all passing for this release
  ([how to run them](developer-guide.md#tests)). They cover:
  - every reader and writer, including format detection and the units and conventions of each program;
  - the chemistry routines;
  - the selection language;
  - the launcher, including its local-server security checks;
  - the packaging;
  - the interface itself, in a headless browser.
- **Large files:** synthetic 600 MB Gaussian and ORCA logs, a 1.2 GB log and a 1.1 GB DCD trajectory, with
  `tests/run_all.py --big`.
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
- **The macOS app** has never been compiled or run on a Mac. It was written and checked on Linux.
- **The Linux Qt window** has been tested against a simulated Qt in the test suite, and the app-window fallbacks in
  a headless browser. It still needs confirmation on real desktops.
- **The macromolecule views** were tested on generated protein, membrane, glycan and DNA structures, not yet on
  structures from the PDB or from simulation set-up tools.
- **Windows** is not supported by an installer. The single HTML file should work in a modern browser, but this has
  not been tested.

## Known limitations

| Area | Limitation |
|---|---|
| Orbitals | Orbitals, densities and potentials are **not computed** from fchk or Molden basis sets. Isosurfaces need a cube or XSF grid from your program (`cubegen`, `orca_plot`) |
| Surfaces | No molecular surfaces (solvent-accessible or solvent-excluded) and no electrostatic potential mapped on a surface |
| QM/MM energies | The layer and component energies are read and stored, but only the total energy is plotted |
| Secondary structure | A simplified DSSP: α-helices and β-ladders only (no 3₁₀ or π helices, no turns), computed once per structure rather than per trajectory frame. Records in PDB or mmCIF files are used when present |
| Cartoons | Rebuilt during trajectory playback only for systems up to 20,000 atoms |
| Glycans | Identification by shape works for six-membered (pyranose) rings; unknown five-membered (furanose) rings get white symbols |
| Trajectories | Long binary trajectories are thinned automatically to fit memory; the stride and frame range cannot yet be chosen in the interface. Amber NetCDF (`.nc`) trajectories are not read |
| Readers | No readers for Molpro or Dalton outputs (Molpro inputs can be written) |
| Gaussian inputs | ONIOM layer columns in `.gjf` files are skipped; layers are read from outputs |
| Clean-up | A generic force field for tidying sketches; not a published force field and not a physical method |
| Building | No conformer search, crystal builder or download from PubChem or the PDB |
| Browser-only use | The single HTML file has no Open Babel, saves through downloads and is limited by the browser's memory |
| Language | The interface and documentation are in English only |

## Ideas for future versions

Not promised, and in no particular order:
- orbitals and densities computed from fchk and Molden files;
- molecular surfaces and electrostatic potentials;
- a stride and frame-range control for trajectories;
- Amber NetCDF trajectories;
- Molpro and Dalton outputs;
- conformer search and crystal tools.

Validating the readers against real outputs comes before all of these. [Reports with example files](../CONTRIBUTING.md)
help most.
