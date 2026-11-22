# Status and limitations

BlueMolecule is [vibe-coded](vibe-coding.md). It was written and tested with synthetic example files and, from 2.0.1,
also checked against public collections of real output files. This page lists, as honestly as possible, what has
been checked, what has not, and what the program does not do yet.

## What has been checked

- **About 34,500 automated checks in 35 test suites** (23,700 of them cut example outputs as a running job would
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
- **Real output files (2.0.1):** every reader was run on the public test files of the cclib project (437 outputs of
  Gaussian, ORCA, GAMESS-US and Firefly, NWChem, Q-Chem, Psi4, MOPAC and xtb, and Gaussian and Q-Chem fchk files)
  and of the MDAnalysis project (structures, topologies and trajectories), none of them included here. Geometries,
  energies, frequencies and IR intensities, excited states, dipoles, thermochemistry, orbital energies and charges
  were compared with cclib's own reading of the same files; every file was also read cut at 10 places and damaged
  on purpose; and the interface was driven with 23 of them and 31 example files in a headless browser (every
  control, then 1,900 random actions). The differences found are fixed in 2.0.1 (see the [changelog](../CHANGELOG.md));
  `tests/test_real_layouts.js` keeps their layouts.
- **Gaussian frequencies:** real Gaussian 09 and 16 frequency outputs (public test files of the cclib project, not
  included here, with up to 18 blocks of normal modes and the hpmodes layout) read the same in full and condensed.
  `tests/test_freq_sections.js` adds layouts seen in big jobs: values too wide for their column (asterisks), unknown
  lines among the labels, atom numbers beyond the column width.
- **The desktop program:** every change is built on Linux, Windows and macOS by the cloud builds and started once
  there with a self-test: the 3D view (WebGL 2) in each system's web engine, a 900 MB log and a molecule opened
  through the native bridge, a save, and a picture of the view.
- **Script injection:** crafted files with HTML in titles, atom and residue names, labels and file names. None of it
  may become markup in the page (`tests/browser_xss_test.js`).
- **Screenshots** of the interface, produced by the browser tests, were inspected during development.

## Not yet verified

- **Real output files beyond those collections.** Most example files were generated to imitate each program's
  output. Real outputs from other versions, print levels, keywords and unusual job types may contain layouts the
  readers do not handle. In particular, these have not been checked against real files:
  - ORCA QM/MM outputs other than ORCA 6.1 optimisations (one real ORCA 6.1 output was checked for 2.0.2), as
    their wording differs between versions;
  - Gaussian ONIOM microiterations;
  - CP2K QM/MM;
  - xtb ONIOM.
- **Real huge outputs.** Condensed reading was tested on synthetic logs of up to 10 GB and on every example output.
  A real multi-gigabyte QM/MM output has not been opened with it yet.
- **The desktop installers by hand.** The cloud builds start the program on each system, but nobody has yet
  installed the Windows and macOS installers on a personal computer, clicked through their warnings, used **Open
  with** from Explorer or the Finder, or installed an update with the built-in updater.
- **The 1.x macOS app and Windows launcher** (no longer released) were never compiled or run on a Mac or a PC.
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
| Readers | No readers for Molpro or Dalton outputs (Molpro inputs can be written), GAMESS punch files (`.dat`) or Turbomole's program outputs (`dscf.out`, `escf.out`, `ricc2.out`...; its `coord` and `gradient` files are read) |
| Frequencies by finite differences | Psi4 frequency jobs computed from displaced geometries, and Q-Chem anharmonic jobs, show every displaced geometry as a step |
| fchk files | The last geometry only (not the optimisation history) and no excited states |
| Gaussian inputs | ONIOM layer columns in `.gjf` files are skipped; layers are read from outputs |
| Clean-up | A generic force field for tidying sketches; not a published force field and not a physical method |
| Building | No conformer search, crystal builder or download from PubChem or the PDB |
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
