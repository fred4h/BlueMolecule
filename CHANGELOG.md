# Changelog

BlueMolecule is [vibe-coded](docs/vibe-coding.md). Each version below was written by an AI assistant from the
author's requests and bug reports.

## 1.5.0 (16 November 2026): big files, trajectories and plots

- **Outputs of several gigabytes open** ([#6](https://github.com/fred4h/BlueMolecule/issues/6)). A QM/MM
  optimisation and frequency log of 5 GB used to be refused for lack of memory: the JavaScript engine behind the
  viewer holds about 4 GB at most. Gaussian and ORCA outputs and XYZ trajectories of 64 MB and more are now
  condensed while they are read. Only the lines the readers use are kept, and geometries and normal-mode
  displacements go straight into compact number arrays. A synthetic 5 GB ONIOM log (8,144 atoms, 24,426
  vibrations) opens in about a minute and a half with under 1 GB of memory, and a 10 GB one in about three minutes.
  Memory messages now say whether the engine or the computer ran out, and what to do.
- **Choose what to load** ([#7](https://github.com/fred4h/BlueMolecule/issues/7)). For outputs and trajectories of
  1 GB or more, a dialog offers every Nth step or frame, a range or the last N, and for vibrations all of them, the
  displacements of the first K, or frequencies and intensities only. It estimates the memory as you choose. For
  outputs, the energies and convergence data of every step are still read, and the first and last geometry of each
  job are always loaded. The same choice works for DCD, TRR, XTC, NetCDF and multi-frame XYZ files. The Info tab
  says what was loaded, and **Load other parts…** reads the file again with another choice.
- **Amber NetCDF trajectories and restarts** ([#8](https://github.com/fred4h/BlueMolecule/issues/8)), Amber's
  default formats since Amber 16. They are read without any library, with the box of every frame and the time:
  `.nc`, `.ncdf`, `.ncrst`, and `.mdcrd`, `.rst7` or `.restrt` files that hold NetCDF. Long trajectories are read in
  slices like DCD. The tests use files written by SciPy, and real pmemd, sander and cpptraj trajectories gave the
  same numbers as SciPy. NetCDF 4 files are refused with a command to convert them.
- **Zoomable plots** ([#9](https://github.com/fred4h/BlueMolecule/issues/9)). Drag a box to zoom into any plot (a
  flat box zooms x only); the magnifier buttons zoom in, zoom out and fit, and a double-click shows the whole plot.
  A click still selects. Each plot keeps its zoom while it is redrawn. CSV, PNG and SVG exports follow the zoom:
  the CSV keeps the rows in view and says so.
- **QM/MM and ONIOM component energies are plotted** ([#10](https://github.com/fred4h/BlueMolecule/issues/10)):
  the total and each component (each level on the model and real systems, or QM, MM and QM–MM interaction), each
  relative to its own final value, with a legend and a CSV in hartree. The Info tab lists the components of the
  geometry shown.
- Fixed: the atom labels under the charge plot of a system with more than about 900 atoms overlapped.
- New tests: condensed reading against the full reading of every example output (`tests/test_condense.js`),
  NetCDF (`tests/test_netcdf.js`), the 1.5 interface in a browser (`tests/browser_15_test.js`), and synthetic 5 GB
  and 10 GB logs in `tests/run_all.py --big`.

## 1.4.2 (16 November 2026): quick fixes

- **Jobs that are still running** ([#1](https://github.com/fred4h/BlueMolecule/issues/1)). A last energy that is
  not printed yet, or cut in the middle of its number, no longer becomes 0 or a truncated number, so ΔE and
  the energy and convergence plots stay correct. The incomplete last line of an output is left out, the step
  without an energy is marked *in progress*, and the Info tab says *Not finished: still running?*. A new test cuts
  every example output at thousands of places.
- **Drop-down menus stay open after a click** ([#2](https://github.com/fred4h/BlueMolecule/issues/2)), including
  the plot picker. BlueMolecule now draws its own menus: click to open, click an option; press, drag and release
  also works, and so do the keys (Enter, Space, ↑ ↓, Home, End, letters, Esc). After a choice with the mouse the
  arrow keys and shortcuts work again at once.
- **The Job name of the calculation-input dialog is the file name**
  ([#3](https://github.com/fred4h/BlueMolecule/issues/3)), made safe for file names and with the program's extension
  (for example `My job 1` → `My_job_1.inp`). It starts as the name of the opened file and is also the Gaussian
  checkpoint name.
- **↑ / ↓ move through the Steps, Vibrations and States lists**
  ([#4](https://github.com/fred4h/BlueMolecule/issues/4)), as ← / → do; Page Up / Page Down jump 10 rows and Home /
  End go to the first and last row, which stays in view.
- **Right-click deletes an atom while building** ([#5](https://github.com/fred4h/BlueMolecule/issues/5)), in the
  Draw, Template, Bond-centric and Manipulate tools, and a bond in the Draw tool. Each deletion can be undone.
  Right-drag still pans, and Navigate and Select never delete. Deleting every atom of a file geometry is refused,
  as with the Delete key.
- New browser test of these interface details (`tests/browser_ui_test.js`).

## 1.4.1 (10 November 2026): first public release

- Full documentation in `docs/`, with a clear notice that the project is vibe-coded. The notice also appears in the
  Help dialog and the package READMEs.
- **An experimental Windows package.** A double-click installer (no administrator rights) adds a Start-menu entry
  and "Open with" for chemistry files. BlueMolecule opens in its own Edge window, through a small launcher compiled
  on the PC from C# source. It has not yet been tested on a real Windows PC.
- MIT License; the packages now include the license.
- New test that crafted files cannot inject markup or scripts into the page (`tests/browser_xss_test.js`).
- One command runs every test suite: `python3 tests/run_all.py [--browser] [--big]`.
- The build needs only Python 3: the app icons and the built-in examples are kept in the repository.
- The launcher's `--help` lists the supported programs.

## 1.4.0 (10 November 2026): representations and macromolecules

- **Display tab** with display rows (selection, style, colour, opacity, size). Later rows override earlier ones, and
  overlays add on top.
- **Selection language** similar to VMD's, with live atom counts, suggestions and error messages: classes, QM/MM
  layers, chains, residues, names, elements, `within`, `same residue as`, boolean operators.
- **Presets**, chosen automatically when a file opens: QM/MM layers (with a legend of the methods), Protein–ligand,
  Membrane, Glycan / proteoglycan, Nucleic acid, All atoms.
- **Macromolecule styles:**
  - protein cartoon with secondary structure (from the file, or a simplified DSSP);
  - backbone trace;
  - nucleic-acid base slabs;
  - the lipid style;
  - filled sugar rings and 3D SNFG glycan symbols, with sugars identified by name or by ring stereochemistry.
- Residue classification for protein, nucleic acid, water, ions, lipids, carbohydrates and ligands, shown in the Info
  tab and on hover.

## 1.3.0 (9 November 2026): large files, QM/MM and ONIOM, own window on Linux

- **Files of any size.** Outputs larger than the ~512 MB string limit of the JavaScript engine used to open as
  empty. They are now read in pieces, with progress, cancelling and a memory guard. Long binary trajectories are
  thinned to fit in memory.
- **QM/MM and ONIOM.** One geometry per optimisation step with the total energy only, so the energy graph no longer
  jumps between sub-calculations. Layers, link atoms, frozen atoms, MM types and charges are read for Gaussian
  ONIOM, ORCA QM/MM and QM/XTB, CP2K QM/MM and xtb ONIOM.
- **Its own window on Linux:**
  - a native Qt WebEngine window, with Chromium-family or Firefox app windows as fallbacks;
  - `bluemolecule --diagnose`;
  - KDE and Wayland window matching.

## 1.2.0 (8 November 2026): file compatibility, contents overview, convergence plots

- A code audit for bugs, with regression tests.
- The app's own icon in the dock or taskbar: a native macOS app, plus app windows on Linux.
- **Contents overview** in the Info tab, with a selectable item for everything found in a file, and job selection for
  multi-job files.
- **Optimisation plots:** energy, ΔE, each convergence criterion against its threshold, RMSD to the final geometry,
  largest atomic step.
- **Many more readers:** GAMESS, NWChem, Q-Chem, Psi4, xtb, MOPAC, CP2K, Turbomole, VASP, Quantum ESPRESSO, CASTEP,
  CIF and mmCIF, Mol2, CML, CJSON, XSF, LAMMPS, Amber, CHARMM/NAMD, GROMACS trajectories. Also scans, IRC, NMR and
  ECD.
- **Save as** in 24 formats. Open Babel is used for formats BlueMolecule does not handle itself.

## 1.1.0 (5 November 2026): building and editing

- **Avogadro-style build tools:** navigate, select, draw, template, bond-centric, manipulate.
- SMILES input and output, and a library of molecules, groups, ligands and coordination geometries.
- Hydrogens, bond-order perception, a clean-up force field, alignment, point-group detection and symmetrization,
  and a coordinate editor.
- Calculation inputs for 10 programs.

## 1.0.0 (4 November 2026): first version

- A quantum-chemistry viewer for Gaussian and ORCA outputs, XYZ, cube, Molden, fchk, MOL/SDF, PDB and GRO files.
- Optimisation steps, vibrations with animation, IR, Raman and UV-Vis spectra, and cube isosurfaces.
- Labels and measurements, image export, and Gaussian and ORCA inputs.
- Launchers for macOS and Linux.
