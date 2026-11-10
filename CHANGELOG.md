# Changelog

BlueMolecule is [vibe-coded](docs/vibe-coding.md). Each version below was written by an AI assistant from the
author's requests and bug reports.

## 1.4.1 (10 November 2026): first public release

- Full documentation in `docs/`, with a clear notice that the project is vibe-coded. The notice also appears in the
  Help dialog and the package READMEs.
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
