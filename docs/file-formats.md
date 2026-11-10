# File formats

BlueMolecule reads **47 file formats** and saves **24**, plus calculation inputs for **10 programs**. With
[Open Babel](installation.md#open-babel-optional) installed, it can also open and save the many formats Open Babel
supports.

## How a file's format is recognised

The format is detected **from the content first and the extension second**. Many programs write `.out` or `.log`
files, so the first part of the file is checked against the signature of every reader: the program banner,
characteristic headers, block names.
- When the content is unambiguous, it decides.
- Otherwise the extension breaks the tie.
- Gzip-compressed files (`.gz`) are decompressed while reading.
- If no reader matches and Open Babel is available, the file is converted with Open Babel (to SDF, else CML) and
  opened. The Info tab says so.

The readers were written from public descriptions of each format and tested on example files written to imitate
each program's output ([why that matters](vibe-coding.md#how-the-program-was-checked)).

## Quantum-chemistry outputs

| Program | Files | What is read (all of it exercised by the test files) |
|---|---|---|
| **Gaussian** 09/16 | `.log` `.out` | Every geometry (input and standard orientation), optimisation criteria with thresholds, energies, several jobs (`--Link1--`), vibrations with IR and Raman intensities (also from ROA jobs), thermochemistry, excited states (TD-DFT, CIS) with ECD, relaxed scans, IRC, NMR shieldings, Mulliken/Hirshfeld/CM5/NPA/ESP/APT charges, orbital energies, dipole, ONIOM layers and energies |
| **ORCA** 4/5/6 | `.out` `.log` | Geometries, optimisation criteria, energies, several jobs (`$new_job`, compound jobs), vibrations with IR and Raman, thermochemistry, excited states (TD-DFT, RPA) with ECD, scans, IRC (including the IRC path summary), NMR, Mulliken/Löwdin/Hirshfeld charges, orbital energies, dipole, QM/MM, QM/XTB and QM/QM2 layers and energies |
| **GAMESS-US / Firefly** | `.log` `.out` `.gamout` | Geometries, optimisation criteria, energies, vibrations with IR, thermochemistry, charges, orbital energies, dipole |
| **NWChem** | `.out` `.log` `.nwout` | Geometries, optimisation criteria, energies, vibrations with IR, thermochemistry, TD-DFT excited states, charges, orbital energies, dipole |
| **Q-Chem** | `.out` `.log` `.qcout` | Geometries, optimisation criteria, energies, several jobs, vibrations with IR, thermochemistry, TD-DFT excited states, charges, orbital energies, dipole |
| **Psi4** | `.out` `.log` | Geometries, optimisation criteria, energies, vibrations with IR, thermochemistry, charges, orbital energies, dipole |
| **xtb** | `.out` `.log` `.xtbout` | Optimisation cycles, energies, gradient norm, vibrations with IR, thermochemistry, charges, orbital energies, dipole, ONIOM regions |
| **MOPAC** | `.out` `.arc` `.mopout` | Geometries, optimisation cycles, heats of formation, vibrations with IR, charges, orbital energies, dipole |
| **CP2K** | `.out` `.log` | Optimisation cycles and energies, charges, the cell, QM/MM regions and link atoms |
| **Turbomole** | `gradient` | Every optimisation cycle with energy and gradient |
| **VASP** | `OUTCAR` | Ionic steps with energies and forces, the cell |
| **Quantum ESPRESSO** | `.pwo` `.out` | Relaxation steps with energies, the cell |

## Program inputs

Opened for their geometry (and charge, multiplicity, route or keywords where present):

| Program | Files | Notes |
|---|---|---|
| Gaussian | `.gjf` `.com` `.gau` `.gjc` | Cartesian and Z-matrix (including variables) geometries; ONIOM layer columns are skipped (layers are read from outputs) |
| ORCA | `.inp` | `* xyz`, `* int` (internal coordinates) and `* xyzfile` blocks |
| GAMESS | `.inp` | `$DATA` with C1 symmetry |
| NWChem | `.nw` `.nwi` | `geometry` block |
| Q-Chem | `.in` `.qcin` `.inp` | `$molecule` |
| Psi4 | `.dat` `.in` `.psi4` | `molecule { }` block |
| MOPAC | `.mop` | Cartesian and Z-matrix |
| CP2K | `.inp` `.restart` | `&COORD` and `&CELL` |
| Quantum ESPRESSO | `.pwi` `.in` | pw.x input: `ibrav` lattices (all except 13) or `CELL_PARAMETERS`; crystal, alat, bohr or angstrom positions |
| CASTEP | `.cell` | Absolute and fractional positions, lattice |
| Turbomole | `coord` `.tmol` | `$coord` in Bohr |

## Wavefunction and volume files

| Format | Files | What is read |
|---|---|---|
| Gaussian formatted checkpoint | `.fchk` `.fch` `.fck` | Structure, energy, vibrations with IR, charges, orbital energies, dipole |
| Molden | `.molden` `.mold` `.molf` | Structure, optimisation geometries and energies, vibrations with IR and Raman |
| Gaussian cube | `.cube` `.cub` | Volumetric data for isosurfaces, including cubes with several orbitals; Bohr and Å grids |
| XCrySDen XSF / AXSF | `.xsf` `.axsf` | Structure, cell, animation steps and `DATAGRID_3D` volumes |

Isosurfaces need a cube or XSF grid. BlueMolecule does not yet compute orbitals on a grid from fchk or Molden
files ([limitations](limitations.md)).

## Structures and crystals

| Format | Files |
|---|---|
| XYZ and extended XYZ, single or multi-frame | `.xyz` `.txt` `.trj` |
| MDL Molfile / SD file (V2000 and V3000) | `.mol` `.sdf` `.sd` `.mdl` |
| Protein Data Bank, with HELIX/SHEET, chains, B-factors, occupancies, CRYST1 | `.pdb` `.ent` |
| mmCIF / PDBx, with secondary structure | `.cif` `.mmcif` `.mcif` |
| CIF (small-molecule crystals, symmetry expanded) | `.cif` |
| Tripos Mol2, with charges and substructures | `.mol2` `.ml2` `.sy2` |
| Chemical Markup Language | `.cml` |
| Avogadro CJSON | `.cjson` |
| VASP POSCAR / CONTCAR | `.vasp` `POSCAR` `CONTCAR` |
| GROMACS structure (multi-frame allowed) | `.gro` |
| LAMMPS data | `.data` `.lmp` `.lammps` |

## Molecular dynamics: topologies and trajectories

| Format | Files | Kind |
|---|---|---|
| Amber topology | `.prmtop` `.parm7` | Topology: atoms, residues, charges, box |
| Amber coordinates / restart | `.inpcrd` `.rst7` `.restrt` | Coordinates for a topology |
| CHARMM / NAMD PSF | `.psf` | Topology: atoms, segments, residues, charges |
| CHARMM / NAMD DCD | `.dcd` | Binary trajectory |
| GROMACS TRR | `.trr` | Binary trajectory |
| GROMACS XTC | `.xtc` | Compressed binary trajectory |
| LAMMPS dump | `.lammpstrj` `.dump` | Trajectory |
| VASP XDATCAR | `XDATCAR` | Trajectory |

Open the topology or structure **first**, then the coordinates or trajectory with the same number of atoms. Long
binary trajectories are read with a stride so that they fit in memory. See
[Large files and trajectories](large-files-and-trajectories.md).

## Saving

**Export → Save as…** saves the current geometry, or all geometries for the trajectory formats.

| Format | Extension | Notes |
|---|---|---|
| XYZ | `.xyz` | Current geometry |
| XYZ trajectory | `.xyz` | All geometries, with energies in the comment lines |
| Extended XYZ | `.xyz` | Current geometry, with the cell when there is one |
| Extended XYZ trajectory | `.xyz` | All geometries |
| MDL Molfile | `.mol` | With bonds and bond orders |
| SD file | `.sdf` | Current geometry |
| SD file, all geometries | `.sdf` | One record per geometry |
| Tripos Mol2 | `.mol2` | With bond orders and charges |
| Protein Data Bank | `.pdb` | Keeps atom names, residues and the cell |
| CIF (P1) | `.cif` | Periodic structures |
| VASP POSCAR | `.vasp` | Periodic structures |
| Quantum ESPRESSO pw.x input | `.pwi` | Periodic structures |
| CP2K input | `.inp` | `&COORD` / `&CELL` |
| XCrySDen XSF | `.xsf` | With the volume when the document has one |
| GROMACS .gro | `.gro` | With residues and the box |
| LAMMPS data (atomic) | `.data` | Atoms and box |
| Molden | `.molden` | Atoms and vibrations |
| Gaussian cube | `.cube` | The volume (Bohr) |
| Gaussian input | `.gjf` | Current geometry; for full control use Calculation input |
| ORCA input | `.inp` | Current geometry |
| Turbomole coord | `coord` | Bohr |
| Chemical Markup Language | `.cml` | With bonds |
| Avogadro CJSON | `.cjson` | With bonds and the cell |
| SMILES | `.smi` | From the perceived bonds and bond orders |

With Open Babel installed, an **Open Babel** group in the same dialog lists every output format of your Open
Babel version.

**Calculation inputs** for Gaussian, ORCA, Psi4, Q-Chem, NWChem, GAMESS, MOPAC, xtb, Molpro and Turbomole are
written by **Export → Calculation input…** ([details](building-and-editing.md#calculation-inputs)).

**Other exports:**
- **Images** of the view (PNG, up to 4× resolution, optionally transparent);
- **plots** as CSV, PNG or SVG;
- **vibration animations** as multi-frame XYZ.

## Units

Internally, and in every export: lengths in Å (Bohr where a format requires it), energies in Hartree, dipoles in
Debye, frequencies in cm⁻¹. The interface converts energy differences to kcal/mol, kJ/mol or eV on request.
