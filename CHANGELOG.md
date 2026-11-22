# Changelog

BlueMolecule is [vibe-coded](docs/vibe-coding.md). Each version below was written by an AI assistant from the
author's requests and bug reports.

## 2.2.0 (November 2026): images and scenes

**Saving images** (issue #16):
- **Formats**: PNG, JPEG, WebP, TIFF and PDF. BlueMolecule writes every format itself, so an image is the same
  file in every browser and in the desktop program; WebKit, the engine of the macOS and Linux desktop programs,
  cannot write WebP otherwise.
- **Size**: a width in pixels, centimetres or inches at a chosen resolution (DPI), written into the file. PDF makes
  one page of that size.
- **Options**: transparent background (PNG, WebP, TIFF, PDF), JPEG quality, compression.
- **Colour profile**: sRGB (recommended) or Display P3, embedded in the file.
- **Large images**: drawn in tiles of 2,048 pixels and put together, up to 250 million pixels (posters larger than
  the graphics card allows in one piece).

**Scenes and effects** (issue #15), in the display settings:
- **Presets**: Default, Publication (white, matte, ambient occlusion, thin outlines), Studio (gradient, glossy,
  shadows on a floor), Dark and Cartoon (flat shading, thick outlines).
- **Each effect on its own**: material (standard, glossy, matte, flat), ambient occlusion, outlines, depth of field,
  shadows, floor and background gradient.
- **Large structures**: above 60,000 atoms, shadows, ambient occlusion and depth of field wait until they are asked
  for.
- Saved images keep the effects; tiles are drawn with a margin, so they match across their edges.

## 2.1.0 (November 2026): transition states and frequencies by finite differences

**Transition states** (issue #14): the Vibrations tab now helps to check a transition state.
- A line at the top says what the frequencies make of the structure: a minimum, a transition state (exactly one
  imaginary frequency) or a higher-order saddle point. Imaginary frequencies below 50i cm⁻¹, often numerical
  noise, are counted apart.
- Opening a file with an imaginary mode, or a transition-state search, selects its largest imaginary mode, and the
  message on opening gives the number of imaginary frequencies.
- **Arrow length** sets the length of the displacement arrows, apart from the amplitude (the largest atom shift of
  the animation).
- **Both ends** shows the two ends of the motion as translucent structures (blue −, orange +) on a still view; a
  bond that forms is drawn at one end only.
- **Bonds along mode N** lists the pairs whose distance changes most, with the distances at both ends: along an
  imaginary mode, which bonds form and which break. A row measures that distance.
- **Forward and reverse inputs…** writes two optimisation inputs (*name*_minus and *name*_plus) from the geometry
  displaced both ways, in any program of the calculation-input dialog. A Gaussian route or ORCA keyword line taken
  from the file becomes an optimisation (`opt=(ts,calcfc)` becomes `opt`, `OptTS` becomes `Opt`).
- **Negative Hessian eigenvalues**: a new optimisation plot, and a line in the Steps tab, count them at each step of
  Gaussian and Q-Chem optimisations. A good transition-state search stays at 1.

**Frequencies computed by finite differences** (issue #32) showed every displaced geometry as a step of the job:
74 for divinylbenzene in Psi4 (`frequencies(..., dertype='gradient')`), 50 for water in Q-Chem anharmonic (VPT2) jobs.
The energy plots treated them as optimisation steps, and Psi4's normal modes and dipole belonged to the last
displaced geometry. Only the reference geometry is kept now, with its energy, dipole and vibrations.

Also fixed:
- The calculation-input dialog wrote the **keywords of an ORCA output separated by commas** (`! B3LYP,def2-SVP,Opt`),
  which ORCA does not read; they are separated by spaces again.
- With **Use the route from this file**, a Gaussian route with `geom=check` or `guess=read` made Gaussian read the old
  geometry from a checkpoint file instead of the exported one; these keywords are now left out.
- The displacement of **Displace along the selected mode** went back to 0.15 Å whenever another mode was selected.

## 2.0.2 (November 2026): ORCA 6 QM/MM energies

Checked against a real ORCA 6.1 QM/MM optimisation (the synthetic example imitated another layout):

- **The energies of ORCA 6 QM/MM optimisations were wrong.** Each step prints three final energies: the MM energy
  of the full system ("FINAL SINGLE POINT ENERGY (MM)"), the QM region's, and the QM/MM total. All three were read
  as one energy, and the MM line comes before the step's geometry, so every step showed the **next step's MM energy**
  (about +2 Eh) and only the last step the QM/MM total: ΔE versus the final step was millions of kcal/mol and the
  energy plots were flat. Each step now has its QM/MM total as its energy, with its own MM and QM energies as
  components; a job still running leaves out a step whose total is not printed yet.
- **Outputs that print only the QM region** (the usual case with `*pdbfile`) showed all of its atoms as MM. The QM
  region is now shown as QM, with its **link atoms** marked; ORCA 6's QM/MM summary (sizes of the QM/MM system, the
  QM subsystem and the active region) is read.

## 2.0.1 (November 2026): checked against real output files

Every reader was run on the public test files of the cclib project (437 outputs of nine programs) and of the
MDAnalysis project, and compared with cclib's reading of the same files; each file was also read cut at ten places
and damaged on purpose, and written in every format and read back; and the interface was driven with real files in a
headless browser. The test files of earlier versions imitated each program's output; the real files showed these
differences, now fixed:

- **Correlated energies (Gaussian):** MP4 energies were never read, because Gaussian prints them in the middle of a
  line (`E4(SDTQ)= … UMP4(SDTQ)= …`), so an MP4 or MP5 job showed its MP3 energy. MP4 (SDTQ, SDQ, DQ), MP5 and CCD
  are now read, and Gaussian 09's CCSD energy, which ends with "CI/CC converged" instead of "Wavefunction amplitudes
  converged".
- **ORCA TD-DFT, CIS, ADC(2) and EOM:** ORCA's "FINAL SINGLE POINT ENERGY" is then the energy of an excited state, and
  it was shown as the energy of the structure; single points now show the ground state, as for Gaussian (MP2 for
  ADC(2)). With triplets, ORCA 5 listed them twice (once as "singlets" without an oscillator strength), and in
  ORCA 6 the **ECD rotatory strengths of the singlets were attached to the triplets**, so the ECD spectrum was
  wrong. The ground-state CCSD energy of EOM-CCSD runs is read.
- **Jobs still running:** a geometry cut while the program was writing it became a frame with only some of the atoms
  in GAMESS, NWChem, Psi4, Q-Chem, xtb, MOPAC and CP2K outputs (Gaussian and ORCA were right). A charge table cut the
  same way is left out until it is whole.
- **Q-Chem orbital energies:** when Q-Chem prints the orbitals' symmetries (most jobs on molecules with symmetry), only
  the first 8 orbital energies were read, so the HOMO, LUMO and gap were missing and the **MO energy-level plot
  stopped with an error** (36 of the 54 public Q-Chem outputs). All are read now, and orbital energies that stop
  before the HOMO are left out with a note instead of breaking the plot.
- **Q-Chem inputs of several jobs** (`@@@`) merged the molecules of all jobs into one; the first job's molecule is
  now shown, with every job's type.
- **ORCA inputs** whose coordinates use `%paras` values (`{R1}`) could not be opened; they are filled in (a scanned
  value takes its first value).
- **More files open:** MOPAC 2012 archives (`FINAL GEOMETRY OBTAINED … CHARGE`), `.mop` inputs without a method
  keyword, Firefly outputs ("FINAL ENERGY IS": no energies before), GAMESS-US outputs whose population table is
  printed twice (no charges before), Psi4, Q-Chem and NWChem inputs of big molecules, and ORCA and xtb `.engrad`
  files (a new reader). A Turbomole `mpshift.out` is no longer mistaken for a PDB file, nor an `.engrad` file for a
  Gaussian input.
- **fchk files:** a TD or CIS job showed the excited state's energy as the structure's; Q-Chem's fchk files get
  their dipole and program name.
- **NWChem** Mulliken charges are read to 5 decimals when the output has them (2 before), and ab initio MD runs
  (`task dft qmd`) show every step with its energy and time (only the first geometry was shown).
- **CP2K inputs** show their charge, multiplicity, functional and job type. Two readers handled CP2K inputs, and
  the one chosen for `.inp` files skipped these; a `.restart` file's cell given in bohr or nm was read as Å.
- **MOPAC** IR intensities are squared transition dipoles (D²), but the table and the spectrum called them km/mol
  and offered ε; they now say D², and the IR spectrum of a MOPAC job is relative. Open-shell (UHF) MOPAC jobs show
  their alpha and beta orbital energies (none were read).
- **SD files** whose records hold the same atoms (conformers, or BlueMolecule's own "SDF (all geometries)") open
  with all their geometries; only the first was shown. Molfiles of more than 999 atoms are written as V3000, as SD
  files were.
- **Calculation inputs** for Psi4, Q-Chem, NWChem, GAMESS, MOPAC, xtb and Turbomole ran two coordinates together
  when one was 1000 Å or more.
- A damaged element number shows as a dummy atom with a note, instead of breaking the view; messages about
  unreadable coordinates name the file and the reader; GROMACS files keep atom names that tell the element; an xtb
  output without coordinates names its input structure.
- The **Info** tab names the level of a geometry's energy when the job computed several (SCF, MP2, CCSD(T)…) and
  lists the others below it, such as a TD-DFT job's excited state, so the energy used can be checked.
- The update dialog shows each release's notes (the release builds now pass them to the update file), as plain text
  without Markdown marks.
- Desktop program: a link to a web server on the computer itself (`http://localhost:…`) opens in the web browser,
  like every other link, instead of replacing BlueMolecule in its window.
- Linux: the MIME types also cover `.cub`, `.ent`, `.mmcif` and `.axsf`, like the program's Open-with entries.

## 2.0.0 (November 2026): a real program with updates

- **BlueMolecule is a desktop program of its own** ([#11](https://github.com/fred4h/BlueMolecule/issues/11)), built
  with Tauri 2: a native window with its own name and icon on Linux, Windows and macOS, using each system's web
  engine. Installers: `.deb` and AppImage for Linux, a per-user `.exe` for Windows (no administrator rights), a
  `.dmg` with a universal app for macOS. Chemistry files get **Open with → BlueMolecule** (never as the default
  program for `.log` or `.out`), files open from the command line and from a second launch into the running
  window, and big files are read on demand through the native bridge, as in the 1.x Qt window. Open Babel, the
  system Open and Save dialogs and the clipboard work on all three systems; web links open in the default browser.
- **Built-in updater** ([#12](https://github.com/fred4h/BlueMolecule/issues/12)). **Help → Check for updates…**
  shows the release notes of a newer version and, if you agree, downloads it, checks its signature against
  BlueMolecule's own key, installs it and restarts. Checking is BlueMolecule's only use of the network: it happens
  when you ask, or at start-up only if you tick that box. `.deb` installations are pointed to the new `.deb`.
- **Cloud builds** ([#13](https://github.com/fred4h/BlueMolecule/issues/13)). Every change is built on Linux,
  Windows and macOS and started once there with a self-test: the 3D view in each engine, a 900 MB log and a
  molecule through the native bridge, a save and a picture of the view. Publishing a release builds, signs and
  attaches the installers and the update files. The installers are not code-signed (macOS: signed ad hoc), so
  Windows and macOS warn once before the first start.
- **Big files in WebKit** (Safari, and the desktop program on macOS and Linux): reading a big output kept the whole
  file in memory as well, because WebKit's JavaScript engine kept every short line that BlueMolecule kept as a
  view of the 32 MB piece of text around it. Lines are now copied: opening the 900 MB log of the cloud self-test
  takes at most 0.8 GB of memory on macOS and 1.1 GB on Linux (before: 1.7 GB and rising), as against 0.6 GB on
  Windows.
- The 1.x packages (the Python launcher and Qt window, and the experimental Windows and macOS packages that
  compiled a launcher during installation) are no longer released; `BlueMolecule.html` still is.
- New tests: `tests/browser_tauri_bridge_test.js` (the page side of the desktop program, with a stand-in for
  Tauri, and the self-test script).

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
- Fixed: a Gaussian frequency job showed only its first three vibrations, which could not be animated, when a value
  of the first block was too wide for its column. Gaussian then prints asterisks right after the dashes
  (`IR Inten    --************`), as for the huge IR intensities of charged QM/MM systems. Such lines, unknown lines
  among a block's labels and atom numbers printed as asterisks no longer end the section. A section that still ends
  early says so in the Reader notes, and a vibration without displacements in the file says so when clicked.
- New tests: condensed reading against the full reading of every example output (`tests/test_condense.js`),
  NetCDF (`tests/test_netcdf.js`), Gaussian frequency sections of many blocks (`tests/test_freq_sections.js`), the
  1.5 interface in a browser (`tests/browser_15_test.js`), and synthetic 5 GB and 10 GB logs in
  `tests/run_all.py --big`.

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
