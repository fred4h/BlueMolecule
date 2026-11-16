# Large files and trajectories

## Very large outputs

Long optimisations, QM/MM runs on whole proteins and long MD outputs easily produce log files of hundreds of
megabytes. Older versions of BlueMolecule showed such files as **empty**. The JavaScript engine behind the viewer
cannot hold more than about 512 MB of text in a single string, so the file could not be read in one piece.

Since version 1.3, files of **64 MB and more are read in pieces**:
- The file is read in 32 MB slices and decoded into lines as it goes. Characters and line breaks that fall across
  two slices are handled correctly.
- Every reader works on the list of lines rather than on one huge text.
- A progress bar shows the megabytes read so far, with a **Cancel** button.
- Before reading, BlueMolecule estimates the memory the file will need, about 1.7 times its size. If that clearly
  does not fit, it stops with a message instead of crashing the window. In the Chromium-based windows it also
  watches the memory while reading, and stops early with an explanation if memory runs out.
- In the Linux Qt window and the macOS app, large files are handed to the viewer piece by piece. The file is never
  copied into the window in one go.

On the test machine, synthetic 600 MB Gaussian and ORCA logs open in a few seconds, and a 1.2 GB log in about 6
seconds. Your times depend on your computer and on how much of the file is geometry.

A few formats are read as one whole text and are limited to about 300 million characters: JSON-based CJSON and
XML-based CML. This limit hardly matters in practice.

**Compressed files** (`.gz`) are decompressed while reading, so a 600 MB log compressed to 100 MB opens without
being unpacked first. Because the uncompressed size is not known in advance, the memory check happens during
reading.

## Outputs of several gigabytes

A QM/MM optimisation and frequency job on a whole protein prints the geometry of every step and the displacements
of tens of thousands of normal modes. As text, that makes a log of many gigabytes, more than the JavaScript engine
behind the viewer can hold: about 4 GB for all of the page's data in the Chromium-based windows, however much memory
the computer has. Before version 1.5 such files were refused for lack of memory.

Since version 1.5, **Gaussian and ORCA outputs and XYZ trajectories of 64 MB and more are condensed while they are
read**. Geometries and normal-mode displacements go straight into compact number arrays, and the many lines that no
reader uses are dropped as they go by. What you see is the same as when the whole text is read: the test suite
checks this on every example output, and on hundreds of cut versions of them.

On the test machine, a synthetic 5 GB Gaussian ONIOM optimisation and frequency log (8,144 atoms, 20 steps and
24,426 vibrations) opens in about a minute and a half with under 1 GB of memory, and a 10 GB one (11,531 atoms) in
about three minutes with about 1.5 GB.

### Choosing what to load

For an output or trajectory of **1 GB or more**, BlueMolecule first asks what to load. The dialog shows what the
start and the end of the file say (the atoms, the jobs and the number of steps or frames) and offers:

| Choice | What is loaded |
|---|---|
| All steps (frames) | Every geometry; the default when it fits in memory |
| Every Nth step (frame) | Steps 1, N + 1, 2N + 1 and so on |
| Only the last N steps (frames) | The end of the run |
| Steps (frames) from … to … | A range, counted from 1 |
| Vibrations: all, with their displacements | Every normal mode, ready to animate |
| Displacements of the first K only | The first K modes can be animated; the others keep their frequencies and intensities |
| Frequencies and intensities only | The spectra and the Vibrations list, without animation |

The dialog proposes a choice that fits in about half of the computer's free memory. As you change it, it estimates
the memory of what will be kept.

For outputs, the **energies and convergence data of every step are read whatever you choose**, so the convergence
plots stay complete. A point whose geometry was not loaded says *(geometry not loaded)* when you hover over it. The
first and last geometry of every job are always loaded. Vibrations without displacements are listed dimmed.

The **Info** tab then has a **Loaded part of the file** section. It shows how many geometries or frames were loaded
and how, how many vibrations have displacements, and the memory they use. **Load other parts…** opens the dialog
again with the previous choice and reads the file again with the new one.

The dialog does not appear for compressed (`.gz`) files, which cannot be read in pieces.

### When memory runs out anyway

The message says which limit was reached, and what to do:
- **The JavaScript engine's limit** (about 4 GB in Chromium-based windows, whatever the computer has): load fewer
  steps or vibrations, or close other documents.
- **The computer's memory**: close other applications or documents, or load less.

## Molecular-dynamics trajectories

MD data usually comes in two files: a **topology or structure** that names the atoms, and **coordinates or a
trajectory** that only has positions.

1. Open the topology or structure: Amber `prmtop`, CHARMM/NAMD `psf`, `pdb` or `gro`.
2. Open the coordinates or trajectory: `inpcrd`/`rst7`, Amber NetCDF (`nc`, `ncrst`), `dcd`, `trr`, `xtc`,
   LAMMPS dump. It must have the same number of atoms. BlueMolecule attaches it to the matching structure and plays
   it with the playbar.

A topology opened on its own shows a **Coordinates pending** card with an **Open coordinates…** button. A
trajectory opened without a matching topology is refused, with a message naming what to open first.

**Amber NetCDF** is Amber's default output format since Amber 16. BlueMolecule reads it without any extra library:
- trajectories (`.nc`, `.ncdf`, or `.mdcrd` files written with `ioutfm = 1`);
- restarts (`.ncrst`, or `.rst7` and `.restrt` files that hold NetCDF; text `.rst7` files are still read as text).

The box of every frame (truncated octahedra included) and the time are read; velocities and forces are skipped.
NetCDF 4 files, which have HDF5 inside and which some tools write on request, are recognised and refused with a
command to convert them: `nccopy -k 64-bit-offset in.nc out.nc`. Text (ASCII) `mdcrd` trajectories are not read.

**Long binary trajectories** (DCD, TRR, XTC, NetCDF) are read in slices:
- From **1 GB**, the loading dialog above offers every Nth frame, a frame range or the last N frames.
- Below that, when all the frames would take more than about 512 MB of memory, BlueMolecule keeps **every Nth
  frame** automatically. It says so in a message and under **Reader notes**, e.g. *Loaded 2,000 of 20,000 frames
  (every 10th)*.

Either way, the Info tab shows the frames loaded, e.g. *Frames of run.nc: 2,000 of 20,000 (every 10th)*. **Load
other parts…** reads the trajectory again with another choice.

Multi-frame XYZ files of 1 GB or more get the same dialog. LAMMPS dumps, XDATCAR and multi-frame GRO files are read
completely.

## Very many atoms

Systems of hundreds of thousands of atoms can be displayed. A 290,000-atom synthetic membrane system opens in about
a second in the tests. A few things keep large systems fast and readable:
- The automatic presets hide water and draw ions as small spheres, lipids as lines and proteins as cartoons
  ([representations](representations.md)).
- **Lines** and **Hidden** are the cheapest styles; **Ghost** gives context without clutter.
- Some building tools are switched off above certain sizes, such as bond-order perception above 2,000 atoms and
  clean-up above 3,000 ([limits](building-and-editing.md#size-limits)).
- Cartoons are rebuilt during trajectory playback only up to 20,000 atoms.

## Tips

- **Working on a cluster?** BlueMolecule runs on your own computer. Copy the output over, or compress it there
  first with `gzip` and open the `.gz` directly.
- **Strip solvent** you do not need before saving a trajectory for viewing. That is often the single biggest saving.
- If a file will not open for lack of memory, close other documents and applications, or open the file in the
  Linux Qt window or a Chromium-based browser. Firefox and Safari are not memory-checked in advance.
