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

On the test machine, synthetic 600 MB Gaussian and ORCA logs open in about 3 seconds, and a 1.2 GB log in about 6
seconds. Your times depend on your computer and on how much of the file is geometry.

A few formats are read as one whole text and are limited to about 300 million characters: JSON-based CJSON and
XML-based CML. This limit hardly matters in practice.

**Compressed files** (`.gz`) are decompressed while reading, so a 600 MB log compressed to 100 MB opens without
being unpacked first. Because the uncompressed size is not known in advance, the memory check happens during
reading.

## Molecular-dynamics trajectories

MD data usually comes in two files: a **topology or structure** that names the atoms, and **coordinates or a
trajectory** that only has positions.

1. Open the topology or structure: Amber `prmtop`, CHARMM/NAMD `psf`, `pdb` or `gro`.
2. Open the coordinates or trajectory: `inpcrd`/`rst7`, `dcd`, `trr`, `xtc`, LAMMPS dump. It must have the same
   number of atoms. BlueMolecule attaches it to the matching structure and plays it with the playbar.

A topology opened on its own shows a **Coordinates pending** card with an **Open coordinates…** button. A
trajectory opened without a matching topology is refused, with a message naming what to open first.

**Long binary trajectories** (DCD, TRR, XTC) are read in slices too. When all the frames' coordinates would take
more than about 512 MB of memory, BlueMolecule automatically keeps **every Nth frame**. It says so in a message and
under **Reader notes**, e.g. *Loaded 2,000 of 20,000 frames (every 10th)*. The interface does not yet let you
choose the stride or a frame range yourself ([limitations](limitations.md)). Text trajectories (multi-frame XYZ,
LAMMPS dump, XDATCAR, multi-frame GRO) are read completely.

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
