# Troubleshooting

## First step: `bluemolecule --diagnose`

On Linux and macOS, this prints:
- the Pythons and Qt bindings that were found;
- the browsers that were found;
- which window kind BlueMolecule would use, and why;
- where the log file is.

Include its output in bug reports about starting the program.

- Log on Linux: `~/.local/state/bluemolecule/bluemolecule.log`
- Logs on macOS: `~/Library/Logs/BlueMolecule.log` and `~/Library/Logs/BlueMolecule-install.log`
- In the Linux Qt window, **F12** opens the developer tools: the console shows JavaScript errors.

## The window opens as a browser tab or shows the browser's icon

BlueMolecule only gets its own taskbar entry and icon in the **native Qt window** (Linux) or the **native app**
(macOS). It falls back to a browser app window when those are not available.

- **Linux:** install Qt WebEngine for Python, e.g. on Ubuntu or Kubuntu:
  `sudo apt install python3-pyqt6.qtwebengine`. Then check with `bluemolecule --diagnose`. The launcher also checks
  `/usr/bin/python3`, so a conda Python earlier on your `PATH` does not get in the way.
- **Forcing a window kind:** `bluemolecule --tier qt`, `--tier chromium` or `--tier firefox`, or set
  `BLUEMOLECULE_TIER`.
- **A browser tab instead of an app window** happens with `--tab`, or when no Chromium-family browser or Firefox is
  found. Install one, or point to it with `--browser /path/to/browser`.
- **KDE on Wayland, browser app window:** the browser runs under XWayland so that the window matches
  BlueMolecule's desktop entry. If the taskbar still shows the browser's icon, run `./install.sh` again. It
  refreshes KDE's desktop-file cache. Then log out and in.

## "Open With BlueMolecule" is missing for a file type

The installer registers `.xyz .cube .fchk .molden .mol .sdf .mol2 .pdb .cif .gro .gjf .cml .cjson .xsf`. Outputs
such as `.log` and `.out` are plain text to the desktop. Open them from inside BlueMolecule, by drag and drop, or
from the terminal with `bluemolecule file.log`. You can also choose BlueMolecule once in your file manager's
**Open With → Other application** dialog.

## A file does not open, or opens wrongly

- **"The format of … was not recognised"**: the content did not match any reader. Check that the file is complete
  and is one of the [supported formats](file-formats.md). With Open Babel installed, unknown formats are converted
  automatically.
- **A trajectory "has coordinates but no atom types"**: open the topology or structure (prmtop, psf, pdb, gro)
  first, then the trajectory ([details](large-files-and-trajectories.md#molecular-dynamics-trajectories)).
- **Missing data**, such as no vibrations, no excited states or wrong layers: look at **Reader notes** in the Info
  tab. Then please [report it](../CONTRIBUTING.md#reporting-a-file-that-is-read-wrongly) with a trimmed example
  file. This is the most valuable kind of bug report.
- **"Not enough memory"**: see [Large files](large-files-and-trajectories.md). Close other documents, or use the
  Qt window or a Chromium-based browser.

## The QM/MM energy graph looks wrong

BlueMolecule plots one total energy per optimisation step ([details](qmmm-oniom.md#energies-one-point-per-optimisation-step)).
If the curve still jumps between unrelated values, the output contains a layout the reader does not know. Please
report it with the output.

## The 3D view is blank or slow

BlueMolecule needs **WebGL 2**.
- On Linux, check that your graphics driver is installed. On old or virtual machines, the browser may fall back to
  slow software rendering.
- For very large systems, use the Lines or Hidden styles for the environment
  ([representations](representations.md)).

## Saving does not show a dialog

In the single HTML file, saving goes through the browser's downloads. Check your downloads folder. The installed
Linux and macOS versions show the system's Save dialog.

## Open Babel is not found

Install it ([how](installation.md#open-babel-optional)), make sure `obabel` runs in a terminal, and restart
BlueMolecule. The launcher looks on your `PATH` and in the usual Homebrew and conda folders. The Open Babel group
in Save as… shows "checking…" until `obabel` answers.
