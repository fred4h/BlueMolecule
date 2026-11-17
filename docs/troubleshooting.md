# Troubleshooting

## First step: the log and the version

- **Version:** Help (the **?** button) shows it in its title.
- **Log** of the desktop program:
  - Linux: `~/.local/share/org.bluemolecule.BlueMolecule/logs/bluemolecule.log`
  - Windows: `%LOCALAPPDATA%\org.bluemolecule.BlueMolecule\logs\bluemolecule.log`
  - macOS: `~/Library/Logs/org.bluemolecule.BlueMolecule/bluemolecule.log`
- Starting it from a terminal shows its messages: `bluemolecule calc.log` on Linux,
  `/Applications/BlueMolecule.app/Contents/MacOS/bluemolecule` on macOS.

Include the version, your system and what the log says in bug reports.

## The window does not start, or stays white (Linux)

The desktop program uses WebKitGTK. With some graphics drivers, especially NVIDIA on Wayland, its window stays
white or the program closes. Start it with one of these, and tell us which one helped:

```bash
WEBKIT_DISABLE_DMABUF_RENDERER=1 bluemolecule
WEBKIT_DISABLE_COMPOSITING_MODE=1 bluemolecule
```

If one works, you can put it before `bluemolecule` in the `Exec=` line of a copy of
`/usr/share/applications/BlueMolecule.desktop` in `~/.local/share/applications/`.

## Two BlueMolecule entries, or the old window opens (Linux)

A BlueMolecule 1.x installation is still there: its command in `~/.local/bin` comes before the new one. Remove it
with `./uninstall.sh` from the unpacked 1.x folder ([details](installation.md#coming-from-bluemolecule-1x)).

## Windows and macOS warnings

- **"Windows protected your PC"** when you start the installer: click **More info → Run anyway**. The installer is
  not code-signed.
- **macOS says it cannot check BlueMolecule for malicious software:** open it once, then click **Open Anyway** in
  **System Settings → Privacy & Security** (on older macOS: right-click the app, **Open**). See
  [installation](installation.md#macos).
- **macOS says BlueMolecule "is damaged":** this is the download quarantine. Run
  `xattr -dr com.apple.quarantine /Applications/BlueMolecule.app` in Terminal.

## Updates

- **"This is a test build of BlueMolecule, which does not look for updates":** only the installers from the
  Releases page update themselves; builds made elsewhere (for example the test builds of the cloud) do not.
- **The check did not succeed:** the computer may be offline, or GitHub unreachable. Try again later, or download
  the new version from the Releases page.

## "Open With BlueMolecule" is missing for a file type

On Linux and macOS, the installer registers `.xyz .cube .fchk .molden .mol .sdf .mol2 .pdb .cif .gro .gjf .cml
.cjson .xsf`; on Windows, also `.log`, `.out`, `.inp` and the MD formats. On Linux and macOS, outputs such as `.log`
and `.out` are plain text to the desktop. Open them from inside BlueMolecule, by drag and drop, or
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
- **"Not enough memory"**: see [Large files](large-files-and-trajectories.md#when-memory-runs-out-anyway). Open
  the file again and load less of it (every Nth step, the last steps, vibrations without displacements), close
  other documents, or use the Qt window or a Chromium-based browser.
- **A NetCDF file is "a NetCDF 4 file"**: convert it to NetCDF 3 with `nccopy -k 64-bit-offset in.nc out.nc`
  ([details](large-files-and-trajectories.md#molecular-dynamics-trajectories)).

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

In Firefox and Safari, the single HTML file saves to the downloads folder. Check that folder. Chrome, Edge and the
installed versions show a Save dialog.

## Open Babel is not found

Install it ([how](installation.md#open-babel-optional)), make sure `obabel` runs in a terminal, and restart
BlueMolecule. It looks on your `PATH` and in the usual Homebrew, conda and `Program Files` folders. The Open Babel group
in Save as… shows "checking…" until `obabel` answers.
