# Installation

Since version 2.0, BlueMolecule is a desktop program of its own on Linux, Windows and macOS: a native window with
BlueMolecule's name and icon, an installer, **Open with** entries for chemistry files and a built-in updater. All
downloads are on the [Releases page](https://github.com/fred4h/BlueMolecule/releases):

| Download | For |
|---|---|
| `BlueMolecule_VERSION_amd64.deb` | Ubuntu, Kubuntu, Debian, Linux Mint and other Debian-based Linux (64-bit Intel/AMD) |
| `BlueMolecule_VERSION_amd64.AppImage` | Any other 64-bit Linux: one file that runs without installing |
| `BlueMolecule_VERSION_x64-setup.exe` | Windows 10 and 11 |
| `BlueMolecule_VERSION_universal.dmg` | macOS 10.15 or later, Apple silicon and Intel Macs |
| `BlueMolecule.html` | Any computer with a modern browser: the whole viewer in one offline file |

They all contain the same viewer. Everything runs on your computer and files are never uploaded.

> [!NOTE]
> The installers are **not code-signed** (that needs paid certificates). Windows and macOS therefore warn once,
> before the first start, that the publisher is unknown; the steps below say how to go on. Updates are checked
> against BlueMolecule's own signing key either way ([Updates](#updates)).

## Linux

### The .deb package (Ubuntu, Kubuntu, Debian…)

```bash
sudo apt install ./BlueMolecule_2.0.1_amd64.deb
```

This installs the `bluemolecule` command, a menu entry under Science / Education, and the chemistry file types:
chemistry files then offer **Open With → BlueMolecule**. This covers `.xyz .cube .fchk .molden .mol .sdf .mol2 .pdb
.cif .gro .gjf .cml .cjson .xsf`. Gaussian `.log` and ORCA `.out` files are plain text: open them from inside
BlueMolecule, by dragging them onto its window, with **Open With → Other application**, or from a terminal:
`bluemolecule calc.log`.

The package needs WebKitGTK 4.1 (`libwebkit2gtk-4.1-0`), which apt installs with it; it is part of Ubuntu 22.04,
Debian 12 and later. Uninstall with `sudo apt remove blue-molecule` (the package name; the command is `bluemolecule`).

### The AppImage (any distribution)

```bash
chmod +x BlueMolecule_2.0.1_amd64.AppImage
./BlueMolecule_2.0.1_amd64.AppImage calc.log
```

It brings its own libraries and needs nothing installed. Tools such as AppImageLauncher or Gear Lever add it to the
menu. The built-in updater updates the AppImage in place.

### Coming from BlueMolecule 1.x

Remove the 1.x installation first: run `./uninstall.sh` from the unpacked 1.x folder (for example
`BlueMolecule-1.5.0-linux`). Otherwise its `~/.local/bin/bluemolecule` command comes first on your `PATH`, and the
menu shows two BlueMolecule entries. View settings and calculation defaults made in 1.x do not carry over to the
new window.

## Windows

1. Run `BlueMolecule_2.0.1_x64-setup.exe`. If Windows says it protected your PC, click **More info → Run anyway**.
2. The installer needs **no administrator rights**. It installs for your user, in `%LOCALAPPDATA%\BlueMolecule`,
   adds BlueMolecule to the **Start menu**, and adds **Open with → BlueMolecule** for chemistry files: Gaussian
   `.log`/`.gjf`, ORCA `.out`/`.inp`, `.xyz`, `.fchk`, `.cube`, `.molden`, `.pdb`, `.cif`, `.sdf`, `.mol2`, `.gro`
   and MD files. It never changes the program a file type opens with by default.
3. The window uses Microsoft Edge WebView2, which Windows 10 and 11 already have; the installer fetches it if not.

Uninstall in **Settings → Apps → Installed apps → BlueMolecule**. This also removes the Open-with entries.

## macOS

1. Open `BlueMolecule_2.0.1_universal.dmg` and drag **BlueMolecule** to **Applications**.
2. The first time, macOS says it cannot check the app for malicious software (it is not notarized):
   - **macOS 15 Sequoia and later:** try to open it once, then go to **System Settings → Privacy & Security**,
     scroll down to the message about BlueMolecule and click **Open Anyway**.
   - **Earlier versions:** right-click (or Control-click) BlueMolecule in Applications and choose **Open**.
   - Or, in Terminal: `xattr -dr com.apple.quarantine /Applications/BlueMolecule.app`
3. After that, BlueMolecule opens normally, also from **Open With** in the Finder for chemistry files.

## Updates

**Help → Check for updates…** asks GitHub for the latest release. If there is a newer one, BlueMolecule shows its
release notes and, when you agree, downloads it, checks its signature, installs it and restarts. The signature is
made with a key that only BlueMolecule's maintainer holds, so a damaged or altered download is never installed.

Checking is **BlueMolecule's only use of the network**, and it happens only when you ask. In the same dialog you can
tick **Check for updates when BlueMolecule starts**; it is off until you do.

On Linux, the AppImage updates itself. A `.deb` installation is updated by installing the new `.deb`: the update
dialog then offers to open the Releases page instead.

## Any computer: the single HTML file

Download `BlueMolecule.html` and open it in Chrome, Edge, Firefox or Safari. It needs no installation and no
network. The differences from the desktop program:
- saving shows a Save dialog in Chrome and Edge, and goes to the downloads folder in other browsers;
- Open Babel conversions are not available;
- very large files are limited by the browser's memory;
- there are no Open-with entries and no updater.

## Open Babel (optional)

Install Open Babel and BlueMolecule uses it automatically, through the `obabel` program, to:
- open formats it cannot read itself;
- save in many more formats (the **Open Babel** group in **Save as…**).

| System | Command |
|---|---|
| Ubuntu / Debian | `sudo apt install openbabel` |
| Fedora | `sudo dnf install openbabel` |
| conda (any OS) | `conda install -c conda-forge openbabel` |
| macOS (Homebrew) | `brew install open-babel` |
| Windows | the installer from the [Open Babel releases](https://github.com/openbabel/openbabel/releases), or conda |

BlueMolecule looks for `obabel` on your `PATH` and in the usual Homebrew, conda and `Program Files` locations, so
it is found also when BlueMolecule is started from the menu or the Dock.

## Building from source

```bash
git clone https://github.com/fred4h/BlueMolecule.git
cd BlueMolecule
python3 tools/build.py              # dist/BlueMolecule.html (and the 1.x-style packages)
```

The desktop program is built with Rust and Node.js: see [desktop/README.md](../desktop/README.md). The cloud builds
of the repository build and test it on all three systems. For running the tests, see the
[developer guide](developer-guide.md#tests).

## Older versions (1.x)

Versions 1.0 to 1.5 came as a Linux package with a Python launcher and a Qt window, and as experimental Windows and
macOS packages that compiled a small launcher during installation. They are still on the Releases page. Their
guide is [the 1.5 version of this page](https://github.com/fred4h/BlueMolecule/blob/v1.5.0/docs/installation.md).
