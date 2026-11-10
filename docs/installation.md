# Installation

BlueMolecule comes in three forms, all from the [Releases page](https://github.com/fred4h/BlueMolecule/releases):

| Download | For | What you get |
|---|---|---|
| `BlueMolecule-VERSION-linux.tar.gz` | Linux (main platform) | The `bluemolecule` command, a menu entry, file associations and its own window |
| `BlueMolecule-VERSION-windows.zip` | Windows 10 and 11 (experimental) | A Start-menu entry, "Open with" for chemistry files and its own Edge window |
| `BlueMolecule-VERSION-macOS.zip` | macOS (experimental) | A native app in Applications and the `bluemolecule` command |
| `BlueMolecule.html` | Any computer with a modern browser | The whole viewer in one offline file |

All of them contain the same viewer. Everything runs on your computer and files are never uploaded.

## Linux

### Requirements

- `python3`: only the standard library is needed.
- For its own window, any **one** of:
  - Python's Qt WebEngine (recommended);
  - Chrome, Chromium, Edge, Brave or Vivaldi;
  - Firefox.
- Optional: [Open Babel](#open-babel-optional) for extra file formats.

It works on X11 and Wayland desktops (KDE Plasma, GNOME and others).

### Install

```bash
tar xzf BlueMolecule-1.4.1-linux.tar.gz
cd BlueMolecule-1.4.1-linux
./install.sh
```

`install.sh` installs for your user, under `~/.local`:
- `~/.local/share/bluemolecule/`: the viewer, the launcher and the Qt window script.
- `~/.local/bin/bluemolecule`: the command. If `~/.local/bin` is not on your `PATH`, the installer says how to add
  it.
- A desktop entry, so BlueMolecule appears in the applications menu with its own icon.
- MIME types, so chemistry files offer **Open With → BlueMolecule**. This covers `.xyz .cube .fchk .molden .mol
  .sdf .mol2 .pdb .cif .gro .gjf .cml .cjson .xsf`. Gaussian `.log` and ORCA `.out` files are plain text, so open
  them from inside BlueMolecule or from the command line.

Running `./install.sh` again updates or repairs an installation, including one from an older version.

- **System-wide:** `sudo PREFIX=/usr/local ./install.sh`
- **Without installing:** `./bluemolecule calc.log` from the unpacked folder.
- **Uninstall:** `./uninstall.sh`. Your settings in `~/.local/share/bluemolecule/` are kept unless you delete that
  folder.

### Its own window: the native Qt window

BlueMolecule tries three kinds of window, in this order:

1. **A native Qt WebEngine window.** This is used when a Python on your computer has PyQt6, PySide6 or PyQt5 with
   QtWebEngine. It is the best option, especially on KDE: a real application window with BlueMolecule's name and
   icon in the taskbar.
2. **An app window of Chrome, Chromium, Edge, Brave or Vivaldi**, with no tabs or address bar, in a separate profile.
   On Wayland it runs under XWayland so that the desktop matches the window to BlueMolecule's icon.
3. **Firefox** with a separate BlueMolecule profile.

At the end, `install.sh` tells you which window you will get. If Qt WebEngine is missing, it shows the command that
installs it and asks before running it; nothing is ever installed with `sudo` without asking. On Ubuntu and Kubuntu:

```bash
sudo apt install python3-pyqt6.qtwebengine
```

On Fedora, `sudo dnf install python3-pyqt6-webengine`. On Arch, `sudo pacman -S python-pyqt6-webengine`.
`BLUEMOLECULE_NO_PROMPT=1 ./install.sh` skips the question.

The launcher checks several Pythons: the one running it, `/usr/bin/python3`, `/usr/local/bin/python3` and `python3` on
your `PATH`. A conda or virtual-environment Python that comes first on `PATH` therefore does not hide the system's
Qt. To see which window will be used and why:

```bash
bluemolecule --diagnose
```

### Command line

```text
bluemolecule [options] [files ...]

  files                  any supported files; .gz-compressed files are fine
  --tier auto|qt|chromium|firefox
                         force a window kind (also: environment variable BLUEMOLECULE_TIER)
  --browser PATH         use this browser executable
  --tab                  open in a normal browser tab instead of an app window
  --port N               preferred local port (default 47591; the next 19 are tried when it is busy;
                         0 = any free port)
  --no-open              print the address instead of opening a window
  --diagnose             show the environment and which window kind would be used, and why
  --version
```

With the browser window types, the launcher serves the viewer on `http://127.0.0.1:47591/<random token>/`. The
fixed port lets the browser keep your settings between sessions. The launcher quits when the window is closed.
It writes a log to `~/.local/state/bluemolecule/bluemolecule.log`. In the Qt window, **F12** opens the developer
tools.

## Windows (experimental)

> [!WARNING]
> The Windows version has **not yet been tested on a real Windows PC**. It was written and checked on Linux, and
> its small launcher is compiled on your PC during installation. Please report what happens.

### Install

1. Extract `BlueMolecule-1.4.1-windows.zip`: right-click it and choose **Extract All**.
2. Double-click **Install BlueMolecule.cmd**. If Windows says it protected your PC, click **More info → Run anyway**.

The installer needs **no administrator rights and nothing installed first**. It uses only Windows PowerShell, the
C# compiler of the .NET Framework and Microsoft Edge, all of which come with Windows 10 and 11. It:
- copies BlueMolecule into `%LOCALAPPDATA%\Programs\BlueMolecule`;
- compiles the small **BlueMolecule.exe** launcher from `BlueMolecule.cs`, so the launcher has BlueMolecule's name
  and icon;
- adds BlueMolecule to the **Start menu**;
- adds **Open with → BlueMolecule** for chemistry files, including Gaussian `.log`/`.gjf`, ORCA `.out`/`.inp`,
  `.xyz`, `.fchk`, `.cube`, `.molden`, `.pdb`, `.cif`, `.sdf`, `.mol2`, `.gro` and MD files. It never changes the
  program that a file type opens with by default;
- registers an uninstaller under **Settings → Apps → Installed apps → BlueMolecule**.

Running the installer again updates the installation. If the launcher cannot be compiled, the installer says why
and makes the Start-menu entry open BlueMolecule directly. In that case **Open with** is not available.

### How it runs

BlueMolecule opens in an **app window of Microsoft Edge**, with no tabs or address bar and a separate BlueMolecule
profile in `%LOCALAPPDATA%\BlueMolecule\browser`. Chrome, Brave or Vivaldi are used when Edge is missing. To
choose the browser, set the environment variable `BLUEMOLECULE_BROWSER` to the browser's `.exe`.

Files opened with **Open with**, or dropped on the program's icon, are handed to the window inside a temporary copy
of the page. That copy is kept for a day in `%LOCALAPPDATA%\BlueMolecule\launch`. This needs no local server and
no special browser settings, but it works for files up to **64 MB**. Open larger files from the window with
**Open…** or by drag and drop, which read files of any size. Saving shows the Windows Save dialog.

- Diagnostics: `"%LOCALAPPDATA%\Programs\BlueMolecule\BlueMolecule.exe" --diagnose`
- Log: `%LOCALAPPDATA%\BlueMolecule\bluemolecule.log`
- Open Babel conversions are not available in the Windows version yet.
- **Uninstall:** use **Settings → Apps → Installed apps → BlueMolecule → Uninstall**. This removes the program, the
  Start-menu entry and the Open-with entries. Your settings in `%LOCALAPPDATA%\BlueMolecule` stay until you delete
  that folder.

## macOS (experimental)

> [!WARNING]
> The macOS version has **not yet been tested on a real Mac**. It was written and checked on Linux, and the Swift
> code is compiled on your Mac during installation. Please report what happens.

1. Unzip `BlueMolecule-1.4.1-macOS.zip`.
2. Double-click **Install BlueMolecule.command**. The installer:
   - compiles the small native app (a WKWebView window, Swift source in `native/macos/`) with Apple's free Command
     Line Tools; if they are missing, it starts their installer, and you run the BlueMolecule installer again
     afterwards;
   - signs the app ad hoc and installs it in `~/Applications` (or `/Applications`);
   - adds the `bluemolecule` terminal command.
3. If macOS says the installer "cannot be opened" because it was downloaded, do one of these:
   - right-click it and choose **Open**;
   - allow it under **System Settings → Privacy & Security**;
   - run `bash "Install BlueMolecule.command"` in Terminal.

If the native app cannot be built, the installer installs a simpler browser-window app instead and tells you why.
That app opens BlueMolecule in a Chromium-family browser's app window, or in your default browser. Logs are in
`~/Library/Logs/BlueMolecule.log` and `~/Library/Logs/BlueMolecule-install.log`.

## Any computer: the single HTML file

Download `BlueMolecule.html` and open it in Chrome, Edge, Firefox or Safari. It needs no installation and no
network. The differences from the installed versions:
- saving shows a Save dialog in Chrome and Edge, and goes to the downloads folder in other browsers;
- Open Babel conversions are not available;
- very large files are limited by the browser's memory;
- your browser may not offer to open chemistry files with it.

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

The launcher looks for `obabel` on your `PATH` and in the usual Homebrew and conda locations.

## Building from source

```bash
git clone https://github.com/fred4h/BlueMolecule.git
cd BlueMolecule
python3 tools/build.py              # the package for this computer's OS, plus dist/BlueMolecule.html
python3 tools/build.py --target all # the Linux, Windows and macOS packages
```

The build needs only Python 3 and writes everything into `dist/`. For running the tests, see the
[developer guide](developer-guide.md#tests).
