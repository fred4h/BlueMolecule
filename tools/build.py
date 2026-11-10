#!/usr/bin/env python3
"""Build BlueMolecule: a single self-contained HTML file plus Linux, Windows and macOS packages.

    python3 tools/build.py [--target linux|windows|macos|all]     # writes dist/; default: this computer's OS family

Needs only Python 3. src/examples.js and native/icons/*.png are generated files kept in the repository: delete
them to regenerate them (tools/make_examples.py needs numpy; tools/render_icons.js needs Node.js + Playwright;
tools/make_ico.py, for native/icons/BlueMolecule.ico, needs Pillow).

dist/BlueMolecule-VERSION-macOS.zip   "Install BlueMolecule.command" builds the native app (native/macos)
                                      on the user's Mac; fallback/BlueMolecule.app is the browser-window app
dist/BlueMolecule-VERSION-linux.tar.gz  launcher, Qt WebEngine window (native/linux/bluemolecule_qt.py),
                                      install.sh/uninstall.sh, desktop file, icons, MIME types
dist/BlueMolecule-VERSION-windows.zip "Install BlueMolecule.cmd" compiles the launcher (native/windows/BlueMolecule.cs)
                                      on the user's PC and adds the Start-menu and Open-with entries
dist/BlueMolecule.html is always rebuilt; the package of a target that is not built is left as it was.
"""
import argparse
import os
import time
import re
import shutil
import struct
import subprocess
import sys
import tarfile
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC, DIST, LAUNCHER = ROOT / 'src', ROOT / 'dist', ROOT / 'launcher'
NATIVE = ROOT / 'native'
VERSION = '1.4.1'
APP_ID = 'org.bluemolecule.BlueMolecule'
LINUX_ICON_SIZES = (48, 64, 128, 256, 512)
ICON_SIZES = (16, 32, 48, 64, 128, 256, 512, 1024)   # native/icons/icon_N.png, rendered from tools/icon.svg


def inline_html():
    """Inline every local <link rel="stylesheet" href="x.css"> and <script src="x.js"> from src/index.html."""
    html = (SRC / 'index.html').read_text(encoding='utf-8')

    def css(m):
        return '<style>\n' + (SRC / m.group(1)).read_text(encoding='utf-8') + '\n</style>'
    html, nc = re.subn(r'<link rel="stylesheet" href="([\w./-]+\.css)">(?:<!--BM:CSS[\w:-]*-->)?', css, html)

    def js(m):
        code = (SRC / m.group(1)).read_text(encoding='utf-8').replace('</script', '<\\/script')
        return '<script>\n' + code + '\n</script>'
    html, n = re.subn(r'<script src="([\w./-]+\.js)"></script>(?:<!--BM:JS:[\w-]+-->)?', js, html)
    if n < 10 or nc < 1 or '<script src=' in html or 'rel="stylesheet" href=' in html:
        raise SystemExit('build: inlining failed (%d scripts, %d stylesheets inlined)' % (n, nc))
    print('inlined %d scripts and %d stylesheets' % (n, nc))
    return html


def icns(png_dir, out):
    entries = [('icp4', 16), ('icp5', 32), ('icp6', 64), ('ic07', 128), ('ic08', 256), ('ic09', 512), ('ic10', 1024),
               ('ic11', 32), ('ic12', 64), ('ic13', 256), ('ic14', 512)]
    body = b''
    for code, size in entries:
        data = (png_dir / ('icon_%d.png' % size)).read_bytes()
        body += code.encode('ascii') + struct.pack('>I', len(data) + 8) + data
    out.write_bytes(b'icns' + struct.pack('>I', len(body) + 8) + body)


def add_native_bridge(html):
    """Inline src/native-bridge.js (BM.native) just before app.js, unless index.html already includes it."""
    if 'BM:native-bridge' in html:
        return html
    code = (SRC / 'native-bridge.js').read_text(encoding='utf-8').replace('</script', '<\\/script')
    app = '<script>\n' + (SRC / 'app.js').read_text(encoding='utf-8').replace('</script', '<\\/script')
    i = html.find(app)
    if i < 0:
        i = html.rfind('</body>')
    if i < 0:
        raise SystemExit('build: found no place for native-bridge.js')
    print('inlined native-bridge.js')
    return html[:i] + '<script>\n' + code + '\n</script>\n' + html[i:]


FALLBACK_PLIST = '''<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
  <key>CFBundleName</key><string>BlueMolecule</string>
  <key>CFBundleDisplayName</key><string>BlueMolecule</string>
  <key>CFBundleIdentifier</key><string>org.bluemolecule.viewer.browserwindow</string>
  <key>CFBundleVersion</key><string>{v}</string>
  <key>CFBundleShortVersionString</key><string>{v}</string>
  <key>CFBundlePackageType</key><string>APPL</string>
  <key>CFBundleExecutable</key><string>BlueMolecule</string>
  <key>CFBundleIconFile</key><string>AppIcon</string>
  <key>LSMinimumSystemVersion</key><string>11.0</string>
  <key>NSHighResolutionCapable</key><true/>
  <key>LSApplicationCategoryType</key><string>public.app-category.education</string>
</dict>
</plist>
'''

MAC_EXEC = r'''#!/bin/bash
# BlueMolecule for macOS: opens the viewer in an app-style window (Chrome, Edge, Brave, Chromium or
# Vivaldi when installed), otherwise in your default browser (Safari works too).
RES="$(cd "$(dirname "$0")/../Resources" && pwd)"
HTML="$RES/BlueMolecule.html"
URL="file://${HTML// /%20}"
for app in "Google Chrome" "Microsoft Edge" "Brave Browser" "Chromium" "Vivaldi"; do
  for base in "/Applications" "$HOME/Applications"; do
    exe="$base/$app.app/Contents/MacOS/$app"
    if [ -x "$exe" ]; then
      "$exe" --app="$URL" >/dev/null 2>&1 &
      exit 0
    fi
  done
done
exec open "$HTML"
'''

MAC_README = '''BlueMolecule {v} for macOS
==========================

Install: double-click "Install BlueMolecule.command" once.
  It builds the native BlueMolecule app on your Mac (about a minute), puts it in your Applications
  folder and adds the "bluemolecule" terminal command. It needs Apple's free Command Line Tools; if
  they are missing, it starts their installer: let that finish, then double-click the installer
  again. Running it again later is safe (it updates the app).

  If macOS says the installer "cannot be opened" because it was downloaded from the internet:
  right-click it and choose Open, or allow it under System Settings > Privacy & Security, or run
      bash "Install BlueMolecule.command"
  in Terminal from this folder.

Then:
  - open chemistry files (.xyz .cube .fchk .molden .pdb .mol .sdf .gjf ...) with Open With >
    BlueMolecule in Finder, drop any file (also Gaussian .log / ORCA .out) on the Dock icon or the
    window, or use File > Open (Cmd+O);
  - from Terminal:  bluemolecule optimisation.log freq.out orbital.cube
  - saving and exporting shows the normal macOS Save dialog;
  - formats BlueMolecule cannot read or write itself can go through Open Babel when it is installed
    (brew install open-babel, or conda install -c conda-forge openbabel).

If the native app cannot be built, the installer puts in the simpler browser-window version
(fallback/BlueMolecule.app: opens BlueMolecule in Chrome, Edge, Brave, Chromium or Vivaldi, else your
default browser) and shows why. Logs: ~/Library/Logs/BlueMolecule.log and BlueMolecule-install.log.
'''

LINUX_README = '''BlueMolecule {v} for Linux
==========================

Install for your user:      ./install.sh        (system-wide: sudo PREFIX=/usr/local ./install.sh)
Run:                        bluemolecule calc.log freq.out orbital.cube
                            (or start BlueMolecule from the applications menu, or use Open With >
                            BlueMolecule on .xyz .cube .fchk .molden .pdb .mol .sdf .gjf ... files)
Run without installing:     ./bluemolecule file.log
Uninstall:                  ./uninstall.sh

BlueMolecule opens in a window of its own, with its name and icon in the taskbar. It picks, in order:
  1. a native Qt WebEngine window, when Python has PyQt6, PySide6 or PyQt5 with QtWebEngine
     (recommended, especially on KDE; install.sh offers to install it, e.g. on Ubuntu/Kubuntu
         sudo apt install python3-pyqt6.qtwebengine );
  2. an app window of Chrome, Chromium, Edge, Brave or Vivaldi (on Wayland it runs under XWayland so
     that the desktop recognises the window);
  3. Firefox with a separate BlueMolecule profile (no tabs or toolbars).
The launcher needs python3 (standard library only); for windows 2 and 3 it serves the viewer on
http://127.0.0.1 (a fixed port, so settings persist) and quits when the window is closed. Open Babel
(obabel), when installed, handles extra file formats. Everything runs locally; on a remote cluster,
run it on your own machine and copy the files over (or use X forwarding).

Which window and why:       bluemolecule --diagnose   (log: ~/.local/state/bluemolecule/bluemolecule.log)
Options: --tier auto|qt|chromium|firefox (or BLUEMOLECULE_TIER=...), --tab (normal browser tab),
         --browser /path/to/browser, --port N, --no-open (print the URL). F12 in the Qt window opens
         the developer tools.
'''


ABOUT = '''
About: BlueMolecule is a vibe-coded project -- its code was written by an AI assistant (Claude, by
Anthropic) from its author's instructions and tested on synthetic files, so check numbers that matter
against the original output files. Documentation, source code and bug reports:
https://github.com/fred4h/BlueMolecule  (MIT License, see LICENSE).
'''


def copy_license(folder):
    if (ROOT / 'LICENSE').is_file():
        shutil.copy(ROOT / 'LICENSE', folder / 'LICENSE')


def add_tree_zip(zf, base, arcbase, exec_names):
    for p in sorted(base.rglob('*')):
        arc = str(Path(arcbase) / p.relative_to(base))
        if p.is_dir():
            zi = zipfile.ZipInfo(arc + '/', time.localtime()[:6])
            zi.external_attr = (0o40755 << 16) | 0x10
            zi.create_system = 3
            zf.writestr(zi, b'')
        else:
            zi = zipfile.ZipInfo(arc, time.localtime()[:6])
            mode = 0o100755 if (p.name in exec_names or os.access(p, os.X_OK)) else 0o100644
            zi.external_attr = mode << 16
            zi.create_system = 3
            zi.compress_type = zipfile.ZIP_DEFLATED
            zf.writestr(zi, p.read_bytes())


def build_macos(icons):
    name = 'BlueMolecule-%s-macOS' % VERSION
    mac = DIST / 'build' / name
    src = mac / 'native' / 'macos'
    src.mkdir(parents=True)
    shutil.copy(NATIVE / 'macos' / 'BlueMolecule.swift', src / 'BlueMolecule.swift')
    shutil.copy(NATIVE / 'macos' / 'Info.plist', src / 'Info.plist')     # template: the installer fills in @VERSION@
    installer = (NATIVE / 'macos' / 'Install BlueMolecule.command').read_text(encoding='utf-8')
    if 'VERSION="@VERSION@"' not in installer:
        raise SystemExit('build: VERSION placeholder missing in Install BlueMolecule.command')
    cmd = mac / 'Install BlueMolecule.command'
    cmd.write_text(installer.replace('VERSION="@VERSION@"', 'VERSION="%s"' % VERSION, 1), encoding='utf-8')
    cmd.chmod(0o755)
    shutil.copy(DIST / 'BlueMolecule.html', mac / 'BlueMolecule.html')
    icns(icons, mac / 'AppIcon.icns')
    shutil.copy(LAUNCHER / 'bluemolecule.py', mac / 'bluemolecule.py')
    (mac / 'bluemolecule.py').chmod(0o755)

    # fallback: the browser-window app (a shell script), installed when the native build is not possible
    app = mac / 'fallback' / 'BlueMolecule.app' / 'Contents'
    (app / 'MacOS').mkdir(parents=True)
    (app / 'Resources').mkdir(parents=True)
    (app / 'Info.plist').write_text(FALLBACK_PLIST.format(v=VERSION))
    (app / 'PkgInfo').write_text('APPL????')
    exe = app / 'MacOS' / 'BlueMolecule'
    exe.write_text(MAC_EXEC)
    exe.chmod(0o755)
    shutil.copy(DIST / 'BlueMolecule.html', app / 'Resources' / 'BlueMolecule.html')
    shutil.copy(LAUNCHER / 'bluemolecule.py', app / 'Resources' / 'bluemolecule.py')
    (app / 'Resources' / 'bluemolecule.py').chmod(0o755)
    shutil.copy(mac / 'AppIcon.icns', app / 'Resources' / 'AppIcon.icns')

    (mac / 'README.txt').write_text(MAC_README.format(v=VERSION) + ABOUT)
    copy_license(mac)
    with zipfile.ZipFile(DIST / (name + '.zip'), 'w') as zf:
        add_tree_zip(zf, mac, name, {'BlueMolecule', 'bluemolecule.py', 'Install BlueMolecule.command'})


def build_linux(icons):
    name = 'BlueMolecule-%s-linux' % VERSION
    lin = DIST / 'build' / name
    (lin / 'icons').mkdir(parents=True)
    shutil.copy(DIST / 'BlueMolecule.html', lin / 'BlueMolecule.html')
    shutil.copy(LAUNCHER / 'bluemolecule.py', lin / 'bluemolecule')
    (lin / 'bluemolecule').chmod(0o755)
    shutil.copy(NATIVE / 'linux' / 'bluemolecule_qt.py', lin / 'bluemolecule_qt.py')
    (lin / 'bluemolecule_qt.py').chmod(0o755)
    for s in LINUX_ICON_SIZES:
        shutil.copy(icons / ('icon_%d.png' % s), lin / 'icons' / ('%s-%d.png' % (APP_ID, s)))
    for f in ('install.sh', 'uninstall.sh'):
        shutil.copy(NATIVE / 'linux' / f, lin / f)
        (lin / f).chmod(0o755)
    for f in (APP_ID + '.desktop', APP_ID + '.xml'):
        shutil.copy(NATIVE / 'linux' / f, lin / f)
    (lin / 'README.txt').write_text(LINUX_README.format(v=VERSION) + ABOUT)
    copy_license(lin)
    with tarfile.open(DIST / (name + '.tar.gz'), 'w:gz') as tf:
        tf.add(lin, arcname=name)


WINDOWS_README = '''BlueMolecule {v} for Windows (experimental)
===========================================

This Windows version has not yet been tested on a real Windows PC. Please report what happens:
https://github.com/fred4h/BlueMolecule/issues

Install: extract this zip file first (right-click it > Extract All), then double-click
"Install BlueMolecule.cmd". If Windows says it protected your PC, click "More info" > "Run anyway".
  It installs for your user only (no administrator rights) in %LOCALAPPDATA%\\Programs\\BlueMolecule:
  - it compiles the small BlueMolecule.exe launcher with the C# compiler that is part of Windows
    (the .NET Framework), so nothing has to be installed first;
  - it adds BlueMolecule to the Start menu and to "Open with" for chemistry files (.log .out .xyz .gjf
    .fchk .cube .molden .pdb .cif .sdf .mol2 .gro ...), without changing their default programs;
  - it registers an uninstaller: Settings > Apps > Installed apps > BlueMolecule.
  Running it again updates the installation.

BlueMolecule opens in an app window of Microsoft Edge (or Chrome, Brave or Vivaldi) with its own
browser profile. Files given to it with "Open with" are handed to the window inside a temporary copy of
the page; files larger than 64 MB must be opened from the window instead (Open... or drag and drop).
Open Babel conversions are not available in the Windows version yet.

Without installing: open BlueMolecule.html in Edge or Chrome.
Diagnostics:  "%LOCALAPPDATA%\\Programs\\BlueMolecule\\BlueMolecule.exe" --diagnose
Log:          %LOCALAPPDATA%\\BlueMolecule\\bluemolecule.log
'''


def windows_text(path):
    """A Windows source file (cmd, ps1, cs): ASCII only, version filled in, CRLF line ends."""
    text = path.read_text(encoding='utf-8')
    bad = sorted({c for c in text if ord(c) > 126 or (ord(c) < 32 and c not in '\n\t')})
    if bad:
        raise SystemExit('build: %s must be ASCII (Windows PowerShell 5.1 and cmd read it as ANSI): %r' % (path.name, bad))
    v4 = '.'.join((VERSION.split('.') + ['0', '0', '0'])[:4])
    text = text.replace('@VERSION4@', v4).replace('@VERSION@', VERSION)
    return text.replace('\r\n', '\n').replace('\n', '\r\n').encode('ascii')


def build_windows():
    name = 'BlueMolecule-%s-windows' % VERSION
    win = DIST / 'build' / name
    sub = win / 'windows'
    sub.mkdir(parents=True)
    for f in ('install.ps1', 'BlueMolecule.cs'):
        if '@VERSION@' not in (NATIVE / 'windows' / f).read_text(encoding='utf-8'):
            raise SystemExit('build: VERSION placeholder missing in native/windows/' + f)
    for f in ('install.ps1', 'uninstall.ps1', 'BlueMolecule.cs'):
        (sub / f).write_bytes(windows_text(NATIVE / 'windows' / f))
    (win / 'Install BlueMolecule.cmd').write_bytes(windows_text(NATIVE / 'windows' / 'Install BlueMolecule.cmd'))
    shutil.copy(NATIVE / 'icons' / 'BlueMolecule.ico', sub / 'BlueMolecule.ico')
    shutil.copy(DIST / 'BlueMolecule.html', win / 'BlueMolecule.html')
    (win / 'README.txt').write_bytes((WINDOWS_README.format(v=VERSION) + ABOUT).replace('\n', '\r\n').encode('ascii'))
    copy_license(win)
    with zipfile.ZipFile(DIST / (name + '.zip'), 'w') as zf:
        add_tree_zip(zf, win, name, set())


def host_target():
    """The package for this computer's OS family: macos on a Mac, windows on Windows, linux on Linux, else all."""
    if sys.platform == 'darwin':
        return 'macos'
    if sys.platform == 'win32':
        return 'windows'
    return 'linux' if sys.platform.startswith('linux') else 'all'


PACKAGE_SUFFIX = {'macos': ('macOS', '.zip'), 'linux': ('linux', '.tar.gz'), 'windows': ('windows', '.zip')}


def outputs(target):
    """What one target writes in dist/ (the package and its build folder)."""
    label, ext = PACKAGE_SUFFIX[target]
    name = 'BlueMolecule-%s-%s' % (VERSION, label)
    return [DIST / (name + ext), DIST / 'build' / name]


def clean(targets):
    """Remove what this build writes again; other targets' packages and other files in dist/ stay."""
    for p in [DIST / 'BlueMolecule.html', DIST / 'build' / 'icons'] + [q for t in targets for q in outputs(t)]:
        if p.is_dir():
            shutil.rmtree(p)
        elif p.exists():
            p.unlink()


def check_python_sources():
    """The launcher and the Qt window must at least compile (they are not run during the build)."""
    for p in (LAUNCHER / 'bluemolecule.py', NATIVE / 'linux' / 'bluemolecule_qt.py'):
        compile(p.read_text(encoding='utf-8'), str(p), 'exec')


def main(argv=None):
    ap = argparse.ArgumentParser(description='Build dist/BlueMolecule.html and the installable packages.')
    ap.add_argument('--target', choices=('linux', 'macos', 'windows', 'all'), default=host_target(),
                    help="package(s) to build (default for this computer: %s)" % host_target())
    args = ap.parse_args(argv)
    targets = ('macos', 'linux', 'windows') if args.target == 'all' else (args.target,)
    if 'linux' in targets:
        check_python_sources()
    if not (SRC / 'examples.js').is_file():   # generated and committed; delete it to regenerate (needs numpy)
        subprocess.run([sys.executable, str(ROOT / 'tools' / 'make_examples.py')], check=True)
    clean(targets)
    DIST.mkdir(exist_ok=True)
    html = add_native_bridge(inline_html())
    (DIST / 'BlueMolecule.html').write_text(html, encoding='utf-8')
    print('BlueMolecule.html: %.0f KB' % (len(html.encode()) / 1024))

    icons = NATIVE / 'icons'   # committed PNGs; delete them after editing tools/icon.svg to render them again
    if not all((icons / ('icon_%d.png' % s)).is_file() for s in ICON_SIZES):
        subprocess.run(['node', str(ROOT / 'tools' / 'render_icons.js'), str(icons)], check=True)   # Node.js + Playwright
    if 'macos' in targets:
        build_macos(icons)
    if 'linux' in targets:
        build_linux(icons)
    if 'windows' in targets:
        build_windows()

    built = {outputs(t)[0].name for t in targets} | {'BlueMolecule.html'}
    for p in sorted(DIST.iterdir()):
        if p.is_file():
            print('%-40s %8.0f KB%s' % (p.name, p.stat().st_size / 1024, '' if p.name in built else '   (not rebuilt)'))


if __name__ == '__main__':
    main()
