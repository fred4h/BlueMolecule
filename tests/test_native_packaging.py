#!/usr/bin/env python3
"""Checks of the native packaging (no Mac or Swift compiler needed):

  * Swift sanity: (), [] and {} balance outside strings and comments; every message-handler name
    used by src/native-bridge.js is registered and handled in BlueMolecule.swift; the JS embedded in
    the Swift source parses; a few API/structure rules for a single-file AppKit program.
  * Info.plist template: document types, imported UTIs, no generic .log/.out/.com/.txt.
  * "Install BlueMolecule.command": bash syntax, no bash-4-only features (macOS has bash 3.2), and a
    dry run on a simulated Mac (fake xcrun/swiftc/codesign/...) for the native, compiler-error and
    missing-Command-Line-Tools paths.
  * Linux: desktop file, shared-mime-info XML, the Qt window compiles; install.sh/uninstall.sh into a
    temporary prefix: 1.0/1.1 leftovers removed, run twice, KDE/desktop caches refreshed (fake tools),
    the Qt WebEngine hint shown without running sudo, profiles kept by uninstall.sh.
  * The built packages in dist/ when they exist.
Usage: python3 tests/test_native_packaging.py
"""
import os
import plistlib
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MAC, LIN = ROOT / 'native' / 'macos', ROOT / 'native' / 'linux'
SWIFT = (MAC / 'BlueMolecule.swift').read_text()
BRIDGE = (ROOT / 'src' / 'native-bridge.js').read_text()
INSTALLER = (MAC / 'Install BlueMolecule.command').read_text()
VERSION = re.search(r"^VERSION = '([\d.]+)'", (ROOT / 'tools' / 'build.py').read_text(), re.M).group(1)
EXTS = {'cube', 'fchk', 'molden', 'xyz', 'mol', 'sdf', 'mol2', 'pdb', 'cif', 'gro', 'gjf', 'cml', 'cjson', 'xsf'}
GENERIC = {'log', 'out', 'com', 'txt'}
TMP = Path(tempfile.mkdtemp(prefix='bm-native-test-'))
passes = failures = 0


def check(cond, msg):
    global passes, failures
    if cond:
        passes += 1
    else:
        failures += 1
        print('  FAIL:', msg)


def swift_balance(src):
    """Problems with (), [], {} nesting, skipping comments and string literals (incl. raw/multi-line
    strings and \\( ) interpolation)."""
    problems, stack, mode = [], [], None      # mode None = code, else (hashes, multiline)
    pairs = {')': '(', ']': '[', '}': '{'}
    i, n, line = 0, len(src), 1
    while i < n:
        c = src[i]
        if c == '\n':
            line += 1
            if mode is not None and not mode[1]:
                problems.append('line %d: newline in a string literal' % line)
                mode = None
            i += 1
            continue
        if mode is None:
            if src.startswith('//', i):
                j = src.find('\n', i)
                i = n if j < 0 else j
            elif src.startswith('/*', i):
                depth, i = 1, i + 2
                while i < n and depth:
                    if src.startswith('/*', i):
                        depth, i = depth + 1, i + 2
                    elif src.startswith('*/', i):
                        depth, i = depth - 1, i + 2
                    else:
                        line += src[i] == '\n'
                        i += 1
                if depth:
                    problems.append('unterminated block comment')
            elif c in '#"':
                j = i
                while j < n and src[j] == '#':
                    j += 1
                if j < n and src[j] == '"':
                    multi = src.startswith('"""', j)
                    mode = (j - i, multi)
                    i = j + (3 if multi else 1)
                else:
                    i = max(j, i + 1)      # #if, #selector, #available ...
            elif c in '([{':
                stack.append((c, line, None))
                i += 1
            elif c in ')]}':
                if not stack or stack[-1][0] != pairs[c]:
                    problems.append('line %d: unexpected %s' % (line, c))
                else:
                    opener = stack.pop()
                    if opener[2] is not None:
                        mode = opener[2]   # back inside the string after \( ... )
                i += 1
            else:
                i += 1
            continue
        hashes, multi = mode
        esc = '\\' + '#' * hashes
        term = ('"""' if multi else '"') + '#' * hashes
        if src.startswith(esc, i):
            k = i + len(esc)
            if k < n and src[k] == '(':
                stack.append(('(', line, mode))
                mode = None
                i = k + 1
            else:
                line += k < n and src[k] == '\n'
                i = k + 1
        elif src.startswith(term, i):
            mode = None
            i += len(term)
        else:
            i += 1
    if mode is not None:
        problems.append('unterminated string literal')
    problems += ['line %d: unclosed %s' % (ln, o) for o, ln, _ in stack]
    return problems


def bash_ok(path_or_text):
    if isinstance(path_or_text, Path):
        r = subprocess.run(['bash', '-n', str(path_or_text)], capture_output=True, text=True)
    else:
        r = subprocess.run(['bash', '-n'], input=path_or_text, capture_output=True, text=True)
    return r.returncode == 0, r.stderr


try:
    # ---------------- Swift sanity ----------------
    print('Swift sanity')
    check(swift_balance('let a = "x(\\(f(1)))" /* ( */ // [\nlet b = #"""\n{"\n"""#\n') == [], 'balance checker self-test (good input)')
    check(swift_balance('func f() { g(1 }') != [] and swift_balance('let s = "abc\n"') != [], 'balance checker self-test (bad input)')
    problems = swift_balance(SWIFT)
    check(problems == [], 'Swift brackets/strings balance: %s' % problems[:5])
    m = re.search(r'let kHandlerNames = \[([^\]]*)\]', SWIFT)
    swift_names = set(re.findall(r'"(bm\w+)"', m.group(1))) if m else set()
    js_names = set(re.findall(r"(?:post|call)\(\s*'(bm\w+)'", BRIDGE)) | set(re.findall(r'\.(bm[A-Z]\w*)\b', BRIDGE))
    # bmRead is sent only for large files of the Linux Qt window; bmHost is its QWebChannel object, not a message
    linux_only = {'bmRead', 'bmHost'}
    check(js_names and js_names - linux_only <= swift_names, 'JS handler names %s registered in Swift %s' % (sorted(js_names), sorted(swift_names)))
    qt_names = set(re.findall(r"'(bm\w+)': self\.on_", (LIN / 'bluemolecule_qt.py').read_text()))
    check(js_names - {'bmHost'} <= qt_names and 'bmHost' in js_names and "registerObject('bmHost'" in (LIN / 'bluemolecule_qt.py').read_text(),
          'JS message names %s handled by the Linux Qt window %s' % (sorted(js_names), sorted(qt_names)))
    check(js_names >= {'bmReady', 'bmSave', 'bmObabel', 'bmCopy'}, 'bridge uses bmReady/bmSave/bmObabel/bmCopy: %s' % sorted(js_names))
    for name in swift_names:
        check('case "%s":' % name in SWIFT, 'Swift handles %s' % name)
    check('messageHandlers.bmLog' in SWIFT and 'bmLog' in swift_names, 'error forwarding posts to a registered bmLog')
    for fn in ('_resolve', '_receiveFile'):
        check('BM.native.%s(' % fn in SWIFT and re.search(r'function %s\(' % fn, BRIDGE) and re.search(r'\b%s,' % fn, BRIDGE),
              'Swift calls BM.native.%s, which the bridge defines and exports' % fn)
    js = re.search(r'let kErrorForwardingJS = #"""\n(.*?)\n"""#', SWIFT, re.S)
    r = subprocess.run(['node', '-e', 'new Function(require("fs").readFileSync(0, "utf8"))'], input=js.group(1) if js else '!',
                       capture_output=True, text=True)
    check(js and r.returncode == 0, 'embedded error-forwarding JS parses: %s' % r.stderr.strip()[:200])
    for snippet in re.findall(r'callJS\("([^"]*)"', SWIFT) + re.findall(r'evaluateJavaScript\("([^"]*)"\)', SWIFT):
        r = subprocess.run(['node', '-e', 'new Function("m", require("fs").readFileSync(0, "utf8"))'], input=snippet, capture_output=True, text=True)
        check(r.returncode == 0, 'JS snippet parses: %s' % snippet)
    start = SWIFT.index('bmApplication.run()')
    check(SWIFT.rstrip().endswith('bmApplication.run()'), 'top-level start-up code is at the end of the file')
    for g in ('let kHandlerNames', 'let kErrorForwardingJS', 'let kWindowAutosaveName'):
        check(0 <= SWIFT.find(g) < start, '%s is defined before the top-level code runs' % g)
    check('@main' not in SWIFT and 'NSApplicationMain' not in SWIFT, 'no @main next to top-level code')
    check(re.search(r'#if compiler\(>=5\.8\)\s*if #available\(macOS 13\.3, \*\) \{\s*webView\.isInspectable = true', SWIFT), 'isInspectable guarded by #available(macOS 13.3, *)')
    for needle in ('loadFileURL(url, allowingReadAccessTo:', 'runOpenPanelWith parameters: WKOpenPanelParameters',
                   'func application(_ application: NSApplication, open urls: [URL])', 'NSSavePanel()', 'NSPasteboard.general',
                   'setFrameAutosaveName', 'observe(\\.title', 'Library/Logs/BlueMolecule.log', 'process.executableURL',
                   'process.standardInput', 'process.standardOutput', 'ObabelRunner.timeout', 'static let timeout: TimeInterval = 60'):
        check(needle in SWIFT, 'Swift contains %r' % needle)
    for d in ('/opt/homebrew/bin', '/usr/local/bin', '/miniconda3/bin', '/anaconda3/bin', '/miniforge3/bin', '/opt/miniconda3/bin'):
        check('"%s"' % d in SWIFT or '+ "%s"' % d in SWIFT, 'obabel search path %s' % d)
    check(SWIFT.count('kHandlerNames') >= 2, 'handler names registered from kHandlerNames')

    # ---------------- Info.plist ----------------
    print('Info.plist')
    plist = plistlib.loads((MAC / 'Info.plist').read_text().replace('@VERSION@', VERSION).encode())
    check(plist['CFBundleIdentifier'] == 'org.bluemolecule.viewer' and plist['CFBundleExecutable'] == 'BlueMolecule'
          and plist['CFBundleIconFile'] == 'AppIcon' and plist['LSMinimumSystemVersion'] == '11.0', 'bundle basics')
    check(plist['CFBundleShortVersionString'] == VERSION, 'version filled in')
    utis = {d['UTTypeIdentifier']: d for d in plist['UTImportedTypeDeclarations']}
    exts = {e for d in utis.values() for e in d['UTTypeTagSpecification']['public.filename-extension']}
    check(EXTS <= exts and not (exts & GENERIC), 'imported UTIs cover %s and no generic extensions (have %s)' % (sorted(EXTS), sorted(exts)))
    viewer = [d for d in plist['CFBundleDocumentTypes'] if d['LSHandlerRank'] == 'Alternate']
    check(len(viewer) == 1 and viewer[0]['CFBundleTypeRole'] == 'Viewer' and set(viewer[0]['LSItemContentTypes']) == set(utis),
          'one Viewer/Alternate document type listing every imported UTI')
    others = [d for d in plist['CFBundleDocumentTypes'] if d['LSHandlerRank'] != 'Alternate']
    check(all(d['LSHandlerRank'] == 'None' for d in others), 'catch-all types never claim files (rank None)')

    # ---------------- installer script ----------------
    print('macOS installer')
    ok, err = bash_ok(MAC / 'Install BlueMolecule.command')
    check(ok, 'installer bash syntax: %s' % err)
    for bad in (r'declare -A', r'\bmapfile\b', r'\breadarray\b', r'\$\{\w+,,\}', r'\$\{\w+\^\^\}', r'\s\|&\s', r'&>>', r'\bset -u\b', r'local \w+=\('):
        check(not re.search(bad, INSTALLER), 'installer avoids bash-4-only %s' % bad)
    for needle in ('VERSION="@VERSION@"', 'xcode-select --install', 'xcrun --find swiftc', '-swift-version 5', '"$arch-apple-macos11.0"',
                   'for arch in arm64 x86_64', 'lipo -create', 'codesign --force --deep -s -', '"$LSREGISTER" -f', '$HOME/Applications',
                   '/Applications', 'exec open -a "$APP"', 'FALLBACK_APP', 'rm -f "$cli"'):
        check(needle in INSTALLER, 'installer contains %r' % needle)
    for body in re.findall(r"<<'EOF'\n(.*?)\nEOF\n", INSTALLER, re.S):
        ok, err = bash_ok(body)
        check(ok, 'embedded CLI script syntax: %s' % err)

    # dry run on a simulated Mac
    pkg = TMP / 'pkg'
    (pkg / 'native' / 'macos').mkdir(parents=True)
    shutil.copy(MAC / 'BlueMolecule.swift', pkg / 'native' / 'macos')
    shutil.copy(MAC / 'Info.plist', pkg / 'native' / 'macos')
    (pkg / 'Install BlueMolecule.command').write_text(INSTALLER.replace('VERSION="@VERSION@"', 'VERSION="%s"' % VERSION, 1))
    for f, data in (('BlueMolecule.html', '<html></html>'), ('AppIcon.icns', 'icns'), ('bluemolecule.py', '#!/usr/bin/env python3\n')):
        (pkg / f).write_text(data)
    fb = pkg / 'fallback' / 'BlueMolecule.app' / 'Contents'
    (fb / 'MacOS').mkdir(parents=True)
    (fb / 'Resources').mkdir()
    (fb / 'MacOS' / 'BlueMolecule').write_text('#!/bin/bash\n# browser-window fallback\n')
    (fb / 'Resources' / 'bluemolecule.py').write_text('#!/usr/bin/env python3\n')
    fake = TMP / 'fakebin'
    fake.mkdir()
    tools = {
        'uname': 'case "$1" in -s) echo Darwin;; -m) echo arm64;; *) echo Darwin;; esac',
        'xcode-select': 'case "$1" in -p) [ -n "$FAKE_NO_CLT" ] && exit 2; echo /Library/Developer/CommandLineTools;; '
                        '--install) echo install >> "$FAKE_LOG";; esac',
        'xcrun': 'if [ "$1" = --find ]; then echo /fake/swiftc; exit 0; fi\n'
                 'if [ "$1" = --sdk ]; then echo /fake/MacOSX.sdk; exit 0; fi\n'
                 'shift; if [ "$1" = --version ]; then echo "Fake Swift version 6.0"; exit 0; fi\n'
                 'echo "swiftc $*" >> "$FAKE_LOG"\n'
                 'if [ -n "$FAKE_SWIFT_FAIL" ]; then echo "BlueMolecule.swift:42:7: error: cannot find foo in scope"; exit 1; fi\n'
                 'out=""; while [ $# -gt 0 ]; do [ "$1" = -o ] && out="$2"; shift; done; printf "MACHO" > "$out"',
        'lipo': 'out=""; ins=""; while [ $# -gt 0 ]; do case "$1" in -output) out="$2"; shift;; -create) ;; *) ins="$ins $1";; esac; shift; done; '
                'cat $ins > "$out"; echo "lipo$ins" >> "$FAKE_LOG"',
        'plutil': 'exec python3 -c "import plistlib,sys; plistlib.load(open(sys.argv[1],\'rb\'))" "$2"',
        'codesign': 'echo "codesign $*" >> "$FAKE_LOG"',
        'xattr': 'exit 0', 'pgrep': 'exit 1', 'osascript': 'exit 0',
        'ditto': 'cp -R "$1" "$2"',
        'open': 'echo "open $*" >> "$FAKE_LOG"',
    }
    for name, body in tools.items():
        (fake / name).write_text('#!/bin/bash\n' + body + '\n')
        (fake / name).chmod(0o755)

    def dry_run(label, **extra):
        home = TMP / ('home-' + label)
        home.mkdir()
        env = dict(os.environ, HOME=str(home), PATH=str(fake) + os.pathsep + os.environ['PATH'], FAKE_LOG=str(home / 'fake.log'),
                   BLUEMOLECULE_INSTALL_DIR=str(home / 'Applications'), BLUEMOLECULE_BIN_DIR=str(home / 'bin'), TMPDIR=str(TMP))
        env.update(extra)
        r = subprocess.run(['bash', str(pkg / 'Install BlueMolecule.command')], stdin=subprocess.DEVNULL, capture_output=True,
                           text=True, env=env, cwd=str(TMP), timeout=60)
        log = (home / 'fake.log').read_text() if (home / 'fake.log').exists() else ''
        return r, home, home / 'Applications' / 'BlueMolecule.app', log

    r, home, app, log = dry_run('native')
    check(r.returncode == 0, 'native dry run exits 0: %s' % (r.stdout + r.stderr)[-600:])
    check((app / 'Contents' / 'MacOS' / 'BlueMolecule').read_text() == 'MACHOMACHO', 'universal binary assembled (lipo of arm64 + x86_64)')
    info = plistlib.loads((app / 'Contents' / 'Info.plist').read_bytes())
    check(info['CFBundleVersion'] == VERSION and (app / 'Contents' / 'PkgInfo').read_text() == 'APPL????', 'Info.plist version and PkgInfo')
    check(all((app / 'Contents' / 'Resources' / f).exists() for f in ('AppIcon.icns', 'BlueMolecule.html', 'bluemolecule.py')), 'resources copied')
    check('-target arm64-apple-macos11.0' in log and '-target x86_64-apple-macos11.0' in log and '-swift-version 5' in log, 'compiled for both architectures')
    check('codesign --force --deep -s - ' in log and log.strip().endswith('open %s' % app), 'signed ad hoc, then opened')
    cli = home / 'bin' / 'bluemolecule'
    check(cli.exists() and os.access(cli, os.X_OK) and 'exec open -a "$APP"' in cli.read_text() and str(app) in cli.read_text(), 'native CLI installed')
    xyz = TMP / 'my water.xyz'
    xyz.write_text('1\n\nO 0 0 0\n')
    rr = subprocess.run([str(cli), 'my water.xyz', 'missing.log'], cwd=str(TMP), capture_output=True, text=True,
                        env=dict(os.environ, PATH=str(fake) + os.pathsep + os.environ['PATH'], FAKE_LOG=str(home / 'cli.log')))
    cli_log = (home / 'cli.log').read_text().strip()
    check(cli_log == 'open -a %s %s' % (app, xyz) and 'missing.log: no such file' in rr.stderr, 'CLI runs open -a with absolute paths: %r' % cli_log)
    rr = subprocess.run([str(cli), '--version'], capture_output=True, text=True)
    check(rr.stdout.strip() == 'BlueMolecule ' + VERSION, 'CLI --version: %r' % rr.stdout)
    # idempotent: a second run over the same installation
    env2 = dict(os.environ, HOME=str(home), PATH=str(fake) + os.pathsep + os.environ['PATH'], FAKE_LOG=str(home / 'fake2.log'),
                BLUEMOLECULE_INSTALL_DIR=str(home / 'Applications'), BLUEMOLECULE_BIN_DIR=str(home / 'bin'), TMPDIR=str(TMP))
    r2 = subprocess.run(['bash', str(pkg / 'Install BlueMolecule.command')], stdin=subprocess.DEVNULL, capture_output=True, text=True, env=env2, timeout=60)
    check(r2.returncode == 0 and (app / 'Contents' / 'MacOS' / 'BlueMolecule').exists() and not (app / 'BlueMolecule.app').exists(),
          'second run replaces the installation cleanly')

    r, home, app, log = dry_run('swifterror', FAKE_SWIFT_FAIL='1')
    out = r.stdout + r.stderr
    check(r.returncode == 0 and 'cannot find foo in scope' in out, 'compiler error is shown: %s' % out[-400:])
    check((app / 'Contents' / 'MacOS' / 'BlueMolecule').read_text().startswith('#!/bin/bash'), 'fallback app installed after a compiler error')
    check('python3' in (home / 'bin' / 'bluemolecule').read_text() and str(app) in (home / 'bin' / 'bluemolecule').read_text(), 'fallback CLI uses the Python launcher')

    r, home, app, log = dry_run('noclt', FAKE_NO_CLT='1')
    out = r.stdout + r.stderr
    check('install' in log and 'swiftc' not in log, 'missing Command Line Tools: xcode-select --install, no compile')
    check('double-click' in out and (app / 'Contents' / 'MacOS' / 'BlueMolecule').read_text().startswith('#!/bin/bash'),
          'missing tools: explains the next step and installs the fallback meanwhile')

    # relative folders (with spaces) in the environment: the terminal command must hold absolute paths
    r, home, _, log = dry_run('relative', BLUEMOLECULE_INSTALL_DIR='rel apps', BLUEMOLECULE_BIN_DIR='rel bin')
    rel_app, rel_cli = TMP / 'rel apps' / 'BlueMolecule.app', TMP / 'rel bin' / 'bluemolecule'
    check(r.returncode == 0 and (rel_app / 'Contents' / 'MacOS' / 'BlueMolecule').exists() and rel_cli.exists(),
          'relative install folders work: %s' % (r.stdout + r.stderr)[-300:])
    check(rel_cli.exists() and ('APP="%s"' % rel_app) in rel_cli.read_text(), 'CLI stores the absolute app path for a relative install folder')

    # ---------------- Linux ----------------
    print('Linux packaging')
    desktop = (LIN / 'org.bluemolecule.BlueMolecule.desktop').read_text()
    keys = dict(l.split('=', 1) for l in desktop.splitlines() if '=' in l)
    check(keys.get('StartupWMClass') == 'org.bluemolecule.BlueMolecule' and keys.get('Icon') == 'org.bluemolecule.BlueMolecule', 'desktop StartupWMClass/Icon')
    check(keys.get('Exec', '').endswith('%F') and keys.get('Type') == 'Application', 'desktop Exec takes files')
    ns = {'m': 'http://www.freedesktop.org/standards/shared-mime-info'}
    root = ET.parse(str(LIN / 'org.bluemolecule.BlueMolecule.xml')).getroot()
    types = {t.get('type'): [g.get('pattern') for g in t.findall('m:glob', ns)] for t in root.findall('m:mime-type', ns)}
    globs = {g[2:] for gl in types.values() for g in gl}
    check(EXTS <= globs and not (globs & GENERIC), 'MIME globs cover the chemistry extensions only: %s' % sorted(globs))
    check(globs == exts, 'Linux MIME globs and macOS UTIs list the same extensions')
    check(set(filter(None, keys.get('MimeType', '').split(';'))) == set(types), 'desktop MimeType matches the MIME package')
    for f in ('install.sh', 'uninstall.sh'):
        ok, err = bash_ok(LIN / f)
        check(ok, '%s syntax: %s' % (f, err))
    check('StartupWMClass=org.bluemolecule.BlueMolecule' in desktop and 'Icon=org.bluemolecule.BlueMolecule' in desktop, 'desktop file names the window class and icon')
    QT_SHELL = (LIN / 'bluemolecule_qt.py').read_text()
    check("setDesktopFileName(APP_ID)" in QT_SHELL and "APP_ID = 'org.bluemolecule.BlueMolecule'" in QT_SHELL, 'Qt window uses the desktop-file name as its app_id')
    try:
        compile(QT_SHELL, 'bluemolecule_qt.py', 'exec')
        check(True, '')
    except SyntaxError as e:
        check(False, 'bluemolecule_qt.py compiles: %s' % e)
    ISIZES = (48, 64, 128, 256, 512)
    lpkg = TMP / 'linuxpkg'
    (lpkg / 'icons').mkdir(parents=True)
    for f in ('install.sh', 'uninstall.sh', 'org.bluemolecule.BlueMolecule.desktop', 'org.bluemolecule.BlueMolecule.xml', 'bluemolecule_qt.py'):
        shutil.copy(LIN / f, lpkg / f)
    (lpkg / 'BlueMolecule.html').write_text('<html></html>')
    NO_QT = '#!/bin/sh\nif [ "$1" = --check-qt ]; then echo "sudo apt install python3-pyqt6.qtwebengine"; exit 1; fi\n'
    (lpkg / 'bluemolecule').write_text(NO_QT)            # a launcher stand-in: Qt WebEngine is missing
    for s in ISIZES:
        (lpkg / 'icons' / ('org.bluemolecule.BlueMolecule-%d.png' % s)).write_bytes(b'png')
    lbin, tools_log = TMP / 'lbin', TMP / 'tools.log'
    lbin.mkdir()
    for tool in ('kbuildsycoca6', 'update-desktop-database', 'sudo'):    # record calls; sudo must never run unasked
        (lbin / tool).write_text('#!/bin/sh\necho "%s $*" >> "%s"\n' % (tool, tools_log))
        (lbin / tool).chmod(0o755)
    prefix = TMP / 'prefix with space'
    share = prefix / 'share'
    env = dict(os.environ, PREFIX=str(prefix), HOME=str(TMP / 'lhome'), XDG_DATA_HOME=str(TMP / 'lhome' / 'data'),
               PATH=str(lbin) + os.pathsep + os.environ.get('PATH', ''))
    env.pop('BLUEMOLECULE_NO_PROMPT', None)
    old_files = [share / 'applications' / 'bluemolecule.desktop', share / 'icons' / 'hicolor' / '48x48' / 'apps' / 'bluemolecule.png',
                 share / 'icons' / 'hicolor' / '32x32' / 'apps' / 'bluemolecule.png', share / 'icons' / 'hicolor' / '256x256' / 'apps' / 'bluemolecule.png',
                 share / 'icons' / 'hicolor' / 'scalable' / 'apps' / 'bluemolecule.svg', share / 'pixmaps' / 'bluemolecule.png',
                 share / 'mime' / 'packages' / 'bluemolecule.xml', share / 'BlueMolecule' / 'BlueMolecule.html', share / 'bluemolecule' / 'bluemolecule.py']
    for p in old_files:                                   # what BlueMolecule 1.0/1.1 left behind
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text('old')
    (prefix / 'bin').mkdir(parents=True)
    (prefix / 'bin' / 'bluemolecule').write_text('#!/bin/sh\n# 1.0 copied the launcher here\n')
    r = subprocess.run(['bash', str(lpkg / 'install.sh')], capture_output=True, text=True, env=env, stdin=subprocess.DEVNULL)
    inst = share / 'applications' / 'org.bluemolecule.BlueMolecule.desktop'
    data = share / 'bluemolecule'
    check(r.returncode == 0 and inst.exists(), 'install.sh into a temporary prefix: %s' % r.stderr[-300:])
    check(('Exec="%s/bin/bluemolecule" %%F' % prefix) in inst.read_text() and 'StartupWMClass=org.bluemolecule.BlueMolecule' in inst.read_text(),
          'installed desktop file runs the absolute launcher path')
    check(all((share / 'icons' / 'hicolor' / ('%dx%d' % (s, s)) / 'apps' / 'org.bluemolecule.BlueMolecule.png').exists() for s in ISIZES), 'hicolor icons 48-512')
    check((share / 'mime' / 'packages' / 'org.bluemolecule.BlueMolecule.xml').exists(), 'MIME package installed')
    check(os.path.islink(prefix / 'bin' / 'bluemolecule') and os.readlink(str(prefix / 'bin' / 'bluemolecule')) == str(data / 'bluemolecule')
          and (data / 'BlueMolecule.html').exists() and os.access(str(data / 'bluemolecule_qt.py'), os.X_OK), 'launcher (linked), page and Qt window installed')
    check([p for p in old_files if p.exists()] == [] and not (share / 'BlueMolecule').exists(), 'BlueMolecule 1.0/1.1 leftovers removed')
    tools = tools_log.read_text() if tools_log.exists() else ''
    check('kbuildsycoca6' in tools and ('update-desktop-database %s' % (share / 'applications')) in tools, 'KDE (ksycoca) and desktop-file caches refreshed: %r' % tools)
    check('    sudo apt install python3-pyqt6.qtwebengine' in r.stdout and 'sudo' not in tools, 'Qt WebEngine hint shown, sudo not run without asking:\n' + r.stdout)
    listing = sorted(str(p.relative_to(prefix)) for p in prefix.rglob('*'))
    (data / 'browser').mkdir()                            # a browser profile with the user's settings
    (data / 'qtwebengine').mkdir()                        # the Qt window's profile
    (lpkg / 'bluemolecule').write_text('#!/bin/sh\n[ "$1" = --check-qt ] && echo "PyQt6 /usr/bin/python3"\n')
    r = subprocess.run(['bash', str(lpkg / 'install.sh')], capture_output=True, text=True, env=env, stdin=subprocess.DEVNULL)
    again = sorted(str(p.relative_to(prefix)) for p in prefix.rglob('*') if p.name not in ('browser', 'qtwebengine'))
    check(r.returncode == 0 and again == listing and (data / 'browser').is_dir(), 'install.sh again: same files, settings kept')
    check('Window: native Qt WebEngine window (PyQt6 /usr/bin/python3)' in r.stdout, 'install.sh reports the native window when Qt is there: ' + r.stdout[-300:])
    tools_log.write_text('')
    r = subprocess.run(['bash', str(lpkg / 'uninstall.sh')], capture_output=True, text=True, env=env)
    check(r.returncode == 0 and not inst.exists() and not (prefix / 'bin' / 'bluemolecule').exists() and not (data / 'bluemolecule_qt.py').exists()
          and not any((share / 'icons' / 'hicolor' / ('%dx%d' % (s, s)) / 'apps' / 'org.bluemolecule.BlueMolecule.png').exists() for s in ISIZES),
          'uninstall.sh removes the files')
    check((data / 'browser').is_dir() and (data / 'qtwebengine').is_dir(), 'uninstall.sh keeps the profiles (settings)')
    check('kbuildsycoca6' in tools_log.read_text(), 'uninstall.sh refreshes the KDE cache')
    r = subprocess.run(['bash', str(lpkg / 'uninstall.sh')], capture_output=True, text=True, env=env)
    check(r.returncode == 0, 'uninstall.sh again is harmless: %s' % r.stderr[-200:])

    # ---------------- built packages (each when present: a build may make only one of them) ----------------
    zpath = ROOT / 'dist' / ('BlueMolecule-%s-macOS.zip' % VERSION)
    tpath = ROOT / 'dist' / ('BlueMolecule-%s-linux.tar.gz' % VERSION)
    if zpath.exists() or tpath.exists():
        print('Built packages')
    if zpath.exists():
        base = 'BlueMolecule-%s-macOS/' % VERSION
        with zipfile.ZipFile(zpath) as zf:
            names = set(zf.namelist())
            for f in ('Install BlueMolecule.command', 'README.txt', 'BlueMolecule.html', 'AppIcon.icns', 'bluemolecule.py',
                      'native/macos/BlueMolecule.swift', 'native/macos/Info.plist', 'fallback/BlueMolecule.app/Contents/MacOS/BlueMolecule',
                      'fallback/BlueMolecule.app/Contents/Info.plist', 'fallback/BlueMolecule.app/Contents/Resources/BlueMolecule.html'):
                check(base + f in names, 'zip has %s' % f)
            info = zf.getinfo(base + 'Install BlueMolecule.command')
            check((info.external_attr >> 16) & 0o111, 'installer is executable in the zip')
            check(('VERSION="%s"' % VERSION) in zf.read(base + 'Install BlueMolecule.command').decode(), 'installer carries the version')
            html = zf.read(base + 'BlueMolecule.html').decode()
            check('BM:native-bridge' in html and 'rel="apple-touch-icon"' in html, 'packaged page has the native bridge and icons')
    if tpath.exists():
        lbase = 'BlueMolecule-%s-linux/' % VERSION
        with tarfile.open(tpath) as tf:
            names = set(tf.getnames())
            for f in ('bluemolecule', 'bluemolecule_qt.py', 'install.sh', 'uninstall.sh', 'README.txt', 'BlueMolecule.html',
                      'org.bluemolecule.BlueMolecule.desktop', 'org.bluemolecule.BlueMolecule.xml') + tuple(
                          'icons/org.bluemolecule.BlueMolecule-%d.png' % s for s in ISIZES):
                check(lbase + f in names, 'tar.gz has %s' % f)
            check(tf.getmember(lbase + 'install.sh').mode & 0o111 and tf.getmember(lbase + 'bluemolecule').mode & 0o111, 'scripts executable in the tar.gz')
            check(tf.extractfile(lbase + 'bluemolecule_qt.py').read().decode() == QT_SHELL, 'tar.gz has the current Qt window')
            check('bluemolecule --diagnose' in tf.extractfile(lbase + 'README.txt').read().decode(), 'README explains --diagnose')
            html = tf.extractfile(lbase + 'BlueMolecule.html').read().decode()
            check('BM:native-bridge' in html and "'linux-app'" in html, 'packaged page has the native bridge with the Linux app')
    if not (zpath.exists() or tpath.exists()):
        print('Built packages: dist/ not built yet, skipped')
finally:
    shutil.rmtree(TMP, ignore_errors=True)

print('\nNative packaging tests: %d passed, %d failed' % (passes, failures))
sys.exit(1 if failures else 0)
