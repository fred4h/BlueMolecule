#!/usr/bin/env python3
"""Checks of the Windows package that run without Windows (there is no Windows machine in the test setup):

  - the package built by tools/build.py --target windows: file list, ASCII and CRLF, version numbers, the icon;
  - native/windows/BlueMolecule.cs: balanced delimiters and only C# 5 syntax (the .NET Framework's csc.exe);
  - native/windows/*.ps1: balanced delimiters, no PowerShell 7-only syntax, no downloads or elevation, registry
    writes only under HKEY_CURRENT_USER, and uninstall.ps1 removing exactly what install.ps1 adds.

The page contract of the launcher is tested in a browser by tests/browser_windows_launch_test.js.
Run:  python3 tests/test_windows_package.py
"""
import re
import struct
import subprocess
import sys
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
WIN = ROOT / 'native' / 'windows'
VERSION = re.search(r"^VERSION = '([\d.]+)'", (ROOT / 'tools' / 'build.py').read_text(), re.M).group(1)
passes = fails = 0


def check(cond, msg):
    global passes, fails
    if cond:
        passes += 1
    else:
        fails += 1
        print('  FAIL:', msg)


# ------------------------------------------------------------------ lexers: blank out strings and comments
def ps_code(text):
    """PowerShell with strings, comments and ${...} variable names blanked out."""
    out, i, n = [], 0, len(text)
    while i < n:
        c = text[i]
        if text.startswith('<#', i):
            j = text.find('#>', i + 2)
            i = n if j < 0 else j + 2
        elif c == '#':
            j = text.find('\n', i)
            i = n if j < 0 else j
        elif c == "'":
            j = i + 1
            while j < n and not (text[j] == "'" and text[j + 1:j + 2] != "'"):
                j += 2 if text[j] == "'" else 1
            out.append("''")
            i = j + 1
        elif c == '"':
            j = i + 1
            while j < n and text[j] != '"':
                if text.startswith('$(', j):
                    raise ValueError('a $( ) subexpression inside a double-quoted string (line %d)' % (text.count('\n', 0, j) + 1))
                j += 2 if text[j] == '`' else 1
            out.append('""')
            i = j + 1
        elif text.startswith('${', i):
            i = text.find('}', i) + 1
            out.append('$v')
        else:
            out.append(c)
            i += 1
    return ''.join(out)


def cs_code(text):
    """C# with strings, chars and comments blanked out."""
    out, i, n = [], 0, len(text)
    while i < n:
        if text.startswith('//', i):
            j = text.find('\n', i)
            i = n if j < 0 else j
        elif text.startswith('/*', i):
            i = text.find('*/', i + 2) + 2
        elif text.startswith('@"', i):
            j = i + 2
            while j < n and not (text[j] == '"' and text[j + 1:j + 2] != '"'):
                j += 2 if text[j] == '"' else 1
            out.append('""')
            i = j + 1
        elif text[i] in '"\'':
            q, j = text[i], i + 1
            while j < n and text[j] != q:
                if text[j] == '\n':
                    raise ValueError('line break inside a literal (line %d)' % (text.count('\n', 0, j) + 1))
                j += 2 if text[j] == '\\' else 1
            out.append(q + q)
            i = j + 1
        else:
            out.append(text[i])
            i += 1
    return ''.join(out)


def balance(code):
    """None when (), {} and [] are balanced, else the line of the first problem."""
    stack = []
    for k, ch in enumerate(code):
        if ch in '({[':
            stack.append((')}]'['({['.index(ch)], k))
        elif ch in ')}]':
            if not stack or stack[-1][0] != ch:
                return code.count('\n', 0, k) + 1
            stack.pop()
    return code.count('\n', 0, stack[-1][1]) + 1 if stack else None


def assignments(code, names):
    """Values of the top-level `$Name = ...` assignments (continuation lines of a list included)."""
    out = {}
    for name in names:
        m = re.search(r'^\$' + name + r' = (.+?)(?=\n\$|\n\n|\nfunction )', code, re.M | re.S)
        lines = [ln for ln in m.group(1).split('\n') if not ln.strip().startswith('#')] if m else None
        out[name] = re.sub(r'\s+', ' ', ' '.join(lines)).strip() if m else None
    return out


# ------------------------------------------------------------------ the C# launcher
cs = (WIN / 'BlueMolecule.cs').read_text(encoding='utf-8')
check(all(ord(c) < 127 for c in cs), 'BlueMolecule.cs is ASCII')
code = cs_code(cs)
check(balance(code) is None, 'BlueMolecule.cs: balanced delimiters (problem at line %s)' % balance(code))
for pattern, what in [(r'\$"', 'string interpolation'), (r'\?\.', 'null-conditional ?.'), (r'=>', 'lambdas or expression bodies'),
                      (r'\bnameof\b', 'nameof'), (r'\bwhen\s*\(', 'exception filters'), (r'\bout\s+var\b', 'out var'),
                      (r'\?\?=', '??='), (r'\bawait\b|\basync\b', 'async/await'), (r'\busing static\b', 'using static')]:
    check(not re.search(pattern, code), 'BlueMolecule.cs uses no %s (C# 5 only)' % what)
defined = set(re.findall(r'\bstatic\s+[\w<>\[\]]+\s+(\w+)\s*\(', code))
for name in ('Main', 'Show', 'Log', 'FindBrowser', 'OpenWindow', 'Quote', 'WriteLaunchPage', 'JsonString', 'RemoveOldPages', 'Diagnostics'):
    check(name in defined, 'BlueMolecule.cs defines %s()' % name)
called = set(re.findall(r'(?<![\w.])([A-Z]\w+)\s*\(', code)) - {'Main'}
local = {c for c in called if c in defined}
check(local == defined - {'Main'}, 'every helper of BlueMolecule.cs is used: unused %s' % sorted(defined - {'Main'} - local))
check('@VERSION@' in cs and '@VERSION4@' in cs, 'BlueMolecule.cs has the version placeholders')
check('System.Windows.Forms' in cs and 'Microsoft.Win32' in cs, 'BlueMolecule.cs: namespaces for the message box and the registry')
check('UseShellExecute = false' in cs and '--user-data-dir=' in cs and '--app=' in cs, 'BlueMolecule.cs opens an app window with its own profile')
check('HKEY_LOCAL_MACHINE' not in cs and 'SetValue(' not in cs and 'CreateSubKey' not in cs, 'BlueMolecule.cs only reads the registry')

# ------------------------------------------------------------------ the PowerShell scripts
scripts = {}
for name in ('install.ps1', 'uninstall.ps1'):
    text = (WIN / name).read_text(encoding='utf-8')
    check(all(ord(c) < 127 for c in text), name + ' is ASCII (Windows PowerShell 5.1 reads it as ANSI)')
    try:
        code = ps_code(text)
    except ValueError as e:
        check(False, name + ': ' + str(e))
        continue
    scripts[name] = code
    check(balance(code) is None, name + ': balanced delimiters (problem at line %s)' % balance(code))
    check('Set-StrictMode -Version 2' in code, name + ' runs in strict mode')
    for pattern, what in [(r'&&|\|\|', 'pipeline chains && ||'), (r'\?\?|\?\.', 'null operators ?? ?.'), (r'::new\(', '::new()'),
                          (r'-AsHashtable|-Parallel\b|\$PSStyle', 'PowerShell 7 parameters'),
                          (r'Set-ExecutionPolicy|-Verb\s+RunAs|Start-Process[^\n]*RunAs', 'changing the execution policy or elevating'),
                          (r'Invoke-Expression|\biex\b|Invoke-WebRequest|Invoke-RestMethod|DownloadString|DownloadFile|Start-BitsTransfer',
                           'downloads or evaluated strings')]:
        check(not re.search(pattern, code, re.I), name + ' uses no ' + what)

inst, unin = scripts.get('install.ps1', ''), scripts.get('uninstall.ps1', '')
inst_text = (WIN / 'install.ps1').read_text(encoding='utf-8')
unin_text = (WIN / 'uninstall.ps1').read_text(encoding='utf-8')
names = ('Dest', 'Shortcut', 'Classes', 'ProgId', 'AppKey', 'UninstallKey', 'Extensions')
a, b = assignments(inst_text, names), assignments(unin_text, names)
for n in names:
    check(a[n] is not None and a[n] == b[n], 'install.ps1 and uninstall.ps1 agree on $%s: %s | %s' % (n, a[n], b[n]))
for n in ('Classes', 'AppKey', 'UninstallKey'):
    check((a[n] or '').startswith("'HKCU:") or (a[n] or '').startswith('$Classes'), '$%s is under HKEY_CURRENT_USER: %s' % (n, a[n]))
exts = re.findall(r"'(\w+)'", a['Extensions'] or '')
check(len(exts) >= 30 and len(set(exts)) == len(exts) and 'com' not in exts and 'exe' not in exts,
      'extension list: %d unique, without .com/.exe' % len(exts))
# HKLM only inside Find-Browser, and only read there
fb = re.search(r'function Find-Browser \{.*?\n\}', inst_text, re.S)
check(fb is not None and 'HKLM' not in inst_text.replace(fb.group(0), ''), 'install.ps1 mentions HKLM only in Find-Browser (to read App Paths)')
for writer in ('New-Item -Path $Key', 'Set-Item -LiteralPath $Key', 'New-ItemProperty -LiteralPath $Key'):
    check(writer in inst_text, 'install.ps1 writes the registry only through its helpers (' + writer + ')')
check(not re.search(r'(New-Item|Set-Item|New-ItemProperty|Set-ItemProperty|Remove-Item)[^\n]*HK(LM|CR|U):', inst_text),
      'install.ps1 never writes a literal registry root (only the HKCU variables)')
funcs = re.findall(r'^function ([\w-]+)', inst_text, re.M)
for f in funcs:
    check(len(re.findall(r'(?<![\w-])' + re.escape(f) + r'(?![\w-])', inst)) >= 2, 'install.ps1 uses its function ' + f)
for f in ('Say', 'New-KeyIfMissing', 'Set-DefaultValue', 'Set-StringValue', 'Set-DwordValue', 'Find-Browser', 'Invoke-Compiler'):
    check(f in funcs, 'install.ps1 defines ' + f)
# what install adds, uninstall removes
for piece in ["'\\.' + $e + '\\OpenWithProgids'", '$AppKey', '$UninstallKey', "$Classes + '\\' + $ProgId", '$Shortcut', '$Dest']:
    check(piece in inst_text and piece in unin_text, 'install adds and uninstall removes ' + piece)
check('Remove-ItemProperty -LiteralPath ($Classes + \'\\.\' + $e + \'\\OpenWithProgids\') -Name $ProgId' in unin_text,
      'uninstall.ps1 removes only BlueMolecule\'s own Open with entry for each file type')
check(re.search(r"'/target:winexe'", inst_text) and "'/r:System.Windows.Forms.dll'" in inst_text and "'/win32icon:'" in inst_text,
      'install.ps1 compiles a windowless exe with the icon')
check('Invoke-Compiler $csc ($common + @($source))' in inst_text, 'install.ps1 retries the compile without the icon')
check(re.search(r"Framework64', 'Framework'", inst_text) and 'v4.0.30319\\csc.exe' in inst_text, 'install.ps1 looks for the .NET Framework compiler')
for f in ('BlueMolecule.html', 'LICENSE', 'README.txt'):
    check("'" + f + "'" in inst_text, 'install.ps1 copies ' + f + ' from the package root')
for f in ('BlueMolecule.cs', 'BlueMolecule.ico', 'uninstall.ps1'):
    check("'" + f + "'" in inst_text, 'install.ps1 copies ' + f + ' from windows/')

cmd = (WIN / 'Install BlueMolecule.cmd').read_text(encoding='utf-8')
check(all(ord(c) < 127 for c in cmd), 'Install BlueMolecule.cmd is ASCII')
check('powershell.exe -NoProfile -ExecutionPolicy Bypass -File "%~dp0windows\\install.ps1"' in cmd, 'the .cmd runs install.ps1 for this run only (-ExecutionPolicy Bypass)')
check(cmd.count('pause') >= 2, 'the .cmd keeps its window open at the end')
block = re.search(r'if not exist [^\n]*\(\n(.*?)\n\)', cmd, re.S)
check(block is not None and not re.search(r'(?<!\^)[()]', block.group(1)), 'parentheses inside the .cmd if-block are escaped')

# ------------------------------------------------------------------ the icon
ico = (ROOT / 'native' / 'icons' / 'BlueMolecule.ico').read_bytes()
res, typ, count = struct.unpack_from('<HHH', ico, 0)
check(res == 0 and typ == 1 and count == 7, 'BlueMolecule.ico header: %d images' % count)
end, sizes = 6 + 16 * count, []
for k in range(count):
    w, h, colors, _, planes, bpp, size, off = struct.unpack_from('<BBBBHHII', ico, 6 + 16 * k)
    s = w or 256
    sizes.append(s)
    hdr = struct.unpack_from('<IiiHH', ico, off)
    mask = ((s + 31) // 32) * 4 * s
    check(off == end and hdr == (40, s, 2 * s, 1, 32) and size == 40 + s * s * 4 + mask and bpp == 32 and planes == 1,
          'icon image %d px: 32-bit BMP with its AND mask, in place' % s)
    end = off + size
check(end == len(ico) and sizes == [16, 24, 32, 48, 64, 128, 256], 'icon sizes %s fill the file' % sizes)

# ------------------------------------------------------------------ the built package
r = subprocess.run([sys.executable, str(ROOT / 'tools' / 'build.py'), '--target', 'windows'], stdout=subprocess.PIPE,
                   stderr=subprocess.STDOUT, universal_newlines=True)
check(r.returncode == 0, 'tools/build.py --target windows: ' + r.stdout[-300:])
name = 'BlueMolecule-%s-windows' % VERSION
zpath = ROOT / 'dist' / (name + '.zip')
if zpath.exists():
    with zipfile.ZipFile(zpath) as zf:
        files = sorted(n for n in zf.namelist() if not n.endswith('/'))
        expected = sorted(name + '/' + f for f in ('BlueMolecule.html', 'Install BlueMolecule.cmd', 'LICENSE', 'README.txt',
                                                   'windows/BlueMolecule.cs', 'windows/BlueMolecule.ico', 'windows/install.ps1', 'windows/uninstall.ps1'))
        check(files == expected, 'package contents: %s' % files)
        for f in ('Install BlueMolecule.cmd', 'README.txt', 'windows/BlueMolecule.cs', 'windows/install.ps1', 'windows/uninstall.ps1'):
            data = zf.read(name + '/' + f)
            check(b'\r\n' in data and data.count(b'\n') == data.count(b'\r\n'), f + ': CRLF line ends')
            check(b'@VERSION' not in data and all(c < 127 for c in data), f + ': ASCII, version filled in')
        check(('const string Version = "%s";' % VERSION).encode() in zf.read(name + '/windows/BlueMolecule.cs'), 'launcher version ' + VERSION)
        v4 = '.'.join((VERSION.split('.') + ['0', '0', '0'])[:4])
        check(('[assembly: AssemblyVersion("%s")]' % v4).encode() in zf.read(name + '/windows/BlueMolecule.cs'), 'assembly version ' + v4)
        check(("$BMVersion = '%s'" % VERSION).encode() in zf.read(name + '/windows/install.ps1'), 'installer version ' + VERSION)
        check(zf.read(name + '/windows/BlueMolecule.ico') == ico, 'the package carries the icon')
        check(zf.read(name + '/BlueMolecule.html') == (ROOT / 'dist' / 'BlueMolecule.html').read_bytes(), 'the package carries the built page')
        readme = zf.read(name + '/README.txt').decode('ascii')
        check('experimental' in readme and 'vibe-coded' in readme and 'github.com/fred4h/BlueMolecule' in readme, 'README: experimental, vibe-coded, link')
        check('Copyright (c) 2026' in zf.read(name + '/LICENSE').decode('utf-8'), 'the package carries the license')
else:
    check(False, 'the Windows package was not built')

print('Windows package tests: %d passed, %d failed' % (passes, fails))
sys.exit(1 if fails else 0)
