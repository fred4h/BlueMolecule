#!/usr/bin/env python3
"""BlueMolecule launcher.

    bluemolecule [FILE ...] [--tier auto|qt|chromium|firefox] [--browser PATH] [--tab] [--port N] [--no-open]
    bluemolecule --diagnose      # the environment, the window kind it would use and why

Opens BlueMolecule in a window of its own, with BlueMolecule's name and icon in the taskbar. The window
kind ("tier") is chosen at start-up and logged (see --diagnose for the log file):
  1. Linux: a native Qt WebEngine window (bluemolecule_qt.py next to this file), when PyQt6, PySide6 or
     PyQt5 with QtWebEngine is found in this Python, /usr/bin/python3 or python3 on PATH; the launcher
     then runs that Python (exec) and the window identifies itself to the desktop as
     org.bluemolecule.BlueMolecule (Wayland app_id), which KDE matches to the desktop file.
  2. an app window of Chrome, Chromium, Edge, Brave or Vivaldi with BlueMolecule's window class and a
     separate browser profile (on Linux under XWayland, so that the window class is honoured);
  3. Firefox with a dedicated profile that hides tabs and toolbars;
  then the default browser. For tiers 2 and 3 BlueMolecule.html is served from a small server on
http://127.0.0.1. The server prefers a fixed port so settings saved by the viewer persist between
launches. It hands the viewer exactly the files named on the command line (nothing else on disk is
reachable), runs Open Babel conversions on request and quits soon after the window is closed.
Only the Python standard library is used (Python 3.7+).
"""
import argparse
import base64
import binascii
import http.server
import json
import os
import re
import secrets
import shlex
import shutil
import socketserver
import subprocess
import sys
import threading
import time
import urllib.parse
import webbrowser
from pathlib import Path

VERSION = '1.4.1'
APP_ID = 'org.bluemolecule.BlueMolecule'   # window class, Wayland app_id, desktop file and icon name
QT_SHELL = 'bluemolecule_qt.py'
QT_BINDINGS = ('PyQt6', 'PySide6', 'PyQt5')   # in order of preference
SYSTEM_PYTHONS = ('/usr/bin/python3', '/usr/local/bin/python3')
TIERS = ('auto', 'qt', 'chromium', 'firefox')
DEFAULT_PORT = 47591       # preferred port; the next PORT_SPAN - 1 ports are tried when it is busy
PORT_SPAN = 20
PING_TIMEOUT = 25          # s without a heartbeat from a visible window before quitting
HIDDEN_TIMEOUT = 180       # s while the window is hidden (browsers throttle hidden timers to ~1/min)
CLOSE_GRACE = 4            # s after the page reported that it closed (a reload pings again sooner)
CONNECT_TIMEOUT = 120      # s for the first connection
OBABEL_TIMEOUT = 60        # s for one Open Babel run
MAX_BODY = 256 * 1024 * 1024
MAC_BROWSERS = ['Google Chrome', 'Microsoft Edge', 'Brave Browser', 'Chromium', 'Vivaldi']
LINUX_BROWSERS = ['google-chrome', 'google-chrome-stable', 'chromium', 'chromium-browser', 'microsoft-edge',
                  'microsoft-edge-stable', 'brave-browser', 'brave', 'vivaldi', 'vivaldi-stable']
OBABEL_DIRS = ['/opt/homebrew/bin', '/usr/local/bin', '~/miniconda3/bin', '~/anaconda3/bin', '~/miniforge3/bin',
               '/opt/miniconda3/bin', '~/mambaforge/bin', '/opt/anaconda3/bin', '/opt/conda/bin']
FORMAT_RE = re.compile(r'^[A-Za-z0-9_-]{1,16}$')
# obabel options that would break the stdin -> stdout contract or write files (-O, -m)
BLOCKED_OPTIONS = ('-O', '-o', '-i', '-m', '-L', '-H', '-V')

FIREFOX_USER_JS = '''// BlueMolecule window profile, written by the bluemolecule launcher (changes are overwritten).
user_pref("toolkit.legacyUserProfileCustomizations.stylesheets", true);
user_pref("browser.shell.checkDefaultBrowser", false);
user_pref("browser.startup.homepage_override.mstone", "ignore");
user_pref("startup.homepage_welcome_url", "");
user_pref("startup.homepage_welcome_url.additional", "");
user_pref("browser.aboutwelcome.enabled", false);
user_pref("datareporting.policy.dataSubmissionPolicyBypassNotification", true);
user_pref("datareporting.policy.firstRunURL", "");
user_pref("browser.startup.page", 0);
user_pref("browser.sessionstore.resume_from_crash", false);
user_pref("browser.tabs.warnOnClose", false);
user_pref("browser.tabs.inTitlebar", 0);
user_pref("browser.toolbars.bookmarks.visibility", "never");
user_pref("browser.link.open_newwindow", 2);
user_pref("browser.download.useDownloadDir", false);
user_pref("browser.translations.automaticallyPopup", false);
'''

FIREFOX_USER_CHROME = '''/* BlueMolecule window profile: no tab strip, toolbars or sidebar (written by the bluemolecule launcher). */
#TabsToolbar, #nav-bar, #PersonalToolbar, #toolbar-menubar, #titlebar,
#sidebar-box, #sidebar-main, #sidebar-splitter { visibility: collapse !important; }
#navigator-toolbox { border-bottom: none !important; }
'''


def find_html(required=True):
    here = Path(__file__).resolve().parent
    for p in (here / 'BlueMolecule.html', here.parent / 'share' / 'bluemolecule' / 'BlueMolecule.html',
              here.parent / 'Resources' / 'BlueMolecule.html', here.parent / 'dist' / 'BlueMolecule.html'):
        if p.is_file():
            return p
    if not required:
        return None
    sys.exit('bluemolecule: BlueMolecule.html was not found next to the launcher (%s)' % here)


def find_qt_shell():
    """bluemolecule_qt.py next to the launcher (package, install) or in native/linux (source tree)."""
    here = Path(__file__).resolve().parent
    for p in (here / QT_SHELL, here.parent / 'native' / 'linux' / QT_SHELL):
        if p.is_file():
            return p
    return None


# ---------------------------------------------------------------- log

def log_path():
    """$XDG_STATE_HOME/bluemolecule/bluemolecule.log (~/Library/Logs/BlueMolecule-launcher.log on macOS)."""
    if sys.platform == 'darwin':
        return Path.home() / 'Library' / 'Logs' / 'BlueMolecule-launcher.log'
    xdg = os.environ.get('XDG_STATE_HOME')
    return (Path(xdg) if xdg else Path.home() / '.local' / 'state') / 'bluemolecule' / 'bluemolecule.log'


def log_event(text, echo=True):
    """A line on stderr and in the log file (which is kept below about 512 KB)."""
    if echo:
        print('bluemolecule: ' + text, file=sys.stderr, flush=True)
    try:
        p = log_path()
        p.parent.mkdir(parents=True, exist_ok=True)
        if p.exists() and p.stat().st_size > 512 * 1024:
            os.replace(str(p), str(p) + '.1')
        with open(str(p), 'a', encoding='utf-8') as fh:
            fh.write('%s [%d] %s\n' % (time.strftime('%Y-%m-%d %H:%M:%S'), os.getpid(), text))
    except OSError:
        pass


# ---------------------------------------------------------------- tier 1: Qt WebEngine for Python

# Run by each candidate Python: which bindings it has. 'fast' only looks for the modules (a launch);
# 'full' imports one binding (named in argv[2]) and reports its versions (--diagnose, --check-qt).
QT_PROBE = r'''
import importlib, importlib.util, json, sys
mode = sys.argv[1] if len(sys.argv) > 1 else 'fast'
names = sys.argv[2:] or ['PyQt6', 'PySide6', 'PyQt5']
out = {'python': sys.executable, 'version': '%d.%d.%d' % tuple(sys.version_info[:3]), 'bindings': {}}
for b in names:
    info = {'ok': False}
    try:
        if importlib.util.find_spec(b) is None:
            info['error'] = 'not installed'
        elif mode != 'full':
            missing = [m for m in ('QtWebEngineWidgets', 'QtWebChannel') if importlib.util.find_spec(b + '.' + m) is None]
            if missing:
                info['error'] = 'no ' + ' or '.join(missing)
            else:
                info['ok'] = True
        else:
            core = importlib.import_module(b + '.QtCore')
            importlib.import_module(b + '.QtWebEngineWidgets')
            importlib.import_module(b + '.QtWebChannel')
            info.update(ok=True, qt=core.qVersion(),
                        version=getattr(core, 'PYQT_VERSION_STR', '') or getattr(importlib.import_module(b), '__version__', ''))
            try:
                web = importlib.import_module(b + '.QtWebEngineCore')
                for key, fn in (('webengine', 'qWebEngineVersion'), ('chromium', 'qWebEngineChromiumVersion')):
                    if hasattr(web, fn):
                        info[key] = getattr(web, fn)()
            except Exception:
                pass
    except Exception as e:
        info = {'ok': False, 'error': '%s: %s' % (type(e).__name__, e)}
    out['bindings'][b] = info
print(json.dumps(out))
'''

# One-line install hints for Qt WebEngine for Python (by /etc/os-release ID and ID_LIKE).
QT_HINTS = (
    (('debian', 'ubuntu', 'linuxmint', 'pop', 'neon', 'elementary', 'zorin', 'raspbian', 'kali'),
     'sudo apt install python3-pyqt6.qtwebengine'),
    (('fedora', 'rhel', 'centos', 'rocky', 'almalinux', 'nobara', 'ultramarine'), 'sudo dnf install python3-pyqt6-webengine'),
    (('arch', 'manjaro', 'endeavouros', 'cachyos', 'garuda', 'artix'), 'sudo pacman -S --needed python-pyqt6-webengine'),
    (('opensuse', 'opensuse-tumbleweed', 'opensuse-leap', 'opensuse-slowroll', 'suse', 'sles'),
     'sudo zypper install python3-PyQt6-WebEngine'),
)


def parse_os_release(text):
    out = {}
    for line in text.splitlines():
        m = re.match(r'^([A-Z0-9_]+)=(.*)$', line.strip())
        if m:
            v = m.group(2).strip()
            if len(v) >= 2 and v[0] == v[-1] and v[0] in '"\'':
                v = v[1:-1]
            out[m.group(1)] = v
    return out


def read_os_release():
    for p in ('/etc/os-release', '/usr/lib/os-release'):
        try:
            with open(p, encoding='utf-8') as fh:
                return parse_os_release(fh.read())
        except OSError:
            continue
    return {}


def install_hint(os_release):
    """The command that installs PyQt6 with QtWebEngine on this distribution."""
    ids = [os_release.get('ID', '').lower()] + os_release.get('ID_LIKE', '').lower().split()
    for keys, cmd in QT_HINTS:
        if any(i in keys for i in ids if i):
            return cmd
    return 'python3 -m pip install --user PyQt6 PyQt6-WebEngine'


def python_candidates(current=None, which=shutil.which, isfile=os.path.isfile):
    """Pythons that may have Qt: the system Python first (its PyQt matches the desktop's Qt), then this one
    (often conda) and python3 on PATH; each real interpreter once."""
    out, seen = [], set()
    for p in (SYSTEM_PYTHONS[0], current if current is not None else sys.executable, which('python3'), SYSTEM_PYTHONS[1]):
        if not p or not isfile(p):
            continue
        key = os.path.realpath(p)
        if key not in seen:
            seen.add(key)
            out.append(p)
    return out


def parse_probe(python, returncode, stdout, stderr=''):
    lines = [l for l in (stdout or '').strip().splitlines() if l.startswith('{')]
    try:
        info = json.loads(lines[-1])
        if not isinstance(info.get('bindings'), dict):
            raise ValueError('no bindings')
    except (IndexError, ValueError, AttributeError):
        err = (stderr or '').strip().splitlines()
        return {'python': python, 'bindings': {},
                'error': 'probe failed (exit status %s)%s' % (returncode, ': ' + err[-1] if err else '')}
    info['python'] = python     # as found, not the realpath the interpreter reports
    return info


def probe_interpreter(python, full=False, timeout=None):
    """{'python', 'version', 'bindings': {name: {'ok', 'qt', 'version', 'error'}}} for one interpreter."""
    def run(args, t):
        try:
            r = subprocess.run([python, '-c', QT_PROBE] + args, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE, timeout=t, cwd='/', universal_newlines=True)
        except subprocess.TimeoutExpired:
            return {'python': python, 'bindings': {}, 'error': 'probe timed out after %d s' % t}
        except OSError as e:
            return {'python': python, 'bindings': {}, 'error': str(e)}
        return parse_probe(python, r.returncode, r.stdout, r.stderr)
    if not full:
        return run(['fast'], timeout or 15)
    info, errors = {'python': python, 'bindings': {}}, []
    for b in QT_BINDINGS:          # one binding per process: two Qt copies in one process can clash
        part = run(['full', b], timeout or 60)
        if part.get('version'):
            info['version'] = part['version']
        if part.get('error'):      # the probe itself failed (e.g. the Qt import crashed Python)
            errors.append(part['error'])
            info['bindings'][b] = {'ok': False, 'error': part['error']}
        else:
            info['bindings'].update(part['bindings'])
    if len(errors) == len(QT_BINDINGS) and len(set(errors)) == 1:
        info['error'] = errors[0]     # e.g. the interpreter does not run at all
    return info


def qt_skip_list(env=None):
    """Pythons whose Qt window failed to start in this launch (set by bluemolecule_qt.py when it falls back)."""
    env = os.environ if env is None else env
    return {os.path.realpath(p) for p in (env.get('BLUEMOLECULE_QT_SKIP') or '').split(os.pathsep) if p}


def pick_qt(interpreters, skip=()):
    """The first usable binding in order of preference (PyQt6, PySide6, PyQt5), across the interpreters."""
    for b in QT_BINDINGS:
        for info in interpreters:
            if os.path.realpath(info.get('python') or '') in skip:
                continue
            if (info.get('bindings') or {}).get(b, {}).get('ok'):
                return {'python': info['python'], 'binding': b, 'info': info['bindings'][b]}
    return None


def qt_command(python, binding, shell, html, files, launcher=None, launcher_python=None, forced=False):
    """argv for the Qt window (run with os.execv)."""
    cmd = [python, str(shell), '--binding', binding, '--html', str(html)]
    if launcher:
        cmd += ['--launcher', str(launcher), '--launcher-python', launcher_python or sys.executable]
    if forced:
        cmd.append('--no-fallback')
    return cmd + ['--'] + [str(f) for f in files]


# ---------------------------------------------------------------- choosing the window

def use_x11(env=None, platform=None):
    """Whether Chromium-family browsers can be sent to X11 (XWayland or an X11 session)."""
    env = os.environ if env is None else env
    platform = platform or sys.platform
    return platform != 'darwin' and not platform.startswith('win') and bool(env.get('DISPLAY'))


def choose_tier(platform, want='auto', browser=None, tab=False, qt_shell=None, qt=None, chromium=(), firefox=None, x11=False):
    """How to show BlueMolecule (pure). Returns {'tier': 'qt'|'chromium'|'firefox'|'default'|None, 'reason', ...};
    tier None means the tier asked for with --tier is not available ('error' says why).
    'qt_missing' is True when the Qt window would have been used but no Qt was found (then a hint is shown)."""
    linux = platform != 'darwin' and not platform.startswith('win')
    if browser:
        kind = browser_kind(browser)
        return {'tier': kind, 'kind': kind, 'exe': browser, 'tab': tab, 'x11': kind == 'chromium' and x11 and not tab,
                'reason': '--browser %s' % browser, 'qt_missing': False}
    if not linux:
        why = 'the Qt window is for Linux (on a Mac the BlueMolecule app is the native window)'
    elif tab:
        why = '--tab asks for a browser tab'
    elif want not in ('auto', 'qt'):
        why = '--tier %s' % want
    elif not qt_shell:
        why = '%s is missing next to the launcher' % QT_SHELL
    elif not qt:
        why = 'no Python with PyQt6, PySide6 or PyQt5 + QtWebEngine was found'
    else:
        return {'tier': 'qt', 'python': qt['python'], 'binding': qt['binding'], 'tab': False, 'qt_missing': False,
                'reason': '%s with QtWebEngine found in %s' % (qt['binding'], qt['python'])}
    missing = linux and not tab and want in ('auto', 'qt') and bool(qt_shell) and not qt   # then the hint helps
    base = {'tab': tab, 'reason': why, 'qt_missing': missing}
    if want == 'qt':
        return dict(base, tier=None, error='the Qt window is not available: ' + why)
    if want in ('auto', 'chromium') and chromium:
        return dict(base, tier='chromium', kind='chromium', exe=chromium[0], x11=x11 and not tab)
    if want == 'chromium':
        return dict(base, tier=None, error='no Chrome, Chromium, Edge, Brave or Vivaldi was found')
    if want in ('auto', 'firefox') and firefox:
        return dict(base, tier='firefox', kind='firefox', exe=firefox)
    if want == 'firefox':
        return dict(base, tier=None, error='Firefox was not found')
    return dict(base, tier='default', reason=why + '; no Chromium-family browser or Firefox was found')


def describe_plan(plan):
    """One line for the log: the tier and why."""
    t = plan.get('tier')
    if t is None:
        return 'no window: ' + plan.get('error', '?')
    if t == 'qt':
        return 'tier 1, native Qt WebEngine window: %s (Wayland app_id %s)' % (plan['reason'], APP_ID)
    name = os.path.basename(plan.get('exe') or '')
    if plan.get('tab'):
        s = 'a browser tab of %s' % (name or 'the default browser')
    elif t == 'chromium':
        s = 'tier 2, %s app window' % name + (' under XWayland (--ozone-platform=x11, window class %s)' % APP_ID
                                              if plan.get('x11') else ' (no DISPLAY: the window class may not be honoured)')
    elif t == 'firefox':
        s = 'tier 3, %s window with its own profile (app_id %s)' % (name, APP_ID)
    else:
        s = 'fallback, a tab of the default browser'
    return s + ' -- ' + plan.get('reason', '')


def plan_launch(args, full=False, env=None, platform=None):
    """choose_tier() with the facts of this computer; full=True imports Qt in each Python (slower, exact)."""
    env = os.environ if env is None else env
    platform = platform or sys.platform
    want = getattr(args, 'tier', 'auto') or 'auto'
    shell = find_qt_shell()
    linux = platform != 'darwin' and not platform.startswith('win')
    interpreters, qt = [], None
    if linux and (full or (shell and not args.tab and not args.browser and want in ('auto', 'qt'))):
        interpreters = [probe_interpreter(p, full) for p in python_candidates()]
        qt = pick_qt(interpreters, qt_skip_list(env))
    chromium, firefox = chromium_browsers(), firefox_browser()
    plan = choose_tier(platform, want, args.browser, args.tab, shell, qt, chromium, firefox, use_x11(env, platform))
    plan.update(interpreters=interpreters, qt_shell=shell, chromium=chromium, firefox=firefox)
    return plan


# ---------------------------------------------------------------- browsers

def data_dir():
    """Per-user data folder: $XDG_DATA_HOME/bluemolecule (~/Library/Application Support/BlueMolecule on macOS)."""
    xdg = os.environ.get('XDG_DATA_HOME')
    if xdg:
        return Path(xdg) / 'bluemolecule'
    if sys.platform == 'darwin':
        return Path.home() / 'Library' / 'Application Support' / 'BlueMolecule'
    return Path.home() / '.local' / 'share' / 'bluemolecule'


def snap_name(exe):
    """Name of the snap behind a browser command (snaps cannot use profiles in hidden folders), else None."""
    if exe.startswith('/snap/bin/'):
        return os.path.basename(exe)
    real = os.path.realpath(exe)
    if real.startswith('/snap/'):
        parts = real.split('/')
        return parts[2] if len(parts) > 2 else None
    try:   # Ubuntu's /usr/bin/firefox and /usr/bin/chromium-browser are small wrappers around /snap/bin/...
        if os.path.getsize(real) < 65536:
            with open(real, 'rb') as fh:
                head = fh.read(65536)
            if not head.startswith(b'\x7fELF'):
                m = re.search(rb'/snap/bin/([A-Za-z0-9-]+)', head)
                if m:
                    return m.group(1).decode('ascii')
    except OSError:
        pass
    return None


def profile_dir(kind, exe):
    """Dedicated browser profile, so BlueMolecule's window class and settings stay separate."""
    name = 'firefox' if kind == 'firefox' else 'browser'
    snap = snap_name(exe)
    if snap:
        return Path.home() / 'snap' / snap / 'common' / ('bluemolecule-' + name)
    return data_dir() / name


def browser_kind(exe):
    base = os.path.basename(exe).lower()
    return 'firefox' if ('firefox' in base or 'librewolf' in base) else 'chromium'


def chromium_browsers():
    found = []
    if sys.platform == 'darwin':
        for name in MAC_BROWSERS:
            for base in ('/Applications', os.path.expanduser('~/Applications')):
                exe = Path(base) / (name + '.app') / 'Contents' / 'MacOS' / name
                if exe.is_file():
                    found.append(str(exe))
    else:
        for name in LINUX_BROWSERS:
            p = shutil.which(name)
            if p and p not in found:
                found.append(p)
    return found


def firefox_browser():
    if sys.platform == 'darwin':
        for base in ('/Applications', os.path.expanduser('~/Applications')):
            exe = Path(base) / 'Firefox.app' / 'Contents' / 'MacOS' / 'firefox'
            if exe.is_file():
                return str(exe)
        return None
    return shutil.which('firefox') or shutil.which('firefox-esr')


def ensure_firefox_profile(profile):
    """Write the dedicated profile's user.js (enables userChrome.css) and chrome/userChrome.css."""
    chrome = Path(profile) / 'chrome'
    chrome.mkdir(parents=True, exist_ok=True)
    (Path(profile) / 'user.js').write_text(FIREFOX_USER_JS, encoding='utf-8')
    (chrome / 'userChrome.css').write_text(FIREFOX_USER_CHROME, encoding='utf-8')


def firefox_profile_running(profile):
    """True when a Firefox instance holds the profile lock (then the URL goes to that instance)."""
    lock = Path(profile) / '.parentlock'
    if not lock.exists():
        return False
    try:
        import fcntl
        fd = os.open(str(lock), os.O_RDWR)
    except (ImportError, OSError):
        return False
    try:
        fcntl.lockf(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError:
        return True
    else:
        fcntl.lockf(fd, fcntl.LOCK_UN)
        return False
    finally:
        os.close(fd)


def browser_command(kind, exe, url, profile, tab=False, running=False, platform=None, x11=False):
    """argv that opens url in BlueMolecule's own window (or a normal tab with tab=True)."""
    platform = platform or sys.platform
    if tab:
        return [exe, url]
    if kind == 'firefox':
        cmd = [exe] if running else [exe, '--new-instance']
        if platform != 'darwin':
            cmd += ['--class', APP_ID, '--name', APP_ID]
        cmd += ['--profile', str(profile)]
        return cmd + (['--new-window', url] if running else [url])
    cmd = [exe, '--app=' + url, '--class=' + APP_ID]
    if platform != 'darwin':
        cmd.append('--name=' + APP_ID)
    cmd += ['--user-data-dir=' + str(profile), '--no-first-run', '--no-default-browser-check']
    if x11:
        # Why X11 is forced: on a Wayland session Chrome, Brave & co. run as native Wayland clients, and
        # there they derive the window's app_id from the --app URL and ignore --class. KWin then finds no
        # desktop file for the window and groups it with the browser (Brave's name and icon in the taskbar).
        # Under XWayland the window's WM_CLASS is the --class value, which KDE matches against
        # StartupWMClass=org.bluemolecule.BlueMolecule in our desktop file. (The Qt window, tier 1, needs
        # none of this: it sets its Wayland app_id itself.)
        cmd.append('--ozone-platform=x11')
    return cmd


def browser_env(kind, env, platform=None, tab=False):
    """Environment for the browser: the desktop-file name for Chromium, the app_id for Firefox."""
    platform = platform or sys.platform
    env = dict(env)
    if tab:
        return env
    if kind == 'chromium':
        env['CHROME_DESKTOP'] = APP_ID + '.desktop'
    elif platform != 'darwin':
        env['MOZ_APP_REMOTINGNAME'] = APP_ID    # Firefox's Wayland app_id (and its remote-control name)
    return env


def open_window(url, browser=None, tab=False, kinds=('chromium', 'firefox')):
    """Open url in BlueMolecule's own window; returns (what was used, kind, whether X11 was forced)."""
    if browser:
        candidates = [(browser_kind(browser), browser)]
    else:
        candidates = [('chromium', exe) for exe in chromium_browsers()] if 'chromium' in kinds else []
        ff = firefox_browser() if 'firefox' in kinds else None
        if ff:
            candidates.append(('firefox', ff))
    for kind, exe in candidates:
        profile = profile_dir(kind, exe)
        running = x11 = False
        if not tab:
            try:
                if kind == 'firefox':
                    ensure_firefox_profile(profile)
                    running = firefox_profile_running(profile)
                else:
                    profile.mkdir(parents=True, exist_ok=True)
            except OSError as e:
                print('bluemolecule: cannot prepare the browser profile %s: %s' % (profile, e), file=sys.stderr)
            x11 = kind == 'chromium' and use_x11()
        try:
            subprocess.Popen(browser_command(kind, exe, url, profile, tab, running, x11=x11), env=browser_env(kind, os.environ, tab=tab),
                             stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, start_new_session=True)
            return exe, kind, x11
        except OSError:
            continue
    if sys.platform == 'darwin':
        subprocess.Popen(['open', url])
    elif shutil.which('xdg-open'):
        subprocess.Popen(['xdg-open', url], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, start_new_session=True)
    else:
        webbrowser.open(url)
    return 'default browser', 'default', False


# ---------------------------------------------------------------- Open Babel

def find_obabel():
    exe = shutil.which('obabel')
    if exe:
        return exe
    for d in OBABEL_DIRS:
        p = os.path.join(os.path.expanduser(d), 'obabel')
        if os.path.isfile(p) and os.access(p, os.X_OK):
            return p
    return None


def option_allowed(opt):
    if not isinstance(opt, str) or not opt or len(opt) > 200 or any(ord(c) < 32 for c in opt):
        return False
    return opt[:2] not in BLOCKED_OPTIONS


def obabel_request(req, exe=None, timeout=OBABEL_TIMEOUT):
    """Run Open Babel for {op: 'formats'|'convert', input (base64), inFormat, outFormat, options[]}.

    Returns {ok, output (base64), error} -- the same contract as the macOS app's bmObabel handler.
    """
    if not isinstance(req, dict):
        return {'ok': False, 'output': '', 'error': 'Bad Open Babel request.'}
    op = req.get('op')
    if op == 'formats':
        args, data = ['-L', 'formats'], b''
    elif op == 'convert':
        inf, outf = req.get('inFormat'), req.get('outFormat')
        if not (isinstance(inf, str) and FORMAT_RE.match(inf) and isinstance(outf, str) and FORMAT_RE.match(outf)):
            return {'ok': False, 'output': '', 'error': 'Invalid format code.'}
        opts = req.get('options') or []
        if not isinstance(opts, list) or len(opts) > 32:
            return {'ok': False, 'output': '', 'error': 'Invalid Open Babel options.'}
        for o in opts:
            if not option_allowed(o):
                return {'ok': False, 'output': '', 'error': 'Open Babel option not allowed: %r' % (o,)}
        try:
            data = base64.b64decode(req.get('input') or '', validate=True)
        except (binascii.Error, ValueError, TypeError):
            return {'ok': False, 'output': '', 'error': 'The input is not valid base64.'}
        args = ['-i' + inf, '-o' + outf] + opts
    else:
        return {'ok': False, 'output': '', 'error': 'Unknown Open Babel request.'}
    exe = exe or find_obabel()
    if not exe:
        return {'ok': False, 'output': '', 'error': 'Open Babel (obabel) was not found. Install it, e.g. '
                'conda install -c conda-forge openbabel, brew install open-babel or sudo apt install openbabel.'}
    try:
        p = subprocess.run([exe] + args, input=data, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=timeout)
    except subprocess.TimeoutExpired:
        return {'ok': False, 'output': '', 'error': 'Open Babel did not finish within %d s.' % timeout}
    except OSError as e:
        return {'ok': False, 'output': '', 'error': 'Could not run %s: %s' % (exe, e)}
    err = p.stderr.decode('utf-8', 'replace').strip()
    if p.returncode != 0 or not p.stdout:
        return {'ok': False, 'output': '', 'error': err or 'Open Babel failed (exit status %d).' % p.returncode}
    return {'ok': True, 'output': base64.b64encode(p.stdout).decode('ascii'), 'error': ''}


# ---------------------------------------------------------------- server

class Session:
    def __init__(self, html_path, files):
        self.token = secrets.token_urlsafe(18)
        self.html = html_path.read_text(encoding='utf-8')
        self.files = files
        # monotonic: it stops while the computer sleeps, so waking up is not mistaken for a closed window
        self.started = time.monotonic()
        self.last_ping = None
        self.hidden = False
        self.closed_at = None
        self.lock = threading.Lock()

    def page(self):
        launch = {
            'files': [{'name': p.name, 'url': '/%s/file/%d/%s' % (self.token, i, urllib.parse.quote(p.name))}
                      for i, p in enumerate(self.files)],
            'ping': '/%s/ping' % self.token,
            'obabel': '/%s/obabel' % self.token,
            'version': VERSION,
        }
        tag = '<script>window.BM_LAUNCH=%s;</script>' % json.dumps(launch).replace('</', '<\\/')
        return self.html.replace('<head>', '<head>' + tag, 1).encode('utf-8')

    def touch(self, state=''):
        """A heartbeat; state is '', 'visible', 'hidden' or 'closed' (sent by native-bridge.js)."""
        with self.lock:
            now = time.monotonic()
            if state == 'closed':
                self.closed_at = now
                return
            self.last_ping, self.closed_at = now, None
            if state in ('hidden', 'visible'):
                self.hidden = state == 'hidden'

    def idle_reason(self, now=None):
        """Why the launcher should quit now, or None to keep serving."""
        now = time.monotonic() if now is None else now
        with self.lock:
            if self.closed_at is not None and now - self.closed_at > CLOSE_GRACE:
                return 'closed'
            if self.last_ping is None:
                return 'never connected' if now - self.started > CONNECT_TIMEOUT else None
            limit = HIDDEN_TIMEOUT if self.hidden else PING_TIMEOUT
            return 'no heartbeat' if now - self.last_ping > limit else None


class Handler(http.server.BaseHTTPRequestHandler):
    server_version = 'BlueMolecule/' + VERSION

    def log_message(self, fmt, *args):  # keep the terminal quiet
        pass

    def _host_ok(self):
        host = (self.headers.get('Host') or '').split(':')[0]
        return host in ('127.0.0.1', 'localhost')

    def _origin_ok(self):
        origin = self.headers.get('Origin')
        port = self.server.server_address[1]
        return origin is None or origin in ('http://127.0.0.1:%d' % port, 'http://localhost:%d' % port)

    def _send(self, code, body=b'', ctype='text/plain; charset=utf-8'):
        self.send_response(code)
        self.send_header('Content-Type', ctype)
        self.send_header('Content-Length', str(len(body)))
        self.send_header('Cache-Control', 'no-store')
        self.send_header('X-Content-Type-Options', 'nosniff')
        self.end_headers()
        if body and self.command != 'HEAD':
            self.wfile.write(body)

    def _json(self, obj, code=200):
        return self._send(code, json.dumps(obj).encode('utf-8'), 'application/json')

    def do_GET(self):
        S = self.server.session
        if not self._host_ok():
            return self._send(403, b'forbidden')
        parts = urllib.parse.urlparse(self.path).path.strip('/').split('/')
        if not parts or parts[0] != S.token:
            return self._send(404, b'not found')
        if len(parts) == 1 or parts[1] == '':
            S.touch()
            return self._send(200, S.page(), 'text/html; charset=utf-8')
        if len(parts) >= 3 and parts[1] == 'file' and parts[2].isdigit():
            i = int(parts[2])
            if 0 <= i < len(S.files):
                try:
                    data = S.files[i].read_bytes()
                except OSError as e:
                    return self._send(500, str(e).encode())
                return self._send(200, data, 'application/octet-stream')
        return self._send(404, b'not found')

    do_HEAD = do_GET

    def do_POST(self):
        S = self.server.session
        if not self._host_ok() or not self._origin_ok():
            return self._send(403, b'forbidden')
        u = urllib.parse.urlparse(self.path)
        parts = u.path.strip('/').split('/')
        if len(parts) != 2 or parts[0] != S.token:
            return self._send(404, b'not found')
        if parts[1] == 'ping':
            S.touch(urllib.parse.parse_qs(u.query).get('state', [''])[0])
            return self._send(204)
        if parts[1] == 'obabel':
            try:
                n = int(self.headers.get('Content-Length') or 0)
            except ValueError:
                n = 0
            if n <= 0 or n > MAX_BODY:
                return self._json({'ok': False, 'output': '', 'error': 'Bad request size.'}, 413 if n > MAX_BODY else 400)
            try:
                req = json.loads(self.rfile.read(n).decode('utf-8'))
            except (ValueError, UnicodeDecodeError):
                return self._json({'ok': False, 'output': '', 'error': 'The request is not valid JSON.'}, 400)
            return self._json(obabel_request(req))
        return self._send(404, b'not found')


class Server(socketserver.ThreadingMixIn, http.server.HTTPServer):
    daemon_threads = True
    allow_reuse_address = True
    allow_reuse_port = False   # a busy preferred port must be detected, never shared


def bind_server(port):
    """Listen on 127.0.0.1: the preferred port or the next free one of its range (0 = any free port)."""
    tries = [0] if port == 0 else [p for p in range(port, port + PORT_SPAN) if 0 < p < 65536] + [0]
    err = None
    for p in tries:
        try:
            return Server(('127.0.0.1', p), Handler)
        except OSError as e:
            err = e
    raise err


def describe_interpreter(info):
    head = '%s (Python %s)' % (info.get('python'), info.get('version') or '?')
    if info.get('error') and not any(b.get('ok') for b in (info.get('bindings') or {}).values()):
        return head + ': ' + info['error']
    parts = []
    for b in QT_BINDINGS:
        d = (info.get('bindings') or {}).get(b)
        if d is None:
            continue
        if d.get('ok'):
            extra = ', '.join(x for x in ('Qt ' + d['qt'] if d.get('qt') else '',
                                          'QtWebEngine ' + d['webengine'] if d.get('webengine') else '',
                                          'Chromium ' + d['chromium'] if d.get('chromium') else '') if x)
            parts.append('%s %s OK%s' % (b, d.get('version') or '', ' (%s)' % extra if extra else ''))
        else:
            parts.append('%s: %s' % (b, d.get('error') or 'no'))
    return head + ': ' + '; '.join(parts)


def desktop_entry_report(env=None):
    """Where the desktop file is and whether its StartupWMClass matches the window class."""
    env = os.environ if env is None else env
    home = env.get('XDG_DATA_HOME') or str(Path.home() / '.local' / 'share')
    dirs = [home] + [d for d in (env.get('XDG_DATA_DIRS') or '/usr/local/share:/usr/share').split(':') if d]
    lines = []
    for d in dirs:
        p = Path(d) / 'applications' / (APP_ID + '.desktop')
        if p.is_file():
            try:
                keys = dict(l.split('=', 1) for l in p.read_text(encoding='utf-8', errors='replace').splitlines() if '=' in l)
            except OSError:
                keys = {}
            ok = keys.get('StartupWMClass') == APP_ID and keys.get('Icon') == APP_ID
            lines.append('%s (%s)' % (p, 'StartupWMClass and Icon OK' if ok else 'StartupWMClass/Icon do not match %s' % APP_ID))
        old = Path(d) / 'applications' / 'bluemolecule.desktop'
        if old.is_file():
            lines.append('%s is left over from BlueMolecule 1.1: run install.sh again or delete it' % old)
    return lines or ['not installed (run install.sh; KDE needs it to show the name and icon)']


def diagnose(args, out=None):
    """bluemolecule --diagnose: the environment, the tier that would be used, and why."""
    out = out or sys.stdout
    env = os.environ
    osr = read_os_release()
    plan = plan_launch(args, full=True)

    def row(key, value):
        print('  %-13s %s' % (key + ':', value), file=out)
    print('BlueMolecule %s diagnostics' % VERSION, file=out)
    row('system', osr.get('PRETTY_NAME') or sys.platform)
    row('desktop', '%s, %s session (WAYLAND_DISPLAY=%s, DISPLAY=%s%s)' % (
        env.get('XDG_CURRENT_DESKTOP') or '?', env.get('XDG_SESSION_TYPE') or '?', env.get('WAYLAND_DISPLAY') or '-',
        env.get('DISPLAY') or '-', ', QT_QPA_PLATFORM=' + env['QT_QPA_PLATFORM'] if env.get('QT_QPA_PLATFORM') else ''))
    row('launcher', '%s (Python %s, %s)' % (Path(__file__).resolve(), '%d.%d.%d' % sys.version_info[:3], sys.executable))
    row('page', find_html(required=False) or 'BlueMolecule.html not found')
    row('qt shell', plan['qt_shell'] or '%s not found' % QT_SHELL)
    for info in plan['interpreters']:
        row('python', describe_interpreter(info))
    skip = qt_skip_list(env)
    if skip:
        row('qt skipped', ', '.join(sorted(skip)) + ' (BLUEMOLECULE_QT_SKIP)')
    row('chromium', ', '.join(plan['chromium']) or 'none found')
    row('firefox', plan['firefox'] or 'not found')
    for line in desktop_entry_report(env):
        row('desktop file', line)
    row('open babel', find_obabel() or 'not found (optional)')
    row('log file', log_path())
    print('Tier: ' + describe_plan(plan), file=out)
    if plan['tier'] == 'qt':
        cmd = qt_command(plan['python'], plan['binding'], plan['qt_shell'], find_html(required=False) or 'BlueMolecule.html', ['FILE...'])
    elif plan['tier'] in ('chromium', 'firefox'):
        cmd = browser_command(plan['kind'], plan['exe'], 'http://127.0.0.1:%d/TOKEN/' % args.port,
                              profile_dir(plan['kind'], plan['exe']), plan.get('tab'), x11=plan.get('x11', False))
    else:
        cmd = None
    if cmd:
        print('Command: ' + ' '.join(shlex.quote(str(c)) for c in cmd), file=out)
    if plan.get('qt_missing'):
        print('Hint: for the native window install Qt WebEngine for Python:  ' + install_hint(osr), file=out)
    return 0 if plan['tier'] else 1


def check_qt():
    """bluemolecule --check-qt (used by install.sh): 'BINDING PYTHON' and status 0, or the install hint and 1."""
    qt = pick_qt([probe_interpreter(p, full=True) for p in python_candidates()])
    if qt:
        print('%s %s' % (qt['binding'], qt['python']))
        return 0
    print(install_hint(read_os_release()))
    return 1


def main(argv=None):
    ap = argparse.ArgumentParser(prog='bluemolecule', description='Open quantum-chemistry files in BlueMolecule.')
    ap.add_argument('files', nargs='*', help='output, input, structure, volume and MD files of the supported programs '
                    '(Gaussian, ORCA, GAMESS, NWChem, Q-Chem, Psi4, xtb, MOPAC, CP2K, VASP, PDB, mmCIF, GROMACS, Amber, '
                    'CHARMM, LAMMPS and more; .gz allowed)')
    env_tier = os.environ.get('BLUEMOLECULE_TIER', 'auto')
    ap.add_argument('--tier', choices=TIERS, default=env_tier if env_tier in TIERS else 'auto',
                    help='window kind: auto (default; also $BLUEMOLECULE_TIER) tries qt (native Qt WebEngine window, '
                         'Linux), then chromium (Chrome/Chromium/Edge/Brave/Vivaldi app window), then firefox')
    ap.add_argument('--browser', help='path of the browser executable to use')
    ap.add_argument('--tab', action='store_true', help='open in a normal browser tab instead of an app window')
    ap.add_argument('--port', type=int, default=DEFAULT_PORT,
                    help='preferred port of the local server (default %d; the next %d are tried when it is busy; '
                         '0 = any free port)' % (DEFAULT_PORT, PORT_SPAN - 1))
    ap.add_argument('--no-open', action='store_true', help='print the address instead of opening a browser')
    ap.add_argument('--diagnose', action='store_true', help='show the environment and which window kind would be used, and why')
    ap.add_argument('--check-qt', action='store_true', help=argparse.SUPPRESS)
    ap.add_argument('--version', action='version', version='BlueMolecule ' + VERSION)
    args = ap.parse_args(argv)
    if env_tier not in TIERS:
        print('bluemolecule: ignoring BLUEMOLECULE_TIER=%s (use one of %s)' % (env_tier, ', '.join(TIERS)), file=sys.stderr)
    if args.diagnose:
        return diagnose(args)
    if args.check_qt:
        return check_qt()

    html = find_html()
    files = []
    for f in args.files:
        p = Path(f).expanduser()
        if not p.is_file():
            print('bluemolecule: skipping %s (not a file)' % f, file=sys.stderr)
            continue
        files.append(p.resolve())

    plan = None
    if not args.no_open:
        plan = plan_launch(args)
        log_event(describe_plan(plan))
        if plan.get('qt_missing'):
            log_event('hint: for the native BlueMolecule window install Qt WebEngine for Python:  ' + install_hint(read_os_release()))
        if plan['tier'] is None:
            return 2
        if plan['tier'] == 'qt':
            cmd = qt_command(plan['python'], plan['binding'], plan['qt_shell'], html, files, Path(__file__).resolve(),
                             sys.executable, forced=args.tier == 'qt')
            sys.stdout.flush()
            sys.stderr.flush()
            try:
                os.execv(cmd[0], cmd)     # the Qt window replaces this process
            except OSError as e:
                plan = dict(choose_tier(sys.platform, 'auto', None, False, None, None, plan['chromium'], plan['firefox'], use_x11()),
                            reason='%s could not be run (%s)' % (cmd[0], e))
                log_event(describe_plan(plan))

    session = Session(html, files)
    try:
        srv = bind_server(args.port)
    except OSError as e:
        print('bluemolecule: cannot start the local server: %s' % e, file=sys.stderr)
        return 1
    srv.session = session
    port = srv.server_address[1]
    if args.port and not args.port <= port < args.port + PORT_SPAN:
        print('bluemolecule: ports %d-%d are busy, using %d (settings saved in the viewer will not carry over)'
              % (args.port, args.port + PORT_SPAN - 1, port), file=sys.stderr)
    url = 'http://127.0.0.1:%d/%s/' % (port, session.token)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    if args.no_open:
        print(url, flush=True)
    else:
        kinds = (plan['tier'],) if args.tier in ('chromium', 'firefox') and not args.browser else ('chromium', 'firefox')
        how, kind, x11 = open_window(url, args.browser, args.tab, kinds)
        if kind != plan.get('tier') or (kind == 'chromium' and how != plan.get('exe')):
            log_event('opened %s instead' % os.path.basename(how))
        print('BlueMolecule: opened %d file(s) in %s. Close the window (or press Ctrl+C) to quit.' %
              (len(files), os.path.basename(how)), file=sys.stderr)
    try:
        while True:
            time.sleep(0.5)
            why = session.idle_reason()
            if why:
                if why == 'never connected':
                    print('bluemolecule: the viewer never connected; stopping.', file=sys.stderr)
                break
    except KeyboardInterrupt:
        pass
    srv.shutdown()
    return 0


if __name__ == '__main__':
    sys.exit(main())
