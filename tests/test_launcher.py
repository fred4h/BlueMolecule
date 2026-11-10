#!/usr/bin/env python3
"""Tests for launcher/bluemolecule.py and native/linux/bluemolecule_qt.py: --no-open serving,
token/host/origin checks, the fixed port (also when it is busy), heartbeat states, the Open Babel
endpoint (with a fake obabel on PATH), the browser command lines, tier selection, Qt discovery and
install hints, --diagnose/--check-qt, the Qt shell's message handling (Bridge) and its JavaScript as
native-bridge.js sees it, a run of the Qt shell against recording fake PyQt6/PySide6/PyQt5 modules (no
real Qt is needed or checked), the exec into the Qt window and its fall-back to a (fake) Brave window.
Usage: python3 tests/test_launcher.py
"""
import re
import base64
import http.client
import importlib.util
import json
import os
import shutil
import socket
import subprocess
import sys
import tempfile
import threading
import time
from pathlib import Path

sys.dont_write_bytecode = True        # no __pycache__ next to the launcher or the Qt shell
ROOT = Path(__file__).resolve().parent.parent
TMP = Path(tempfile.mkdtemp(prefix='bm-launcher-test-'))
passes = failures = 0


def check(cond, msg):
    global passes, failures
    if cond:
        passes += 1
    else:
        failures += 1
        print('  FAIL:', msg)


# The launcher looks for BlueMolecule.html next to itself: run a copy beside a small stand-in page.
APPDIR = TMP / 'app'
APPDIR.mkdir()
LAUNCHER = APPDIR / 'bluemolecule'
shutil.copy(ROOT / 'launcher' / 'bluemolecule.py', LAUNCHER)
(APPDIR / 'BlueMolecule.html').write_text('<!doctype html><html><head><title>t</title></head><body>viewer</body></html>')
spec = importlib.util.spec_from_file_location('bm_launcher', str(ROOT / 'launcher' / 'bluemolecule.py'))
L = importlib.util.module_from_spec(spec)
spec.loader.exec_module(L)

FAKE_OBABEL = r'''#!/usr/bin/env python3
import sys, time
a = sys.argv[1:]
if a == ['-L', 'formats']:
    print('mol2 -- Sybyl Mol2 format'); print('xyz -- XYZ cartesian coordinates format')
    print('log -- Generic output file format [Read-only]'); sys.exit(0)
if '--sleep' in a:
    time.sleep(30)
data = sys.stdin.buffer.read()
if '-ibad' in a:
    sys.stderr.write('Cannot read input format "bad"\n'); sys.exit(1)
if '-iempty' in a:
    sys.stderr.write('0 molecules converted\n'); sys.exit(0)
sys.stdout.buffer.write(b'CONVERTED ' + ' '.join(a).encode() + b'\n' + data)
sys.stderr.write('1 molecule converted\n')
'''
BIN = TMP / 'bin'
BIN.mkdir()
(BIN / 'obabel').write_text(FAKE_OBABEL)
(BIN / 'obabel').chmod(0o755)
ENV = dict(os.environ, PATH=str(BIN) + os.pathsep + os.environ.get('PATH', ''), PYTHONUNBUFFERED='1', PYTHONDONTWRITEBYTECODE='1')


FAKE_QT = r'''"""Recording stand-ins for the parts of Qt that bluemolecule_qt.py uses (tests only: no Qt API is checked)."""
import base64, os
CALLS, ON_EXEC = [], []


def _arg(a):
    if a is None or isinstance(a, (str, int, float, bool)):
        return a
    return a._n if isinstance(a, Any) else type(a).__name__


def rec(name, args):
    CALLS.append([name, [_arg(a) for a in args]])


class Any(object):
    def __init__(self, name):
        self._n = name
    def __getattr__(self, name):
        if name.startswith('__'):
            raise AttributeError(name)
        return Any(self._n + '.' + name)
    def __call__(self, *a, **k):
        rec(self._n, a)
        return Any(self._n + '()')
    def __eq__(self, other):
        return isinstance(other, Any) and other._n == self._n
    def __ne__(self, other):
        return not self.__eq__(other)
    def __hash__(self):
        return hash(self._n)


class Meta(type):
    def __getattr__(cls, name):
        if name.startswith('__'):
            raise AttributeError(name)
        return Any(cls._qt + '.' + name)


class Base(metaclass=Meta):
    _qt = 'Base'
    def __init__(self, *a, **k):
        rec(self._qt + '()', a)
    def __getattr__(self, name):
        if name.startswith('__'):
            raise AttributeError(name)
        return Any(self._qt + '.' + name)


class _Bound(object):
    def __init__(self):
        self.slots = []
    def connect(self, fn):
        self.slots.append(fn)
    def emit(self, *a):
        for fn in list(self.slots):
            fn(*a)


class Signal(object):
    def __init__(self, *types):
        self.name = None
    def __set_name__(self, owner, name):
        self.name = name
    def __get__(self, obj, owner=None):
        return self if obj is None else obj.__dict__.setdefault('_sig_' + self.name, _Bound())


def Slot(*types):
    return lambda fn: fn


def qVersion():
    return '6.9.0'


def qWebEngineVersion():
    return '6.9.0'


def delete(obj):
    rec('delete', [obj])


class Bytes(object):
    def __init__(self, b):
        self.b = b
    def data(self):
        return self.b
    def toBase64(self):
        return Bytes(base64.b64encode(self.b))


def _cls(name, **extra):
    return Meta(name, (Base,), dict(extra, _qt=name))


QObject, QTimer, QUrl, QIODevice, QByteArray, Qt = (_cls(n) for n in ('QObject', 'QTimer', 'QUrl', 'QIODevice', 'QByteArray', 'Qt'))
QCoreApplication, QKeySequence, QDesktopServices, QGuiApplication, QShortcut, QMessageBox = (_cls(n) for n in (
    'QCoreApplication', 'QKeySequence', 'QDesktopServices', 'QGuiApplication', 'QShortcut', 'QMessageBox'))
QWebEngineProfile, QWebEngineSettings, QWebEngineScript, QWebChannel, QWebEnginePage = (_cls(n) for n in (
    'QWebEngineProfile', 'QWebEngineSettings', 'QWebEngineScript', 'QWebChannel', 'QWebEnginePage'))
QIcon = _cls('QIcon', isNull=lambda self: False)


class QFile(Base):
    _qt = 'QFile'
    def open(self, mode):
        rec('QFile.open', [mode])
        return True
    def readAll(self):
        return Bytes(b'var QWebChannel = function () {}; // fake qwebchannel.js')
    def close(self):
        pass


class QApplication(Base):
    _qt = 'QApplication'
    def exec(self):
        rec('QApplication.exec', [])
        for fn in ON_EXEC:
            fn()
        return 0
    exec_ = exec


class QFileDialog(Base):
    _qt = 'QFileDialog'
    @staticmethod
    def getSaveFileName(parent, caption, start):
        rec('QFileDialog.getSaveFileName', [caption, start])
        return (os.environ.get('FAKE_SAVE_PATH', ''), '')


class QWebEngineView(Base):
    _qt = 'QWebEngineView'
    def saveGeometry(self):
        return Bytes(b'geometry')
    def restoreGeometry(self, g):
        rec('QWebEngineView.restoreGeometry', [g])
        return False
    def isVisible(self):
        return False
'''


def make_fake_qt(root, spec):
    """Fake Qt bindings built on FAKE_QT; spec maps a binding to 'ok' or 'broken' (QtWebEngine does not import)."""
    root.mkdir(parents=True)
    (root / 'bmfakeqt.py').write_text(FAKE_QT)
    core = 'from bmfakeqt import QObject, QTimer, QUrl, QFile, QIODevice, QByteArray, Qt, QCoreApplication, qVersion\n'
    for b, state in spec.items():
        pyside, qt5 = b == 'PySide6', b == 'PyQt5'
        web = 'QWebEnginePage, QWebEngineProfile, QWebEngineSettings, QWebEngineScript'
        files = {
            '__init__.py': "__version__ = '6.9.0'\n" if pyside else '',
            'QtCore.py': core + ('from bmfakeqt import Slot, Signal\n' if pyside else
                                 "from bmfakeqt import Slot as pyqtSlot, Signal as pyqtSignal\nPYQT_VERSION_STR = '%s'\n"
                                 % ('5.15.11' if qt5 else '6.9.0')),
            'QtGui.py': 'from bmfakeqt import QIcon, QKeySequence, QDesktopServices, QGuiApplication%s\n' % ('' if qt5 else ', QShortcut'),
            'QtWidgets.py': 'from bmfakeqt import QApplication, QFileDialog, QMessageBox%s\n' % (', QShortcut' if qt5 else ''),
            'QtWebEngineWidgets.py': 'from bmfakeqt import QWebEngineView%s\n' % ((', ' + web) if qt5 else ''),
            'QtWebChannel.py': 'from bmfakeqt import QWebChannel\n',
        }
        if not qt5:
            files['QtWebEngineCore.py'] = 'from bmfakeqt import %s, qWebEngineVersion\n' % web
        if pyside:
            (root / 'shiboken6.py').write_text('from bmfakeqt import delete\n')
        else:
            files['sip.py'] = 'from bmfakeqt import delete\n'
        if state == 'broken':
            files['QtWebEngineWidgets.py'] = "raise ImportError('libQt6WebEngineCore.so.6: cannot open shared object file')\n"
        (root / b).mkdir()
        for name, text in files.items():
            (root / b / name).write_text(text)


# Runs Shell.run() of bluemolecule_qt.py against the fake Qt; while "exec" runs, it plays the page's part.
SMOKE_DRIVER = r'''
import importlib.util, json, sys, time
sys.path.insert(0, sys.argv[1])
spec = importlib.util.spec_from_file_location('bmqt', sys.argv[2])
S = importlib.util.module_from_spec(spec)
spec.loader.exec_module(S)
import bmfakeqt as F
binding, html, launcher, files = sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6:]
Q = S.import_qt(binding)
shell = S.Shell(Q, html, files, S.load_launcher(launcher))
logs = []
shell.log = logs.append
out = {}

class Url(object):
    def __init__(self, u):
        self.u = u
    def scheme(self):
        return self.u.split(':', 1)[0]
    def toLocalFile(self):
        return self.u[7:] if self.u.startswith('file://') else ''

class Item(object):
    def downloadFileName(self):
        return 'pic.png'
    def __getattr__(self, n):
        return lambda *a: F.rec('download.' + n, a)

def page():
    h = shell.host
    h.post('bmReady', json.dumps({'version': '1.2.0', 'app': True}))
    h.post('bmCopy', json.dumps({'text': 'CCO'}))
    h.post('bmLog', json.dumps({'level': 'error', 'message': 'boom'}))
    h.post('bmObabel', json.dumps({'id': 'bm7', 'op': 'formats'}))
    for _ in range(250):
        if any('"bm7"' in str(c) for c in F.CALLS):
            break
        time.sleep(0.02)
    out['nav'] = [shell.navigation(Url('https://example.org/'), True), shell.navigation(Url('file://' + html), True),
                  shell.navigation(Url('file://' + files[1]), True), shell.navigation(Url('data:text/html,x'), True),
                  shell.navigation(Url('data:text/html,x'), False)]
    out['popup'] = type(shell.popup()).__name__
    shell.toggle_devtools()
    shell.loaded(True)
    shell.crashed('CrashedTerminationStatus', 11)
    shell.download(Item())
    shell.console(Q.QWebEnginePage.JavaScriptConsoleMessageLevel.ErrorMessageLevel, 'oops', 3, 'file:///x/BlueMolecule.html')
    shell.closing()

F.ON_EXEC.append(page)
out['rc'] = shell.run()
out['calls'] = F.CALLS
out['logs'] = logs
print(json.dumps(out))
'''

# Feeds JavaScript made by the Qt shell's Bridge to native-bridge.js (kind linux-app, mocked QWebChannel).
NODE_CONTRACT = r'''
const fs = require('fs'), vm = require('vm');
const SRC = fs.readFileSync(process.argv[1], 'utf8');
const scripts = JSON.parse(fs.readFileSync(process.argv[2], 'utf8'));
const posted = [];
const bmHost = { post(name, text) { posted.push([name, JSON.parse(text)]); } };
const g = { console, setTimeout, clearTimeout, atob, btoa, TextEncoder, TextDecoder, Blob, File, URL,
            BM_QT: { shell: 'linux-app', readLimit: 6291456 }, qt: { webChannelTransport: {} },
            QWebChannel: function (t, cb) { setTimeout(() => cb({ objects: { bmHost } }), 0); } };
g.window = g;
vm.createContext(g);
vm.runInContext(SRC, g);
const N = g.BM.native, got = [];
N.onOpenFile((a, b) => got.push(b === undefined ? { handle: true, name: a.name, size: a.size, read: typeof a.readChunk }
                                                 : { name: a, text: Buffer.from(b).toString('utf8') }));
(async () => {
  const saved = N.save('x.xyz', 'data');      // the first request: id bm1
  await new Promise(r => setTimeout(r, 10));   // ... registered once its bytes are read
  for (const s of scripts) vm.runInContext(s, g);
  const ok = await saved;
  await new Promise(r => setTimeout(r, 20));
  console.log(JSON.stringify({ kind: N.kind, got, saved: ok, lastSave: N.lastSave, posted }));
  process.exit(0);                             // announceReady() would keep polling for BM.app
})();
'''

FAKE_SHELL = r'''import json, os, sys
with open(os.environ['FAKE_SHELL_OUT'], 'w') as fh:
    json.dump({'argv': sys.argv, 'exe': sys.executable}, fh)
'''

FAKE_BROWSER = r'''#!/usr/bin/env python3
import json, os, sys, urllib.request
with open(os.environ['FAKE_BROWSER_OUT'], 'w') as fh:
    json.dump({'argv': sys.argv, 'env': {k: os.environ.get(k) for k in ('CHROME_DESKTOP', 'BLUEMOLECULE_QT_SKIP')}}, fh)
url = [a for a in sys.argv if a.startswith('--app=')][0][6:]
web = urllib.request.build_opener(urllib.request.ProxyHandler({}))
web.open(url, timeout=10).read()
web.open(urllib.request.Request(url + 'ping?state=closed', data=b'', method='POST'), timeout=10).read()
'''


class FakeUI:
    def __init__(self):
        self.save_answer, self.asked, self.clip, self.errors = None, [], [], []

    def ask_save_path(self, suggested, mime, done):
        self.asked.append((suggested, mime))
        done(self.save_answer)

    def set_clipboard(self, text):
        self.clip.append(text)

    def background(self, work, done):
        done(work())

    def error(self, title, text):
        self.errors.append((title, text))


JS_CALL_RE = re.compile(r"^\(function \(m\) \{ try \{ (.*?); \} catch \(e\) \{ console\.error\('BlueMolecule bridge: ' \+ e\); \} \}\)\((.*)\);$", re.S)


def js_args(script):
    m = JS_CALL_RE.match(script)
    return (m.group(1), json.loads(m.group(2))) if m else (None, None)


def free_port_block(n=3):
    """A port p such that p..p+n-1 are currently free."""
    for _ in range(200):
        s = socket.socket()
        s.bind(('127.0.0.1', 0))
        p = s.getsockname()[1]
        s.close()
        if p + n >= 65535:
            continue
        socks, ok = [], True
        for q in range(p, p + n):
            t = socket.socket()
            try:
                t.bind(('127.0.0.1', q))
                socks.append(t)
            except OSError:
                ok = False
                t.close()
                break
        for t in socks:
            t.close()
        if ok:
            return p
    raise RuntimeError('no free port block')


def start(*args):
    p = subprocess.Popen([sys.executable, str(LAUNCHER), '--no-open'] + list(args), stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE, env=ENV, text=True)
    timer = threading.Timer(15, p.kill)
    timer.start()
    url = p.stdout.readline().strip()
    timer.cancel()
    port = int(url.split(':')[2].split('/')[0]) if url.startswith('http://127.0.0.1:') else None
    token = url.rstrip('/').rsplit('/', 1)[-1] if port else None
    return p, url, port, token


def stop(p):
    if p.poll() is None:
        p.terminate()
        try:
            p.wait(5)
        except subprocess.TimeoutExpired:
            p.kill()


def req(port, method, path, body=None, headers=None):
    c = http.client.HTTPConnection('127.0.0.1', port, timeout=30)
    c.request(method, path, body=body, headers=headers or {})
    r = c.getresponse()
    data = r.read()
    c.close()
    return r.status, data


def obabel(port, token, payload, headers=None):
    h = {'Content-Type': 'application/json'}
    h.update(headers or {})
    status, data = req(port, 'POST', '/%s/obabel' % token, json.dumps(payload).encode(), h)
    return status, (json.loads(data) if data.startswith(b'{') else data)


procs = []
try:
    # ---------------- serving with files ----------------
    print('Serving, token and host checks')
    f1, f2 = TMP / 'water.xyz', TMP / 'with space.cube'
    f1.write_text('3\nwater\nO 0 0 0\nH 0 0 0.96\nH 0.93 0 -0.24\n')
    f2.write_bytes(b'\x00\x01binary\xff')
    base = free_port_block(3)
    p, url, port, token = start('--port', str(base), str(f1), str(f2), str(TMP / 'missing.log'))
    procs.append(p)
    check(port == base, 'preferred port used: %s (wanted %d)' % (url, base))
    status, page = req(port, 'GET', '/%s/' % token)
    check(status == 200 and b'viewer' in page, 'page served (%d)' % status)
    head = page.decode()
    launch = json.loads(head.split('window.BM_LAUNCH=', 1)[1].split(';</script>', 1)[0]) if 'window.BM_LAUNCH=' in head else {}
    check(head.index('window.BM_LAUNCH=') < head.index('<title>'), 'BM_LAUNCH injected at the start of <head>')
    check([f['name'] for f in launch.get('files', [])] == ['water.xyz', 'with space.cube'], 'BM_LAUNCH files: %s' % launch.get('files'))
    check(launch.get('ping') == '/%s/ping' % token and launch.get('obabel') == '/%s/obabel' % token, 'ping and obabel endpoints: %s' % launch)
    check(launch.get('version') == L.VERSION, 'version in BM_LAUNCH')
    status, data = req(port, 'GET', launch['files'][1]['url'])
    check(status == 200 and data == f2.read_bytes(), 'binary file served byte for byte')
    check(req(port, 'GET', '/%s/file/0/water.xyz' % token)[1] == f1.read_bytes(), 'first file served')
    check(req(port, 'GET', '/%s/file/7/x' % token)[0] == 404, 'unknown file index -> 404')
    check(req(port, 'GET', '/%s/file/../../etc/passwd' % token)[0] == 404, 'path tricks -> 404')
    check(req(port, 'GET', '/wrongtoken/')[0] == 404, 'wrong token -> 404')
    check(req(port, 'GET', '/')[0] == 404, 'root -> 404')
    check(req(port, 'GET', '/%s/' % token, headers={'Host': 'evil.example:%d' % port})[0] == 403, 'foreign Host -> 403')
    check(req(port, 'GET', '/%s/' % token, headers={'Host': 'localhost:%d' % port})[0] == 200, 'localhost Host accepted')
    check(req(port, 'POST', '/%s/ping' % token)[0] == 204, 'ping -> 204')
    check(req(port, 'POST', '/wrongtoken/ping')[0] == 404, 'ping with wrong token -> 404')
    check(req(port, 'POST', '/%s/ping' % token, headers={'Origin': 'https://evil.example'})[0] == 403, 'foreign Origin -> 403')
    check(req(port, 'POST', '/%s/ping' % token, headers={'Origin': 'http://127.0.0.1:%d' % port})[0] == 204, 'own Origin accepted')

    # ---------------- Open Babel endpoint ----------------
    print('Open Babel endpoint')
    st, r = obabel(port, token, {'op': 'formats'})
    check(st == 200 and r.get('ok') and base64.b64decode(r['output']).decode().startswith('mol2 -- Sybyl'), 'formats: %s' % r)
    xyz = f1.read_bytes()
    st, r = obabel(port, token, {'op': 'convert', 'input': base64.b64encode(xyz).decode(), 'inFormat': 'xyz', 'outFormat': 'sdf', 'options': ['-h', '--title', 'my mol']})
    out = base64.b64decode(r.get('output') or '')
    check(r.get('ok') and out.startswith(b'CONVERTED -ixyz -osdf -h --title my mol\n') and out.endswith(xyz), 'convert passes input, formats and options: %r' % out[:80])
    check(r.get('error') == '', 'no error text on success')
    st, r = obabel(port, token, {'op': 'convert', 'input': '', 'inFormat': 'bad', 'outFormat': 'sdf'})
    check(r.get('ok') is False and 'Cannot read input format' in r.get('error', ''), 'obabel failure reported: %s' % r)
    st, r = obabel(port, token, {'op': 'convert', 'input': 'eA==', 'inFormat': 'empty', 'outFormat': 'sdf'})
    check(r.get('ok') is False and '0 molecules converted' in r.get('error', ''), 'empty output is a failure: %s' % r)
    st, r = obabel(port, token, {'op': 'convert', 'input': 'eA==', 'inFormat': 'xyz', 'outFormat': 'sdf', 'options': ['-O', '/tmp/x']})
    check(r.get('ok') is False and 'not allowed' in r.get('error', ''), 'output-file option rejected: %s' % r)
    st, r = obabel(port, token, {'op': 'convert', 'input': 'eA==', 'inFormat': 'x y', 'outFormat': 'sdf'})
    check(r.get('ok') is False and 'format' in r.get('error', ''), 'bad format code rejected')
    st, r = obabel(port, token, {'op': 'convert', 'input': '!!notbase64', 'inFormat': 'xyz', 'outFormat': 'sdf'})
    check(r.get('ok') is False and 'base64' in r.get('error', ''), 'bad base64 rejected')
    st, r = obabel(port, token, {'op': 'rm -rf'})
    check(r.get('ok') is False, 'unknown op rejected')
    st, _ = req(port, 'POST', '/%s/obabel' % token, b'{not json', {'Content-Type': 'application/json'})
    check(st == 400, 'invalid JSON -> 400')
    st, _ = obabel(port, 'wrongtoken', {'op': 'formats'})
    check(st == 404, 'obabel with wrong token -> 404')
    st, _ = obabel(port, token, {'op': 'formats'}, {'Origin': 'http://evil.example'})
    check(st == 403, 'obabel from a foreign origin -> 403')

    # module-level: timeout and "not found"
    r = L.obabel_request({'op': 'convert', 'input': '', 'inFormat': 'xyz', 'outFormat': 'sdf', 'options': ['--sleep']}, exe=str(BIN / 'obabel'), timeout=1)
    check(r['ok'] is False and 'did not finish' in r['error'], 'timeout reported: %s' % r)
    saved_path, saved_dirs = os.environ.get('PATH', ''), L.OBABEL_DIRS
    os.environ['PATH'], L.OBABEL_DIRS = str(TMP / 'empty'), []
    try:
        r = L.obabel_request({'op': 'formats'})
    finally:
        os.environ['PATH'], L.OBABEL_DIRS = saved_path, saved_dirs
    check(r['ok'] is False and 'not found' in r['error'], 'missing obabel reported: %s' % r)
    check(L.option_allowed('--gen3d') and L.option_allowed('-p') and L.option_allowed('7.4'), 'ordinary options allowed')
    check(not any(L.option_allowed(o) for o in ('-O', '-Oout.sdf', '-m', '-ixyz', '-osdf', '-L', '-H', 'a\nb', '', 5)), 'dangerous options blocked')

    # ---------------- fixed port ----------------
    print('Fixed port')
    p2, url2, port2, _ = start('--port', str(base))     # the first launcher holds `base`
    procs.append(p2)
    check(port2 == base + 1, 'busy preferred port -> next port: %s' % url2)
    blocker = socket.socket()
    blocker.bind(('127.0.0.1', base + 2))
    blocker.listen(1)
    p3, url3, port3, _ = start('--port', str(base))
    procs.append(p3)
    check(port3 not in (base, base + 1, base + 2) and port3 is not None, 'ports held by another program are skipped: %s' % url3)
    blocker.close()
    for q in (p2, p3):
        stop(q)
    p4, url4, port4, tok4 = start()   # no files, default port range
    procs.append(p4)
    check(port4 is not None and L.DEFAULT_PORT <= port4 < L.DEFAULT_PORT + L.PORT_SPAN, 'default port range without files: %s' % url4)
    st, page4 = req(port4, 'GET', '/%s/' % tok4)
    check(st == 200 and b'"files": []' in page4, 'served also without files')

    # ---------------- heartbeat ----------------
    print('Heartbeat')
    req(port4, 'POST', '/%s/ping?state=closed' % tok4)
    t0 = time.time()
    try:
        p4.wait(L.CLOSE_GRACE + 6)
    except subprocess.TimeoutExpired:
        pass
    check(p4.poll() is not None and time.time() - t0 < L.CLOSE_GRACE + 5, 'launcher quits a few seconds after "closed" (%.1f s)' % (time.time() - t0))
    S = L.Session(APPDIR / 'BlueMolecule.html', [])
    t = S.started
    check(S.idle_reason(t + 10) is None and S.idle_reason(t + L.CONNECT_TIMEOUT + 1) == 'never connected', 'first-connection timeout')
    S.touch()
    t = S.last_ping
    check(S.idle_reason(t + L.PING_TIMEOUT - 1) is None and S.idle_reason(t + L.PING_TIMEOUT + 1) == 'no heartbeat', 'visible-window timeout')
    S.touch('hidden')
    t = S.last_ping
    check(S.idle_reason(t + 90) is None and S.idle_reason(t + L.HIDDEN_TIMEOUT + 1) == 'no heartbeat', 'hidden windows get the long timeout')
    S.touch('closed')
    check(S.idle_reason(S.closed_at + 1) is None and S.idle_reason(S.closed_at + L.CLOSE_GRACE + 1) == 'closed', 'closed window -> quit after the grace period')
    S.touch('closed')
    S.touch()
    check(S.closed_at is None and S.idle_reason(S.last_ping + L.CLOSE_GRACE + 1) is None, 'a reload (new heartbeat) cancels "closed"')
    S.touch('visible')
    check(S.hidden is False, 'visible clears hidden')
    # a sleeping laptop: the wall clock jumps ahead, the monotonic clock does not -> no "no heartbeat" on wake-up
    real_time = L.time.time
    try:
        L.time.time = lambda: real_time() + 3600
        check(S.idle_reason() is None, 'the heartbeat ignores wall-clock jumps (suspend/resume): %r' % S.idle_reason())
    finally:
        L.time.time = real_time

    # ---------------- browser command lines ----------------
    print('Browser command lines')
    u, prof = 'http://127.0.0.1:47591/TOKEN/', '/home/u/.local/share/bluemolecule/browser'
    c = L.browser_command('chromium', '/usr/bin/google-chrome', u, prof, platform='linux')
    check(c == ['/usr/bin/google-chrome', '--app=' + u, '--class=org.bluemolecule.BlueMolecule', '--name=org.bluemolecule.BlueMolecule',
                '--user-data-dir=' + prof, '--no-first-run', '--no-default-browser-check'], 'chromium app window: %s' % c)
    c = L.browser_command('chromium', '/usr/bin/brave-browser', u, prof, platform='linux', x11=True)
    check(c[-1] == '--ozone-platform=x11' and '--class=org.bluemolecule.BlueMolecule' in c and c.count('--ozone-platform=x11') == 1,
          'chromium on Wayland with XWayland: forced to X11 so that --class becomes the WM_CLASS: %s' % c)
    c = L.browser_command('chromium', '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome', u, prof, platform='darwin', x11=False)
    check(c == ['/Applications/Google Chrome.app/Contents/MacOS/Google Chrome', '--app=' + u, '--class=org.bluemolecule.BlueMolecule',
                '--user-data-dir=' + prof, '--no-first-run', '--no-default-browser-check'], 'chromium on macOS unchanged: %s' % c)
    check(L.use_x11({'DISPLAY': ':0', 'WAYLAND_DISPLAY': 'wayland-0'}, 'linux') and not L.use_x11({'WAYLAND_DISPLAY': 'wayland-0'}, 'linux')
          and not L.use_x11({'DISPLAY': ':0'}, 'darwin'), 'X11 is forced only when an X server (XWayland) is there')
    e = L.browser_env('chromium', {'A': '1'}, 'linux')
    check(e == {'A': '1', 'CHROME_DESKTOP': 'org.bluemolecule.BlueMolecule.desktop'}, 'chromium env names the desktop file: %s' % e)
    e = L.browser_env('firefox', {}, 'linux')
    check(e == {'MOZ_APP_REMOTINGNAME': 'org.bluemolecule.BlueMolecule'}, 'firefox env sets the Wayland app_id: %s' % e)
    check(L.browser_env('firefox', {}, 'darwin') == {} and L.browser_env('chromium', {}, 'linux', tab=True) == {}, 'no identity env on macOS or for tabs')
    c = L.browser_command('firefox', '/usr/bin/firefox', u, '/p/ff', platform='linux')
    check(c == ['/usr/bin/firefox', '--new-instance', '--class', 'org.bluemolecule.BlueMolecule', '--name', 'org.bluemolecule.BlueMolecule',
                '--profile', '/p/ff', u], 'firefox new instance: %s' % c)
    c = L.browser_command('firefox', '/usr/bin/firefox', u, '/p/ff', running=True, platform='linux')
    check(c == ['/usr/bin/firefox', '--class', 'org.bluemolecule.BlueMolecule', '--name', 'org.bluemolecule.BlueMolecule',
                '--profile', '/p/ff', '--new-window', u], 'firefox already running -> new window: %s' % c)
    c = L.browser_command('firefox', '/Applications/Firefox.app/Contents/MacOS/firefox', u, '/p/ff', platform='darwin')
    check('--class' not in c and c[1] == '--new-instance' and c[-1] == u, 'firefox on macOS has no X11 class: %s' % c)
    check(L.browser_command('chromium', 'chrome', u, prof, tab=True) == ['chrome', u], '--tab opens a normal tab')
    check(L.browser_kind('/usr/bin/firefox-esr') == 'firefox' and L.browser_kind('/usr/bin/chromium') == 'chromium', 'browser kind')
    os.environ['XDG_DATA_HOME'] = str(TMP / 'xdg')
    check(L.profile_dir('chromium', '/usr/bin/google-chrome') == TMP / 'xdg' / 'bluemolecule' / 'browser', 'chromium profile under $XDG_DATA_HOME')
    check(L.snap_name('/snap/bin/chromium') == 'chromium', 'snap browser detected')
    wrapper = TMP / 'firefox-wrapper'
    wrapper.write_text('#!/bin/sh\nexec /snap/bin/firefox "$@"\n')
    check(L.snap_name(str(wrapper)) == 'firefox' and L.profile_dir('firefox', str(wrapper)) == Path.home() / 'snap' / 'firefox' / 'common' / 'bluemolecule-firefox',
          'snap wrapper script -> profile in ~/snap/firefox/common')
    check(L.snap_name(sys.executable) is None, 'ordinary executable is not a snap')
    ff = TMP / 'ffprofile'
    L.ensure_firefox_profile(ff)
    userjs = (ff / 'user.js').read_text()
    css = (ff / 'chrome' / 'userChrome.css').read_text()
    check('toolkit.legacyUserProfileCustomizations.stylesheets", true' in userjs and 'browser.download.useDownloadDir", false' in userjs, 'user.js enables userChrome.css and Save dialogs')
    check('#TabsToolbar' in css and '#nav-bar' in css and 'visibility: collapse' in css, 'userChrome.css hides tabs and toolbars')
    check(L.firefox_profile_running(ff) is False, 'fresh profile is not running')
    (ff / '.parentlock').write_bytes(b'')
    holder = subprocess.Popen([sys.executable, '-c', 'import fcntl,os,sys,time; fd=os.open(sys.argv[1], os.O_RDWR); '
                               'fcntl.lockf(fd, fcntl.LOCK_EX); print("locked", flush=True); time.sleep(20)', str(ff / '.parentlock')],
                              stdout=subprocess.PIPE, text=True)
    holder.stdout.readline()
    check(L.firefox_profile_running(ff) is True, 'locked profile is detected as running')
    holder.kill()
    holder.wait()
    check(L.firefox_profile_running(ff) is False, 'unlocked profile is free again')

    # ---------------- tier selection (pure) ----------------
    print('Tier selection')
    QTF, SH, BRAVE, FF = {'python': '/usr/bin/python3', 'binding': 'PyQt6'}, '/d/bluemolecule_qt.py', ['/usr/bin/brave-browser'], '/usr/bin/firefox'
    p = L.choose_tier('linux', 'auto', None, False, SH, QTF, BRAVE, FF, True)
    check(p['tier'] == 'qt' and p['python'] == '/usr/bin/python3' and p['binding'] == 'PyQt6' and not p['qt_missing'], 'Qt found: tier 1 %s' % p)
    p = L.choose_tier('linux', 'auto', None, False, SH, None, BRAVE, FF, True)
    check(p['tier'] == 'chromium' and p['exe'] == BRAVE[0] and p['x11'] and p['qt_missing'] and 'PyQt6' in p['reason'], 'no Qt: tier 2 under X11 %s' % p)
    p = L.choose_tier('linux', 'auto', None, False, SH, None, BRAVE, FF, False)
    check(p['tier'] == 'chromium' and not p['x11'], 'no DISPLAY: X11 is not forced')
    p = L.choose_tier('linux', 'auto', None, False, None, QTF, [], FF, True)
    check(p['tier'] == 'firefox' and p['exe'] == FF and not p['qt_missing'] and 'missing' in p['reason'], 'no shell, no Chromium: tier 3, no hint %s' % p)
    p = L.choose_tier('linux', 'auto', None, False, SH, None, [], None, True)
    check(p['tier'] == 'default' and p['qt_missing'], 'nothing found: default browser')
    p = L.choose_tier('darwin', 'auto', None, False, SH, QTF, ['/Applications/Brave Browser.app/Contents/MacOS/Brave Browser'], None, False)
    check(p['tier'] == 'chromium' and not p['qt_missing'] and not p['x11'], 'never the Qt window on macOS: %s' % p)
    p = L.choose_tier('linux', 'auto', None, True, SH, QTF, BRAVE, FF, True)
    check(p['tier'] == 'chromium' and p['tab'] and not p['x11'] and not p['qt_missing'], '--tab: a browser tab, no hint')
    p = L.choose_tier('linux', 'auto', '/opt/ff/firefox', False, SH, QTF, BRAVE, None, True)
    check(p['tier'] == 'firefox' and p['exe'] == '/opt/ff/firefox' and not p['qt_missing'], '--browser wins')
    p = L.choose_tier('linux', 'chromium', None, False, SH, QTF, BRAVE, FF, True)
    check(p['tier'] == 'chromium' and not p['qt_missing'], '--tier chromium skips Qt')
    p = L.choose_tier('linux', 'firefox', None, False, SH, QTF, BRAVE, FF, True)
    check(p['tier'] == 'firefox', '--tier firefox')
    p = L.choose_tier('linux', 'qt', None, False, SH, None, BRAVE, FF, True)
    check(p['tier'] is None and p['qt_missing'] and 'not available' in p['error'], '--tier qt without Qt is an error')
    check(L.choose_tier('linux', 'firefox', None, False, SH, QTF, BRAVE, None, True)['tier'] is None, '--tier firefox without Firefox is an error')
    d1 = L.describe_plan(L.choose_tier('linux', 'auto', None, False, SH, QTF, BRAVE, FF, True))
    d2 = L.describe_plan(L.choose_tier('linux', 'auto', None, False, SH, None, BRAVE, FF, True))
    check(d1.startswith('tier 1') and 'PyQt6' in d1 and 'app_id org.bluemolecule.BlueMolecule' in d1, d1)
    check(d2.startswith('tier 2, brave-browser') and 'ozone-platform=x11' in d2 and 'no Python with PyQt6' in d2, d2)

    # ---------------- Qt discovery and hints ----------------
    print('Qt discovery and install hints')
    conda = '/home/u/miniconda3/bin/python3'
    c = L.python_candidates(conda, lambda n: conda, lambda p: p in ('/usr/bin/python3', conda, '/usr/local/bin/python3'))
    check(c == ['/usr/bin/python3', conda, '/usr/local/bin/python3'], 'system Python first, then conda, each once: %s' % c)
    c = L.python_candidates(sys.executable)
    check(len({os.path.realpath(x) for x in c}) == len(c) and any(os.path.realpath(x) == os.path.realpath(sys.executable) for x in c),
          'real interpreters are listed once: %s' % c)
    interps = [{'python': '/a', 'bindings': {'PyQt5': {'ok': True}, 'PyQt6': {'ok': False}}}, {'python': '/b', 'bindings': {'PyQt6': {'ok': True}}}]
    check(L.pick_qt(interps) == {'python': '/b', 'binding': 'PyQt6', 'info': {'ok': True}}, 'PyQt6 in any Python beats PyQt5')
    check(L.pick_qt(interps, skip={os.path.realpath('/b')})['binding'] == 'PyQt5', 'a Python whose Qt window failed is skipped')
    check(L.pick_qt([{'python': '/a', 'bindings': {}}]) is None, 'no bindings: None')
    check(L.qt_skip_list({'BLUEMOLECULE_QT_SKIP': '/a' + os.pathsep + '/b'}) == {'/a', '/b'} and L.qt_skip_list({}) == set(), 'BLUEMOLECULE_QT_SKIP')
    pr = L.parse_probe('/x', 1, '', 'Traceback\nSyntaxError: boom')
    check(pr['bindings'] == {} and pr['error'].endswith('SyntaxError: boom'), 'a failing probe is reported: %s' % pr)
    real = L.probe_interpreter(sys.executable)
    check(real.get('version') == '%d.%d.%d' % sys.version_info[:3] and set(real['bindings']) == {'PyQt6', 'PySide6', 'PyQt5'},
          'probe of this Python: %s' % real)
    cmd = L.qt_command('/usr/bin/python3', 'PyQt6', '/d/bluemolecule_qt.py', '/d/BlueMolecule.html', ['/f/a b.xyz'], '/d/bluemolecule', '/conda/python3')
    check(cmd == ['/usr/bin/python3', '/d/bluemolecule_qt.py', '--binding', 'PyQt6', '--html', '/d/BlueMolecule.html', '--launcher', '/d/bluemolecule',
                  '--launcher-python', '/conda/python3', '--', '/f/a b.xyz'], 'Qt window command: %s' % cmd)
    check('--no-fallback' in L.qt_command('/p', 'PyQt5', '/s', '/h', [], forced=True), '--tier qt: no fall-back')
    osr = L.parse_os_release('NAME="Ubuntu"\nID=ubuntu\nID_LIKE=debian\nPRETTY_NAME="Ubuntu 26.10"\nVERSION_CODENAME=\'x\'\n')
    check(osr['PRETTY_NAME'] == 'Ubuntu 26.10' and osr['VERSION_CODENAME'] == 'x', 'os-release parsed: %s' % osr)
    check(L.install_hint(osr) == 'sudo apt install python3-pyqt6.qtwebengine', 'Kubuntu/Ubuntu hint')
    check(L.install_hint({'ID': 'neon', 'ID_LIKE': 'ubuntu debian'}).startswith('sudo apt install'), 'KDE neon hint')
    check(L.install_hint({'ID': 'fedora'}) == 'sudo dnf install python3-pyqt6-webengine', 'Fedora hint')
    check(L.install_hint({'ID': 'manjaro', 'ID_LIKE': 'arch'}).startswith('sudo pacman -S'), 'Arch-family hint')
    check(L.install_hint({'ID': 'opensuse-tumbleweed', 'ID_LIKE': 'opensuse suse'}).startswith('sudo zypper install'), 'openSUSE hint')
    check('pip install' in L.install_hint({'ID': 'gentoo'}) and 'pip install' in L.install_hint({}), 'other systems: pip hint')

    # ---------------- the Qt shell's message handling (no Qt) ----------------
    print('Qt shell: compile, Bridge, JavaScript contract')
    SHELL_SRC = ROOT / 'native' / 'linux' / 'bluemolecule_qt.py'
    import py_compile
    try:
        py_compile.compile(str(SHELL_SRC), cfile=str(TMP / 'shell.pyc'), doraise=True)
        check(True, '')
    except py_compile.PyCompileError as e:
        check(False, 'bluemolecule_qt.py compiles: %s' % e)
    spec = importlib.util.spec_from_file_location('bm_qt_shell', str(SHELL_SRC))
    S = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(S)
    check(not any(m.split('.')[0] in ('PyQt6', 'PySide6', 'PyQt5', 'shiboken6') for m in sys.modules), 'importing the shell does not import Qt')
    check(S.INLINE_LIMIT == 32 * 1024 * 1024 and len(base64.b64encode(b'x' * S.READ_LIMIT)) == 8 * 1024 * 1024,
          'limits: inline below 32 MB, 8 MiB of base64 per bmRead reply')
    check(S.safe_file_name('../../etc/passwd') == 'passwd' and S.safe_file_name('a\\b\x00c.xyz') == 'bc.xyz'
          and S.safe_file_name('..') == 'untitled' and S.safe_file_name(None) == 'untitled', 'proposed file names are cleaned')
    body, a = js_args(S.js_call('BM.native._resolve(m.id, m.result)', {'id': 'bm1', 'result': {'t': 'é </script>'}}))
    check(body == 'BM.native._resolve(m.id, m.result)' and a['result']['t'] == 'é </script>', 'js_call JSON round trip (ASCII-only)')
    check(all(ord(ch) < 128 for ch in S.js_call('x', {'t': 'é '})), 'js_call output is ASCII')
    small, big = TMP / 'small.xyz', TMP / 'big.cube'
    small.write_text('3\nwater\nO 0 0 0\nH 0 0 0.96\nH 0.93 0 -0.24\n')
    big.write_bytes(bytes(range(256)) * 40)          # 10240 bytes: "large" with the limits below
    S.INLINE_LIMIT, S.READ_LIMIT = 4096, 1000
    try:
        ui, js, logs = FakeUI(), [], []
        B = S.Bridge(js.append, ui, obabel=lambda req: {'ok': True, 'output': '', 'error': '', 'req': req}, log=logs.append, files=[small, big])
        check(js == [] and B.receive('bmReady', '{"version":"1.2.0","app":true}') and len(js) == 2, 'files wait for bmReady')
        body, a = js_args(js[0])
        check(body == 'BM.native._receiveFile(m.name, m.data)' and a['name'] == 'small.xyz' and base64.b64decode(a['data']) == small.read_bytes(),
              'small file inline: %s' % body)
        body, a = js_args(js[1])
        check(body == 'BM.native._receiveLargeFile(m.name, m.size, m.fileId)' and a == {'name': 'big.cube', 'size': 10240, 'fileId': 'f1'},
              'large file as a descriptor: %s %s' % (body, a))
        B.receive('bmRead', json.dumps({'id': 'bm3', 'fileId': 'f1', 'offset': 100, 'length': 5000}))
        a = js_args(js[-1])[1]
        r = a['result']
        check(a['id'] == 'bm3' and r['ok'] and r['offset'] == 100 and r['length'] == 1000 and r['eof'] is False
              and base64.b64decode(r['base64']) == big.read_bytes()[100:1100], 'bmRead: at most READ_LIMIT bytes per reply')
        B.receive('bmRead', json.dumps({'id': 'bm4', 'fileId': 'f1', 'offset': 10000, 'length': 1000}))
        r = js_args(js[-1])[1]['result']
        check(r['ok'] and r['length'] == 240 and r['eof'] and base64.b64decode(r['base64']) == big.read_bytes()[10000:], 'bmRead: the last bytes')
        B.receive('bmRead', json.dumps({'id': 'bm5', 'fileId': 'f1', 'offset': 99999, 'length': 10}))
        r = js_args(js[-1])[1]['result']
        check(r['ok'] and r['length'] == 0 and r['base64'] == '' and r['eof'], 'bmRead past the end: empty')
        for bad in ({'fileId': 'f9', 'offset': 0, 'length': 1}, {'fileId': 'f1', 'offset': -1, 'length': 1}, {'fileId': 'f1', 'offset': 0, 'length': 1.5},
                    {'fileId': 'f1', 'offset': True, 'length': 1}, {'fileId': 'f1', 'offset': '0', 'length': 1}, {'offset': 0, 'length': 1}):
            B.receive('bmRead', json.dumps(dict(bad, id='bmX')))
            a = js_args(js[-1])[1]
            check(a['id'] == 'bmX' and a['result']['ok'] is False and a['result']['error'], 'bad bmRead refused: %s' % bad)
        n = len(js)
        B.receive('bmRead', json.dumps({'fileId': 'f1', 'offset': 0, 'length': 1}))
        check(len(js) == n, 'a request without an id gets no reply')
        check(not B.receive('bmNope', '{}') and not B.receive('bmRead', '{not json') and len(js) == n, 'unknown messages and bad JSON are ignored')
        B.receive('bmCopy', json.dumps({'text': 'C1=CC=CC=C1'}))
        check(ui.clip == ['C1=CC=CC=C1'], 'bmCopy -> clipboard')
        B.receive('bmLog', json.dumps({'level': 'error', 'message': 'x' * 9000}))
        B.receive('log', json.dumps({'message': 'plain'}))
        check(logs[-2].startswith('[page error] xxx') and len(logs[-2]) < 8100 and logs[-1] == '[page log] plain', 'bmLog / log -> log file')
        B.receive('bmObabel', json.dumps({'id': 'bm6', 'op': 'formats'}))
        a = js_args(js[-1])[1]
        check(a['id'] == 'bm6' and a['result']['ok'] and a['result']['req'] == {'op': 'formats'}, 'bmObabel: the request without its id, reply by id')
        S.Bridge(js.append, ui).receive('bmObabel', json.dumps({'id': 'bm7', 'op': 'formats'}))
        check(js_args(js[-1])[1]['result']['ok'] is False, 'bmObabel without the launcher: an error reply')
        out_dir = TMP / 'saved'
        out_dir.mkdir()
        target = out_dir / 'saved.xyz'
        ui.save_answer = str(target)
        payload = {'id': 'bm8', 'name': '../../etc/passwd', 'mime': 'chemical/x-xyz', 'base64': base64.b64encode(b'3\nsaved\n').decode()}
        B.receive('bmSave', json.dumps(payload))
        a = js_args(js[-1])[1]
        check(target.read_bytes() == b'3\nsaved\n' and a == {'id': 'bm8', 'result': {'ok': True, 'path': os.path.realpath(str(target))}},
              'bmSave writes the chosen file: %s' % a)
        check(ui.asked[-1] == (str(TMP / 'passwd'), 'chemical/x-xyz'), 'Save dialog proposes the bare name in the folder of the opened file: %s' % (ui.asked[-1],))
        target.chmod(0o640)
        B.receive('bmSave', json.dumps(dict(payload, id='bm9', base64=base64.b64encode(b'new').decode())))
        check(target.read_bytes() == b'new' and (target.stat().st_mode & 0o777) == 0o640 and os.listdir(str(out_dir)) == ['saved.xyz'],
              'overwrite keeps the permissions and leaves no temporary file')
        check(ui.asked[-1][0] == os.path.join(os.path.realpath(str(out_dir)), 'passwd'), 'the next Save dialog starts in the last folder')
        ui.save_answer = None
        B.receive('bmSave', json.dumps(dict(payload, id='bm10')))
        check(js_args(js[-1])[1] == {'id': 'bm10', 'result': {'ok': False, 'cancelled': True}}, 'cancelled Save dialog')
        ui.save_answer = str(TMP / 'no such dir' / 'x.xyz')
        B.receive('bmSave', json.dumps(dict(payload, id='bm11')))
        r = js_args(js[-1])[1]['result']
        check(r['ok'] is False and r['error'] and ui.errors and ui.errors[-1][0] == 'Could not save x.xyz', 'unwritable target: error reply and a message')
        B.receive('bmSave', json.dumps(dict(payload, id='bm12', base64='abc')))
        check(js_args(js[-1])[1]['result'] == {'ok': False, 'error': 'The data to save could not be decoded.'}, 'bad base64 is refused')
        del js[:]
        B.open_paths([str(small)])
        check(len(js) == 1 and js_args(js[0])[1]['name'] == 'small.xyz', 'a file dropped later is delivered at once')
        B.page_reloaded()
        del js[:]
        B.receive('bmReady', '{}')
        check([js_args(x)[1]['name'] for x in js] == ['small.xyz', 'big.cube'] and js_args(js[1])[1]['fileId'] == 'f2',
              'after a reload or crash the start-up files are delivered again')
        B3 = S.Bridge(js.append, ui, files=[TMP / 'nope.xyz'])
        B3.receive('bmReady', '{}')
        check(ui.errors[-1][0] == 'Could not open nope.xyz', 'a missing file is reported')

        # the same JavaScript, run by native-bridge.js (kind linux-app) in node
        if shutil.which('node'):
            uni = TMP / 'café.xyz'
            uni.write_text('1\nü\nH 0 0 0\n', encoding='utf-8')
            scripts = []
            B4 = S.Bridge(scripts.append, ui, files=[uni, big])
            ui.save_answer = str(TMP / 'contract.xyz')
            B4.receive('bmReady', '{}')
            B4.receive('bmSave', json.dumps({'id': 'bm1', 'name': 'x.xyz', 'mime': 'text/plain', 'base64': base64.b64encode(b'data').decode()}))
            (TMP / 'scripts.json').write_text(json.dumps(scripts))
            r = subprocess.run(['node', '-e', NODE_CONTRACT, str(ROOT / 'src' / 'native-bridge.js'), str(TMP / 'scripts.json')],
                               capture_output=True, text=True, timeout=30)
            try:
                res = json.loads(r.stdout)
            except ValueError:
                res = {}
            check(res.get('kind') == 'linux-app', 'native-bridge.js detects the Qt window: %s %s' % (r.stdout[-300:], r.stderr[-300:]))
            check(res.get('got') == [{'name': 'café.xyz', 'text': '1\nü\nH 0 0 0\n'}, {'handle': True, 'name': 'big.cube', 'size': 10240, 'read': 'function'}],
                  'onOpenFile gets (name, bytes) and the large-file handle: %s' % res.get('got'))
            check(res.get('saved') is True and (res.get('lastSave') or {}).get('path') == os.path.realpath(str(TMP / 'contract.xyz'))
                  and res.get('posted', [[None]])[0][0] == 'bmSave', 'the Python reply resolves the page\'s bmSave: %s' % res.get('lastSave'))
        else:
            print('  (node not found: JavaScript contract skipped)')
    finally:
        S.INLINE_LIMIT, S.READ_LIMIT = 32 * 1024 * 1024, 6 * 1024 * 1024

    # ---------------- the Qt shell against recording fake bindings ----------------
    print('Qt shell with fake PyQt6 / PySide6 / PyQt5 (structure only, no real Qt)')
    FAKE_ALL = TMP / 'fakeqt-all'
    make_fake_qt(FAKE_ALL, {'PyQt6': 'ok', 'PySide6': 'ok', 'PyQt5': 'ok'})
    (TMP / 'smoke.py').write_text(SMOKE_DRIVER)
    page_html = APPDIR / 'BlueMolecule.html'
    for binding in ('PyQt6', 'PySide6', 'PyQt5'):
        home = TMP / ('smoke-' + binding)
        env = dict(ENV, XDG_DATA_HOME=str(home / 'data'), XDG_CACHE_HOME=str(home / 'cache'), FAKE_SAVE_PATH=str(home / 'dl' / 'pic.png'))
        r = subprocess.run([sys.executable, str(TMP / 'smoke.py'), str(FAKE_ALL), str(SHELL_SRC), binding, str(page_html),
                            str(ROOT / 'launcher' / 'bluemolecule.py'), str(small), str(f1)], capture_output=True, text=True, env=env, timeout=60)
        try:
            res = json.loads(r.stdout.strip().splitlines()[-1])
        except (ValueError, IndexError):
            check(False, '%s smoke run: %s %s' % (binding, r.stdout[-400:], r.stderr[-800:]))
            continue
        calls = res['calls']
        names = [c[0] for c in calls]

        def called(name, args=None):
            return any(c[0] == name and (args is None or c[1] == args) for c in calls)

        def args_of(suffix):
            return [c[1] for c in calls if c[0].endswith(suffix)]
        tag = binding + ': '
        check(res['rc'] == 0 and 'QApplication.exec' in names, tag + 'runs the event loop')
        check(called('QGuiApplication.setDesktopFileName', ['org.bluemolecule.BlueMolecule']), tag + 'desktop file name = Wayland app_id')
        check(called('QCoreApplication.setApplicationName', ['BlueMolecule']) and called('QCoreApplication.setOrganizationName', ['BlueMolecule'])
              and called('QGuiApplication.setApplicationDisplayName', ['BlueMolecule']) and 'QGuiApplication.setWindowIcon' in names,
              tag + 'application name, organisation and icon')
        check(called('QWebEngineProfile()', ['BlueMolecule']) and
              called('QWebEngineProfile.setPersistentStoragePath', [str(home / 'data' / 'bluemolecule' / 'qtwebengine')]), tag + 'persistent profile under $XDG_DATA_HOME/bluemolecule')
        attrs = {a[0].rsplit('.', 1)[-1] for a in args_of('.setAttribute') if a and a[-1] is True}
        check({'WebGLEnabled', 'Accelerated2dCanvasEnabled', 'LocalContentCanAccessFileUrls', 'LocalStorageEnabled', 'AA_ShareOpenGLContexts'} <= attrs,
              tag + 'WebGL, 2D canvas, file access, localStorage: %s' % sorted(attrs))
        check(called('QWebChannel.registerObject', ['bmHost', 'Host']) and called('QWebEnginePage.setWebChannel', ['QWebChannel']), tag + 'QWebChannel with bmHost')
        sources = [a[0] for a in args_of('QWebEngineScript.setSourceCode')]
        check(any('fake qwebchannel.js' in s for s in sources) and any(s.startswith('window.BM_QT = Object.freeze(') and '"linux-app"' in s for s in sources)
              and args_of('QWebEngineScript.setInjectionPoint') == [['QWebEngineScript.InjectionPoint.DocumentCreation']] * 2
              and args_of('QWebEngineScript.setWorldId') == [[0]] * 2, tag + 'qwebchannel.js and the BM_QT marker before the page')
        check(called('QUrl.fromLocalFile', [str(page_html)]) and called('QWebEngineView.setUrl', ['QUrl.fromLocalFile()']) and 'QWebEngineView.show' in names,
              tag + 'loads the bundled page')
        check(called('QKeySequence()', ['F12']) and called('QKeySequence()', ['Ctrl+Shift+I']) and 'QWebEnginePage.setDevToolsPage' in names, tag + 'F12 developer tools')
        scripts = [a[0] for a in args_of('QWebEnginePage.runJavaScript')]
        received = [js_args(s)[1]['name'] for s in scripts if s.startswith('(function') and '_receiveFile' in s]
        check(received == ['small.xyz', 'water.xyz', 'water.xyz'], tag + 'files delivered on bmReady and a dropped file: %s' % received)
        check(any('BM.native._resolve' in s and '"bm7"' in s for s in scripts), tag + 'Open Babel reply from the worker thread')
        check(any(s == S.WEBGL_CHECK for s in scripts), tag + 'WebGL 2 check after loading')
        check(called('QGuiApplication.clipboard().setText', ['CCO']), tag + 'bmCopy -> clipboard')
        check(res['nav'] == [False, True, False, False, True] and 'QDesktopServices.openUrl' in names and res['popup'] == 'ExternalPage',
              tag + 'navigation policy: %s' % res['nav'])
        check(called('QTimer.singleShot', [300, 'QWebEngineView.reload']), tag + 'a crashed page is reloaded once')
        check('download.accept' in names and called('download.setDownloadFileName', ['pic.png']), tag + 'downloads go through a Save dialog')
        check('[page error] boom' in res['logs'] and '[page error] oops (BlueMolecule.html:3)' in res['logs'], tag + 'page log and console errors')
        geo = home / 'data' / 'bluemolecule' / 'qt-window.json'
        check(geo.is_file() and json.loads(geo.read_text()) == {'geometry': base64.b64encode(b'geometry').decode()}, tag + 'window size saved on close')
        deleted = [c[1][0] for c in calls if c[0] == 'delete']
        check(deleted == ['QWebEngineView', 'Window', 'QWebEngineProfile'], tag + 'views deleted before the profile: %s' % deleted)

    # ---------------- end to end: --diagnose, --check-qt, exec into the Qt window, fall-back ----------------
    if sys.platform.startswith('linux'):
        print('Diagnose, exec into the Qt window and fall-back to a browser window')
        BIN2 = TMP / 'bin2'                     # PATH with a fake Brave only (keeps real browsers out)
        BIN2.mkdir()
        (BIN2 / 'brave-browser').write_text(FAKE_BROWSER)
        (BIN2 / 'brave-browser').chmod(0o755)
        os.symlink(sys.executable, str(BIN2 / 'python3'))
        FAKE_BROKEN = TMP / 'fakeqt-broken'
        make_fake_qt(FAKE_BROKEN, {'PyQt6': 'broken'})
        base_env = dict(ENV, PATH=str(BIN2), DISPLAY=':99', WAYLAND_DISPLAY='wayland-0', XDG_SESSION_TYPE='wayland', XDG_CURRENT_DESKTOP='KDE',
                        XDG_STATE_HOME=str(TMP / 'state'), XDG_DATA_HOME=str(TMP / 'e2e-data'), FAKE_SHELL_OUT=str(TMP / 'shell.json'),
                        FAKE_BROWSER_OUT=str(TMP / 'browser.json'))
        base_env.pop('PYTHONPATH', None)
        logfile = TMP / 'state' / 'bluemolecule' / 'bluemolecule.log'

        def run_launcher(args, extra=None, timeout=60):
            return subprocess.run([sys.executable, str(LAUNCHER)] + args, capture_output=True, text=True, timeout=timeout,
                                  env=dict(base_env, **(extra or {})))
        (APPDIR / 'bluemolecule_qt.py').write_text(FAKE_SHELL)   # records how the launcher runs it
        r = run_launcher(['--diagnose'])
        check(r.returncode == 0 and 'Tier: tier 2, brave-browser app window under XWayland' in r.stdout and '--ozone-platform=x11' in r.stdout
              and 'Hint: for the native window install Qt WebEngine for Python:' in r.stdout and 'KDE, wayland session' in r.stdout,
              '--diagnose without Qt: tier 2 and the hint\n' + r.stdout[-1500:] + r.stderr[-500:])
        r = run_launcher(['--diagnose'], {'PYTHONPATH': str(FAKE_ALL)})
        check(r.returncode == 0 and 'Tier: tier 1, native Qt WebEngine window: PyQt6' in r.stdout and 'PyQt6 6.9.0 OK (Qt 6.9.0, QtWebEngine 6.9.0)' in r.stdout
              and 'Hint:' not in r.stdout and '--binding PyQt6' in r.stdout, '--diagnose with Qt: tier 1\n' + r.stdout[-1500:] + r.stderr[-500:])
        r = run_launcher(['--diagnose'], {'PYTHONPATH': str(FAKE_BROKEN)})
        check('PyQt6: ImportError: libQt6WebEngineCore.so.6' in r.stdout and 'Tier: tier 2' in r.stdout, '--diagnose shows why a binding fails\n' + r.stdout[-800:])
        r = run_launcher(['--check-qt'])
        check(r.returncode == 1 and r.stdout.strip() == L.install_hint(L.read_os_release()), '--check-qt without Qt: the hint, status 1: %r' % r.stdout)
        r = run_launcher(['--check-qt'], {'PYTHONPATH': str(FAKE_ALL)})
        check(r.returncode == 0 and r.stdout.startswith('PyQt6 /'), '--check-qt with Qt: %r' % r.stdout)

        # tier 1: the launcher runs the Qt window with the Python that has the binding (the stand-in shell records it)
        r = run_launcher([str(f1)], {'PYTHONPATH': str(FAKE_ALL)}, timeout=30)
        try:
            got = json.loads((TMP / 'shell.json').read_text())
        except (OSError, ValueError):
            got = {'argv': []}
        exp = ['--binding', 'PyQt6', '--html', str(APPDIR.resolve() / 'BlueMolecule.html'), '--launcher', str(LAUNCHER.resolve()),
               '--launcher-python', sys.executable, '--', str(f1.resolve())]
        check(r.returncode == 0 and got['argv'][1:] == exp, 'exec into the Qt window: %s / %s' % (got, r.stderr[-400:]))
        check(os.path.realpath(got.get('exe', '')) == os.path.realpath(L.python_candidates()[0]), 'with the first Python that has the binding')
        check(logfile.is_file() and 'tier 1, native Qt WebEngine window: PyQt6' in logfile.read_text(), 'the choice is logged')

        # the real shell with a binding whose QtWebEngine does not import: it falls back (each Python once), then tier 2
        shutil.copy(SHELL_SRC, APPDIR / 'bluemolecule_qt.py')
        logfile.unlink()
        r = run_launcher([str(f1)], {'PYTHONPATH': str(FAKE_BROKEN)}, timeout=90)
        log = logfile.read_text() if logfile.is_file() else ''
        try:
            br = json.loads((TMP / 'browser.json').read_text())
        except (OSError, ValueError):
            br = {'argv': [], 'env': {}}
        npy = len(L.python_candidates(sys.executable, lambda n: str(BIN2 / n)))   # as the launcher sees them with PATH=BIN2
        check(r.returncode == 0 and log.count('the Qt window cannot start') == npy and log.count('tier 1,') == npy and 'tier 2, brave-browser' in log,
              'fall-back after every Python failed: %d Pythons\n%s\n%s' % (npy, log[-1500:], r.stderr[-800:]))
        check('--ozone-platform=x11' in br['argv'] and '--class=org.bluemolecule.BlueMolecule' in br['argv'] and '--name=org.bluemolecule.BlueMolecule' in br['argv']
              and ('--user-data-dir=' + str(TMP / 'e2e-data' / 'bluemolecule' / 'browser')) in br['argv']
              and br['env'].get('CHROME_DESKTOP') == 'org.bluemolecule.BlueMolecule.desktop', 'Brave started as tier 2: %s' % br)
        check(len((br['env'].get('BLUEMOLECULE_QT_SKIP') or '').split(os.pathsep)) == npy, 'each failed Python is skipped once: %s' % br['env'])
finally:
    for q in procs:
        stop(q)
    shutil.rmtree(TMP, ignore_errors=True)

print('\nLauncher tests: %d passed, %d failed' % (passes, failures))
sys.exit(1 if failures else 0)
