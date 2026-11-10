#!/usr/bin/env python3
"""BlueMolecule for Linux: a native Qt WebEngine window around BlueMolecule.html (tier 1 of the launcher).

    python3 bluemolecule_qt.py [--binding PyQt6|PySide6|PyQt5] [--html PAGE] [FILE ...]

Normally started by the bluemolecule launcher, which picks a Python that has PyQt6, PySide6 or PyQt5 with
QtWebEngine. The window is org.bluemolecule.BlueMolecule to the desktop: QGuiApplication.setDesktopFileName
sets the Wayland app_id that KWin matches to org.bluemolecule.BlueMolecule.desktop (name and icon in the
taskbar); under X11 the WM_CLASS instance is the same name. Settings (localStorage) live in a persistent
profile under $XDG_DATA_HOME/bluemolecule/qtwebengine. The page (src/native-bridge.js, kind 'linux-app')
talks to the window through QWebChannel (channel.objects.bmHost.post(name, json)) with the messages of the
macOS app plus bmRead:
  bmReady  {version, app}                    the page is ready: files from the command line are delivered
  bmSave   {id, name, mime, base64}          Save dialog; reply BM.native._resolve(id, {ok, path})
  bmObabel {id, op, input, inFormat, outFormat, options}   Open Babel; reply {ok, output, error}
  bmCopy   {text}                            plain text to the clipboard
  bmRead   {id, fileId, offset, length}      bytes of a large file; reply {ok, base64, offset, length, eof}
  bmLog    {level, message}                  a line in the log file
Files under 32 MB reach the page as BM.native._receiveFile(name, base64), larger ones as
BM.native._receiveLargeFile(name, size, fileId); the page then reads them in pieces of at most 6 MiB
(8 MiB of base64) with bmRead. F12 opens the developer tools. If no Qt binding can be imported, the
launcher is started again (with this Python excluded) and opens a browser window instead.
The message handling (class Bridge) does not need Qt, so it is tested without it.
"""
import argparse
import base64
import binascii
import importlib
import importlib.machinery
import importlib.util
import json
import os
import sys
import tempfile
import threading
import types
from pathlib import Path

APP_ID = 'org.bluemolecule.BlueMolecule'
VERSION = '1.4.1'
BINDINGS = ('PyQt6', 'PySide6', 'PyQt5')
INLINE_LIMIT = 32 * 1024 * 1024     # files below this size go to the page in one piece
READ_LIMIT = 6 * 1024 * 1024        # bytes per bmRead reply (8 MiB of base64)
ICON_SIZES = (48, 64, 128, 256, 512)
EXTERNAL_SCHEMES = ('http', 'https', 'mailto', 'ftp')
HERE = Path(__file__).resolve().parent
WEBGL_CHECK = ("(function () { try { var c = document.createElement('canvas'); return !!c.getContext('webgl2'); }"
               " catch (e) { return false; } })()")


# ---------------------------------------------------------------- helpers (no Qt)

def data_dir():
    xdg = os.environ.get('XDG_DATA_HOME')
    return (Path(xdg) if xdg else Path.home() / '.local' / 'share') / 'bluemolecule'


def cache_dir():
    xdg = os.environ.get('XDG_CACHE_HOME')
    return (Path(xdg) if xdg else Path.home() / '.cache') / 'bluemolecule'


def find_launcher(path=None):
    """Path of the bluemolecule launcher shipped with this file (it runs Open Babel and writes the log)."""
    for p in ([Path(path)] if path else []) + [HERE / 'bluemolecule', HERE / 'bluemolecule.py',
                                                HERE.parent.parent / 'launcher' / 'bluemolecule.py']:
        if p.is_file():
            return p
    return None


def load_launcher(path=None):
    """The launcher as a module (its main() does not run), or None."""
    p = find_launcher(path)
    if p is None:
        return None
    try:
        sys.dont_write_bytecode = True    # no __pycache__ next to an installed launcher
        loader = importlib.machinery.SourceFileLoader('bluemolecule_launcher', str(p))
        spec = importlib.util.spec_from_loader('bluemolecule_launcher', loader)
        mod = importlib.util.module_from_spec(spec)
        loader.exec_module(mod)
        return mod
    except Exception as e:   # a broken launcher must not keep the window from opening
        print('bluemolecule: could not load the launcher %s: %s' % (p, e), file=sys.stderr)
        return None


def find_page(path=None):
    for p in ([Path(path)] if path else []) + [HERE / 'BlueMolecule.html', HERE.parent.parent / 'dist' / 'BlueMolecule.html']:
        if p.is_file():
            return p.resolve()
    return None


def icon_files(here=HERE):
    """PNG icons: icons/ next to this file (package), the hicolor theme of the install prefix, or dist/ (source tree)."""
    out = []
    for s in ICON_SIZES:
        for p in (here / 'icons' / ('%s-%d.png' % (APP_ID, s)),
                  here.parent / 'icons' / 'hicolor' / ('%dx%d' % (s, s)) / 'apps' / (APP_ID + '.png'),
                  here.parent.parent / 'dist' / 'build' / 'icons' / ('icon_%d.png' % s)):
            if p.is_file():
                out.append(p)
                break
    return out


def js_call(body, args):
    """JavaScript that runs `body` with `m` bound to args (JSON, ASCII-only, so any string is safe)."""
    data = json.dumps(args, ensure_ascii=True, separators=(',', ':'))
    return "(function (m) { try { %s; } catch (e) { console.error('BlueMolecule bridge: ' + e); } })(%s);" % (body, data)


def safe_file_name(name):
    """A file name proposed by the page, without folders or control characters."""
    name = str(name or '').replace('\\', '/').split('/')[-1]
    name = ''.join(c for c in name if ord(c) >= 32 and c != '\x7f').strip()
    if name in ('', '.', '..'):
        name = 'untitled'
    return name[:200]


def write_atomic(path, data):
    """Write bytes so that an existing file is replaced only once the new one is complete."""
    target = os.path.realpath(path)
    folder = os.path.dirname(target) or '.'
    fd, tmp = tempfile.mkstemp(prefix='.' + os.path.basename(target) + '.', suffix='.part', dir=folder)
    try:
        with os.fdopen(fd, 'wb') as fh:
            fh.write(data)
        try:
            os.chmod(tmp, os.stat(target).st_mode & 0o7777 if os.path.exists(target) else 0o666 & ~current_umask())
        except OSError:
            pass
        os.replace(tmp, target)
    except BaseException:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise
    return target


def current_umask():
    mask = os.umask(0o22)
    os.umask(mask)
    return mask


def is_count(x):
    return isinstance(x, int) and not isinstance(x, bool) and x >= 0


class Bridge:
    """The page's messages and the replies to them, independent of Qt.

    run_js(script) runs JavaScript in the page. ui provides ask_save_path(suggested, mime, done(path|None)),
    set_clipboard(text), background(work, done(result)) and error(title, text). obabel(request) -> reply is
    the launcher's obabel_request; log(text) writes a log line.
    """

    def __init__(self, run_js, ui, obabel=None, log=None, files=()):
        self.run_js, self.ui, self.obabel = run_js, ui, obabel
        self.log = log or (lambda text: None)
        self.initial = [Path(f) for f in files]
        self.pending = list(self.initial)
        self.ready = False
        self.large = {}            # fileId -> (path, size): the only files bmRead can read
        self.last_dir = self.initial[0].parent if self.initial else None
        self.handlers = {'bmReady': self.on_ready, 'bmSave': self.on_save, 'bmObabel': self.on_obabel,
                         'bmCopy': self.on_copy, 'bmRead': self.on_read, 'bmLog': self.on_log, 'log': self.on_log}

    def receive(self, name, payload):
        """One message from the page: its name and its JSON text."""
        handler = self.handlers.get(name)
        if handler is None:
            self.log('bridge: unknown message %r' % (name,))
            return False
        try:
            msg = json.loads(payload) if payload else {}
        except ValueError:
            self.log('bridge: %s carried invalid JSON' % name)
            return False
        handler(msg if isinstance(msg, dict) else {})
        return True

    def resolve(self, id, result):
        self.run_js(js_call('BM.native._resolve(m.id, m.result)', {'id': id, 'result': result}))

    # ---- documents
    def on_ready(self, msg):
        self.ready = True
        paths, self.pending = self.pending, []
        for p in paths:
            self.deliver(p)

    def open_paths(self, paths):
        """More documents (e.g. a file dropped outside the page's drop area): delivered once the page is ready."""
        self.pending += [Path(p) for p in paths]
        if self.ready:
            self.on_ready({})

    def page_reloaded(self):
        """After a reload (or a crashed page) the page announces itself again and gets the start-up files again."""
        self.ready = False
        self.pending = list(self.initial) + [p for p in self.pending if p not in self.initial]

    def deliver(self, path):
        try:
            st = path.stat()
            if not path.is_file():
                raise OSError('not a regular file')
            if st.st_size < INLINE_LIMIT:
                data = path.read_bytes()
        except OSError as e:
            self.log('could not open %s: %s' % (path, e))
            self.ui.error('Could not open %s' % path.name, str(e))
            return
        self.last_dir = path.parent
        if st.st_size < INLINE_LIMIT:
            self.run_js(js_call('BM.native._receiveFile(m.name, m.data)',
                                {'name': path.name, 'data': base64.b64encode(data).decode('ascii')}))
        else:
            fid = 'f%d' % (len(self.large) + 1)
            self.large[fid] = (path, st.st_size)
            self.run_js(js_call('BM.native._receiveLargeFile(m.name, m.size, m.fileId)',
                                {'name': path.name, 'size': st.st_size, 'fileId': fid}))

    def on_read(self, msg):
        id = msg.get('id')
        if not isinstance(id, str):
            return
        entry = self.large.get(msg.get('fileId'))
        if entry is None:
            return self.resolve(id, {'ok': False, 'error': 'Unknown file.'})
        offset, length = msg.get('offset'), msg.get('length')
        if not (is_count(offset) and is_count(length)):
            return self.resolve(id, {'ok': False, 'error': 'Bad read request (offset and length must be whole numbers >= 0).'})
        path, size = entry
        try:
            with open(path, 'rb') as fh:
                fh.seek(offset)
                data = fh.read(min(length, READ_LIMIT))
        except OSError as e:
            return self.resolve(id, {'ok': False, 'error': 'Could not read %s: %s' % (path.name, e)})
        self.resolve(id, {'ok': True, 'offset': offset, 'length': len(data), 'eof': offset + len(data) >= size,
                          'base64': base64.b64encode(data).decode('ascii')})

    # ---- save, copy, Open Babel, log
    def on_save(self, msg):
        id = msg.get('id')
        if not isinstance(id, str):
            return
        try:
            data = base64.b64decode(str(msg.get('base64') or ''))
        except (binascii.Error, ValueError):
            return self.resolve(id, {'ok': False, 'error': 'The data to save could not be decoded.'})
        name = safe_file_name(msg.get('name'))

        def chosen(path):
            if not path:
                return self.resolve(id, {'ok': False, 'cancelled': True})
            try:
                written = write_atomic(path, data)
            except OSError as e:
                self.log('could not save %s: %s' % (path, e))
                self.ui.error('Could not save %s' % os.path.basename(path), str(e))
                return self.resolve(id, {'ok': False, 'error': str(e)})
            self.last_dir = Path(written).parent
            self.resolve(id, {'ok': True, 'path': written})
        self.ui.ask_save_path(str((self.last_dir or Path.home()) / name), str(msg.get('mime') or ''), chosen)

    def on_copy(self, msg):
        self.ui.set_clipboard(str(msg.get('text') or ''))

    def on_obabel(self, msg):
        id = msg.get('id')
        if not isinstance(id, str):
            return
        req = {k: v for k, v in msg.items() if k != 'id'}
        if self.obabel is None:
            return self.resolve(id, {'ok': False, 'output': '', 'error': 'Open Babel support (the bluemolecule launcher) was not found.'})
        self.ui.background(lambda: self.obabel(req), lambda result: self.resolve(id, result))

    def on_log(self, msg):
        self.log('[page %s] %s' % (str(msg.get('level') or 'log')[:20], str(msg.get('message') or '')[:8000]))


# ---------------------------------------------------------------- Qt

def E(obj, scope, name):
    """A Qt enum member, scoped (PyQt6, PySide6: Class.Scope.Name) or not (PyQt5: Class.Name); None if missing."""
    group = getattr(obj, scope, None)
    member = getattr(group, name, None) if group is not None else None
    return member if member is not None else getattr(obj, name, None)


def import_qt(binding):
    """The classes the window needs from one binding (raises ImportError and the like)."""
    def mod(sub):
        return importlib.import_module(binding + '.' + sub)
    core, gui, widgets = mod('QtCore'), mod('QtGui'), mod('QtWidgets')
    web_widgets = mod('QtWebEngineWidgets')       # QtWebEngine must be imported before the QApplication exists
    web = web_widgets if binding == 'PyQt5' else mod('QtWebEngineCore')
    channel = mod('QtWebChannel')
    delete = None
    try:
        if binding == 'PySide6':
            delete = importlib.import_module('shiboken6').delete
        else:
            try:
                delete = mod('sip').delete
            except ImportError:
                delete = importlib.import_module('sip').delete
    except (ImportError, AttributeError):
        pass
    version = getattr(core, 'PYQT_VERSION_STR', None) or getattr(importlib.import_module(binding), '__version__', '?')
    return types.SimpleNamespace(
        binding=binding, binding_version=version, qt_version=core.qVersion(), delete=delete,
        QObject=core.QObject, QTimer=core.QTimer, QUrl=core.QUrl, QFile=core.QFile, QIODevice=core.QIODevice,
        QByteArray=core.QByteArray, Qt=core.Qt, QCoreApplication=core.QCoreApplication,
        Slot=getattr(core, 'pyqtSlot', None) or core.Slot, Signal=getattr(core, 'pyqtSignal', None) or core.Signal,
        QIcon=gui.QIcon, QKeySequence=gui.QKeySequence, QDesktopServices=gui.QDesktopServices,
        QGuiApplication=gui.QGuiApplication, QShortcut=getattr(gui, 'QShortcut', None) or widgets.QShortcut,
        QApplication=widgets.QApplication, QFileDialog=widgets.QFileDialog, QMessageBox=widgets.QMessageBox,
        QWebEngineView=web_widgets.QWebEngineView, QWebEnginePage=web.QWebEnginePage,
        QWebEngineProfile=web.QWebEngineProfile, QWebEngineSettings=web.QWebEngineSettings,
        QWebEngineScript=web.QWebEngineScript, QWebChannel=channel.QWebChannel)


def build_classes(Q):
    """Qt subclasses, made once the binding is known."""

    class Relay(Q.QObject):
        """Runs callables on the GUI thread; worker threads (Open Babel) emit them."""
        call = Q.Signal(object)

        def __init__(self):
            super().__init__()
            self.call.connect(self._run)

        @Q.Slot(object)
        def _run(self, fn):
            if callable(fn):
                fn()

    class Host(Q.QObject):
        """channel.objects.bmHost in the page; post(name, json) is its only entry point."""

        def __init__(self, receive):
            super().__init__()
            self._receive = receive

        @Q.Slot(str, str)
        def post(self, name, payload):
            self._receive(str(name), str(payload))

    class Page(Q.QWebEnginePage):
        def __init__(self, profile, parent, shell):
            super().__init__(profile, parent)
            self._shell = shell

        def javaScriptConsoleMessage(self, level, message, line, source):
            self._shell.console(level, message, line, source)

        def acceptNavigationRequest(self, url, nav_type, is_main_frame):
            return self._shell.navigation(url, is_main_frame)

        def createWindow(self, window_type):
            return self._shell.popup()

    class ExternalPage(Q.QWebEnginePage):
        """Target of window.open() and target=_blank links: web links open in the user's browser."""

        def acceptNavigationRequest(self, url, nav_type, is_main_frame):
            if url.scheme().lower() in EXTERNAL_SCHEMES:
                Q.QDesktopServices.openUrl(url)
                self.deleteLater()
                return False
            return url.scheme().lower() == 'about'

    class Window(Q.QWebEngineView):
        def __init__(self, shell):
            super().__init__()
            self._shell = shell

        def closeEvent(self, event):
            self._shell.closing()
            super().closeEvent(event)

    return types.SimpleNamespace(Relay=Relay, Host=Host, Page=Page, ExternalPage=ExternalPage, Window=Window)


class QtUI:
    """What the Bridge needs from Qt: dialogs, the clipboard and a worker thread."""

    def __init__(self, shell):
        self.shell = shell

    def ask_save_path(self, suggested, mime, done):
        Q, shell = self.shell.Q, self.shell

        def run():   # later, from the event loop, not inside the QWebChannel call
            try:
                path = Q.QFileDialog.getSaveFileName(shell.window, 'Save', suggested)[0]
            except Exception as e:
                shell.log('Save dialog failed: %s' % e)
                path = ''
            done(path or None)
        Q.QTimer.singleShot(0, run)

    def set_clipboard(self, text):
        self.shell.Q.QGuiApplication.clipboard().setText(text)

    def background(self, work, done):
        relay = self.shell.relay

        def worker():
            try:
                result = work()
            except Exception as e:
                result = {'ok': False, 'output': '', 'error': '%s: %s' % (type(e).__name__, e)}
            relay.call.emit(lambda: done(result))
        threading.Thread(target=worker, name='bluemolecule-worker', daemon=True).start()

    def error(self, title, text):
        Q, shell = self.shell.Q, self.shell
        Q.QTimer.singleShot(0, lambda: Q.QMessageBox.warning(shell.window, title, text))


class Shell:
    """The window: application identity, profile, page, channel and Bridge."""

    def __init__(self, Q, html, files=(), launcher=None):
        self.Q, self.html, self.files, self.L = Q, Path(html), list(files), launcher
        self.C = build_classes(Q)
        self.app = self.profile = self.window = self.page = self.dev = None
        self.crashes = 0
        self.shortcuts = []

    def log(self, text):
        if self.L is not None and hasattr(self.L, 'log_event'):
            self.L.log_event(text)
        else:
            print('bluemolecule: ' + text, file=sys.stderr)

    def run(self):
        Q = self.Q
        os.environ.setdefault('RESOURCE_NAME', APP_ID)        # X11/XWayland: the WM_CLASS instance name
        for name in ('AA_ShareOpenGLContexts',) + (('AA_EnableHighDpiScaling', 'AA_UseHighDpiPixmaps') if Q.binding == 'PyQt5' else ()):
            attr = E(Q.Qt, 'ApplicationAttribute', name)
            if attr is not None:
                Q.QCoreApplication.setAttribute(attr, True)
        Q.QCoreApplication.setApplicationName('BlueMolecule')
        Q.QCoreApplication.setOrganizationName('BlueMolecule')
        Q.QCoreApplication.setOrganizationDomain('bluemolecule.org')
        Q.QCoreApplication.setApplicationVersion(VERSION)
        self.app = app = Q.QApplication(['bluemolecule'])
        Q.QGuiApplication.setApplicationDisplayName('BlueMolecule')
        Q.QGuiApplication.setDesktopFileName(APP_ID)       # the Wayland app_id: KWin matches it to APP_ID.desktop
        self.icon = Q.QIcon()
        for p in icon_files():
            self.icon.addFile(str(p))
        if self.icon.isNull():
            self.icon = Q.QIcon.fromTheme(APP_ID)
        Q.QGuiApplication.setWindowIcon(self.icon)

        # a persistent profile of its own: localStorage (the viewer's settings) survives restarts
        self.profile = profile = Q.QWebEngineProfile('BlueMolecule')
        profile.setPersistentStoragePath(str(data_dir() / 'qtwebengine'))
        profile.setCachePath(str(cache_dir() / 'qtwebengine'))
        settings = profile.settings()
        for name in ('WebGLEnabled', 'Accelerated2dCanvasEnabled', 'LocalContentCanAccessFileUrls',
                     'LocalStorageEnabled', 'JavascriptCanAccessClipboard'):
            attr = E(Q.QWebEngineSettings, 'WebAttribute', name)
            if attr is not None:
                settings.setAttribute(attr, True)
        profile.downloadRequested.connect(self.download)

        self.window = win = self.C.Window(self)
        win.setWindowTitle('BlueMolecule')
        win.setWindowIcon(self.icon)
        self.page = page = self.C.Page(profile, win, self)
        self.relay = self.C.Relay()
        self.bridge = Bridge(lambda script: self.page.runJavaScript(script), QtUI(self),
                             obabel=getattr(self.L, 'obabel_request', None), log=self.log, files=self.files)
        self.host = self.C.Host(self.bridge.receive)
        self.channel = Q.QWebChannel(page)
        self.channel.registerObject('bmHost', self.host)
        page.setWebChannel(self.channel)
        self.inject('qwebchannel', self.qwebchannel_js())
        marker = {'shell': 'linux-app', 'version': VERSION, 'binding': Q.binding, 'qt': Q.qt_version,
                  'readLimit': READ_LIMIT, 'inlineLimit': INLINE_LIMIT}
        self.inject('bm-qt', 'window.BM_QT = Object.freeze(%s);' % json.dumps(marker))
        page.titleChanged.connect(lambda title: win.setWindowTitle(title or 'BlueMolecule'))
        page.loadFinished.connect(self.loaded)
        page.renderProcessTerminated.connect(self.crashed)
        for name in ('Back', 'Forward', 'Reload', 'ReloadAndBypassCache', 'ViewSource', 'SavePage'):
            action = E(Q.QWebEnginePage, 'WebAction', name)    # not in the right-click menu of an app window
            if action is not None:
                try:
                    page.action(action).setVisible(False)
                except Exception:
                    pass
        win.setPage(page)
        for keys in ('F12', 'Ctrl+Shift+I'):
            shortcut = Q.QShortcut(Q.QKeySequence(keys), win)
            shortcut.activated.connect(self.toggle_devtools)
            self.shortcuts.append(shortcut)
        self.restore_geometry()
        win.setUrl(Q.QUrl.fromLocalFile(str(self.html)))
        win.show()
        run = getattr(app, 'exec', None) or getattr(app, 'exec_')
        rc = run()
        self.teardown()
        return rc if isinstance(rc, int) else 0

    def inject(self, name, source):
        """A script that runs before the page's own scripts, in the page's JavaScript world."""
        Q = self.Q
        s = Q.QWebEngineScript()
        s.setName(name)
        s.setSourceCode(source)
        s.setInjectionPoint(E(Q.QWebEngineScript, 'InjectionPoint', 'DocumentCreation'))
        s.setWorldId(0)                  # MainWorld, where native-bridge.js runs
        s.setRunsOnSubFrames(False)
        self.page.scripts().insert(s)

    def qwebchannel_js(self):
        Q = self.Q
        f = Q.QFile(':/qtwebchannel/qwebchannel.js')
        if f.open(E(Q.QIODevice, 'OpenModeFlag', 'ReadOnly')):
            try:
                return bytes(f.readAll().data()).decode('utf-8')
            finally:
                f.close()
        self.log('qwebchannel.js is missing from Qt WebChannel: saving, Open Babel and command-line files will not work')
        return ''

    # ---- page events
    def console(self, level, message, line, source):
        P = self.Q.QWebEnginePage
        if level == E(P, 'JavaScriptConsoleMessageLevel', 'InfoMessageLevel') and not os.environ.get('BLUEMOLECULE_DEBUG'):
            return
        kind = 'error' if level == E(P, 'JavaScriptConsoleMessageLevel', 'ErrorMessageLevel') else 'warning'
        self.log('[page %s] %s (%s:%s)' % (kind, message, str(source).rsplit('/', 1)[-1], line))

    def navigation(self, url, is_main_frame):
        scheme = url.scheme().lower()
        if scheme in EXTERNAL_SCHEMES:
            self.Q.QDesktopServices.openUrl(url)
            return False
        if not is_main_frame:
            return True
        if scheme == 'file':
            path = url.toLocalFile()
            if path and os.path.realpath(path) != os.path.realpath(str(self.html)):
                # a file dropped outside the page's drop area: open it in the viewer instead of leaving the page
                if os.path.isfile(path):
                    self.bridge.open_paths([path])
                return False
            return True
        return False          # data:, blob:, about: ... would replace the viewer

    def popup(self):
        return self.C.ExternalPage(self.profile, self.page)

    def loaded(self, ok):
        if not ok:
            self.log('the page did not load: %s' % self.html)
            return
        self.page.runJavaScript(WEBGL_CHECK, 0, lambda r: self.log(
            'WebGL 2 available' if r else 'WebGL 2 is NOT available in Qt WebEngine (try: bluemolecule --tier chromium)'))

    def crashed(self, status, code):
        if status == E(self.Q.QWebEnginePage, 'RenderProcessTerminationStatus', 'NormalTerminationStatus'):
            return
        self.crashes += 1
        self.log('the page process stopped (status %s, exit code %s)' % (status, code))
        self.bridge.page_reloaded()
        if self.crashes == 1:
            self.Q.QTimer.singleShot(300, self.window.reload)
        else:
            self.Q.QMessageBox.warning(self.window, 'BlueMolecule stopped',
                                       'The viewer stopped again. Reload it from a terminal with\n\n'
                                       '    bluemolecule --tier chromium\n\nto use a browser window instead; '
                                       'bluemolecule --diagnose shows details.')

    def download(self, item):
        """Safety net for downloads the page starts itself (<a download> outside BM.native.save)."""
        Q = self.Q
        name = item.downloadFileName() if hasattr(item, 'downloadFileName') else os.path.basename(item.path())
        start = str((self.bridge.last_dir or Path.home()) / safe_file_name(name))
        path = Q.QFileDialog.getSaveFileName(self.window, 'Save', start)[0]
        if not path:
            item.cancel()
            return
        if hasattr(item, 'setDownloadFileName'):
            item.setDownloadDirectory(os.path.dirname(path))
            item.setDownloadFileName(os.path.basename(path))
        else:
            item.setPath(path)
        item.accept()
        self.bridge.last_dir = Path(path).parent

    def toggle_devtools(self):
        Q = self.Q
        if self.dev is None:
            self.dev = Q.QWebEngineView()
            self.dev.setWindowTitle('BlueMolecule developer tools')
            self.dev.setWindowIcon(self.icon)
            self.dev.resize(1100, 720)
            self.page.setDevToolsPage(self.dev.page())
        if self.dev.isVisible():
            self.dev.hide()
        else:
            self.dev.show()
            self.dev.raise_()
            self.dev.activateWindow()

    # ---- window size, shutdown
    def restore_geometry(self):
        ok = False
        try:
            state = json.loads((data_dir() / 'qt-window.json').read_text(encoding='utf-8'))
            raw = self.Q.QByteArray(str(state['geometry']).encode('ascii'))
            ok = bool(self.window.restoreGeometry(self.Q.QByteArray.fromBase64(raw)))
        except Exception:
            pass
        if not ok:
            self.window.resize(1280, 840)

    def closing(self):
        try:
            geometry = bytes(self.window.saveGeometry().toBase64().data()).decode('ascii')
            p = data_dir() / 'qt-window.json'
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text(json.dumps({'geometry': geometry}), encoding='utf-8')
        except Exception as e:
            self.log('could not save the window size: %s' % e)
        if self.dev is not None:
            self.dev.close()
        self.Q.QTimer.singleShot(0, self.app.quit)

    def teardown(self):
        """Delete the views (and so their pages) before the profile, so that it is flushed to disk cleanly."""
        delete = self.Q.delete
        for obj in (self.dev, self.window, self.profile):
            if obj is not None and delete is not None:
                try:
                    delete(obj)
                except Exception:
                    pass
        self.dev = self.window = self.page = self.profile = None


# ---------------------------------------------------------------- start-up

def fall_back(opts, reason):
    """Start the launcher again without this Python, so that it chooses a browser window (or another Python)."""
    msg = 'the Qt window cannot start with %s: %s' % (sys.executable, reason)
    launcher = find_launcher(opts.launcher)
    L = load_launcher(opts.launcher)
    if L is not None and hasattr(L, 'log_event'):
        L.log_event(msg)
    else:
        print('bluemolecule: ' + msg, file=sys.stderr)
    if opts.no_fallback or launcher is None:
        if launcher is None:
            print('bluemolecule: the launcher was not found next to %s' % __file__, file=sys.stderr)
        return 1
    env = dict(os.environ)
    skip = [p for p in env.get('BLUEMOLECULE_QT_SKIP', '').split(os.pathsep) if p] + [sys.executable]
    env['BLUEMOLECULE_QT_SKIP'] = os.pathsep.join(skip)
    py = opts.launcher_python or sys.executable
    try:
        os.execve(py, [py, str(launcher)] + opts.files, env)
    except OSError as e:
        print('bluemolecule: cannot start the launcher %s: %s' % (launcher, e), file=sys.stderr)
    return 1


def parse_args(argv=None):
    ap = argparse.ArgumentParser(prog='bluemolecule_qt', description='BlueMolecule in a Qt WebEngine window.')
    ap.add_argument('files', nargs='*')
    ap.add_argument('--binding', choices=BINDINGS, help='Qt binding to try first (default: PyQt6, PySide6, PyQt5)')
    ap.add_argument('--html', help='BlueMolecule.html (default: next to this file)')
    ap.add_argument('--launcher', help='the bluemolecule launcher (Open Babel, log, fallback)')
    ap.add_argument('--launcher-python', help='Python that runs the launcher when falling back')
    ap.add_argument('--no-fallback', action='store_true', help='exit instead of falling back to a browser window')
    opts = ap.parse_args(argv)
    opts.files = [os.path.abspath(f) for f in opts.files]
    return opts


def main(argv=None):
    opts = parse_args(argv)
    order = [opts.binding] + [b for b in BINDINGS if b != opts.binding] if opts.binding else list(BINDINGS)
    Q, errors = None, []
    for b in order:
        try:
            Q = import_qt(b)
            break
        except Exception as e:       # ImportError, or a Qt library that does not load
            errors.append('%s: %s' % (b, e))
            if b + '.QtCore' in sys.modules:
                break                # one Qt per process: a second binding could clash with the first
    if Q is None:
        return fall_back(opts, '; '.join(errors))
    html = find_page(opts.html)
    if html is None:
        print('bluemolecule: BlueMolecule.html was not found next to %s' % __file__, file=sys.stderr)
        return 1
    L = load_launcher(opts.launcher)
    shell = Shell(Q, html, [f for f in opts.files], L)
    shell.log('Qt window: %s %s, Qt %s, Python %s (%s), %d file(s)' % (
        Q.binding, Q.binding_version, Q.qt_version, '%d.%d.%d' % sys.version_info[:3], sys.executable, len(opts.files)))
    return shell.run()


if __name__ == '__main__':
    sys.exit(main())
