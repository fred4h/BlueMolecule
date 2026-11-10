// BlueMolecule for macOS: a small native window around BlueMolecule.html (WKWebView).
//
// Built on the user's Mac by "Install BlueMolecule.command":
//     swiftc -O -swift-version 5 -target arm64-apple-macos11.0 BlueMolecule.swift -o BlueMolecule
// One file, AppKit + WebKit only, macOS 11 or later. The page (src/native-bridge.js) talks to it
// through these WKScriptMessage handlers:
//   bmReady  {version, app}                        the page is ready; queued documents are delivered
//   bmSave   {id, name, mime, base64}              Save dialog; reply BM.native._resolve(id, {ok, path})
//   bmObabel {id, op, input, inFormat, outFormat, options}   Open Babel; reply {ok, output, error}
//   bmCopy   {text}                                plain text to the pasteboard
//   bmLog    {level, message}                      appended to ~/Library/Logs/BlueMolecule.log
// Documents opened from Finder, the Dock or `open -a` reach the page as
// BM.native._receiveFile(name, base64).

import Cocoa
import WebKit

let kHandlerNames = ["bmReady", "bmSave", "bmObabel", "bmCopy", "bmLog"]
let kWindowAutosaveName = "BlueMoleculeMainWindow"

/// Runs before the page's own scripts: forwards uncaught errors and console.error to bmLog.
let kErrorForwardingJS = #"""
(function () {
  var h = window.webkit && window.webkit.messageHandlers && window.webkit.messageHandlers.bmLog;
  if (!h) return;
  function str(a) {
    try {
      if (a instanceof Error) return a.name + ': ' + a.message + (a.stack ? '\n' + a.stack : '');
      if (a && typeof a === 'object') return JSON.stringify(a);
      return String(a);
    } catch (e) { return String(a); }
  }
  function post(level, text) {
    try { h.postMessage({ level: level, message: String(text).slice(0, 8000) }); } catch (e) {}
  }
  window.addEventListener('error', function (e) {
    var where = e.filename ? ' (' + String(e.filename).split('/').pop() + ':' + e.lineno + ':' + e.colno + ')' : '';
    post('error', (e.message || 'Script error') + where + (e.error && e.error.stack ? '\n' + e.error.stack : ''));
  });
  window.addEventListener('unhandledrejection', function (e) { post('error', 'Unhandled rejection: ' + str(e.reason)); });
  var original = console.error;
  console.error = function () {
    post('console.error', Array.prototype.map.call(arguments, str).join(' '));
    return original.apply(console, arguments);
  };
})();
"""#

final class AppDelegate: NSObject, NSApplicationDelegate, WKNavigationDelegate, WKUIDelegate, WKScriptMessageHandler {
    private var window: NSWindow!
    private var webView: WKWebView!
    private var titleObserver: NSKeyValueObservation?
    private var pageReady = false
    private var pendingURLs: [URL] = []
    private let pageURL: URL? = Bundle.main.url(forResource: "BlueMolecule", withExtension: "html")

    // MARK: - Application life cycle

    func applicationWillFinishLaunching(_ notification: Notification) {
        NSWindow.allowsAutomaticWindowTabbing = false
        buildMainMenu()
    }

    func applicationDidFinishLaunching(_ notification: Notification) {
        createWindow()
        loadPage()
        window.makeKeyAndOrderFront(nil)
        _ = window.makeFirstResponder(webView)
        NSApp.activate(ignoringOtherApps: true)
    }

    func applicationShouldTerminateAfterLastWindowClosed(_ sender: NSApplication) -> Bool {
        return true
    }

    func applicationSupportsSecureRestorableState(_ app: NSApplication) -> Bool {
        return true
    }

    /// Finder, the Dock and `open -a BlueMolecule file` deliver documents here (also before launch finishes).
    func application(_ application: NSApplication, open urls: [URL]) {
        openURLs(urls)
    }

    // MARK: - Window and page

    private func createWindow() {
        let content = WKUserContentController()
        for name in kHandlerNames {
            content.add(self, name: name)
        }
        content.addUserScript(WKUserScript(source: kErrorForwardingJS, injectionTime: .atDocumentStart, forMainFrameOnly: true))
        let config = WKWebViewConfiguration()
        config.userContentController = content

        let frame = NSRect(x: 0, y: 0, width: 1280, height: 820)
        webView = WKWebView(frame: frame, configuration: config)
        webView.navigationDelegate = self
        webView.uiDelegate = self
        webView.autoresizingMask = [.width, .height]
        #if compiler(>=5.8)
        if #available(macOS 13.3, *) {
            webView.isInspectable = true
        }
        #endif

        window = NSWindow(contentRect: frame,
                          styleMask: [.titled, .closable, .miniaturizable, .resizable],
                          backing: .buffered,
                          defer: false)
        window.title = "BlueMolecule"
        window.minSize = NSSize(width: 720, height: 480)
        window.isReleasedWhenClosed = false
        window.tabbingMode = .disallowed
        window.contentView = webView
        window.center()
        _ = window.setFrameAutosaveName(kWindowAutosaveName)
        _ = window.setFrameUsingName(kWindowAutosaveName)

        titleObserver = webView.observe(\.title, options: [.new]) { [weak self] view, _ in
            let title = view.title ?? ""
            self?.window?.title = title.isEmpty ? "BlueMolecule" : title
        }
    }

    private func loadPage() {
        pageReady = false
        if let url = pageURL {
            _ = webView.loadFileURL(url, allowingReadAccessTo: url.deletingLastPathComponent())
        } else {
            let html = "<body style='font:15px -apple-system;padding:2em'>BlueMolecule.html is missing from the app. "
                + "Please run “Install BlueMolecule.command” again.</body>"
            _ = webView.loadHTMLString(html, baseURL: nil)
        }
    }

    // MARK: - Menus

    private func buildMainMenu() {
        let bar = NSMenu()

        let appMenu = NSMenu(title: "BlueMolecule")
        appMenu.addItem(menuItem("About BlueMolecule", #selector(NSApplication.orderFrontStandardAboutPanel(_:)), ""))
        appMenu.addItem(NSMenuItem.separator())
        let servicesMenu = NSMenu(title: "Services")
        let services = NSMenuItem(title: "Services", action: nil, keyEquivalent: "")
        services.submenu = servicesMenu
        NSApp.servicesMenu = servicesMenu
        appMenu.addItem(services)
        appMenu.addItem(NSMenuItem.separator())
        appMenu.addItem(menuItem("Hide BlueMolecule", #selector(NSApplication.hide(_:)), "h"))
        appMenu.addItem(menuItem("Hide Others", #selector(NSApplication.hideOtherApplications(_:)), "h", [.command, .option]))
        appMenu.addItem(menuItem("Show All", #selector(NSApplication.unhideAllApplications(_:)), ""))
        appMenu.addItem(NSMenuItem.separator())
        appMenu.addItem(menuItem("Quit BlueMolecule", #selector(NSApplication.terminate(_:)), "q"))
        addSubmenu(appMenu, to: bar)

        let fileMenu = NSMenu(title: "File")
        fileMenu.addItem(menuItem("Open…", #selector(showOpenPanel(_:)), "o", target: self))
        fileMenu.addItem(NSMenuItem.separator())
        fileMenu.addItem(menuItem("Close Window", #selector(NSWindow.performClose(_:)), "w"))
        addSubmenu(fileMenu, to: bar)

        let editMenu = NSMenu(title: "Edit")
        editMenu.addItem(menuItem("Undo", Selector(("undo:")), "z"))
        editMenu.addItem(menuItem("Redo", Selector(("redo:")), "z", [.command, .shift]))
        editMenu.addItem(NSMenuItem.separator())
        editMenu.addItem(menuItem("Cut", #selector(NSText.cut(_:)), "x"))
        editMenu.addItem(menuItem("Copy", #selector(NSText.copy(_:)), "c"))
        editMenu.addItem(menuItem("Paste", #selector(NSText.paste(_:)), "v"))
        editMenu.addItem(menuItem("Delete", #selector(NSText.delete(_:)), ""))
        editMenu.addItem(menuItem("Select All", #selector(NSText.selectAll(_:)), "a"))
        addSubmenu(editMenu, to: bar)

        let viewMenu = NSMenu(title: "View")
        viewMenu.addItem(menuItem("Actual Size", #selector(zoomReset(_:)), "0", target: self))
        viewMenu.addItem(menuItem("Zoom In", #selector(zoomIn(_:)), "+", target: self))
        viewMenu.addItem(menuItem("Zoom Out", #selector(zoomOut(_:)), "-", target: self))
        viewMenu.addItem(NSMenuItem.separator())
        viewMenu.addItem(menuItem("Enter Full Screen", #selector(NSWindow.toggleFullScreen(_:)), "f", [.command, .control]))
        addSubmenu(viewMenu, to: bar)

        let windowMenu = NSMenu(title: "Window")
        windowMenu.addItem(menuItem("Minimize", #selector(NSWindow.performMiniaturize(_:)), "m"))
        windowMenu.addItem(menuItem("Zoom", #selector(NSWindow.performZoom(_:)), ""))
        windowMenu.addItem(NSMenuItem.separator())
        windowMenu.addItem(menuItem("Bring All to Front", #selector(NSApplication.arrangeInFront(_:)), ""))
        addSubmenu(windowMenu, to: bar)
        NSApp.windowsMenu = windowMenu

        let helpMenu = NSMenu(title: "Help")
        helpMenu.addItem(menuItem("BlueMolecule Help", #selector(openHelp(_:)), "?", target: self))
        addSubmenu(helpMenu, to: bar)
        NSApp.helpMenu = helpMenu

        NSApp.mainMenu = bar
    }

    private func menuItem(_ title: String, _ action: Selector?, _ key: String,
                          _ modifiers: NSEvent.ModifierFlags = [.command], target: AnyObject? = nil) -> NSMenuItem {
        let item = NSMenuItem(title: title, action: action, keyEquivalent: key)
        item.keyEquivalentModifierMask = modifiers
        if let target = target {
            item.target = target
        }
        return item
    }

    private func addSubmenu(_ menu: NSMenu, to bar: NSMenu) {
        let item = NSMenuItem(title: menu.title, action: nil, keyEquivalent: "")
        item.submenu = menu
        bar.addItem(item)
    }

    @objc func showOpenPanel(_ sender: Any?) {
        let panel = NSOpenPanel()
        panel.allowsMultipleSelection = true
        panel.canChooseDirectories = false
        panel.canChooseFiles = true
        panel.beginSheetModal(for: window) { [weak self] response in
            if response == .OK {
                self?.openURLs(panel.urls)
            }
        }
    }

    @objc func zoomIn(_ sender: Any?) { webView.pageZoom = min(webView.pageZoom + 0.1, 2.0) }
    @objc func zoomOut(_ sender: Any?) { webView.pageZoom = max(webView.pageZoom - 0.1, 0.5) }
    @objc func zoomReset(_ sender: Any?) { webView.pageZoom = 1.0 }

    @objc func openHelp(_ sender: Any?) {
        webView.evaluateJavaScript("(function () { var b = document.getElementById('btn-help'); if (b) b.click(); })();") { _, _ in }
    }

    // MARK: - Documents

    private func openURLs(_ urls: [URL]) {
        for url in urls where url.isFileURL {
            pendingURLs.append(url)
            NSDocumentController.shared.noteNewRecentDocumentURL(url)
        }
        deliverPending()
    }

    /// Hands queued documents to the page once it has posted bmReady.
    private func deliverPending() {
        guard pageReady, webView != nil else { return }
        let urls = pendingURLs
        pendingURLs.removeAll()
        for url in urls {
            do {
                let data = try Data(contentsOf: url)
                callJS("BM.native._receiveFile(m.name, m.data)", ["name": url.lastPathComponent, "data": data.base64EncodedString()])
            } catch {
                appendLog("could not read \(url.path): \(error.localizedDescription)")
                showAlert("Could not open \(url.lastPathComponent)", error.localizedDescription)
            }
        }
    }

    // MARK: - Messages from the page

    func userContentController(_ userContentController: WKUserContentController, didReceive message: WKScriptMessage) {
        let body = (message.body as? [String: Any]) ?? [:]
        switch message.name {
        case "bmReady":
            pageReady = true
            deliverPending()
        case "bmSave":
            handleSave(body)
        case "bmObabel":
            handleObabel(body)
        case "bmCopy":
            let pasteboard = NSPasteboard.general
            _ = pasteboard.clearContents()
            _ = pasteboard.setString((body["text"] as? String) ?? "", forType: .string)
        case "bmLog":
            let level = (body["level"] as? String) ?? "log"
            let text = (body["message"] as? String) ?? String(describing: message.body)
            appendLog("[\(level)] \(text)")
        default:
            break
        }
    }

    private func handleSave(_ body: [String: Any]) {
        guard let id = body["id"] as? String else { return }
        guard let data = Data(base64Encoded: (body["base64"] as? String) ?? "", options: .ignoreUnknownCharacters) else {
            resolve(id, ["ok": false, "error": "The data to save could not be decoded."])
            return
        }
        let panel = NSSavePanel()
        panel.nameFieldStringValue = safeFileName((body["name"] as? String) ?? "untitled")
        panel.canCreateDirectories = true
        panel.isExtensionHidden = false
        panel.beginSheetModal(for: window) { [weak self] response in
            guard let self = self else { return }
            guard response == .OK, let url = panel.url else {
                self.resolve(id, ["ok": false, "cancelled": true])
                return
            }
            do {
                try data.write(to: url, options: .atomic)
                self.resolve(id, ["ok": true, "path": url.path])
            } catch {
                self.appendLog("could not save \(url.path): \(error.localizedDescription)")
                self.showAlert("Could not save \(url.lastPathComponent)", error.localizedDescription)
                self.resolve(id, ["ok": false, "error": error.localizedDescription])
            }
        }
    }

    private func handleObabel(_ body: [String: Any]) {
        guard let id = body["id"] as? String else { return }
        DispatchQueue.global(qos: .userInitiated).async {
            let result = ObabelRunner.run(body)
            DispatchQueue.main.async { [weak self] in
                self?.resolve(id, result)
            }
        }
    }

    private func resolve(_ id: String, _ result: [String: Any]) {
        callJS("BM.native._resolve(m.id, m.result)", ["id": id, "result": result])
    }

    /// Runs `body` in the page with `m` bound to `args` (JSON-encoded, so any string is safe).
    private func callJS(_ body: String, _ args: [String: Any]) {
        guard JSONSerialization.isValidJSONObject(args),
              let json = try? JSONSerialization.data(withJSONObject: args, options: []),
              let text = String(data: json, encoding: .utf8) else {
            appendLog("could not encode a reply for the page")
            return
        }
        let script = "(function (m) { try { " + body + "; } catch (e) { console.error('BlueMolecule bridge: ' + e); } })(" + text + ");"
        webView.evaluateJavaScript(script) { [weak self] _, error in
            if let error = error {
                self?.appendLog("evaluateJavaScript failed: \(error.localizedDescription)")
            }
        }
    }

    // MARK: - WKUIDelegate

    /// <input type=file> in the page.
    func webView(_ webView: WKWebView, runOpenPanelWith parameters: WKOpenPanelParameters,
                 initiatedByFrame frame: WKFrameInfo, completionHandler: @escaping ([URL]?) -> Void) {
        let panel = NSOpenPanel()
        panel.allowsMultipleSelection = parameters.allowsMultipleSelection
        panel.canChooseDirectories = false
        panel.canChooseFiles = true
        panel.beginSheetModal(for: window) { response in
            completionHandler(response == .OK ? panel.urls : nil)
        }
    }

    func webView(_ webView: WKWebView, runJavaScriptAlertPanelWithMessage message: String,
                 initiatedByFrame frame: WKFrameInfo, completionHandler: @escaping () -> Void) {
        let alert = NSAlert()
        alert.messageText = message
        _ = alert.addButton(withTitle: "OK")
        alert.beginSheetModal(for: window) { _ in completionHandler() }
    }

    func webView(_ webView: WKWebView, runJavaScriptConfirmPanelWithMessage message: String,
                 initiatedByFrame frame: WKFrameInfo, completionHandler: @escaping (Bool) -> Void) {
        let alert = NSAlert()
        alert.messageText = message
        _ = alert.addButton(withTitle: "OK")
        _ = alert.addButton(withTitle: "Cancel")
        alert.beginSheetModal(for: window) { response in
            completionHandler(response == .alertFirstButtonReturn)
        }
    }

    /// window.open / target=_blank: web links go to the default browser; no second web view.
    func webView(_ webView: WKWebView, createWebViewWith configuration: WKWebViewConfiguration,
                 for navigationAction: WKNavigationAction, windowFeatures: WKWindowFeatures) -> WKWebView? {
        if let url = navigationAction.request.url, isExternal(url) {
            _ = NSWorkspace.shared.open(url)
        }
        return nil
    }

    // MARK: - WKNavigationDelegate

    /// Only BlueMolecule.html may be shown; web links open in the default browser and a file dropped
    /// on the window that the page did not take is opened as a document instead of replacing the page.
    func webView(_ webView: WKWebView, decidePolicyFor navigationAction: WKNavigationAction,
                 decisionHandler: @escaping (WKNavigationActionPolicy) -> Void) {
        guard let url = navigationAction.request.url else {
            decisionHandler(.cancel)
            return
        }
        let scheme = (url.scheme ?? "").lowercased()
        if let page = pageURL, url.isFileURL, url.standardizedFileURL.path == page.standardizedFileURL.path {
            decisionHandler(.allow)
            return
        }
        if scheme == "about" {
            decisionHandler(.allow)
            return
        }
        if (scheme == "blob" || scheme == "data") && navigationAction.targetFrame?.isMainFrame == false {
            decisionHandler(.allow)
            return
        }
        decisionHandler(.cancel)
        if url.isFileURL {
            openURLs([url])
        } else if isExternal(url) {
            _ = NSWorkspace.shared.open(url)
        }
    }

    func webView(_ webView: WKWebView, didStartProvisionalNavigation navigation: WKNavigation!) {
        pageReady = false
    }

    func webView(_ webView: WKWebView, didFailProvisionalNavigation navigation: WKNavigation!, withError error: Error) {
        appendLog("page failed to load: \(error.localizedDescription)")
    }

    func webViewWebContentProcessDidTerminate(_ webView: WKWebView) {
        appendLog("the web content process ended; reloading")
        loadPage()
    }

    // MARK: - Helpers

    private func isExternal(_ url: URL) -> Bool {
        let scheme = (url.scheme ?? "").lowercased()
        return scheme == "http" || scheme == "https" || scheme == "mailto"
    }

    private func safeFileName(_ name: String) -> String {
        let cleaned = name.replacingOccurrences(of: "/", with: "-")
            .replacingOccurrences(of: ":", with: "-")
            .trimmingCharacters(in: .whitespacesAndNewlines)
        return cleaned.isEmpty ? "untitled" : cleaned
    }

    private func showAlert(_ title: String, _ text: String) {
        let alert = NSAlert()
        alert.messageText = title
        alert.informativeText = text
        alert.alertStyle = .warning
        if let window = window, window.isVisible, window.attachedSheet == nil {
            alert.beginSheetModal(for: window, completionHandler: nil)
        } else {
            _ = alert.runModal()
        }
    }

    private func appendLog(_ text: String) {
        let url = URL(fileURLWithPath: NSHomeDirectory()).appendingPathComponent("Library/Logs/BlueMolecule.log")
        let fm = FileManager.default
        if let attributes = try? fm.attributesOfItem(atPath: url.path),
           let size = attributes[.size] as? NSNumber, size.intValue > 2_000_000 {
            try? fm.removeItem(at: url)
        }
        let line = Data("\(ISO8601DateFormatter().string(from: Date())) \(text)\n".utf8)
        if let handle = try? FileHandle(forWritingTo: url) {
            _ = try? handle.seekToEnd()
            try? handle.write(contentsOf: line)
            try? handle.close()
        } else {
            try? line.write(to: url)
        }
    }
}

/// Finds and runs Open Babel (`obabel`) for bmObabel requests. Called on a background queue.
enum ObabelRunner {
    static let timeout: TimeInterval = 60

    static func searchDirectories() -> [String] {
        let home: String = NSHomeDirectory()
        let path: String = ProcessInfo.processInfo.environment["PATH"] ?? ""
        var dirs: [String] = path.split(separator: ":").map { String($0) }
        let fixed: [String] = ["/opt/homebrew/bin", "/usr/local/bin", "/opt/miniconda3/bin", "/opt/anaconda3/bin", "/opt/local/bin"]
        let inHome: [String] = ["/miniconda3/bin", "/anaconda3/bin", "/miniforge3/bin", "/mambaforge/bin"]
        dirs.append(contentsOf: fixed)
        for sub in inHome {
            dirs.append(home + sub)
        }
        return dirs
    }

    static func find() -> String? {
        for dir in searchDirectories() where !dir.isEmpty {
            let candidate = (dir as NSString).appendingPathComponent("obabel")
            if FileManager.default.isExecutableFile(atPath: candidate) {
                return candidate
            }
        }
        return nil
    }

    static func isFormatCode(_ text: String) -> Bool {
        guard !text.isEmpty, text.count <= 16 else { return false }
        for scalar in text.unicodeScalars {
            let v: UInt32 = scalar.value
            let digit = v >= 48 && v <= 57
            let upper = v >= 65 && v <= 90
            let lower = v >= 97 && v <= 122
            if !(digit || upper || lower || v == 95 || v == 45) {
                return false
            }
        }
        return true
    }

    /// Options that would break the stdin -> stdout contract or write files (-O, -m) are refused.
    static func isAllowedOption(_ option: String) -> Bool {
        if option.isEmpty || option.count > 200 || option.unicodeScalars.contains(where: { $0.value < 32 }) {
            return false
        }
        let blocked = ["-O", "-o", "-i", "-m", "-L", "-H", "-V"]
        return !blocked.contains(where: { option.hasPrefix($0) })
    }

    static func failure(_ message: String) -> [String: Any] {
        return ["ok": false, "output": "", "error": message]
    }

    static func run(_ request: [String: Any]) -> [String: Any] {
        let op = (request["op"] as? String) ?? ""
        var arguments: [String] = []
        var input = Data()
        if op == "formats" {
            arguments = ["-L", "formats"]
        } else if op == "convert" {
            let inFormat = (request["inFormat"] as? String) ?? ""
            let outFormat = (request["outFormat"] as? String) ?? ""
            guard isFormatCode(inFormat), isFormatCode(outFormat) else {
                return failure("Invalid format code.")
            }
            let options = (request["options"] as? [Any]) ?? []
            guard options.count <= 32 else {
                return failure("Too many Open Babel options.")
            }
            var extra: [String] = []
            for option in options {
                guard let text = option as? String, isAllowedOption(text) else {
                    return failure("Open Babel option not allowed: \(option)")
                }
                extra.append(text)
            }
            guard let data = Data(base64Encoded: (request["input"] as? String) ?? "", options: .ignoreUnknownCharacters) else {
                return failure("The input could not be decoded.")
            }
            input = data
            arguments = ["-i" + inFormat, "-o" + outFormat] + extra
        } else {
            return failure("Unknown Open Babel request.")
        }
        guard let executable = find() else {
            return failure("Open Babel (obabel) was not found. Install it with Homebrew (brew install open-babel) "
                + "or conda (conda install -c conda-forge openbabel).")
        }
        return execute(executable, arguments, input)
    }

    /// stdin, stdout and stderr go through temporary files: no pipe buffers to fill, no SIGPIPE.
    static func execute(_ executable: String, _ arguments: [String], _ input: Data) -> [String: Any] {
        let fm = FileManager.default
        let dir = fm.temporaryDirectory.appendingPathComponent("BlueMolecule-obabel-" + UUID().uuidString, isDirectory: true)
        defer { try? fm.removeItem(at: dir) }
        let inURL = dir.appendingPathComponent("stdin")
        let outURL = dir.appendingPathComponent("stdout")
        let errURL = dir.appendingPathComponent("stderr")
        do {
            try fm.createDirectory(at: dir, withIntermediateDirectories: true, attributes: nil)
            try input.write(to: inURL)
            try Data().write(to: outURL)
            try Data().write(to: errURL)
        } catch {
            return failure("Could not create temporary files: \(error.localizedDescription)")
        }

        let process = Process()
        process.executableURL = URL(fileURLWithPath: executable)
        process.arguments = arguments
        let finished = DispatchSemaphore(value: 0)
        process.terminationHandler = { _ in _ = finished.signal() }
        var timedOut = false
        do {
            let stdinHandle = try FileHandle(forReadingFrom: inURL)
            let stdoutHandle = try FileHandle(forWritingTo: outURL)
            let stderrHandle = try FileHandle(forWritingTo: errURL)
            defer {
                try? stdinHandle.close()
                try? stdoutHandle.close()
                try? stderrHandle.close()
            }
            process.standardInput = stdinHandle
            process.standardOutput = stdoutHandle
            process.standardError = stderrHandle
            try process.run()
            if finished.wait(timeout: .now() + ObabelRunner.timeout) == .timedOut {
                timedOut = true
                process.terminate()
                if finished.wait(timeout: .now() + 3) == .timedOut {
                    _ = kill(process.processIdentifier, SIGKILL)
                    _ = finished.wait(timeout: .now() + 3)
                }
            }
        } catch {
            return failure("Could not run \(executable): \(error.localizedDescription)")
        }
        if timedOut || process.isRunning {
            return failure("Open Babel did not finish within \(Int(ObabelRunner.timeout)) seconds.")
        }
        let output = (try? Data(contentsOf: outURL)) ?? Data()
        let errors = String(decoding: (try? Data(contentsOf: errURL)) ?? Data(), as: UTF8.self)
            .trimmingCharacters(in: .whitespacesAndNewlines)
        if process.terminationStatus != 0 || output.isEmpty {
            return failure(errors.isEmpty ? "Open Babel failed (exit status \(process.terminationStatus))." : errors)
        }
        return ["ok": true, "output": output.base64EncodedString(), "error": ""]
    }
}

// MARK: - Start-up (top-level code: this file is the program's main file)

let bmApplication = NSApplication.shared
let bmDelegate = AppDelegate()
bmApplication.delegate = bmDelegate
_ = bmApplication.setActivationPolicy(.regular)
bmApplication.run()
