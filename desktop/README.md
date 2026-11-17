# BlueMolecule desktop app

BlueMolecule 2 runs as a program of its own, built with [Tauri 2](https://tauri.app): a native window on Linux,
Windows and macOS around the same page as `BlueMolecule.html`, with installers and "Open with" entries. The window
uses each system's web engine: WebKitGTK on Linux, Edge WebView2 on Windows, WebKit on macOS.

| Part | What it does |
|---|---|
| `src-tauri/src/main.rs` | The host side of `BM.native` (`src/native-bridge.js`, kind `tauri-app`): documents from the command line, "Open with" and a second launch, read on demand; the Open and Save dialogs; the clipboard; Open Babel; web links in the default browser; the updater |
| `src-tauri/src/selftest.js` | The start-up check run by the cloud builds (`--self-test=report.json`) |
| `src-tauri/tauri.conf.json` | Name, icons, installers, file types; `tauri.<os>.conf.json` adjust it per system |
| `src-tauri/windows/hooks.nsh` | Windows installer: "Open with" entries for chemistry files (never the default program) |
| `src-tauri/linux/` | `.deb` scripts that register the chemistry file types |
| `www/` | The page, copied in for a build (not in git) |

## Building

The cloud builds (`.github/workflows/desktop.yml`) build, start and test the app on all three systems for every
push that touches it, and keep the installers as run artifacts. To build on your own computer you need Rust, Node.js
and, on Linux, `libwebkit2gtk-4.1-dev libgtk-3-dev librsvg2-dev`:

```bash
python3 tools/build.py                       # writes dist/BlueMolecule.html
mkdir -p desktop/www && cp dist/BlueMolecule.html desktop/www/index.html
cd desktop && npx @tauri-apps/cli@2 build    # installers in src-tauri/target/release/bundle/
```

The self-test opens the files it is given, checks WebGL 2, saving and the native bridge, takes a picture of the view
and exits with 0 when everything passed:

```bash
src-tauri/target/release/bluemolecule --self-test=$PWD/report.json ../tests/fixtures/cisplatin.xyz
```

## Updates

Release builds check GitHub for a newer release only when you ask (Help → Check for updates) or tick the start-up
check. An update is installed only if its signature matches BlueMolecule's public key, which release builds get as
`BM_UPDATER_PUBKEY`; the private key is held by the maintainer alone. Other builds cannot update themselves.
