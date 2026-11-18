// BlueMolecule desktop app (Tauri 2): a native window around the same page as BlueMolecule.html, and the host
// side of BM.native (src/native-bridge.js, kind "tauri-app"). The page asks; this answers:
//   bm_ready          documents named on the command line, by "Open with" or by a second launch (then "bm-open" events)
//   bm_read           bytes of such a document, on demand (big files are never copied whole)
//   bm_open_dialog    the system's Open dialog, giving documents read the same way
//   bm_save           the system's Save dialog
//   bm_copy, bm_log   the clipboard and a log file
//   bm_obabel         Open Babel, run as a program (stdin to stdout), as the other BlueMolecule hosts do
//   bm_open_url       web links in the default browser, not in this window
//   bm_update_check, bm_update_install   the built-in updater (signed with BlueMolecule's own key)
//   bm_self_test_done the report of `--self-test=report.json` (used by the cloud builds)
#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

use base64::{engine::general_purpose::STANDARD as B64, Engine as _};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::collections::HashMap;
use std::fs::{File, OpenOptions};
use std::io::{Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::Mutex;
use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};
use tauri::{AppHandle, Emitter, Manager, State, WebviewUrl, WebviewWindowBuilder};

/// Bytes per bm_read reply.
const READ_MAX: u64 = 64 * 1024 * 1024;
/// Seconds for one Open Babel run.
const OBABEL_TIMEOUT: u64 = 60;
/// Where releases announce themselves to the updater.
const UPDATE_URL: &str = "https://github.com/fred4h/BlueMolecule/releases/latest/download/latest.json";
/// The self-test of the cloud builds: what the page checks and reports.
const SELF_TEST_JS: &str = include_str!("selftest.js");

// ---------------------------------------------------------------- documents

#[derive(Clone, Serialize)]
#[serde(rename_all = "camelCase")]
struct DocInfo {
    name: String,
    size: u64,
    file_id: String,
}

#[derive(Default)]
struct DocsInner {
    files: HashMap<String, PathBuf>,
    pending: Vec<DocInfo>,
    ready: bool,
    next: u64,
}

/// Documents the page may read: only files the user named (command line, Open with, the Open dialog).
#[derive(Default)]
struct Docs(Mutex<DocsInner>);

impl Docs {
    fn register(&self, path: &Path) -> Option<DocInfo> {
        let meta = std::fs::metadata(path).ok()?;
        if !meta.is_file() {
            return None;
        }
        let name = path.file_name()?.to_string_lossy().into_owned();
        let mut d = self.0.lock().unwrap();
        d.next += 1;
        let id = format!("f{}", d.next);
        d.files.insert(id.clone(), path.to_path_buf());
        Some(DocInfo { name, size: meta.len(), file_id: id })
    }
    fn path(&self, id: &str) -> Option<PathBuf> {
        self.0.lock().unwrap().files.get(id).cloned()
    }
}

/// Documents for the page: queued until it is ready, then sent as "bm-open" events.
fn deliver(app: &AppHandle, paths: Vec<PathBuf>) {
    let docs = app.state::<Docs>();
    for p in paths {
        if let Some(info) = docs.register(&p) {
            let send = {
                let mut d = docs.0.lock().unwrap();
                if d.ready {
                    Some(info)
                } else {
                    d.pending.push(info);
                    None
                }
            };
            if let Some(info) = send {
                let _ = app.emit("bm-open", info);
            }
        }
    }
}

struct Launch {
    self_test: Option<PathBuf>,
    files: Vec<PathBuf>,
}

/// File names from a command line (flags are skipped; relative names are taken from `cwd`).
fn parse_args(args: &[String], cwd: &Path) -> Launch {
    let mut l = Launch { self_test: None, files: Vec::new() };
    for a in args.iter().skip(1) {
        if let Some(p) = a.strip_prefix("--self-test=") {
            l.self_test = Some(PathBuf::from(p));
            continue;
        }
        if a.starts_with('-') {
            continue; // other flags, e.g. macOS -psn_…
        }
        let p = PathBuf::from(a);
        l.files.push(if p.is_absolute() { p } else { cwd.join(p) });
    }
    l
}

// ---------------------------------------------------------------- self-test

struct SelfTest {
    report: Option<PathBuf>,
}

fn write_report(path: &Path, report: &Value) {
    let _ = std::fs::write(path, serde_json::to_vec_pretty(report).unwrap_or_default());
}

// ---------------------------------------------------------------- small helpers

fn memory() -> (u64, u64) {
    let mut s = sysinfo::System::new();
    s.refresh_memory();
    (s.total_memory(), s.available_memory())
}

fn now() -> u64 {
    SystemTime::now().duration_since(UNIX_EPOCH).map(|d| d.as_secs()).unwrap_or(0)
}

/// A file name the Save dialog can propose: no folders, no control characters.
fn safe_name(name: &str) -> String {
    let base = name.rsplit(|c| c == '/' || c == '\\').next().unwrap_or("");
    let s: String = base.chars().filter(|c| !c.is_control() && !"<>:\"|?*".contains(*c)).take(200).collect();
    if s.trim().is_empty() || s == "." || s == ".." {
        "untitled".to_string()
    } else {
        s
    }
}

// ---------------------------------------------------------------- commands

#[tauri::command]
fn bm_ready(docs: State<'_, Docs>) -> Vec<DocInfo> {
    let mut d = docs.0.lock().unwrap();
    d.ready = true;
    std::mem::take(&mut d.pending)
}

#[tauri::command]
fn bm_info(app: AppHandle) -> Value {
    let (total, avail) = memory();
    json!({ "version": app.package_info().version.to_string(), "os": std::env::consts::OS,
            "memTotal": total, "memAvailable": avail, "updates": updater_key().is_some() })
}

#[tauri::command]
async fn bm_read(docs: State<'_, Docs>, file_id: String, offset: u64, length: u64) -> Result<tauri::ipc::Response, String> {
    let path = docs.path(&file_id).ok_or_else(|| "This file is not open in BlueMolecule.".to_string())?;
    let mut f = File::open(&path).map_err(|e| format!("Could not open {}: {}", path.display(), e))?;
    f.seek(SeekFrom::Start(offset)).map_err(|e| format!("Could not read {}: {}", path.display(), e))?;
    let n = length.min(READ_MAX);
    let mut buf = Vec::with_capacity(n as usize);
    f.take(n).read_to_end(&mut buf).map_err(|e| format!("Could not read {}: {}", path.display(), e))?;
    Ok(tauri::ipc::Response::new(buf))
}

#[tauri::command]
async fn bm_open_dialog(app: AppHandle, multiple: bool) -> Result<Option<Vec<DocInfo>>, String> {
    use tauri_plugin_dialog::DialogExt;
    let picked: Vec<PathBuf> = if multiple {
        match app.dialog().file().set_title("Open").blocking_pick_files() {
            Some(v) => v.into_iter().filter_map(|p| p.into_path().ok()).collect(),
            None => return Ok(None),
        }
    } else {
        match app.dialog().file().set_title("Open").blocking_pick_file() {
            Some(p) => p.into_path().ok().into_iter().collect(),
            None => return Ok(None),
        }
    };
    let docs = app.state::<Docs>();
    Ok(Some(picked.iter().filter_map(|p| docs.register(p)).collect()))
}

#[derive(Serialize)]
struct SaveReply {
    ok: bool,
    cancelled: bool,
    path: Option<String>,
    error: Option<String>,
}

#[tauri::command]
async fn bm_save(app: AppHandle, test: State<'_, SelfTest>, name: String, data: String) -> Result<SaveReply, String> {
    let bytes = match B64.decode(data.as_bytes()) {
        Ok(b) => b,
        Err(e) => return Ok(SaveReply { ok: false, cancelled: false, path: None, error: Some(format!("Bad data: {}", e)) }),
    };
    let target = match &test.report {
        Some(r) => Some(r.with_extension("saved")), // the self-test saves without a dialog
        None => {
            use tauri_plugin_dialog::DialogExt;
            app.dialog().file().set_file_name(safe_name(&name)).blocking_save_file().and_then(|p| p.into_path().ok())
        }
    };
    let Some(path) = target else {
        return Ok(SaveReply { ok: false, cancelled: true, path: None, error: None });
    };
    Ok(match std::fs::write(&path, &bytes) {
        Ok(()) => SaveReply { ok: true, cancelled: false, path: Some(path.display().to_string()), error: None },
        Err(e) => SaveReply { ok: false, cancelled: false, path: None, error: Some(format!("{}: {}", path.display(), e)) },
    })
}

#[tauri::command]
fn bm_copy(app: AppHandle, text: String) -> bool {
    use tauri_plugin_clipboard_manager::ClipboardExt;
    app.clipboard().write_text(text).is_ok()
}

#[tauri::command]
fn bm_log(app: AppHandle, level: String, message: String) {
    if let Ok(dir) = app.path().app_log_dir() {
        let _ = std::fs::create_dir_all(&dir);
        if let Ok(mut f) = OpenOptions::new().create(true).append(true).open(dir.join("bluemolecule.log")) {
            let msg: String = message.chars().take(8000).collect();
            let _ = writeln!(f, "{} [{}] {}", now(), level.chars().take(16).collect::<String>(), msg);
        }
    }
}

#[tauri::command]
fn bm_open_url(app: AppHandle, url: String) -> Result<(), String> {
    if !(url.starts_with("https://") || url.starts_with("http://") || url.starts_with("mailto:")) {
        return Err("Only web and e-mail links open outside BlueMolecule.".to_string());
    }
    use tauri_plugin_opener::OpenerExt;
    app.opener().open_url(url, None::<&str>).map_err(|e| e.to_string())
}

#[tauri::command]
fn bm_self_test_done(app: AppHandle, test: State<'_, SelfTest>, report: Value, image: Option<String>) {
    let Some(path) = &test.report else { return };
    write_report(path, &report);
    if let Some(img) = image {
        if let Ok(b) = B64.decode(img.as_bytes()) {
            let _ = std::fs::write(path.with_extension("png"), b);
        }
    }
    let ok = report.get("ok").and_then(|v| v.as_bool()).unwrap_or(false);
    app.exit(if ok { 0 } else { 1 });
}

// ---------------------------------------------------------------- Open Babel

#[derive(Deserialize)]
#[serde(rename_all = "camelCase")]
struct ObReq {
    op: String,
    input: Option<String>,
    in_format: Option<String>,
    out_format: Option<String>,
    options: Option<Vec<String>>,
}

#[derive(Serialize)]
struct ObReply {
    ok: bool,
    output: String,
    error: String,
}

fn ob_fail(msg: &str) -> ObReply {
    ObReply { ok: false, output: String::new(), error: msg.to_string() }
}

fn format_ok(s: &str) -> bool {
    !s.is_empty() && s.len() <= 16 && s.chars().all(|c| c.is_ascii_alphanumeric() || c == '_' || c == '-')
}

/// Options that would break the stdin -> stdout contract or write files are refused (as in the other hosts).
fn option_ok(o: &str) -> bool {
    !o.is_empty() && o.len() <= 200 && !o.chars().any(|c| (c as u32) < 32)
        && !["-O", "-o", "-i", "-m", "-L", "-H", "-V"].iter().any(|b| o.starts_with(b))
}

fn executable(p: &Path) -> bool {
    p.is_file()
}

fn find_obabel() -> Option<PathBuf> {
    let exe = if cfg!(windows) { "obabel.exe" } else { "obabel" };
    if let Some(paths) = std::env::var_os("PATH") {
        for dir in std::env::split_paths(&paths) {
            let p = dir.join(exe);
            if executable(&p) {
                return Some(p);
            }
        }
    }
    // a program started from the Dock or the Start menu does not see the shell's PATH
    let home = std::env::var_os(if cfg!(windows) { "USERPROFILE" } else { "HOME" }).map(PathBuf::from).unwrap_or_default();
    let mut dirs: Vec<PathBuf> = Vec::new();
    if cfg!(windows) {
        for base in ["C:\\Program Files", "C:\\Program Files (x86)"] {
            if let Ok(rd) = std::fs::read_dir(base) {
                for e in rd.flatten() {
                    if e.file_name().to_string_lossy().to_lowercase().starts_with("openbabel") {
                        dirs.push(e.path());
                    }
                }
            }
        }
        for c in ["miniconda3", "anaconda3", "miniforge3", "mambaforge"] {
            dirs.push(home.join(c).join("Library").join("bin"));
        }
    } else {
        for d in ["/opt/homebrew/bin", "/usr/local/bin", "/usr/bin", "/opt/miniconda3/bin", "/opt/anaconda3/bin", "/opt/conda/bin"] {
            dirs.push(PathBuf::from(d));
        }
        for c in ["miniconda3", "anaconda3", "miniforge3", "mambaforge"] {
            dirs.push(home.join(c).join("bin"));
        }
    }
    dirs.into_iter().map(|d| d.join(exe)).find(|p| executable(p))
}

fn obabel(req: ObReq) -> ObReply {
    let (args, data): (Vec<String>, Vec<u8>) = match req.op.as_str() {
        "formats" => (vec!["-L".into(), "formats".into()], Vec::new()),
        "convert" => {
            let (Some(inf), Some(outf)) = (req.in_format.as_deref(), req.out_format.as_deref()) else {
                return ob_fail("Invalid format code.");
            };
            if !format_ok(inf) || !format_ok(outf) {
                return ob_fail("Invalid format code.");
            }
            let opts = req.options.unwrap_or_default();
            if opts.len() > 32 {
                return ob_fail("Invalid Open Babel options.");
            }
            if let Some(o) = opts.iter().find(|o| !option_ok(o)) {
                return ob_fail(&format!("Open Babel option not allowed: {:?}", o));
            }
            let Ok(data) = B64.decode(req.input.unwrap_or_default().as_bytes()) else {
                return ob_fail("The input is not valid base64.");
            };
            let mut a = vec![format!("-i{}", inf), format!("-o{}", outf)];
            a.extend(opts);
            (a, data)
        }
        _ => return ob_fail("Unknown Open Babel request."),
    };
    let Some(exe) = find_obabel() else {
        return ob_fail("Open Babel (obabel) was not found. Install it, e.g. conda install -c conda-forge openbabel, brew install open-babel or sudo apt install openbabel.");
    };
    let mut cmd = Command::new(&exe);
    cmd.args(&args).stdin(Stdio::piped()).stdout(Stdio::piped()).stderr(Stdio::piped());
    #[cfg(windows)]
    {
        use std::os::windows::process::CommandExt;
        cmd.creation_flags(0x0800_0000); // CREATE_NO_WINDOW: no console window flashes up
    }
    let mut child = match cmd.spawn() {
        Ok(c) => c,
        Err(e) => return ob_fail(&format!("Could not run {}: {}", exe.display(), e)),
    };
    let mut stdin = child.stdin.take();
    let writer = std::thread::spawn(move || {
        if let Some(s) = stdin.as_mut() {
            let _ = s.write_all(&data);
        }
        drop(stdin); // closes obabel's input
    });
    let mut out = child.stdout.take();
    let mut err = child.stderr.take();
    let reader = std::thread::spawn(move || {
        let mut v = Vec::new();
        if let Some(o) = out.as_mut() {
            let _ = o.read_to_end(&mut v);
        }
        v
    });
    let errs = std::thread::spawn(move || {
        let mut v = Vec::new();
        if let Some(e) = err.as_mut() {
            let _ = e.read_to_end(&mut v);
        }
        v
    });
    let t0 = Instant::now();
    let status = loop {
        match child.try_wait() {
            Ok(Some(s)) => break Some(s),
            Ok(None) if t0.elapsed() < Duration::from_secs(OBABEL_TIMEOUT) => std::thread::sleep(Duration::from_millis(20)),
            _ => {
                let _ = child.kill();
                let _ = child.wait();
                break None;
            }
        }
    };
    let _ = writer.join();
    let stdout = reader.join().unwrap_or_default();
    let stderr = String::from_utf8_lossy(&errs.join().unwrap_or_default()).trim().to_string();
    match status {
        None => ob_fail(&format!("Open Babel did not finish within {} s.", OBABEL_TIMEOUT)),
        Some(s) if !s.success() || stdout.is_empty() => {
            ob_fail(if stderr.is_empty() { "Open Babel failed." } else { &stderr })
        }
        Some(_) => ObReply { ok: true, output: B64.encode(&stdout), error: String::new() },
    }
}

#[tauri::command]
async fn bm_obabel(request: ObReq) -> Result<ObReply, String> {
    Ok(tauri::async_runtime::spawn_blocking(move || obabel(request))
        .await
        .unwrap_or_else(|e| ob_fail(&format!("Open Babel failed: {}", e))))
}

// ---------------------------------------------------------------- updates

/// BlueMolecule's public update key, given to release builds (BM_UPDATER_PUBKEY); other builds do not update.
fn updater_key() -> Option<&'static str> {
    option_env!("BM_UPDATER_PUBKEY").filter(|k| !k.trim().is_empty())
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct UpdateInfo {
    supported: bool,
    available: bool,
    current: String,
    version: Option<String>,
    notes: Option<String>,
    date: Option<String>,
    error: Option<String>,
}

fn updater(app: &AppHandle) -> Result<tauri_plugin_updater::Updater, String> {
    use tauri_plugin_updater::UpdaterExt;
    let key = updater_key().ok_or_else(|| "This is a test build of BlueMolecule, which does not look for updates. The installers on the Releases page do.".to_string())?;
    let url = tauri::Url::parse(UPDATE_URL).map_err(|e| e.to_string())?;
    app.updater_builder()
        .pubkey(key)
        .endpoints(vec![url])
        .map_err(|e| e.to_string())?
        .build()
        .map_err(|e| e.to_string())
}

#[tauri::command]
async fn bm_update_check(app: AppHandle) -> Result<UpdateInfo, String> {
    let current = app.package_info().version.to_string();
    let none = |current: String, supported: bool, error: Option<String>| UpdateInfo {
        supported, available: false, current, version: None, notes: None, date: None, error,
    };
    let up = match updater(&app) {
        Ok(u) => u,
        Err(e) => return Ok(none(current, false, Some(e))),
    };
    Ok(match up.check().await {
        Ok(Some(u)) => UpdateInfo {
            supported: true, available: true, current, version: Some(u.version.clone()), notes: u.body.clone(),
            date: u.date.map(|d| d.to_string()), error: None,
        },
        Ok(None) => none(current, true, None),
        Err(e) => none(current, true, Some(e.to_string())),
    })
}

#[tauri::command]
async fn bm_update_install(app: AppHandle) -> Result<(), String> {
    let up = updater(&app)?;
    let Some(update) = up.check().await.map_err(|e| e.to_string())? else {
        return Err("No newer version is available.".to_string());
    };
    update.download_and_install(|_chunk, _total| {}, || {}).await.map_err(|e| e.to_string())?;
    app.restart();
}

// ---------------------------------------------------------------- start

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let cwd = std::env::current_dir().unwrap_or_default();
    let launch = parse_args(&args, &cwd);
    let self_test = launch.self_test.clone();
    let files = launch.files;

    let app = tauri::Builder::default()
        .plugin(tauri_plugin_single_instance::init(|app, argv, cwd| {
            // a second launch (Open with while BlueMolecule runs): its documents go to this window
            let l = parse_args(&argv, Path::new(&cwd));
            deliver(app, l.files);
            if let Some(w) = app.get_webview_window("main") {
                let _ = w.unminimize();
                let _ = w.show();
                let _ = w.set_focus();
            }
        }))
        .plugin(tauri_plugin_dialog::init())
        .plugin(tauri_plugin_clipboard_manager::init())
        .plugin(tauri_plugin_opener::init())
        .plugin(tauri_plugin_updater::Builder::new().build())
        .manage(Docs::default())
        .manage(SelfTest { report: self_test.clone() })
        .invoke_handler(tauri::generate_handler![
            bm_ready, bm_info, bm_read, bm_open_dialog, bm_save, bm_copy, bm_log, bm_open_url, bm_obabel,
            bm_update_check, bm_update_install, bm_self_test_done
        ])
        .setup(move |app| {
            let handle = app.handle().clone();
            let (total, avail) = memory();
            let names: Vec<String> = if self_test.is_some() {
                files.iter().filter_map(|p| p.file_name().map(|n| n.to_string_lossy().into_owned())).collect()
            } else {
                Vec::new()
            };
            let marker = json!({
                "version": app.package_info().version.to_string(), "os": std::env::consts::OS,
                "appImage": std::env::var_os("APPIMAGE").is_some(),
                "memTotal": total, "memAvailable": avail, "updates": updater_key().is_some(),
                "selfTest": self_test.is_some(), "selfTestFiles": names
            });
            let init = format!("window.BM_TAURI = Object.freeze({});\n{}", marker, SELF_TEST_JS);
            let nav = handle.clone();
            WebviewWindowBuilder::new(app, "main", WebviewUrl::App("index.html".into()))
                .title("BlueMolecule")
                .inner_size(1440.0, 900.0)
                .min_inner_size(820.0, 560.0)
                .disable_drag_drop_handler() // files dropped on the window reach the page, as in a browser
                .initialization_script(&init)
                .on_navigation(move |url| {
                    // the page itself: tauri://localhost (macOS, Linux) or http://tauri.localhost (Windows)
                    let local = matches!(url.scheme(), "tauri" | "asset" | "about" | "blob" | "data")
                        || (matches!(url.scheme(), "http" | "https") && url.host_str() == Some("tauri.localhost"));
                    if !local && matches!(url.scheme(), "http" | "https" | "mailto") {
                        use tauri_plugin_opener::OpenerExt;
                        let _ = nav.opener().open_url(url.as_str(), None::<&str>);
                    }
                    local
                })
                .build()?;
            deliver(&handle, files.clone());
            if let Some(report) = self_test.clone() {
                // the cloud builds wait for a report: one is written even if the page never sends its own
                std::thread::spawn(move || {
                    std::thread::sleep(Duration::from_secs(900));
                    write_report(&report, &json!({ "ok": false, "error": "The self-test did not finish within 15 minutes." }));
                    std::process::exit(3);
                });
            }
            Ok(())
        })
        .build(tauri::generate_context!())
        .expect("BlueMolecule could not start");

    app.run(|handle, event| {
        // macOS gives documents opened from the Finder or the Dock as events, not on the command line
        #[cfg(target_os = "macos")]
        if let tauri::RunEvent::Opened { urls } = &event {
            let paths: Vec<PathBuf> = urls.iter().filter_map(|u| u.to_file_path().ok()).collect();
            deliver(handle, paths);
        }
        let _ = (handle, event);
    });
}
