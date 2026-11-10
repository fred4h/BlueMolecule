// BlueMolecule launcher for Windows (experimental).
//
// install.ps1 compiles this file on the user's computer with the C# compiler of the .NET Framework
// (csc.exe, C# 5), so it needs nothing that Windows 10 and 11 do not already have:
//
//   BlueMolecule.exe [files...]
//       Opens BlueMolecule.html (next to this program) in an app window of Microsoft Edge, Chrome, Brave or
//       Vivaldi, with its own browser profile in %LOCALAPPDATA%\BlueMolecule\browser. Files named on the
//       command line ("Open with", drag and drop on the icon) are embedded as data: URLs in window.BM_LAUNCH
//       in a temporary copy of the page in %LOCALAPPDATA%\BlueMolecule\launch, so no local server and no
//       browser flags are needed. Files larger than 64 MB are not embedded: open them from the window.
//   BlueMolecule.exe --diagnose    shows what it found
//   BlueMolecule.exe --tab         opens the default browser instead of an app window
//   BLUEMOLECULE_BROWSER=C:\path\to\browser.exe chooses the browser.
//
// The page contract (window.BM_LAUNCH = {files: [{name, url}], version}) is the one of the Python launcher
// (launcher/bluemolecule.py) without its ping and Open Babel endpoints; tests/test_windows_launch.py checks it.
// Keep this file ASCII and C# 5 (no string interpolation, no ?. operator, no expression-bodied members).
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Reflection;
using System.Text;
using System.Windows.Forms;
using Microsoft.Win32;

[assembly: AssemblyTitle("BlueMolecule")]
[assembly: AssemblyDescription("BlueMolecule molecule viewer")]
[assembly: AssemblyProduct("BlueMolecule")]
[assembly: AssemblyCopyright("Copyright (c) 2026 Frederico Henrique Do. C. R. Ferreira. MIT License.")]
[assembly: AssemblyVersion("@VERSION4@")]
[assembly: AssemblyFileVersion("@VERSION4@")]

static class BlueMoleculeLauncher
{
    const string Version = "@VERSION@";
    const long InlineLimit = 64L * 1024 * 1024;
    static readonly string[] BrowserExes = { "msedge.exe", "chrome.exe", "brave.exe", "vivaldi.exe" };

    static string Here
    {
        get { return Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location); }
    }

    static string DataDir
    {
        get { return Path.Combine(Environment.GetFolderPath(Environment.SpecialFolder.LocalApplicationData), "BlueMolecule"); }
    }

    [STAThread]
    static int Main(string[] args)
    {
        List<string> files = new List<string>();
        bool diagnose = false, tab = false;
        foreach (string a in args)
        {
            if (a == "--diagnose") diagnose = true;
            else if (a == "--tab") tab = true;
            else if (a == "--version") { Show("BlueMolecule " + Version, MessageBoxIcon.Information); return 0; }
            else files.Add(a);
        }
        string page = Path.Combine(Here, "BlueMolecule.html");
        string browser = FindBrowser();
        if (diagnose)
        {
            Show(Diagnostics(page, browser), MessageBoxIcon.Information);
            return 0;
        }
        try
        {
            if (!File.Exists(page)) throw new FileNotFoundException("BlueMolecule.html was not found next to BlueMolecule.exe: " + page);
            List<FileInfo> open = new List<FileInfo>();
            List<string> skipped = new List<string>();
            foreach (string f in files)
            {
                FileInfo fi;
                try { fi = new FileInfo(Path.GetFullPath(f)); }
                catch (Exception) { skipped.Add(f + " (not a valid path)"); continue; }
                if (!fi.Exists) { skipped.Add(f + " (not found)"); continue; }
                if (fi.Length > InlineLimit) { skipped.Add(fi.Name + " (" + (fi.Length / (1024 * 1024)) + " MB)"); continue; }
                open.Add(fi);
            }
            string target = open.Count > 0 ? WriteLaunchPage(page, open) : page;
            OpenWindow(browser, target, tab);
            Log("opened " + target + " with " + (browser ?? "the default browser") + "; files: " + open.Count +
                (skipped.Count > 0 ? "; skipped: " + string.Join(", ", skipped) : ""));
            if (skipped.Count > 0)
            {
                Show("These files were not opened automatically:\n\n" + string.Join("\n", skipped) +
                     "\n\nOpen files larger than 64 MB with Open... in the BlueMolecule window, or drag them into it.",
                     MessageBoxIcon.Warning);
            }
            return 0;
        }
        catch (Exception e)
        {
            Log("error: " + e);
            Show("BlueMolecule could not start:\n\n" + e.Message, MessageBoxIcon.Error);
            return 1;
        }
    }

    static void Show(string text, MessageBoxIcon icon)
    {
        MessageBox.Show(text, "BlueMolecule", MessageBoxButtons.OK, icon);
    }

    static void Log(string line)
    {
        try
        {
            Directory.CreateDirectory(DataDir);
            File.AppendAllText(Path.Combine(DataDir, "bluemolecule.log"),
                DateTime.Now.ToString("yyyy-MM-dd HH:mm:ss") + " " + line + Environment.NewLine, Encoding.UTF8);
        }
        catch (Exception) { }
    }

    // ------------------------------------------------------------------ the browser
    static string FindBrowser()
    {
        string env = Environment.GetEnvironmentVariable("BLUEMOLECULE_BROWSER");
        if (!string.IsNullOrEmpty(env) && File.Exists(env)) return env;
        foreach (string exe in BrowserExes)
        {
            foreach (RegistryKey root in new RegistryKey[] { Registry.CurrentUser, Registry.LocalMachine })
            {
                try
                {
                    using (RegistryKey k = root.OpenSubKey(@"SOFTWARE\Microsoft\Windows\CurrentVersion\App Paths\" + exe))
                    {
                        string p = k == null ? null : k.GetValue("") as string;
                        if (!string.IsNullOrEmpty(p))
                        {
                            p = p.Trim().Trim('"');
                            if (File.Exists(p)) return p;
                        }
                    }
                }
                catch (Exception) { }
            }
        }
        string pf = Environment.GetFolderPath(Environment.SpecialFolder.ProgramFiles);
        string pf86 = Environment.GetFolderPath(Environment.SpecialFolder.ProgramFilesX86);
        string local = Environment.GetFolderPath(Environment.SpecialFolder.LocalApplicationData);
        string[] candidates = {
            Path.Combine(pf86, @"Microsoft\Edge\Application\msedge.exe"),
            Path.Combine(pf, @"Microsoft\Edge\Application\msedge.exe"),
            Path.Combine(pf, @"Google\Chrome\Application\chrome.exe"),
            Path.Combine(pf86, @"Google\Chrome\Application\chrome.exe"),
            Path.Combine(local, @"Google\Chrome\Application\chrome.exe"),
            Path.Combine(pf, @"BraveSoftware\Brave-Browser\Application\brave.exe"),
            Path.Combine(local, @"BraveSoftware\Brave-Browser\Application\brave.exe")
        };
        foreach (string c in candidates)
        {
            if (Path.IsPathRooted(c) && File.Exists(c)) return c;
        }
        return null;
    }

    static void OpenWindow(string browser, string target, bool tab)
    {
        if (browser == null || tab)
        {
            ProcessStartInfo shell = new ProcessStartInfo(target);
            shell.UseShellExecute = true;
            Process.Start(shell);
            return;
        }
        string url = new Uri(target).AbsoluteUri;
        string profile = Path.Combine(DataDir, "browser");
        Directory.CreateDirectory(profile);
        ProcessStartInfo psi = new ProcessStartInfo(browser,
            "--app=" + Quote(url) + " --user-data-dir=" + Quote(profile) + " --no-first-run --no-default-browser-check");
        psi.UseShellExecute = false;
        Process.Start(psi);
    }

    /// <summary>One command-line argument, quoted for CommandLineToArgvW when it needs it.</summary>
    static string Quote(string a)
    {
        if (a.Length > 0 && a.IndexOfAny(new char[] { ' ', '\t', '"' }) < 0) return a;
        StringBuilder sb = new StringBuilder("\"");
        int slashes = 0;
        foreach (char c in a)
        {
            if (c == '\\') { slashes++; continue; }
            if (c == '"') { sb.Append('\\', slashes * 2 + 1); sb.Append('"'); slashes = 0; continue; }
            if (slashes > 0) { sb.Append('\\', slashes); slashes = 0; }
            sb.Append(c);
        }
        sb.Append('\\', slashes * 2);
        sb.Append('"');
        return sb.ToString();
    }

    // ------------------------------------------------------------------ the launch page
    /// <summary>A copy of the page with the files embedded in window.BM_LAUNCH, written just before the first
    /// &lt;/head&gt; (after the charset meta tag) with a UTF-8 byte-order mark.</summary>
    static string WriteLaunchPage(string page, List<FileInfo> files)
    {
        string html = File.ReadAllText(page, Encoding.UTF8);
        int at = html.IndexOf("</head>", StringComparison.Ordinal);
        if (at < 0) throw new InvalidDataException("BlueMolecule.html has no </head> tag.");
        string dir = Path.Combine(DataDir, "launch");
        Directory.CreateDirectory(dir);
        RemoveOldPages(dir);
        string path = Path.Combine(dir, "launch-" + DateTime.Now.ToString("yyyyMMdd-HHmmss") + "-" +
            Guid.NewGuid().ToString("N").Substring(0, 8) + ".html");
        using (StreamWriter w = new StreamWriter(path, false, new UTF8Encoding(true)))
        {
            w.Write(html.Substring(0, at));
            w.Write("<script>window.BM_LAUNCH={\"files\":[");
            for (int i = 0; i < files.Count; i++)
            {
                if (i > 0) w.Write(",");
                w.Write("{\"name\":" + JsonString(files[i].Name) + ",\"url\":\"data:application/octet-stream;base64,");
                w.Write(Convert.ToBase64String(File.ReadAllBytes(files[i].FullName)));
                w.Write("\"}");
            }
            w.Write("],\"version\":" + JsonString(Version) + "};</script>");
            w.Write(html.Substring(at));
        }
        return path;
    }

    /// <summary>A JSON string literal that is also safe inside an HTML script element.</summary>
    static string JsonString(string s)
    {
        StringBuilder sb = new StringBuilder("\"");
        foreach (char c in s)
        {
            if (c == '"') sb.Append("\\\"");
            else if (c == '\\') sb.Append("\\\\");
            else if (c < ' ' || c == '<' || c == '>' || c == '&' || c == '\u2028' || c == '\u2029')
                sb.Append("\\u").Append(((int)c).ToString("x4"));
            else sb.Append(c);
        }
        sb.Append('"');
        return sb.ToString();
    }

    /// <summary>Launch pages older than a day are no longer needed (a reload within a day still works).</summary>
    static void RemoveOldPages(string dir)
    {
        try
        {
            foreach (string f in Directory.GetFiles(dir, "launch-*.html"))
            {
                try { if (File.GetLastWriteTime(f) < DateTime.Now.AddDays(-1)) File.Delete(f); }
                catch (Exception) { }
            }
        }
        catch (Exception) { }
    }

    static string Diagnostics(string page, string browser)
    {
        StringBuilder sb = new StringBuilder();
        sb.AppendLine("BlueMolecule " + Version + " for Windows (experimental)");
        sb.AppendLine();
        sb.AppendLine("Windows: " + Environment.OSVersion.VersionString + (Environment.Is64BitOperatingSystem ? ", 64-bit" : ", 32-bit"));
        sb.AppendLine(".NET runtime: " + Environment.Version);
        sb.AppendLine("Program: " + Assembly.GetExecutingAssembly().Location);
        sb.AppendLine("Page: " + page + (File.Exists(page) ? "" : "  (MISSING)"));
        sb.AppendLine("Browser: " + (browser ?? "none found; the default browser is used, in a normal tab"));
        sb.AppendLine("Profile: " + Path.Combine(DataDir, "browser"));
        sb.AppendLine("Log: " + Path.Combine(DataDir, "bluemolecule.log"));
        return sb.ToString();
    }
}
