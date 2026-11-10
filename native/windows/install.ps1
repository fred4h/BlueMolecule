# Install BlueMolecule for the current user (Windows, experimental).
#
# Run by double-clicking "Install BlueMolecule.cmd", or in a console:
#     powershell -NoProfile -ExecutionPolicy Bypass -File windows\install.ps1
# It needs no administrator rights: it writes only to your profile and to HKEY_CURRENT_USER.
#   - copies BlueMolecule into %LOCALAPPDATA%\Programs\BlueMolecule;
#   - compiles the small BlueMolecule.exe launcher (BlueMolecule.cs) with the C# compiler of the .NET Framework,
#     which is part of Windows;
#   - adds BlueMolecule to the Start menu and to "Open with" for chemistry files (it never changes the program a
#     file type opens with by default);
#   - registers an uninstaller under Settings > Apps.
# Running it again updates the installation. Keep this file ASCII (Windows PowerShell 5.1 reads it as ANSI).

Set-StrictMode -Version 2
$ErrorActionPreference = 'Stop'

$BMVersion = '@VERSION@'
$Package = Split-Path -Parent $PSScriptRoot
$Dest = Join-Path $env:LOCALAPPDATA 'Programs\BlueMolecule'
$Exe = Join-Path $Dest 'BlueMolecule.exe'
$Icon = Join-Path $Dest 'BlueMolecule.ico'
$Shortcut = Join-Path ([Environment]::GetFolderPath('Programs')) 'BlueMolecule.lnk'
$Classes = 'HKCU:\Software\Classes'
$ProgId = 'BlueMolecule.Document'
$AppKey = $Classes + '\Applications\BlueMolecule.exe'
$UninstallKey = 'HKCU:\Software\Microsoft\Windows\CurrentVersion\Uninstall\BlueMolecule'
# File types that get an "Open with > BlueMolecule" entry (same list in uninstall.ps1)
$Extensions = @('log', 'out', 'xyz', 'gjf', 'inp', 'fchk', 'fch', 'cube', 'cub', 'molden', 'mol', 'sdf', 'mol2', 'pdb', 'ent',
    'cif', 'mmcif', 'gro', 'cml', 'cjson', 'xsf', 'axsf', 'vasp', 'prmtop', 'parm7', 'inpcrd', 'rst7', 'psf', 'dcd', 'xtc',
    'trr', 'lammpstrj', 'pwi', 'pwo', 'cell', 'nw', 'mop', 'arc', 'gz')

function Say([string] $Text) { Write-Host $Text }

function New-KeyIfMissing([string] $Key) {
    if (Test-Path -LiteralPath $Key) { return }
    New-KeyIfMissing (Split-Path -Parent $Key)
    [void](New-Item -Path $Key)
}

function Set-DefaultValue([string] $Key, [string] $Value) {
    New-KeyIfMissing $Key
    Set-Item -LiteralPath $Key -Value $Value
}

function Set-StringValue([string] $Key, [string] $Name, [string] $Value) {
    New-KeyIfMissing $Key
    [void](New-ItemProperty -LiteralPath $Key -Name $Name -Value $Value -PropertyType String -Force)
}

function Set-DwordValue([string] $Key, [string] $Name, [int] $Value) {
    New-KeyIfMissing $Key
    [void](New-ItemProperty -LiteralPath $Key -Name $Name -Value $Value -PropertyType DWord -Force)
}

function Find-Browser {
    foreach ($name in @('msedge.exe', 'chrome.exe', 'brave.exe', 'vivaldi.exe')) {
        foreach ($root in @('HKCU:', 'HKLM:')) {
            $key = $root + '\SOFTWARE\Microsoft\Windows\CurrentVersion\App Paths\' + $name
            try {
                $value = (Get-Item -LiteralPath $key -ErrorAction Stop).GetValue('')
                if ($value) {
                    $value = ([string]$value).Trim().Trim('"')
                    if (Test-Path -LiteralPath $value) { return $value }
                }
            } catch { }
        }
    }
    $candidates = @()
    $pf86 = [Environment]::GetEnvironmentVariable('ProgramFiles(x86)')   # missing on 32-bit Windows
    if ($pf86) { $candidates += (Join-Path $pf86 'Microsoft\Edge\Application\msedge.exe') }
    if ($env:ProgramFiles) {
        $candidates += (Join-Path $env:ProgramFiles 'Microsoft\Edge\Application\msedge.exe')
        $candidates += (Join-Path $env:ProgramFiles 'Google\Chrome\Application\chrome.exe')
    }
    if ($env:LOCALAPPDATA) { $candidates += (Join-Path $env:LOCALAPPDATA 'Google\Chrome\Application\chrome.exe') }
    foreach ($c in $candidates) { if (Test-Path -LiteralPath $c) { return $c } }
    return $null
}

function Invoke-Compiler([string] $Csc, [string[]] $Arguments) {
    # Windows PowerShell turns a native program's error output into errors; collect everything as text instead.
    $saved = $ErrorActionPreference
    $ErrorActionPreference = 'Continue'
    try {
        $output = & $Csc @Arguments 2>&1 | ForEach-Object { [string]$_ }
        return @{ ok = ($LASTEXITCODE -eq 0); output = $output }
    } finally {
        $ErrorActionPreference = $saved
    }
}

if ($PSVersionTable.PSVersion.Major -lt 3) {
    Say 'BlueMolecule needs Windows PowerShell 3 or newer (Windows 8 or later).'
    exit 1
}
if (-not (Test-Path -LiteralPath (Join-Path $Package 'BlueMolecule.html'))) {
    Say ('BlueMolecule.html was not found in ' + $Package + '.')
    Say 'Extract the whole zip file first (right-click it > Extract All), then run Install BlueMolecule.cmd again.'
    exit 1
}

Say ('Installing BlueMolecule ' + $BMVersion + ' (experimental Windows version) in ' + $Dest)

# ---- files
[void](New-Item -ItemType Directory -Path $Dest -Force)
if (Test-Path -LiteralPath $Exe) {
    try { Remove-Item -LiteralPath $Exe -Force }
    catch { Say 'BlueMolecule.exe is in use: close BlueMolecule and run the installer again.'; exit 1 }
}
foreach ($f in @('BlueMolecule.html', 'LICENSE', 'README.txt')) {
    Copy-Item -LiteralPath (Join-Path $Package $f) -Destination $Dest -Force
}
foreach ($f in @('BlueMolecule.cs', 'BlueMolecule.ico', 'uninstall.ps1')) {
    Copy-Item -LiteralPath (Join-Path $PSScriptRoot $f) -Destination $Dest -Force
}
# files from a downloaded zip carry the "downloaded from the internet" mark; the installed copies do not need it
Get-ChildItem -LiteralPath $Dest -File | Unblock-File -ErrorAction SilentlyContinue

# ---- the launcher: BlueMolecule.exe, compiled here (with the icon if the compiler accepts it, else without)
$csc = $null
foreach ($fw in @('Framework64', 'Framework')) {
    $c = Join-Path $env:WINDIR ('Microsoft.NET\' + $fw + '\v4.0.30319\csc.exe')
    if (Test-Path -LiteralPath $c) { $csc = $c; break }
}
$built = $false
if ($csc) {
    $source = Join-Path $Dest 'BlueMolecule.cs'
    $common = @('/nologo', '/target:winexe', '/optimize+', '/platform:anycpu', ('/out:' + $Exe), '/r:System.Windows.Forms.dll')
    $r = Invoke-Compiler $csc ($common + @(('/win32icon:' + $Icon), $source))
    if (-not $r.ok) { $r = Invoke-Compiler $csc ($common + @($source)) }
    if ($r.ok -and (Test-Path -LiteralPath $Exe)) {
        $built = $true
        Say 'Compiled the BlueMolecule.exe launcher.'
    } else {
        Say 'Could not compile BlueMolecule.exe (please report this, with the lines below):'
        foreach ($line in @($r.output)) { Say ('    ' + $line) }
    }
} else {
    Say 'The .NET Framework C# compiler (csc.exe) was not found, so BlueMolecule.exe could not be made.'
}

# ---- Start menu
$browser = Find-Browser
$shell = New-Object -ComObject WScript.Shell
$lnk = $shell.CreateShortcut($Shortcut)
if ($built) {
    $lnk.TargetPath = $Exe
    $lnk.Arguments = ''
} elseif ($browser) {
    # without the launcher, the shortcut opens the page in an app window directly (no "Open with")
    $profileDir = Join-Path $env:LOCALAPPDATA 'BlueMolecule\browser'
    [void](New-Item -ItemType Directory -Path $profileDir -Force)
    $pageUrl = (New-Object -TypeName System.Uri -ArgumentList (Join-Path $Dest 'BlueMolecule.html')).AbsoluteUri
    $lnk.TargetPath = $browser
    $lnk.Arguments = '--app=' + $pageUrl + ' --user-data-dir="' + $profileDir + '" --no-first-run --no-default-browser-check'
} else {
    $lnk.TargetPath = Join-Path $Dest 'BlueMolecule.html'
    $lnk.Arguments = ''
}
$lnk.WorkingDirectory = $Dest
$lnk.IconLocation = $Icon + ',0'
$lnk.Description = 'BlueMolecule: viewer for quantum-chemistry and molecular-modelling files'
$lnk.Save()
Say ('Added BlueMolecule to the Start menu (' + $Shortcut + ').')

# ---- Open with (only with the launcher: it is what receives the file)
if ($built) {
    $command = '"' + $Exe + '" "%1"'
    Set-DefaultValue ($Classes + '\' + $ProgId) 'BlueMolecule document'
    Set-DefaultValue ($Classes + '\' + $ProgId + '\DefaultIcon') ($Exe + ',0')
    Set-DefaultValue ($Classes + '\' + $ProgId + '\shell\open\command') $command
    Set-StringValue $AppKey 'FriendlyAppName' 'BlueMolecule'
    Set-DefaultValue ($AppKey + '\DefaultIcon') ($Exe + ',0')
    Set-DefaultValue ($AppKey + '\shell\open\command') $command
    foreach ($e in $Extensions) {
        Set-StringValue ($Classes + '\.' + $e + '\OpenWithProgids') $ProgId ''
        Set-StringValue ($AppKey + '\SupportedTypes') ('.' + $e) ''
    }
    Say ('Added "Open with > BlueMolecule" for .' + ($Extensions -join ' .'))
}

# ---- uninstaller (Settings > Apps)
$powershell = Join-Path $env:SystemRoot 'System32\WindowsPowerShell\v1.0\powershell.exe'
Set-StringValue $UninstallKey 'DisplayName' 'BlueMolecule'
Set-StringValue $UninstallKey 'DisplayVersion' $BMVersion
Set-StringValue $UninstallKey 'Publisher' 'Frederico Henrique Do. C. R. Ferreira'
Set-StringValue $UninstallKey 'DisplayIcon' $Icon
Set-StringValue $UninstallKey 'InstallLocation' $Dest
Set-StringValue $UninstallKey 'URLInfoAbout' 'https://github.com/fred4h/BlueMolecule'
Set-StringValue $UninstallKey 'UninstallString' ('"' + $powershell + '" -NoProfile -ExecutionPolicy Bypass -File "' + (Join-Path $Dest 'uninstall.ps1') + '"')
Set-DwordValue $UninstallKey 'NoModify' 1
Set-DwordValue $UninstallKey 'NoRepair' 1

# ---- let Explorer pick up the new icons and Open with entries
try {
    Add-Type -Namespace BlueMoleculeSetup -Name Shell -MemberDefinition '[DllImport("shell32.dll")] public static extern void SHChangeNotify(int eventId, int flags, IntPtr item1, IntPtr item2);'
    [BlueMoleculeSetup.Shell]::SHChangeNotify(0x08000000, 0, [IntPtr]::Zero, [IntPtr]::Zero)
} catch { }

Say ''
Say 'BlueMolecule is installed.'
Say '  Start it from the Start menu, or right-click a chemistry file and choose Open with > BlueMolecule.'
if ($browser) { Say ('  It opens in an app window of ' + $browser) }
else { Say '  No Edge, Chrome, Brave or Vivaldi was found: it opens in your default browser.' }
if (-not $built) { Say '  "Open with" is not available, because BlueMolecule.exe could not be made (see above).' }
Say '  Uninstall: Settings > Apps > Installed apps > BlueMolecule.'
Say ''
Say 'This Windows version is experimental. Please report problems: https://github.com/fred4h/BlueMolecule/issues'
