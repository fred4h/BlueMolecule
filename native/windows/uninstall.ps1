# Remove BlueMolecule for the current user (Windows, experimental).
#
# Run from Settings > Apps > Installed apps > BlueMolecule > Uninstall, or in a console:
#     powershell -NoProfile -ExecutionPolicy Bypass -File "%LOCALAPPDATA%\Programs\BlueMolecule\uninstall.ps1"
# It removes exactly what install.ps1 added. Your settings and the browser profile in %LOCALAPPDATA%\BlueMolecule
# are kept; delete that folder to remove them too. Keep this file ASCII.

Set-StrictMode -Version 2
$ErrorActionPreference = 'Continue'

$Dest = Join-Path $env:LOCALAPPDATA 'Programs\BlueMolecule'
$Shortcut = Join-Path ([Environment]::GetFolderPath('Programs')) 'BlueMolecule.lnk'
$Classes = 'HKCU:\Software\Classes'
$ProgId = 'BlueMolecule.Document'
$AppKey = $Classes + '\Applications\BlueMolecule.exe'
$UninstallKey = 'HKCU:\Software\Microsoft\Windows\CurrentVersion\Uninstall\BlueMolecule'
# the same list as in install.ps1
$Extensions = @('log', 'out', 'xyz', 'gjf', 'inp', 'fchk', 'fch', 'cube', 'cub', 'molden', 'mol', 'sdf', 'mol2', 'pdb', 'ent',
    'cif', 'mmcif', 'gro', 'cml', 'cjson', 'xsf', 'axsf', 'vasp', 'prmtop', 'parm7', 'inpcrd', 'rst7', 'psf', 'dcd', 'xtc',
    'trr', 'lammpstrj', 'pwi', 'pwo', 'cell', 'nw', 'mop', 'arc', 'gz')

Write-Host 'Removing BlueMolecule...'
Remove-Item -LiteralPath $Shortcut -Force -ErrorAction SilentlyContinue
foreach ($e in $Extensions) {
    # only BlueMolecule's own entry: other programs' Open with entries and the default program stay
    Remove-ItemProperty -LiteralPath ($Classes + '\.' + $e + '\OpenWithProgids') -Name $ProgId -ErrorAction SilentlyContinue
}
Remove-Item -LiteralPath ($Classes + '\' + $ProgId) -Recurse -Force -ErrorAction SilentlyContinue
Remove-Item -LiteralPath $AppKey -Recurse -Force -ErrorAction SilentlyContinue
Remove-Item -LiteralPath $UninstallKey -Recurse -Force -ErrorAction SilentlyContinue
Remove-Item -LiteralPath $Dest -Recurse -Force -ErrorAction SilentlyContinue

try {
    Add-Type -Namespace BlueMoleculeSetup -Name Shell -MemberDefinition '[DllImport("shell32.dll")] public static extern void SHChangeNotify(int eventId, int flags, IntPtr item1, IntPtr item2);'
    [BlueMoleculeSetup.Shell]::SHChangeNotify(0x08000000, 0, [IntPtr]::Zero, [IntPtr]::Zero)
} catch { }

if (Test-Path -LiteralPath $Dest) { Write-Host ('Some files could not be removed (is BlueMolecule open?): ' + $Dest) }
else { Write-Host 'BlueMolecule was removed.' }
Write-Host ('Your settings remain in ' + (Join-Path $env:LOCALAPPDATA 'BlueMolecule') + '; delete that folder to remove them too.')
if ([Environment]::UserInteractive) { [void](Read-Host 'Press Enter to close') }
