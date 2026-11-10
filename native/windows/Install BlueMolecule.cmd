@echo off
rem Install BlueMolecule for the current user (Windows, experimental): double-click this file.
rem It runs windows\install.ps1 with Windows PowerShell; no administrator rights are needed.
setlocal
if not exist "%~dp0windows\install.ps1" (
  echo The installer files were not found next to this file.
  echo Extract the whole zip file first ^(right-click it, Extract All^), then run Install BlueMolecule.cmd again.
  echo.
  pause
  exit /b 1
)
powershell.exe -NoProfile -ExecutionPolicy Bypass -File "%~dp0windows\install.ps1"
echo.
pause
