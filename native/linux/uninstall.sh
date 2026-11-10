#!/usr/bin/env bash
# Remove BlueMolecule (use the same PREFIX as for install.sh; default ~/.local). Safe to run again.
# Settings kept by the BlueMolecule window profiles (browser, firefox, qtwebengine) stay unless you delete them yourself.
set -euo pipefail
PREFIX="${PREFIX:-$HOME/.local}"
APP_ID=org.bluemolecule.BlueMolecule
SHARE="$PREFIX/share"
DATA="$SHARE/bluemolecule"
rm -f "$PREFIX/bin/bluemolecule" "$SHARE/applications/$APP_ID.desktop" "$SHARE/mime/packages/$APP_ID.xml"
for d in "$SHARE"/icons/hicolor/*/apps; do rm -f "$d/$APP_ID.png" "$d/bluemolecule.png" "$d/BlueMolecule.png"; done
# leftovers of BlueMolecule 1.0/1.1 (see install.sh)
rm -f "$SHARE/applications/bluemolecule.desktop" "$SHARE/applications/BlueMolecule.desktop" \
      "$SHARE/mime/packages/bluemolecule.xml" "$SHARE/pixmaps/bluemolecule.png" "$SHARE/pixmaps/BlueMolecule.png" \
      "$SHARE/icons/hicolor/scalable/apps/bluemolecule.svg"
if ! [ "$SHARE/BlueMolecule" -ef "$DATA" ]; then
  rm -f "$SHARE/BlueMolecule/BlueMolecule.html" "$SHARE/BlueMolecule/bluemolecule" "$SHARE/BlueMolecule/bluemolecule.py"
  rmdir "$SHARE/BlueMolecule" 2>/dev/null || true
fi
rm -f "$DATA/BlueMolecule.html" "$DATA/bluemolecule" "$DATA/bluemolecule_qt.py" "$DATA/bluemolecule.py"
rm -rf "$DATA/__pycache__"
rmdir "$DATA" 2>/dev/null || true
if command -v update-mime-database >/dev/null 2>&1; then update-mime-database "$SHARE/mime" >/dev/null 2>&1 || true; fi
if command -v update-desktop-database >/dev/null 2>&1; then update-desktop-database "$SHARE/applications" >/dev/null 2>&1 || true; fi
touch "$SHARE/icons/hicolor" 2>/dev/null || true
if command -v kbuildsycoca6 >/dev/null 2>&1; then kbuildsycoca6 >/dev/null 2>&1 || true
elif command -v kbuildsycoca5 >/dev/null 2>&1; then kbuildsycoca5 >/dev/null 2>&1 || true; fi
echo "BlueMolecule removed from $PREFIX"
DATA_HOME="${XDG_DATA_HOME:-$HOME/.local/share}"
if [ -d "$DATA_HOME/bluemolecule" ]; then
  echo "Your BlueMolecule settings remain in $DATA_HOME/bluemolecule (delete that folder to remove them)."
fi
