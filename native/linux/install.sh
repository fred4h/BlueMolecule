#!/usr/bin/env bash
# Install BlueMolecule for the current user (default prefix ~/.local). Running it again updates or
# repairs an installation (also one of BlueMolecule 1.0/1.1).
# System-wide:  sudo PREFIX=/usr/local ./install.sh
# BLUEMOLECULE_NO_PROMPT=1 skips the question about installing Qt WebEngine for Python.
set -euo pipefail
SRC="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PREFIX="${PREFIX:-$HOME/.local}"
APP_ID=org.bluemolecule.BlueMolecule
SHARE="$PREFIX/share"
DATA="$SHARE/bluemolecule"
BIN="$PREFIX/bin/bluemolecule"
SIZES="48 64 128 256 512"

# ---- leftovers of BlueMolecule 1.0/1.1: other desktop-file, icon and MIME names, and the old share layout
# (a capitalised share/BlueMolecule folder, a launcher called bluemolecule.py). Only those files go; the
# folders themselves only when empty, so settings kept in $DATA (browser profiles) stay.
rm -f "$SHARE/applications/bluemolecule.desktop" "$SHARE/applications/BlueMolecule.desktop" \
      "$SHARE/mime/packages/bluemolecule.xml" "$SHARE/pixmaps/bluemolecule.png" "$SHARE/pixmaps/BlueMolecule.png" \
      "$SHARE/icons/hicolor/scalable/apps/bluemolecule.svg" "$DATA/bluemolecule.py"
for d in "$SHARE"/icons/hicolor/*/apps; do rm -f "$d/bluemolecule.png" "$d/BlueMolecule.png"; done
if [ "$SHARE/BlueMolecule" -ef "$DATA" ]; then :   # a case-insensitive file system: that is $DATA itself
else
  rm -f "$SHARE/BlueMolecule/BlueMolecule.html" "$SHARE/BlueMolecule/bluemolecule" "$SHARE/BlueMolecule/bluemolecule.py"
  rmdir "$SHARE/BlueMolecule" 2>/dev/null || true
fi

# ---- this version
mkdir -p "$DATA" "$PREFIX/bin" "$SHARE/applications" "$SHARE/mime/packages"
install -m 644 "$SRC/BlueMolecule.html" "$DATA/BlueMolecule.html"
install -m 755 "$SRC/bluemolecule" "$DATA/bluemolecule"
install -m 755 "$SRC/bluemolecule_qt.py" "$DATA/bluemolecule_qt.py"
rm -f "$BIN"            # 1.0 put a copy of the launcher here; now it is a link
ln -s "$DATA/bluemolecule" "$BIN"
for s in $SIZES; do
  mkdir -p "$SHARE/icons/hicolor/${s}x${s}/apps"
  install -m 644 "$SRC/icons/$APP_ID-$s.png" "$SHARE/icons/hicolor/${s}x${s}/apps/$APP_ID.png"
done
install -m 644 "$SRC/$APP_ID.xml" "$SHARE/mime/packages/$APP_ID.xml"
esc="$(printf '%s' "$BIN" | sed -e 's/[\\|&]/\\&/g')"
sed -e "s|^Exec=.*|Exec=\"$esc\" %F|" "$SRC/$APP_ID.desktop" > "$SHARE/applications/$APP_ID.desktop"
chmod 644 "$SHARE/applications/$APP_ID.desktop"

# ---- let the desktop see the changes now (KDE keeps its own cache of desktop files: ksycoca)
if command -v update-mime-database >/dev/null 2>&1; then update-mime-database "$SHARE/mime" >/dev/null 2>&1 || true; fi
if command -v update-desktop-database >/dev/null 2>&1; then update-desktop-database "$SHARE/applications" >/dev/null 2>&1 || true; fi
if [ -f "$SHARE/icons/hicolor/index.theme" ] && command -v gtk-update-icon-cache >/dev/null 2>&1; then
  gtk-update-icon-cache -q "$SHARE/icons/hicolor" >/dev/null 2>&1 || true
else
  touch "$SHARE/icons/hicolor" 2>/dev/null || true
fi
if command -v kbuildsycoca6 >/dev/null 2>&1; then kbuildsycoca6 >/dev/null 2>&1 || true
elif command -v kbuildsycoca5 >/dev/null 2>&1; then kbuildsycoca5 >/dev/null 2>&1 || true; fi

echo "BlueMolecule installed in $DATA"

# ---- the native Qt window needs Qt WebEngine for Python; offer the install command, never run it unasked
if command -v python3 >/dev/null 2>&1; then
  if qt="$("$DATA/bluemolecule" --check-qt 2>/dev/null)"; then
    echo "Window: native Qt WebEngine window ($qt)"
  else
    hint="$(printf '%s\n' "$qt" | tail -n 1)"
    if [ -n "$hint" ]; then
      echo "Window: a browser app window, until Qt WebEngine for Python is installed. For the native window"
      echo "(its own taskbar entry, recommended on KDE) run:"
      echo "    $hint"
      if [ -t 0 ] && [ -t 1 ] && [ -z "${BLUEMOLECULE_NO_PROMPT:-}" ]; then
        answer=
        read -r -p "Run it now? [y/N] " answer || answer=
        case "$answer" in
          [yY]*)
            read -r -a cmd <<< "$hint"
            if "${cmd[@]}" && qt="$("$DATA/bluemolecule" --check-qt 2>/dev/null)"; then
              echo "Window: native Qt WebEngine window ($qt)"
            else
              echo "Qt WebEngine is still missing; BlueMolecule uses a browser window meanwhile."
            fi ;;
        esac
      fi
    fi
  fi
fi

echo "Run:  bluemolecule [files...]   (or start BlueMolecule from the applications menu; chemistry"
echo "      files such as .xyz, .cube, .fchk, .molden, .pdb or .gjf offer 'Open With BlueMolecule')"
echo "      bluemolecule --diagnose    shows which window it uses and why"
case ":$PATH:" in *":$PREFIX/bin:"*) ;; *) echo "Note: $PREFIX/bin is not on your PATH; add  export PATH=\"$PREFIX/bin:\$PATH\"  to ~/.bashrc";; esac
