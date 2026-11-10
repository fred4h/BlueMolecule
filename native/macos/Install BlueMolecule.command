#!/bin/bash
# Install BlueMolecule on this Mac: double-click this file once
# (or run it in Terminal:  bash "Install BlueMolecule.command").
#
# It compiles the small native BlueMolecule app from native/macos/BlueMolecule.swift with Apple's
# Command Line Tools, signs it ad hoc, installs it in ~/Applications (or /Applications) and adds the
# "bluemolecule" terminal command. Running it again simply updates the installation. If the native
# build is not possible, the simpler browser-window version (fallback/BlueMolecule.app) is installed
# instead and the reason is shown.
# Optional: BLUEMOLECULE_INSTALL_DIR=/some/Applications and BLUEMOLECULE_BIN_DIR=/some/bin choose the
# folders without asking. Written for the bash 3.2 that ships with macOS.

VERSION="@VERSION@"
HERE="$(cd "$(dirname "$0")" && pwd)"
SWIFT_SRC="$HERE/native/macos/BlueMolecule.swift"
PLIST_SRC="$HERE/native/macos/Info.plist"
FALLBACK_APP="$HERE/fallback/BlueMolecule.app"
LSREGISTER="/System/Library/Frameworks/CoreServices.framework/Frameworks/LaunchServices.framework/Support/lsregister"
LOG="$HOME/Library/Logs/BlueMolecule-install.log"
INTERACTIVE=0
[ -t 0 ] && INTERACTIVE=1

step() { printf '\n\033[1m==> %s\033[0m\n' "$*"; echo "==> $*" >> "$LOG"; }
info() { printf '    %s\n' "$*"; echo "    $*" >> "$LOG"; }
warn() { printf '    \033[33m%s\033[0m\n' "$*"; echo "!!  $*" >> "$LOG"; }
ask() {  # ask "question" y|n  -> status 0 for yes; without a terminal the default is taken
  local answer
  if [ "$INTERACTIVE" != 1 ]; then [ "$2" = y ]; return; fi
  if [ "$2" = y ]; then read -r -p "    $1 [Y/n] " answer; else read -r -p "    $1 [y/N] " answer; fi
  case "$answer" in
    [yY]*) return 0 ;;
    [nN]*) return 1 ;;
    *) [ "$2" = y ] ;;
  esac
}
finish() {
  if [ "$INTERACTIVE" = 1 ]; then echo; read -n 1 -s -r -p "Press any key to close this window."; echo; fi
  exit "$1"
}
sed_escape() { printf '%s' "$1" | sed -e 's/[\\|&]/\\&/g'; }

mkdir -p "$HOME/Library/Logs"
echo "==== $(date)  BlueMolecule $VERSION installer, from $HERE" >> "$LOG"
printf '\033[1mBlueMolecule %s installer\033[0m   (log: %s)\n' "$VERSION" "$LOG"

if [ "$(uname -s)" != Darwin ]; then
  warn "This installer is for macOS. On Linux use the BlueMolecule Linux package (install.sh)."
  finish 1
fi
for f in "$SWIFT_SRC" "$PLIST_SRC" "$HERE/BlueMolecule.html" "$HERE/AppIcon.icns" "$HERE/bluemolecule.py" \
         "$FALLBACK_APP/Contents/MacOS/BlueMolecule"; do
  if [ ! -f "$f" ]; then
    warn "Missing: $f"
    warn "Unpack the whole BlueMolecule zip and run the installer from inside that folder."
    finish 1
  fi
done

WORK="$(mktemp -d "${TMPDIR:-/tmp}/bluemolecule-build.XXXXXX")" || { warn "Cannot create a temporary folder."; finish 1; }
trap 'rm -rf "$WORK"' EXIT

# ---------------------------------------------------------------- where to install
choose_destination() {
  DEST="$HOME/Applications"
  if [ -n "${BLUEMOLECULE_INSTALL_DIR:-}" ]; then
    DEST="$BLUEMOLECULE_INSTALL_DIR"
  elif [ -w /Applications ]; then
    if [ -d /Applications/BlueMolecule.app ]; then
      if ask "BlueMolecule is already in /Applications. Update it there?" y; then DEST=/Applications; fi
    elif ask "Install in /Applications (all users of this Mac) instead of your own ~/Applications?" n; then
      DEST=/Applications
    fi
  fi
  mkdir -p "$DEST" || { warn "Cannot create $DEST"; finish 1; }
  DEST="$(cd "$DEST" && pwd)"   # absolute: the terminal command stores this path
  APP="$DEST/BlueMolecule.app"
}

quit_running_copy() {
  if pgrep -f "BlueMolecule.app/Contents/MacOS/BlueMolecule" >/dev/null 2>&1; then
    if ask "BlueMolecule is running and will be closed for the update. Continue?" y; then
      osascript -e 'tell application id "org.bluemolecule.viewer" to quit' >/dev/null 2>&1 || true
      sleep 2
    else
      finish 1
    fi
  fi
}

# ---------------------------------------------------------------- native build
FAIL_REASON=""
FAIL_LOG=""
BUILT_APP=""

build_native() {
  step "1/5  Checking for Apple's Command Line Tools (Swift compiler)"
  if ! xcode-select -p >/dev/null 2>&1 || ! xcrun --find swiftc >/dev/null 2>&1; then
    warn "The Swift compiler is not installed yet; it comes with Apple's free Command Line Tools."
    xcode-select --install >/dev/null 2>&1 || true
    info "A macOS dialog should now offer to install the Command Line Tools. Click \"Install\" and wait"
    info "until it has finished (a few minutes; no Apple ID needed). Then double-click"
    info "\"Install BlueMolecule.command\" again to get the native app."
    info "(If no dialog appeared, run  xcode-select --install  in Terminal.)"
    FAIL_REASON="the Command Line Tools are not installed yet"
    return 1
  fi
  info "swiftc: $(xcrun --find swiftc)"
  xcrun swiftc --version 2>&1 | head -n 1 | while read -r line; do info "$line"; done

  step "2/5  Compiling BlueMolecule (about a minute)"
  local sdk arch sdk_args ok_arm64=0 ok_x86_64=0
  sdk="$(xcrun --sdk macosx --show-sdk-path 2>/dev/null)"
  sdk_args=()
  if [ -n "$sdk" ]; then sdk_args=(-sdk "$sdk"); fi
  for arch in arm64 x86_64; do
    info "compiling for $arch ..."
    rm -f "$WORK/BlueMolecule-$arch"
    if xcrun swiftc -O -swift-version 5 -target "$arch-apple-macos11.0" "${sdk_args[@]}" \
         -framework Cocoa -framework WebKit -o "$WORK/BlueMolecule-$arch" "$SWIFT_SRC" \
         > "$WORK/swiftc-$arch.log" 2>&1 && [ -f "$WORK/BlueMolecule-$arch" ]; then
      if [ "$arch" = arm64 ]; then ok_arm64=1; else ok_x86_64=1; fi
    else
      warn "could not compile for $arch"
    fi
    cat "$WORK/swiftc-$arch.log" >> "$LOG"
  done
  local machine
  machine="$(uname -m)"
  if [ "$ok_arm64" = 1 ] && [ "$ok_x86_64" = 1 ] &&
     lipo -create -output "$WORK/BlueMolecule" "$WORK/BlueMolecule-arm64" "$WORK/BlueMolecule-x86_64" >> "$LOG" 2>&1; then
    info "universal app (Apple silicon and Intel)"
  elif [ "$ok_arm64" = 1 ] && [ "$machine" = arm64 ]; then
    cp "$WORK/BlueMolecule-arm64" "$WORK/BlueMolecule"; info "app for Apple silicon"
  elif [ "$ok_x86_64" = 1 ] && [ "$machine" = x86_64 ]; then
    cp "$WORK/BlueMolecule-x86_64" "$WORK/BlueMolecule"; info "app for Intel"
  elif [ "$ok_arm64" = 1 ]; then
    cp "$WORK/BlueMolecule-arm64" "$WORK/BlueMolecule"; info "app for Apple silicon"
  else
    FAIL_REASON="the Swift compiler reported errors"
    FAIL_LOG="$WORK/swiftc-arm64.log"
    return 1
  fi

  step "3/5  Assembling BlueMolecule.app"
  local app="$WORK/BlueMolecule.app"
  if ! { mkdir -p "$app/Contents/MacOS" "$app/Contents/Resources" &&
         cp "$WORK/BlueMolecule" "$app/Contents/MacOS/BlueMolecule" &&
         chmod 755 "$app/Contents/MacOS/BlueMolecule" &&
         sed "s/@VERSION@/$VERSION/g" "$PLIST_SRC" > "$app/Contents/Info.plist" &&
         printf 'APPL????' > "$app/Contents/PkgInfo" &&
         cp "$HERE/AppIcon.icns" "$HERE/BlueMolecule.html" "$app/Contents/Resources/" &&
         cp "$HERE/bluemolecule.py" "$app/Contents/Resources/bluemolecule.py" &&
         chmod 755 "$app/Contents/Resources/bluemolecule.py"; }; then
    FAIL_REASON="the app bundle could not be assembled"
    return 1
  fi
  if ! plutil -lint "$app/Contents/Info.plist" >> "$LOG" 2>&1; then
    FAIL_REASON="Info.plist is not valid (see $LOG)"
    return 1
  fi
  xattr -cr "$app" 2>/dev/null || true

  step "4/5  Signing the app (ad hoc, for this Mac)"
  if ! codesign --force --deep -s - "$app" >> "$LOG" 2>&1; then
    FAIL_REASON="codesign failed (see $LOG)"
    return 1
  fi
  BUILT_APP="$app"
}

install_bundle() {  # $1 = app bundle to install as $APP
  quit_running_copy
  if [ -e "$APP" ] && ! rm -rf "$APP"; then
    warn "Cannot remove the old $APP"
    return 1
  fi
  ditto "$1" "$APP" || return 1
  xattr -cr "$APP" 2>/dev/null || true
  "$LSREGISTER" -f "$APP" >/dev/null 2>&1 || true
  touch "$APP"
}

# ---------------------------------------------------------------- command-line tool
install_cli() {  # $1 = native | fallback
  local dir="${BLUEMOLECULE_BIN_DIR:-}" d
  if [ -n "$dir" ]; then mkdir -p "$dir" && dir="$(cd "$dir" && pwd)" || return 1; fi
  for d in /opt/homebrew/bin /usr/local/bin "$HOME/.local/bin" "$HOME/bin"; do
    if [ -z "$dir" ] && [ -d "$d" ] && [ -w "$d" ]; then dir="$d"; fi
  done
  if [ -z "$dir" ]; then dir="$HOME/.local/bin"; mkdir -p "$dir" || return 1; fi
  local cli="$dir/bluemolecule"
  if [ "$1" = native ]; then
    cat > "$WORK/cli" <<'EOF'
#!/bin/bash
# bluemolecule: open files in BlueMolecule.app (installed by "Install BlueMolecule.command").
APP="@APP@"
[ -d "$APP" ] || APP="BlueMolecule"
case "$1" in
  -h|--help) echo "usage: bluemolecule [FILE ...]    opens the files in BlueMolecule"; exit 0 ;;
  -V|--version) echo "BlueMolecule @VERSION@"; exit 0 ;;
esac
files=()
for f in "$@"; do
  if [ -e "$f" ]; then
    case "$f" in /*) files+=("$f") ;; *) files+=("$PWD/$f") ;; esac
  else
    echo "bluemolecule: $f: no such file" >&2
  fi
done
if [ ${#files[@]} -eq 0 ]; then exec open -a "$APP"; fi
exec open -a "$APP" "${files[@]}"
EOF
  else
    cat > "$WORK/cli" <<'EOF'
#!/bin/bash
# bluemolecule: open files in BlueMolecule (browser-window version; needs python3).
exec python3 "@APP@/Contents/Resources/bluemolecule.py" "$@"
EOF
  fi
  rm -f "$cli"   # BlueMolecule 1.1 made this a symlink into the app bundle: never write through it
  if ! sed -e "s|@APP@|$(sed_escape "$APP")|g" -e "s|@VERSION@|$VERSION|g" "$WORK/cli" > "$cli"; then
    warn "Could not write $cli"
    return 1
  fi
  chmod 755 "$cli"
  info "terminal command: $cli"
  case ":$PATH:" in
    *":$dir:"*) ;;
    *) warn "$dir is not on your PATH. Add it with:  echo 'export PATH=\"$dir:\$PATH\"' >> ~/.zprofile" ;;
  esac
}

# ---------------------------------------------------------------- main
choose_destination
info "installing into $DEST"
MODE=""
if build_native; then
  step "5/5  Installing $APP"
  if install_bundle "$BUILT_APP"; then MODE=native; else FAIL_REASON="the app could not be copied to $DEST"; fi
fi

if [ "$MODE" != native ]; then
  warn "The native app was not installed: $FAIL_REASON."
  if [ -n "$FAIL_LOG" ] && [ -s "$FAIL_LOG" ]; then
    echo "    ---- compiler output (also in $LOG) ----"
    sed -n '1,60p' "$FAIL_LOG" | sed 's/^/    /'
    echo "    ----------------------------------------"
  fi
  step "Installing the browser-window version of BlueMolecule instead"
  if ! install_bundle "$FALLBACK_APP"; then
    warn "Could not install $APP either. Please send $LOG to whoever gave you BlueMolecule."
    finish 1
  fi
  MODE=fallback
  info "It opens BlueMolecule in Chrome, Edge, Brave, Chromium or Vivaldi (else your default browser)."
fi

install_cli "$MODE" || warn "The terminal command could not be installed."
for other in "/Applications/BlueMolecule.app" "$HOME/Applications/BlueMolecule.app"; do
  if [ "$other" != "$APP" ] && [ -d "$other" ]; then
    warn "An older copy is still at $other; you can move it to the Trash."
  fi
done

step "Done"
if [ "$MODE" = native ]; then
  info "BlueMolecule $VERSION is installed: $APP"
  info "Open chemistry files from Finder (Open With > BlueMolecule), drop them on the Dock icon,"
  info "or run:  bluemolecule calc.log freq.out orbital.cube"
else
  info "BlueMolecule $VERSION (browser-window version) is installed: $APP"
  info "Run this installer again after fixing the problem above to get the native app."
fi
echo "==== finished ($MODE)" >> "$LOG"
if ask "Open BlueMolecule now?" y; then open "$APP"; fi
finish 0
