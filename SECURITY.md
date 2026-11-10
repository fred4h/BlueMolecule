# Security

## How BlueMolecule handles your files

- **Everything runs on your computer.** BlueMolecule never uploads files and makes no network connections beyond
  its own local launcher.
- **The launcher** (`launcher/bluemolecule.py`), used for the browser-window types, runs a small web server that:
  - listens on `127.0.0.1` only, so other computers cannot reach it;
  - puts every address under a random token that changes each session;
  - checks the `Host` header (against DNS rebinding) and the `Origin` of every request that changes something;
  - serves only the viewer and the files you named on the command line;
  - stops when the window closes.
- **Open Babel** is run as a separate program with an argument list, never through a shell. Format codes are
  validated, and options that could read or write other files are refused.
- **The native windows** (the Linux Qt window and the macOS app) let the page read only the files you opened, and
  save only where you choose in the Save dialog.
- **On Windows**, there is no server at all. The launcher copies the files you open into a temporary page in your
  profile, and only reads the registry, to find the browser. The installer writes only to your user profile and
  to `HKEY_CURRENT_USER`.
- **Text from files is never treated as markup.** Titles, atom and residue names, labels and file names are escaped
  before they are displayed. `tests/browser_xss_test.js` checks this with crafted files.
- **The installers** never ask for administrator rights on their own. The Linux installer offers the command that
  installs Qt WebEngine and asks before running it with `sudo`. The Windows installer needs no administrator rights
  at all.

## Reporting a vulnerability

Please **do not open a public issue** for a security problem. Use GitHub's private vulnerability reporting instead:
on the repository page, choose **Security → Report a vulnerability**. If that is not available, open an issue that
only asks for a private contact, without details.

Please include the version, your platform, and the steps or a file that shows the problem.

BlueMolecule is maintained by one person in their spare time, with AI assistance. Reports are handled as quickly as
possible, but there is no guaranteed response time.
