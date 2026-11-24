# Security

## How BlueMolecule handles your files

- **Everything runs on your computer.** BlueMolecule never uploads files. It uses the network only when you ask:
  - the desktop program's update check (**Help → Check for updates…**, or at start-up if you tick that option): it
    asks GitHub for the latest release, and installs an update only after checking its signature against
    BlueMolecule's own key ([Updates](docs/installation.md#updates));
  - **Download** (PubChem and the PDB), which is off until you tick **Allow downloads** in its dialog. A download
    sends only the name, CID or PDB ID you typed, as part of the address, to `pubchem.ncbi.nlm.nih.gov` or
    `files.rcsb.org`, without cookies or a referrer, and only when you press **Download**. What comes back is read like a file you
    opened ([Downloading](docs/user-guide.md#downloading-from-pubchem-and-the-pdb)).
- **The desktop program** (2.0 and later):
  - reads only the files you opened (with the Open dialog, **Open with**, the command line or by dropping them on
    the window), and saves only where you choose in the Save dialog;
  - shows only BlueMolecule's own page in its window; links open in your web browser;
  - installs for your user only on Windows, without administrator rights. On Linux the AppImage needs no
    installation, and the `.deb` package is installed with your system's package tools, which ask for your password.
- **Text from files is never treated as markup.** Titles, atom and residue names, labels and file names are escaped
  before they are displayed. `tests/browser_xss_test.js` checks this with crafted files.
- **Open Babel** is run as a separate program with an argument list, never through a shell. Format codes are
  validated, and options that could read or write other files are refused.
- **The single HTML file** needs no network, except for downloads you ask for. **The launcher** of the Linux package (`launcher/bluemolecule.py`),
  which opens it in a browser window from the command line, runs a small web server that:
  - listens on `127.0.0.1` only, so other computers cannot reach it;
  - puts every address under a random token that changes each session;
  - checks the `Host` header (against DNS rebinding) and the `Origin` of every request that changes something;
  - serves only the viewer and the files you named on the command line;
  - stops when the window closes.

  Its Qt window lets the page read only the files you opened, and save only where you choose. Its install script
  never runs `sudo` on its own: it offers the command that installs Qt WebEngine and asks first.

## Reporting a vulnerability

Please **do not open a public issue** for a security problem. Use GitHub's private vulnerability reporting instead:
on the repository page, choose **Security → Report a vulnerability**. If that is not available, open an issue that
only asks for a private contact, without details.

Please include the version, your platform, and the steps or a file that shows the problem.

BlueMolecule is maintained by one person in their spare time, with AI assistance. Reports are handled as quickly as
possible, but there is no guaranteed response time.
