# Contributing to BlueMolecule

Thank you for helping. BlueMolecule is a [vibe-coded](docs/vibe-coding.md) project that was tested almost
entirely on synthetic files, so **reports from real use are the most valuable contribution of all**.

## Reporting a file that is read wrongly

[Open an issue](https://github.com/fred4h/BlueMolecule/issues/new) and include:

1. The **program and version** that wrote the file, e.g. Gaussian 16 C.01 or ORCA 6.0.1, and the job type.
2. **What you expected and what BlueMolecule showed.** Screenshots help.
3. **An example file**, if you may share it. A trimmed file is best: one small molecule, or the first and last few
   optimisation cycles.
4. The BlueMolecule version (shown in the Help dialog).

> [!CAUTION]
> Issues are **public**. Do not attach unpublished results, confidential structures or files you are not allowed
> to share. Remove hostnames, user names and paths if they matter to you. Gaussian and ORCA outputs contain them in
> their headers and archive entries. When in doubt, make a small test calculation that reproduces the problem.

## Reporting other bugs

Include:
- your operating system and desktop, e.g. Kubuntu 26.10 with KDE Plasma on Wayland;
- the output of `bluemolecule --diagnose`;
- the steps to reproduce;
- any error from the log (`~/.local/state/bluemolecule/bluemolecule.log`) or from the developer tools (**F12** in the
  Qt window).

## Contributing code

Contributions written by hand and contributions written with AI tools are both welcome, on the same terms:

1. **Tests first.** Every fix comes with a test that fails without it. Every new reader comes with a fixture and
   checks. See the [developer guide](docs/developer-guide.md#tests).
2. **Run the whole suite** before opening a pull request: `python3 tests/run_all.py --browser`.
3. **Clean room.** Do not copy code, data or text from other programs, especially not from closed-source viewers.
   Write readers from public format descriptions and real example files.
4. **Say how it was made.** Mention in the pull request whether the change was written with an AI assistant, and
   which one.
5. **Keep it dependency-free.** The viewer has no runtime dependencies, and the launcher uses only the Python
   standard library.
6. **Keep the changes focused.** One topic per pull request, with the reason in the description.

## Documentation

The documentation is in `docs/`. Corrections are welcome, especially where the documentation and the program
disagree. Screenshots in `docs/images/` are produced by `node tools/docs_screenshots.js` after a build.
