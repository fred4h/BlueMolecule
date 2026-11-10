# How BlueMolecule was made: a vibe-coded project

**BlueMolecule is vibe-coded.** Its author described in plain language what they wanted from a molecule
viewer:
- the features of the desktop programs they use for quantum-chemistry work;
- the file formats of their calculations;
- the problems they ran into.

An AI assistant, **Claude** (made by Anthropic), wrote the program, its tests and this documentation. The author
did not write or review the code line by line. They ran each version, looked at what it did and asked for changes.

This page explains what that means for you as a user or contributor.

## What "vibe coding" means here

"Vibe coding" is a name, coined in 2025, for building software by describing it to an AI and judging the result by
using it, rather than by reading the code. BlueMolecule is vibe-coded in that full sense, not "a person's code
with some AI help":

| The author | The AI assistant |
|---|---|
| Set the goals: a free viewer with the everyday functions of desktop quantum-chemistry viewers, plus building tools | Wrote all the code: about 18,000 lines of JavaScript, Python, Swift, HTML and CSS, tests included |
| Chose the features, using public feature lists of programs such as Chemcraft and Avogadro | Designed the architecture, the file readers and writers, the renderer and the interface |
| Tried each version on their own computers and reported problems, such as "a 600 MB ORCA log opens as empty" or "the ONIOM energy graph jumps back and forth" | Investigated each report, fixed the cause and added a test for it |
| Decided to publish the project and chose the license | Generated the synthetic example files, wrote about 8,000 automated checks, wrote the documentation and checked the code before publication |

Claude worked in a sandboxed workspace with an editor, a terminal and a headless browser, so it could run the tests
and look at screenshots of the interface. The versions from 1.0 to 1.4.1 were made between 4 and 10 November 2026
(see the [changelog](../CHANGELOG.md)).

## How the program was checked

- **Automated tests.** Every feature came with tests: about 8,000 checks in 25 suites, run with
  `python3 tests/run_all.py --browser` ([developer guide](developer-guide.md#tests)). They cover:
  - the readers and writers of each format;
  - format detection;
  - the chemistry routines (SMILES, bond orders, hydrogens, symmetry, the clean-up force field);
  - the selection language;
  - the launcher and the packaging;
  - the interface itself, driven in a headless browser.
- **Synthetic example files.** No real outputs were available while BlueMolecule was written. So each test file
  was generated to imitate a program's output layout, as that layout is publicly documented and commonly seen.
  This is the **main limitation**. Real files vary in ways the examples may not cover: program versions, print
  levels, keywords, warnings and unusual job types.
- **Screenshots.** The browser tests save screenshots of the interface. These were inspected during development.
- **Use by the author.** Problems found while the author used the program led to versions 1.1 to 1.4.
- **A publication audit** before the first public release. It checked for:
  - personal data, credentials and environment paths;
  - code or data copied from other programs;
  - the safety of the local launcher;
  - script injection from crafted files, which is now part of the test suite.

What has and has not been verified is listed in [Status and limitations](limitations.md).

## Clean-room statement

BlueMolecule contains **no code, data, icons or text from Chemcraft, GaussView, Avogadro, VMD, Jmol, Open Babel or
any other program**. No program was decompiled or disassembled to make it.
- The file readers were written from public descriptions of the file layouts. The reader modules say so in their
  headers.
- Feature ideas came from what such programs are publicly known to do.
- Open Babel is optional. When it is installed, BlueMolecule runs it as a separate program for extra conversions.
- The reference data comes from the scientific literature: radii, masses, colours and glycan symbols. The sources are
  listed in the [README](../README.md#license-and-credits).

## What this means for you

- **Check numbers that matter against the output file.** BlueMolecule displays what your quantum-chemistry program
  computed: energies, frequencies, intensities, excitation energies, shieldings. It does not compute them. The
  displayed values were read by code that nobody reviewed by hand. Confirm any value you publish in the original
  output.
- **Spectra are simulated from the computed lines**, using the broadening and scaling you choose
  ([Analysing outputs](analysing-outputs.md#spectra)). They are a visual aid, not an independent calculation.
- **The clean-up force field is not a physical method.** It is a generic model for tidying sketched structures
  before a real calculation ([Building and editing](building-and-editing.md#clean-up)).
- **Report files that are read wrongly.** A small example file is the most useful bug report
  ([CONTRIBUTING.md](../CONTRIBUTING.md)).

## Copyright and license

BlueMolecule is released under the [MIT License](../LICENSE). You may use, change and share it, including
commercially, as long as you keep the copyright and license notice.

Copyright law for AI-generated code is still unsettled. For example, the U.S. Copyright Office has said that
material produced by AI without enough human creative control is not protected by copyright. The license grants
its permissions to the extent the author holds rights in the work; either way, you can use the code under the MIT
terms. This note is not legal advice.

## Contributing to a vibe-coded project

Contributions written by hand and contributions written with AI tools are both welcome, on the same terms:
- the change must come with tests;
- the full test suite must pass;
- the pull request should say how it was made.

See [CONTRIBUTING.md](../CONTRIBUTING.md).
