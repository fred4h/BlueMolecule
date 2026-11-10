#!/usr/bin/env python3
"""Run BlueMolecule's test suites and print one line per suite.

    python3 tests/run_all.py             unit tests: Node.js (readers, writers, chemistry, selections...) and Python (launcher, packaging)
    python3 tests/run_all.py --browser   also the browser tests (needs Playwright with Chromium); builds dist/ first
    python3 tests/run_all.py --big       also the 600 MB / 1.2 GB large-file cases (implies --browser; slow, needs ~4 GB RAM and disk)

A failing suite's full output is printed. The exit status is 0 only when every suite that ran passed.
"""
import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
T = ROOT / 'tests'

NODE_UNIT = ['test_parsers.js', 'test_parsers2.js', 'test_readers_qc.js', 'test_readers_struct.js', 'test_writers.js',
             'test_writers_struct.js', 'test_chem.js', 'test_ff.js', 'test_templates.js', 'test_symmetry.js', 'test_analysis.js',
             'test_audit.js', 'test_qmmm.js', 'test_native_bridge.js', 'test_reps.js', 'test_cartoon.js']
PY_UNIT = ['test_launcher.py', 'test_native_packaging.py']
BROWSER = ['test_renderer.js', 'browser_test.js', 'browser_build_test.js', 'browser_formats_test.js', 'browser_reps_test.js',
           'browser_xss_test.js', 'test_bigfiles.js']


def run(cmd, label):
    t0 = time.time()
    p = subprocess.run(cmd, cwd=str(ROOT), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
    out = p.stdout.strip()
    last = out.splitlines()[-1] if out else ''
    ok = p.returncode == 0
    print('%-4s %-28s %6.1f s  %s' % ('ok' if ok else 'FAIL', label, time.time() - t0, last[:110]))
    if not ok:
        print('\n'.join('     | ' + s for s in out.splitlines()[-60:]))
    return ok


def have_playwright():
    p = subprocess.run(['node', '-e', "require('playwright')"], cwd=str(ROOT), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return p.returncode == 0


def main():
    ap = argparse.ArgumentParser(description='Run the BlueMolecule test suites.')
    ap.add_argument('--browser', action='store_true', help='also run the browser tests (Playwright + Chromium)')
    ap.add_argument('--big', action='store_true', help='also run the large-file cases (implies --browser)')
    args = ap.parse_args()
    results = []
    browser = args.browser or args.big
    if browser:   # first, so that the packaging test also checks the freshly built package
        results.append(run([sys.executable, str(ROOT / 'tools' / 'build.py')], 'build (tools/build.py)'))
    for f in NODE_UNIT:
        results.append(run(['node', str(T / f)], f))
    for f in PY_UNIT:
        results.append(run([sys.executable, str(T / f)], f))
    if browser:
        if not have_playwright():
            print('skip browser tests: Playwright is not installed (npm install playwright && npx playwright install chromium)')
        else:
            for f in BROWSER:
                results.append(run(['node', str(T / f)] + (['--big'] if args.big and f == 'test_bigfiles.js' else []), f))
    failed = results.count(False)
    print('\n%d suite(s) passed, %d failed' % (len(results) - failed, failed))
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
