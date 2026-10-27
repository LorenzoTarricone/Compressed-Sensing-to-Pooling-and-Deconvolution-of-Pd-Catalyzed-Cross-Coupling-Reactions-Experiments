"""Compare a fresh run (output/) with the results reported in the paper (expected/).

Also checks the input files against the SHA-256 hashes in Data/paper/pooled_plates/manifest.json.
Exit status 0 when everything agrees within the stated tolerances.

    python check.py [--skip-simulations]
"""
import argparse
import hashlib
import json
import os
import sys

import numpy as np
import pandas as pd

from common import EXPECTED, OUT, POOLED

# tolerance per file: scores are printed to one or two decimals in the paper, statistics to two or three;
# the simulations count successes out of 120 trials and may shift by a trial or two with another solver build
TOLERANCE = {'scores_aggregate.csv': 1e-3, 'scores_by_condition.csv': 1e-3, 'pooled_readings.csv': 1e-9,
             'reference_334.csv': 1e-9, 'additivity_table.csv': 1e-9, 'titration.csv': 1e-9,
             'si_statistics.json': 5e-3, 'simulation.csv': 0.05}
problems = []


def compare_frames(name, tol):
    a, b = pd.read_csv(os.path.join(OUT, name)), pd.read_csv(os.path.join(EXPECTED, name))
    if list(a.columns) != list(b.columns) or len(a) != len(b):
        problems.append(f'{name}: different shape or columns'); return
    for c in a.columns:
        if pd.api.types.is_numeric_dtype(b[c]):
            diff = np.nanmax(np.abs(a[c].to_numpy(float) - b[c].to_numpy(float))) if len(a) else 0.0
            if not diff <= tol:
                problems.append(f'{name}: column {c} differs by up to {diff:.3g} (tolerance {tol})')
        elif not (a[c].astype(str) == b[c].astype(str)).all():
            n = int((a[c].astype(str) != b[c].astype(str)).sum())
            problems.append(f'{name}: column {c} differs in {n} rows')


def compare_json(a, b, path, tol):
    if isinstance(b, dict):
        for k in b:
            if k not in a:
                problems.append(f'si_statistics.json: {path}{k} missing')
            else:
                compare_json(a[k], b[k], f'{path}{k} / ', tol)
    elif isinstance(b, list):
        if len(a) != len(b):
            problems.append(f'si_statistics.json: {path} has {len(a)} entries, expected {len(b)}')
        for x, y in zip(a, b):
            compare_json(x, y, path, tol)
    elif isinstance(b, (int, float)) and not isinstance(b, bool):
        if abs(a - b) > tol:
            problems.append(f'si_statistics.json: {path} = {a}, expected {b}')
    elif a != b:
        problems.append(f'si_statistics.json: {path} = {a!r}, expected {b!r}')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--skip-simulations', action='store_true')
    args = ap.parse_args()
    manifest = json.load(open(os.path.join(POOLED, 'manifest.json')))
    for f in manifest['files']:
        digest = hashlib.sha256(open(os.path.join(POOLED, f['file']), 'rb').read()).hexdigest()
        if digest != f['sha256']:
            problems.append(f"{f['file']}: SHA-256 does not match the manifest")
    for name, tol in TOLERANCE.items():
        if args.skip_simulations and name == 'simulation.csv':
            continue
        if name.endswith('.json'):
            compare_json(json.load(open(os.path.join(OUT, name))), json.load(open(os.path.join(EXPECTED, name))), '', tol)
        else:
            compare_frames(name, tol)
    if problems:
        print('DIFFERENCES FROM THE PAPER:\n  ' + '\n  '.join(problems))
        sys.exit(1)
    print('All outputs agree with the values reported in the paper.')


if __name__ == '__main__':
    main()
