"""Reproduce every number and SI figure of the paper from the data in Data/paper.

    python run_all.py                     # all steps, then compare with expected/
    python run_all.py --skip-simulations  # leave out step 4 (the slowest)

Outputs go to reproduction/paper/output (or the folder in the CATMIX_OUTPUT environment variable).
"""
import argparse
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
STEPS = [('scores.py', 'deconvolution scores, five campaigns'),
         ('reference_set.py', 'single-catalyst reference set, campaign 334'),
         ('si_statistics.py', 'statistics of SI S1.3, S2.2, S2.4, S2.5, S2.6 and S2.9'),
         ('simulations.py', 'recovery simulations, SI S2.6'),
         ('titration.py', 'XantPhos titration, SI S2.2'),
         ('figures.py', 'SI figures and main-text Figure 5')]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--skip-simulations', action='store_true')
    args = ap.parse_args()
    for script, what in STEPS:
        if args.skip_simulations and script == 'simulations.py':
            continue
        t0 = time.time()
        print(f'== {script}: {what}', flush=True)
        subprocess.run([sys.executable, os.path.join(HERE, script)], cwd=HERE, check=True)
        print(f'   done in {time.time() - t0:.0f} s', flush=True)
    check = [sys.executable, os.path.join(HERE, 'check.py')] + (['--skip-simulations'] if args.skip_simulations else [])
    sys.exit(subprocess.run(check, cwd=HERE).returncode)


if __name__ == '__main__':
    main()
