# Reproducing the paper

Every deconvolution score, statistic and SI figure in the revised paper *Accelerated Pd-catalysed cross-coupling optimisation with compressed sensing* can be regenerated from the data in [`Data/paper`](../Data/paper) with the scripts in [`reproduction/paper`](paper). The release used for the paper is tagged **`v1.0`**.

## Quick start

```bash
python3.12 -m venv .venv
source .venv/bin/activate
pip install -r reproduction/paper/requirements.txt
python reproduction/paper/run_all.py
```

The run takes about a minute and a half on a laptop, most of it in the recovery simulations; `--skip-simulations` leaves those out (about 20 seconds). It ends with

```
All outputs agree with the values reported in the paper.
```

or with a list of the outputs that differ and by how much. Outputs are written to `reproduction/paper/output` (or to the folder named in the `CATMIX_OUTPUT` environment variable); the values reported in the paper are in `reproduction/paper/expected`.

## The steps

| Script | What it does | Writes |
|---|---|---|
| `scores.py` | Runs the Spotfire data function of this repository (`Spotfire/Lorenzo_M_E_deconvolution.sfd`) on the five pooled plates. The script embedded in the `.sfd` file is executed as it is, so the scores are those Spotfire gives for the same inputs. | `scores_aggregate.csv`, `scores_by_condition.csv`, `pooled_readings.csv`, `plates_time_collapsed.csv`, `table_scores_rows.tex`, `run_manifest.json` |
| `reference_set.py` | Builds the single-catalyst reference set of the Buchwald--Hartwig amination (43 precatalysts x 8 solvent/base combinations). | `reference_334.csv` |
| `si_statistics.py` | Condition filter, additivity test, scale of the score, library-wide benchmark, rank stability and validation subsets. | `si_statistics.json`, `additivity_table.csv` |
| `simulations.py` | Recovery simulations on synthetic plates. | `simulation.csv` |
| `titration.py` | XantPhos titration. | `titration.csv` |
| `figures.py` | SI figures and main-text Figure 5. | `figures/*.png` |
| `check.py` | Checks the input hashes and compares every output with `expected/`. | |

## Where each result in the paper comes from

| In the paper | Output |
|---|---|
| Scores and ranks in main-text Figures 2C, 3, 4, 6 and 7, in the text and in SI S5; complete score table of SI S2.9 | `scores_aggregate.csv` (`table_scores_rows.tex` is the LaTeX body of the S2.9 table) |
| Pooled readings of Figure 2B and 2C | `plates_time_collapsed.csv`, `pooled_readings.csv` |
| Accepted tolerance, retained conditions, time weights and input hashes (SI S2.9) | `run_manifest.json` |
| Condition filter (SI S1.3) | `si_statistics.json`: `S1.3 nonzero pools per condition`, `S1.3 filter sweep` |
| Additivity test, its table and figure (SI S2.2; main text) | `si_statistics.json`: `S2.2 additivity`; `additivity_table.csv`; `figures/SI_signal_model.png` |
| XantPhos titration (SI S2.2) | `titration.csv`; `figures/SI_dilution.png` |
| Scale of the score and the conservation identity (SI S2.4) | `si_statistics.json`: `S2.4 scale` |
| Library-wide benchmark, its table (SI S2.5) and main-text Figure 5 | `si_statistics.json`: `S2.5 benchmark`; `figures/SI_benchmark.png` |
| Rank stability (SI S2.6) | `si_statistics.json`: `S2.6 rank stability (400 resamples)`; `figures/SI_rank_stability.png` |
| Recovery simulations (SI S2.6) | `simulation.csv`; `figures/SI_simulation.png` |
| Validation subsets (SI S2.9) | `si_statistics.json`: `validation subsets` |

## Settings

The same settings are used for every campaign:

- all peaks classified as Product in a sample and sampling point are summed;
- zero carry-forward time correction: a reading that drops to exactly zero after a non-zero reading is replaced by the last non-zero one;
- time points are combined with weights proportional to the number of wells in which the limiting starting material is detected (uniform when it is detected in none);
- a solvent/base combination is kept if more than 6 of its 15 pooled readings are positive (if none is kept, the two with the largest total are used; in the five campaigns all six were kept);
- the kept combinations are averaged with equal weight;
- reconstruction: minimise the L1 norm of x subject to ||y - Ax||_2 <= epsilon and x >= 0, with epsilon = 0.01 ||y||_2 x 1, 2, 5, 10, 20, 50, 100, 200, the first feasible solve being accepted;
- score = 7x; scores closer than 1e-4 share a rank;
- no chemical-similarity penalty.

The catalyst list differs at library position 28: (dppdtbpf)PdCl2 in campaigns 334, 424 and 425 (`Data/cats_names_35_dppdtbpf.csv`), rac-Antphos in 810 and 821 (`Data/cats_names_35.csv`). The plate layouts are set per campaign in `scores.py`.

## Notes

- **Benchmark membership.** The pooled catalysts are matched to the reference set in `common.py` (`REFERENCE_NAME`). P(oTol)3 was pooled but never measured on its own, so the benchmark covers 34 of the 35 pooled precatalysts. P(tBu)3 was pooled as its G2 precatalyst and measured on its own as the crotyl chloride, so it is matched by ligand.
- **Validation yields** (`Data/paper/validation_yields.csv`) are externally calibrated LC--MS assay yields of the individually validated catalysts. They enter only the Spearman correlations of the validation subsets.
- **Solver.** CVXPY chooses the solver (Clarabel with the versions above). Another solver or version can move scores in the fourth decimal and the simulated frequencies by a trial or two out of 120; `check.py` allows for this.
- **Seeds.** Rank stability uses seed 0 and 400 resamples; the simulations use seed 7 and 120 trials per setting.

## Earlier versions of the code

The scores in the first submission were produced by several earlier versions of the analysis. [`historical/`](historical) re-implements those versions and documents which of the earlier numbers each one reproduces.
