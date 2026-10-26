# Reproducing historical deconvolution scores

The deconvolution code changed several times between the first campaigns (2024) and the current Spotfire function. Scores reported for different campaigns were produced by different versions, so no single setting of the current function reproduces all of them. This folder re-implements each version so that a reported score can be traced to the code and settings that produced it.

Plate exports are not part of this repository. Pass the path to a Spotfire export on the command line.

## Versions

| Pipeline | Source | Tolerance (aggregate solve) | Sign constraint | Chemical prior |
|---|---|---|---|---|
| `notebook_2024` | `Python Notebooks/3_0_Algorithm.ipynb` (Sept 2024) | fixed, 0.1 | none | none |
| `spotfire_2024` | `.sfd` at `6fc2394` (Sept 2024) | fixed, 0.1 | none | optional: hand-made 2D coordinates, integer truncated, penalty λ‖Dx‖₂ |
| `spotfire_aug2025` | `.sfd` at `de9ec4f` (28 Aug 2025) | fixed, 10 (the `epsilon` input is only passed to the per-condition reconstruction) | x ≥ 0 | optional: α = 5, σ = √3 (bh, sm) or 300 (neutral), penalty λ‖Dx‖₁ |
| `spotfire_current` | `.sfd` from `b2e985f` | 0.01‖y‖₂ × 1, 2, 5, 10, 20, 50, 100, 200, first feasible | x ≥ 0 | as `spotfire_aug2025` |

All versions take the mean (or starting-material-weighted mean) over time points, apply the zero carry-forward correction, keep conditions with more than 6 non-zero pooled readings and average them. Before October 2026 every version kept only the last row of the export when a sample had several peaks tagged as Product. Every `.sfd` committed here before October 2026 assumes the vertical layout with the carbene in H7 to H12; the horizontal and carbene-first mappings first appear in `SpotfireProcessing_new.ipynb` (November 2025).

## What has been checked with the exports available locally

| Campaign | Reported scores | Result |
|---|---|---|
| 425 (barbituric acid) | all 12 values in main-text Figure 3 | reproduced by `notebook_2024` within the rounding of the figure (largest difference 0.04) |
| 424 (Meldrum's acid) | all 12 values in main-text Figure 3 | same ordering except cPhos and the phosphatrane (21.9 vs 21.6), but values differ by up to 9.5. No pipeline, prior or tolerance tested reproduces them. The 425 run used a pre-processed export (`AREA_TOTAL_REDUCED_corrected`); the corresponding 424 file is probably needed |
| 424 | SI section S2.9 (frozen reanalysis) | reproduced by `spotfire_current` with summed Product peaks (largest difference 0.005; accepted ε = 1.7813, as listed in S2.9) |
| 334, 810, 821 | main text and Figures 4, 6, 7 | not tested: the pooled-plate exports are not available locally |

## Running a search

```
python reproduce.py --export 821_plate1.csv --plate 1 \
    --cats ../Data/cats_names_35.csv --layout horizontal --carbene first \
    --embedding-sm M_E_embeddings_sm.csv --embedding-bh M_E_embeddings_bh.csv \
    --embedding-neutral M_E_embeddings_neutral.csv --alias rac-Antphos=dppdtbpf \
    --also-as-run-layout --epsilons "1,5,10" \
    --targets "dppp=99.2,BINAP=89.3,Jackie CPhos=75.0,CPhos=69.8,dtbpf=49.3"
```

For every pipeline the search tries summed and last-row Product peaks, starting-material weights on and off, the chemical prior off and on (each embedding supplied, λ in 0.5, 1, 5, 10), and any extra fixed tolerances given. It scores the aggregate signal and each single solvent/base condition, both as 7x and normalised to a plate maximum of 100 (the `df_complete_data` output). Results are sorted by the largest absolute difference from the targets.

Points specific to the 2025 campaigns:

- The local embedding files have no row for rac-Antphos, so any prior run with rac-Antphos in the catalyst list would have failed unless a row was added. `--alias rac-Antphos=dppdtbpf` gives it the dppdtbpf coordinates. Use the embedding tables exported from the Spotfire analysis if they differ.
- `--also-as-run-layout` adds the vertical/last mapping, the only one in the committed Spotfire functions. If a horizontal or carbene-first plate was analysed with such a function, that is the mapping it received.
- The submitted manuscript says the four Suzuki nominees were selected in dioxane with aqueous KOH, so the single-condition rows of the output are relevant for 821.

## Checks

`historical_pipelines.py` was checked against the saved outputs of `3_0_Algorithm.ipynb` (425) and against the S2.9 table (424) with Python 3.12, NumPy 2.5, pandas 3.1, SciPy and CVXPY 1.9.3 (CLARABEL).
