# Tracing the scores of the first submission to the code that produced them

The deconvolution code changed several times between the first campaigns (2024) and the current Spotfire function. Scores reported for different campaigns were produced by different versions, so no single setting of the current function reproduces all of them. This folder re-implements each version so that a score reported in the first submission can be traced to the code and settings that produced it. The revised paper does not use these versions: every number in it comes from the current function with one set of settings, reproduced by `reproduction/paper` (see `reproduction/README.md`).

The pooled-plate exports are in `Data/paper/pooled_plates`; pass a file from there on the command line.

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
| All five | SI section S2.9 of the revised paper | reproduced by `Spotfire/Lorenzo_M_E_deconvolution.sfd` with `product_aggregation = "sum"`, prior off, campaign layouts and catalyst lists: all 35 scores per campaign within 0.005, accepted ε identical to S2.9 |
| 821 | an independent Spotfire run of 20 October 2026 | reproduced exactly by the same function with `product_aggregation = "last"` |
| 334 | Figure 4, all 35 values | from the processed 6 × 15 table archived at the time: `spotfire_aug2025` with the Suzuki epoch-10 embedding, λ = 10, ε = 0.5 ‖y‖ gives all 35 within 0.17. From today's raw export no setting comes within 6: the export differs from the archived table by up to 25 to 36 area% per cell, and 132 of 190 Product samples carry more than one Product peak |
| 821 | Figure 6, five values | closest: single-condition reconstruction in dioxane / aq. KOH (condition 5), last-row Product peaks, adaptive ε; largest difference 4.1. The prior has almost no effect on this condition |
| 810 | Figure 7, three values | closest: aggregate, no prior, summed peaks, fixed ε = 8; largest difference 0.8. Three values do not identify the settings |

## Running a search

```
python reproduce.py --export 821_plate1.csv --plate 1 \
    --cats ../Data/cats_names_35.csv --layout horizontal --carbene first \
    --embedding-sm M_E_embeddings_sm.csv --embedding-bh M_E_embeddings_bh.csv \
    --embedding-neutral M_E_embeddings_neutral.csv --alias rac-Antphos=dppdtbpf \
    --also-as-run-layout --epsilons "1,5,10,0.5y" \
    --targets "dppp=99.2,BINAP=89.3,Jackie CPhos=75.0,CPhos=69.8,dtbpf=49.3"
```

The raw files used for the table above are the pooled-plate exports in `Data/paper/pooled_plates` (with SHA-256 manifest). Pass only the columns the Spotfire function maps (`ELN_ID`, `SAMPLE_ID`, `PLATENUMBER`, `COMPOUND`, `SAMPLENUMBER`, `AREA_TOTAL_REDUCED`, `ColumnLabel`, `PLATEROW`, `PLATECOLUMN`) when running the `.sfd` script outside Spotfire: it calls `dropna()` on every column, and the full export has gaps in every row.

For every pipeline the search tries summed and last-row Product peaks, starting-material weights on and off, the chemical prior off and on (each embedding supplied, λ in 0.5, 1, 5, 10), and any extra fixed tolerances given. It scores the aggregate signal and each single solvent/base condition, both as 7x and normalised to a plate maximum of 100 (the `df_complete_data` output). Results are sorted by the largest absolute difference from the targets.

Points specific to the 2025 campaigns:

- The local embedding files have no row for rac-Antphos, so any prior run with rac-Antphos in the catalyst list would have failed unless a row was added. `--alias rac-Antphos=dppdtbpf` gives it the dppdtbpf coordinates. Use the embedding tables exported from the Spotfire analysis if they differ.
- `--also-as-run-layout` adds the vertical/last mapping, the only one in the committed Spotfire functions. If a horizontal or carbene-first plate was analysed with such a function, that is the mapping it received.
- The submitted manuscript says the four Suzuki nominees were selected in dioxane with aqueous KOH, so the single-condition rows of the output are relevant for 821.

## Checks

`historical_pipelines.py` was checked against the saved outputs of `3_0_Algorithm.ipynb` (425) and against the S2.9 table (424) with Python 3.12, NumPy 2.5, pandas 3.1, SciPy and CVXPY 1.9.3 (CLARABEL).
