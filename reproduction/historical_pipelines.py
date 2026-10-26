"""Re-implementations of every historical version of the CatMix deconvolution.

Each pipeline reproduces the numerical behaviour of one dated version of the
code, so that a reported score can be traced back to the version and settings
that produced it. Plate exports are not included in this repository; pass the
path to a Spotfire export on the command line (see README.md).

Versions
--------
notebook_2024     Python Notebooks/3_0_Algorithm.ipynb (Sept 2024). Fixed
                  epsilon = 0.1, no sign constraint, no prior. Produced the
                  published scores for campaigns 424 and 425.
spotfire_2024     Spotfire/Lorenzo_M_E_deconvolution.sfd at 6fc2394 (Sept 2024).
                  Fixed epsilon = 0.1, no sign constraint. Optional prior built
                  from the hand-made 2D coordinates (embedding_df.csv), integer
                  truncated, penalising lambda * ||Dx||_2.
spotfire_aug2025  .sfd at de9ec4f (28 Aug 2025). Sign constraint x >= 0. Prior
                  penalises lambda * ||Dx||_1 with D from alpha = 5,
                  sigma = sqrt(3) (bh/sm) or 300 (neutral). The aggregate solve
                  always uses the function default epsilon = 10; the epsilon
                  input only reaches the per-condition reconstruction.
spotfire_current  .sfd from b2e985f onwards. Adaptive epsilon: 0.01 * ||y||_2
                  times 1, 2, 5, 10, 20, 50, 100, 200 until feasible. x >= 0.
                  Same prior as spotfire_aug2025.
"""
import itertools

import cvxpy as cp
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

M_15_35 = np.array([[1,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,1,0],
                    [0,1,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,1,0,0,0,0,1,0,0],
                    [0,0,1,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,0,1,0],
                    [0,0,0,1,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0],
                    [0,0,0,0,1,1,0,0,0,0,0,1,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,0],
                    [1,0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,1],
                    [0,1,0,0,0,0,1,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,1],
                    [0,0,1,0,0,0,0,1,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1,0,0],
                    [0,0,0,1,0,0,0,1,0,0,0,0,1,0,0,1,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,1,0],
                    [0,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,1,0,0,0,0,0,1],
                    [1,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,1,0,0,0,1,0,0,0,0,0,0,1,0,0,1,0,0,0,0],
                    [0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,0,0,0,1,0,0,1,0,0,0],
                    [0,0,1,0,0,0,0,0,0,1,0,0,1,0,0,0,0,1,0,0,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0],
                    [0,0,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,1,0,1,0,0,0,0,1,0,0,0],
                    [0,0,0,0,1,0,0,0,0,1,0,0,0,1,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,1,0,0,1,0,0]])

POOL_SIZE = 7
PIPELINES = ('notebook_2024', 'spotfire_2024', 'spotfire_aug2025', 'spotfire_current')
ROW_INDEX = {c: i for i, c in enumerate('ABCDEFGH')}

# canonical mixture number held by each of the 15 physical rows (see plate maps in the SI)
PLATE_ORDER = {
    ('vertical', 'last'): list(range(1, 16)),
    ('horizontal', 'last'): [1, 3, 5, 7, 9, 11, 13, 15, 2, 4, 6, 8, 10, 12, 14],
    ('vertical', 'first'): list(range(1, 16)),
    ('horizontal', 'first'): [2, 4, 6, 8, 10, 12, 14, 1, 3, 5, 7, 9, 11, 13, 15],
}


def read_export(path):
    """Read a Spotfire CSV export, which may be UTF-16 or UTF-8."""
    for encoding in ('utf-16', 'utf-8-sig'):
        try:
            df = pd.read_csv(path, encoding=encoding)
        except (UnicodeError, pd.errors.ParserError):
            continue
        if df.shape[1] > 1:
            df.columns = df.columns.str.strip()
            return df.rename(columns={'PlateRow': 'PLATEROW', 'PlateColumn': 'PLATECOLUMN'})
    raise ValueError(f"Could not read {path}")


def plate_signal(df, plate_number=None, product_aggregation='last', time_correction=True,
                 sm_weights=True, sm_labels=('lim SM',)):
    """Collapse an export to the 8 x 12 plate of product area percentages."""
    if plate_number is not None and 'PLATENUMBER' in df.columns:
        df = df[df['PLATENUMBER'] == plate_number]
    value = 'AREA_TOTAL_REDUCED' if 'AREA_TOTAL_REDUCED' in df.columns else 'AREA_TOTAL_REDUCED_corrected'
    df = df.dropna(subset=['SAMPLENUMBER', value])
    times = df['SAMPLENUMBER'].astype(int)
    max_time = int(times.max())
    # 2024 exports for plate 425-11 were already filtered to product rows
    product = df[df['COMPOUND'] == 'Product'] if 'COMPOUND' in df.columns else df

    plate = np.zeros((8, 12, max_time))
    keys = [product['PLATEROW'], product['PLATECOLUMN'], product['SAMPLENUMBER'].astype(int)]
    for (row, col, t), group in product.groupby(keys, sort=False):
        areas = group[value].to_numpy(dtype=float)
        plate[ROW_INDEX[row], int(col) - 1, t - 1] = areas.sum() if product_aggregation == 'sum' else areas[-1]

    if time_correction:
        for t in range(1, max_time):
            dropped = (plate[:, :, t] == 0) & (plate[:, :, t - 1] != 0)
            plate[:, :, t][dropped] = plate[:, :, t - 1][dropped]

    weights = np.ones(max_time)
    if sm_weights and 'COMPOUND' in df.columns:
        sm = times[df['COMPOUND'].isin(sm_labels)]
        counts = np.array([(sm == t).sum() for t in range(1, max_time + 1)], dtype=float)
        if counts.sum() > 0:
            weights = counts
    return np.average(plate, axis=2, weights=weights)


def mixtures_by_condition(plate, layout='vertical', carbene='last'):
    """Return the 6 x 15 matrix of pooled readings in canonical mixture order."""
    left, right = plate[:, :6], plate[:, 6:]
    rows = np.concatenate((left, right[:7]) if carbene == 'last' else (left[1:], right), axis=0)
    order = np.argsort(np.array(PLATE_ORDER[(layout, carbene)]) - 1)
    return rows[order].T


def averaged_signal(readings, toll_exp=6):
    keep = [i for i in range(readings.shape[0]) if np.sum(readings[i] > 0) > toll_exp]
    if not keep:
        keep = list(np.argsort(readings.sum(axis=1))[-2:][::-1])
    return readings[keep].mean(axis=0), keep


def prior_matrix(pipeline, cats_names, embedding, rxn_category='sm', sigma=None, alpha=5.0):
    """Build the pairwise similarity matrix D exactly as the given version did."""
    if pipeline == 'spotfire_2024':
        coords = embedding.loc[cats_names, ['x', 'y']].apply(pd.to_numeric).astype('int64')
        sim = np.exp(-cdist(coords, coords) ** 2 / (2 * sigma ** 2)) * 100
        pairs = list(itertools.combinations(range(len(cats_names)), 2))
        D = np.zeros((len(pairs), len(cats_names)), dtype=int)
        for r, (i, k) in enumerate(pairs):
            D[r, i] = -1 * sim[i, k] / 2
            D[r, k] = sim[i, k] / 2
        return D
    if sigma is None:
        sigma = 300.0 if rxn_category == 'neutral' else 1.73205080
    if 'catalyst_name' in embedding.columns:
        coords = embedding.set_index('catalyst_name').loc[cats_names, ['dim_0', 'dim_1']].to_numpy(float)
    else:  # hand-made 2D coordinates indexed by name
        coords = embedding.loc[cats_names].iloc[:, :2].to_numpy(float)
    sim = alpha * np.exp(-cdist(coords, coords) ** 2 / (2 * sigma ** 2))
    pairs = list(itertools.combinations(range(len(cats_names)), 2))
    D = np.zeros((len(pairs), len(cats_names)))
    for r, (i, k) in enumerate(pairs):
        D[r, i], D[r, k] = -sim[i, k], sim[i, k]
    return D


def reconstruct(y, pipeline, D=None, lambdaa=0.0, epsilon=None, solver=None):
    """Solve one reconstruction; return (scores = 7x, epsilon used) or (None, None)."""
    x = cp.Variable(M_15_35.shape[1])
    objective = cp.norm(x, 1)
    if D is not None and lambdaa:
        objective = objective + lambdaa * (cp.norm(D @ x) if pipeline == 'spotfire_2024' else cp.norm(D @ x, 1))
    nonneg = pipeline in ('spotfire_aug2025', 'spotfire_current')

    if epsilon is not None:
        schedule = [epsilon]
    elif pipeline == 'spotfire_aug2025':
        schedule = [10.0]
    elif pipeline == 'spotfire_current':
        schedule = [0.01 * np.linalg.norm(y) * s for s in (1, 2, 5, 10, 20, 50, 100, 200)]
    else:
        schedule = [0.1]

    for eps in schedule:
        constraints = [cp.norm2(y - M_15_35 @ x) <= eps] + ([x >= 0] if nonneg else [])
        problem = cp.Problem(cp.Minimize(objective), constraints)
        problem.solve(solver=solver) if solver else problem.solve()
        if problem.status in ('optimal', 'optimal_inaccurate'):
            return x.value * POOL_SIZE, eps
    return None, None
