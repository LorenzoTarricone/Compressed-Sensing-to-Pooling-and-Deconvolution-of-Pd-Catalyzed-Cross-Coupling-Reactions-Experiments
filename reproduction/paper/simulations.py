"""Step 4: recovery simulations of SI S2.6.

Synthetic plates on the Kirkman design: k of 35 catalysts active (35 to 90 area percentage units), the rest at
background (0.5 to 3); pooled readings from the mean model (gamma = 0), the maximum model (gamma = 1) or a mixture,
with multiplicative log-normal noise of the given coefficient of variation. Each plate is reconstructed with the
solver of the Spotfire function. Seed 7, 120 trials per setting.

Input: pooled_readings.csv and scores_aggregate.csv from scores.py (for the design matrix).
Writes simulation.csv to the output folder.

    python simulations.py
"""
import warnings

import cvxpy as cp
import numpy as np
import pandas as pd

warnings.filterwarnings('ignore')
from common import SCALES, out

members = pd.read_csv(out('pooled_readings.csv'))
members = members[(members.campaign == 424) & (members.condition_index == 1)].sort_values('mixture').members.str.split('; ')
cats = pd.read_csv(out('scores_aggregate.csv')).query('campaign == 424').catalyst.tolist()
A = np.array([[1.0 if c in m else 0.0 for c in cats] for m in members])


def recon(y):
    x = cp.Variable(35)
    for s in SCALES:
        prob = cp.Problem(cp.Minimize(cp.norm(x, 1)), [cp.norm2(y - A @ x) <= 0.01 * np.linalg.norm(y) * s, x >= 0])
        try:
            prob.solve()
        except cp.error.SolverError:
            continue
        if prob.status in ('optimal', 'optimal_inaccurate'):
            return x.value
    return None


def trial(k, cv, gamma, rng):
    x = rng.uniform(0.5, 3.0, 35)                  # inactive background
    act = rng.choice(35, k, replace=False)
    x[act] = rng.uniform(35, 90, k)                # active catalysts
    lin = A @ x / 7.0                              # mean model
    mx = np.array([x[A[i] > 0].max() for i in range(15)])
    y = (1 - gamma) * lin + gamma * mx             # gamma = 0 mean, gamma = 1 max
    if cv > 0:
        y = y * np.exp(rng.normal(0, cv, 15))
    xs = recon(np.maximum(y, 0))
    if xs is None:
        return None
    truth, o = set(np.argsort(-x)[:k]), np.argsort(-xs)
    return dict(top1=int(o[0] in truth), top3=len(set(o[:3]) & set(np.argsort(-x)[:3])), topk=len(set(o[:k]) & truth) / k)


rng = np.random.default_rng(7)
rows = []
for k in [1, 2, 3, 4, 5, 6, 8, 10]:
    for cv in [0.0, 0.15, 0.30, 0.50]:
        for gamma in [0.0, 0.5, 1.0]:
            r = [z for z in (trial(k, cv, gamma, rng) for _ in range(120)) if z]
            rows.append(dict(k=k, cv=cv, gamma=gamma, n=len(r), top1=np.mean([z['top1'] for z in r]),
                             top3_of_3=np.mean([z['top3'] for z in r]), recall_k=np.mean([z['topk'] for z in r])))
S = pd.DataFrame(rows)
S.to_csv(out('simulation.csv'), index=False)
for g in [0.0, 1.0]:
    print(f'gamma = {g}\nP(top-1 correct):\n', S[S.gamma == g].pivot(index='k', columns='cv', values='top1').round(2).to_string())
    print('recall@k:\n', S[S.gamma == g].pivot(index='k', columns='cv', values='recall_k').round(2).to_string())
