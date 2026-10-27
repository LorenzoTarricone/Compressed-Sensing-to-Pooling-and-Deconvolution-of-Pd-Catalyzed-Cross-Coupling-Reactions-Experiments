"""Step 3: the statistics of SI sections S1.3, S2.2, S2.4, S2.5 and S2.6 and the validation subsets of S2.9.

Inputs (output folder): scores_aggregate.csv, pooled_readings.csv and run_manifest.json from scores.py,
reference_334.csv from reference_set.py; validation yields from Data/paper/validation_yields.csv.
Every reconstruction uses the solver of the Spotfire function (x >= 0, adaptive tolerance).

Writes si_statistics.json and additivity_table.csv to the output folder.

    python si_statistics.py
"""
import collections
import json
import os
import warnings

import numpy as np
import pandas as pd
from scipy import stats

from common import DATA, REFERENCE_NAME as NM, out, solve, spearman

warnings.filterwarnings('ignore')


def load_validation():
    """Externally calibrated LC-MS assay yields of the individually validated catalysts (%)."""
    v = pd.read_csv(os.path.join(DATA, 'validation_yields.csv'), dtype={'campaign': str})
    val = {}
    for c, g in v.groupby('campaign', sort=False):
        if set(g.condition) == {'A', 'B'}:
            p = g.pivot(index='catalyst', columns='condition', values='yield_percent')
            val[c] = {k: (int(p.loc[k, 'A']), int(p.loc[k, 'B'])) for k in g.catalyst.unique()}
        else:
            val[c] = dict(zip(g.catalyst, g.yield_percent.astype(int)))
    return val


def main():
    agg = pd.read_csv(out('scores_aggregate.csv'))
    readings = pd.read_csv(out('pooled_readings.csv'))
    manifest = json.load(open(out('run_manifest.json')))
    G = pd.read_csv(out('reference_334.csv'), index_col=0)
    VALIDATION = load_validation()
    res = {}

    # ---- plate 23: pooled readings, design matrix and aggregate scores
    r334 = readings[readings.campaign == 334]
    cats = list(agg[agg.campaign == 334].catalyst)
    members = r334[r334.condition_index == 1].sort_values('mixture').members.str.split('; ').tolist()
    A = np.array([[1.0 if c in m else 0.0 for c in cats] for m in members])
    assert set(A.sum(1)) == {7} and set(A.sum(0)) == {3}
    conds = [r334[r334.condition_index == i].condition.iloc[0] for i in range(1, 7)]
    Y = np.array([r334[r334.condition_index == i].sort_values('mixture').reading.to_numpy() for i in range(1, 7)])
    nz = (Y > 0).sum(1)
    y = Y[nz > 6].mean(0)
    score = dict(zip(agg[agg.campaign == 334].catalyst, agg[agg.campaign == 334].score))
    s_all = np.array([score[c] for c in cats])
    assert np.allclose(solve(y, A), s_all, atol=1e-4), 'solver does not reproduce the scores of scores.py'

    # ---- S1.3 condition filter
    res['S1.3 nonzero pools per condition'] = dict(zip(conds, nz.tolist()))
    sweep = []
    for thr in range(0, 15):
        keep = np.where(nz > thr)[0]
        if len(keep) == 0:
            break
        xs = solve(Y[keep].mean(0), A)
        sweep.append(dict(threshold=thr, n_kept=int(len(keep)), top3=[cats[i] for i in np.argsort(-xs)[:3]]))
    res['S1.3 filter sweep'] = sweep

    # ---- S2.4 scale of the score
    res['S2.4 scale'] = dict(mean_pooled_reading=round(float(y.mean()), 3), max_pooled_reading=round(float(y.max()), 3),
                             norm_y=round(float(np.linalg.norm(y)), 4), sum_y=round(float(y.sum()), 3),
                             seven_thirds_sum_y=round(float(7 / 3 * y.sum()), 3), sum_scores=round(float(s_all.sum()), 3),
                             max_score=round(float(s_all.max()), 3))

    # ---- S2.2 additivity: pooled readings against the single-catalyst reference
    complete = [i for i, m in enumerate(members) if all(n in NM for n in m)]
    rows = []
    for ci, cond in enumerate(conds):
        gcol = cond if cond.endswith('/') else cond + '/'
        for pi in complete:
            v = np.array([G.loc[NM[n], gcol] for n in members[pi]])
            rows.append(dict(cond=cond, pool=pi + 1, obs=Y[ci, pi], SUM=v.sum(), MEAN=v.mean(), MAX=v.max(),
                             SAT=100 * (1 - np.prod(1 - np.clip(v, 0, 99.9) / 100))))
    R = pd.DataFrame(rows)
    R.to_csv(out('additivity_table.csv'), index=False)
    add = {'n_comparisons': len(R), 'pools_with_all_members_measured': [i + 1 for i in complete]}
    for m in ['MEAN', 'MAX', 'SAT', 'SUM']:
        add[m] = dict(pearson=round(float(stats.pearsonr(R.obs, R[m])[0]), 3),
                      spearman=round(float(stats.spearmanr(R.obs, R[m]).correlation), 3),
                      scale_through_origin=round(float(np.sum(R.obs * R[m]) / np.sum(R[m] ** 2)), 3))
    for m in ['MEAN', 'MAX', 'SAT', 'SUM']:
        sl = add[m]['scale_through_origin']
        add[m]['rmse_scaled'] = round(float(np.sqrt(np.mean((R.obs - sl * R[m]) ** 2))), 2)
    sub = R[R.cond.str.contains('NaOtBu')]
    add['pearson_MEAN_NaOtBu_conditions_only'] = round(float(stats.pearsonr(sub.obs, sub.MEAN)[0]), 3)
    add['per_condition_spearman_MEAN'] = {c: round(float(stats.spearmanr(g.obs, g.MEAN).correlation), 3) for c, g in R.groupby('cond')}
    res['S2.2 additivity'] = add

    # ---- S2.5 benchmark of the 334 scores against the reference
    idx = [i for i, c in enumerate(cats) if c in NM]
    names = [cats[i] for i in idx]
    gt = np.array([G.loc[NM[c], [k if k.endswith('/') else k + '/' for k in conds]].max() for c in names])
    ref_cpme = np.array([G.loc[NM[c], 'CPME (10vol)/ NaOtBu (3.0 eq)/'] for c in names])
    s = s_all[idx]
    ordgt, ords = np.argsort(-gt), np.argsort(-s)
    hit = set(np.where(gt >= 25)[0])
    ties = collections.Counter(np.round(s, 4))
    bench = {'n': len(idx), 'spearman_best_of_6': spearman(s, gt), 'spearman_CPME_NaOtBu': spearman(s, ref_cpme),
             'spearman_mean_of_3_pools': spearman((A.T @ y)[idx] / 3, gt),
             'reference_top6': [(names[i], round(float(gt[i]), 1)) for i in ordgt[:6]],
             'score_top6': [(names[i], round(float(s[i]), 2), 'reference rank ' + str(list(ordgt).index(i) + 1)) for i in ords[:6]],
             'n_reference_ge_25': len(hit), 'reference_ge_25': [(names[i], round(float(gt[i]), 1)) for i in ordgt if gt[i] >= 25],
             'recovered_in_top_n_hit': len(hit & set(ords[:len(hit)])),
             'top_k_overlap': {k: len(set(ords[:k]) & set(ordgt[:k])) for k in (1, 2, 3, 5)},
             'tied_groups_of_reference_set': sorted([(float(v), n) for v, n in ties.items() if n > 1], reverse=True),
             'scores_of_reference_ge_25': [(names[i], round(float(s[i]), 2), list(ords).index(i) + 1) for i in ordgt if gt[i] >= 25]}
    def head_metrics(v):
        o = np.argsort(-np.round(v, 4), kind='stable')
        return dict(best_ranked_first=bool(o[0] == ordgt[0]), best_two_in_leading_two=len(set(o[:2]) & set(ordgt[:2])),
                    actives_in_leading_8=len(hit & set(o[:len(hit)])), actives_in_leading_3=len(hit & set(o[:3])))
    bench['head_score'] = head_metrics(s)
    bench['head_mean_of_3_pools'] = head_metrics((A.T @ y)[idx] / 3)
    full = collections.Counter(np.round(s_all, 4))
    bench['tied_groups_of_35'] = sorted([(float(v), n) for v, n in full.items() if n > 1], reverse=True)
    bench['n_nonzero_of_35'] = int((np.round(s_all, 4) > 0).sum())
    res['S2.5 benchmark'] = bench

    # ---- per-condition reconstructions against the reference (S1.3 / limitations)
    per = []
    for ci, cond in enumerate(conds):
        xs = solve(Y[ci], A)
        per.append(dict(condition=cond, top3=[cats[j] for j in np.argsort(-xs)[:3]],
                        spearman_vs_reference=spearman(xs[idx], gt)[0]))
    res['per-condition vs reference'] = per

    # ---- S2.6 rank stability: resample conditions, 30% multiplicative noise
    rng = np.random.default_rng(0)
    first, top3 = collections.Counter(), collections.Counter()
    B = 400
    for _ in range(B):
        yb = Y[rng.choice(6, 6, replace=True)].mean(0) * np.exp(rng.normal(0, 0.30, 15))
        o = np.argsort(-solve(yb, A))
        first[cats[o[0]]] += 1
        for j in o[:3]:
            top3[cats[j]] += 1
    res['S2.6 rank stability (400 resamples)'] = {c: dict(first=round(first[c] / B, 3), top3=round(top3[c] / B, 3))
                                                  for c, _ in top3.most_common(8)}

    # ---- validation subsets: Spearman of score against assay yield
    val = {}
    for c, d in VALIDATION.items():
        sc = dict(zip(agg[agg.campaign == int(c)].catalyst, agg[agg.campaign == int(c)].score))
        rk = dict(zip(agg[agg.campaign == int(c)].catalyst, agg[agg.campaign == int(c)]['rank']))
        lig = list(d)
        if c in ('424', '425'):
            best = [max(d[k]) for k in lig]
            val[c] = dict(n=len(lig), rho_best_of_AB=spearman([sc[k] for k in lig], best)[0],
                          rho_A=spearman([sc[k] for k in lig], [d[k][0] for k in lig])[0],
                          rho_B=spearman([sc[k] for k in lig], [d[k][1] for k in lig])[0],
                          best_validated=max(lig, key=lambda k: max(d[k])))
        else:
            val[c] = dict(n=len(lig), rho=spearman([sc[k] for k in lig], [d[k] for k in lig])[0],
                          best_validated=max(lig, key=lambda k: d[k]))
        val[c]['best_validated_rank'] = rk[val[c]['best_validated']]
        val[c]['scores_and_ranks'] = {k: (round(sc[k], 2), rk[k]) for k in lig}
    res['validation subsets'] = val
    res['accepted epsilon'] = {c: round(v['accepted_epsilon'], 4) for c, v in manifest['campaigns'].items()}

    json.dump(res, open(out('si_statistics.json'), 'w'), indent=1, default=str)
    print('S2.2 comparisons:', res['S2.2 additivity']['n_comparisons'], '| S2.5 benchmark n =', res['S2.5 benchmark']['n'],
          '| Spearman', res['S2.5 benchmark']['spearman_best_of_6'])


if __name__ == '__main__':
    main()
