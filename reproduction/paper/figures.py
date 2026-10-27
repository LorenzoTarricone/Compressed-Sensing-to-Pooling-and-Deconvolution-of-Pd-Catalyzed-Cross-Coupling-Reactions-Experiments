"""Step 6: the figures of the SI (and main-text Figure 5) from the outputs of steps 1 to 5.

Writes to output/figures: SI_signal_model.png, SI_benchmark.png (main-text Figure 5), SI_rank_stability.png,
SI_simulation.png and SI_dilution.png.

    python figures.py
"""
import json
import os
import warnings

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from common import FIGURES, REFERENCE_NAME as NM, out  # noqa: E402

warnings.filterwarnings('ignore')
OUT = FIGURES
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'figure.dpi': 200})
C = {'obs': '#1f4e79', 'fit': '#c0504d', 'g': '#7f7f7f'}
SI = json.load(open(out('si_statistics.json')))

# ---------- signal model (SI S2.2) ----------
R = pd.read_csv(out('additivity_table.csv'))
fig, axes = plt.subplots(1, 3, figsize=(9.2, 3.2))
titles = ['Mean of members\n$(1/7)\\sum_j s_j$', 'Max of members\n$\\max_j s_j$', 'Saturating\n$100[1-\\prod_j(1-s_j/100)]$']
for ax, m, t in zip(axes, ['MEAN', 'MAX', 'SAT'], titles):
    ax.scatter(R[m], R.obs, s=16, c=C['obs'], alpha=.75, edgecolors='none')
    lim = max(R[m].max(), R.obs.max()) * 1.05
    sl = np.sum(R.obs * R[m]) / np.sum(R[m] ** 2)
    xs = np.linspace(0, lim, 50)
    ax.plot(xs, xs, '--', c=C['g'], lw=.9, label='1:1')
    ax.plot(xs, sl * xs, '-', c=C['fit'], lw=1.4, label=f'slope {sl:.2f}')
    r, rho = stats.pearsonr(R.obs, R[m])[0], stats.spearmanr(R.obs, R[m]).correlation
    ax.set_title(t, fontsize=8.5); ax.set_xlim(0, lim); ax.set_ylim(0, lim)
    ax.set_xlabel('predicted pool reading (area %)')
    ax.text(.04, .93, f'$r$ = {r:.2f}\n$\\rho$ = {rho:.2f}', transform=ax.transAxes, va='top', fontsize=8)
    ax.legend(frameon=False, fontsize=7.5, loc='lower right')
axes[0].set_ylabel('observed pool reading (area %)')
fig.tight_layout(); fig.savefig(os.path.join(OUT, 'SI_signal_model.png'), bbox_inches='tight'); plt.close(fig)

# ---------- benchmark (main-text Figure 5, SI S2.5) ----------
agg = pd.read_csv(out('scores_aggregate.csv')).query('campaign == 334').set_index('catalyst')
G = pd.read_csv(out('reference_334.csv'), index_col=0)
conds = [c for c in G.columns if not c.startswith(('MeCN', 'CPME (10vol)/ KOH'))]
names = [c for c in agg.index if c in NM]
score = np.round(agg.loc[names, 'score'].to_numpy(float), 4)
ref = np.array([G.loc[NM[n], conds].max() for n in names])
tied = score <= 1e-4
fig, ax = plt.subplots(figsize=(6.0, 4.2))
ax.scatter(score[~tied], ref[~tied], s=40, c=C['obs'], zorder=4, label='non-zero score')
ax.scatter(score[tied], ref[tied], s=34, facecolors='none', edgecolors=C['fit'], zorder=4, label=f'score 0 ({int(tied.sum())} tied)')
ax.axhline(25, ls=':', c=C['g'], lw=.9)
ax.text(score.max() * 1.02, 25.7, 'active threshold', fontsize=7, c=C['g'], ha='right')
offsets = {'CPhos': (6, 3), 'QPhos': (6, -8), 'cBridP': (6, 2), 'Jackie CPhos': (-6, 6), 'P(tBu)_3': (6, 2),
           'RuPhos': (6, 0), 'SPhos': (6, 2), 'GPhos': (6, 2), 'Triisobutyl-Phosphatrane': (6, 2)}
for n, (dx, dy) in offsets.items():
    i = names.index(n)
    ax.annotate(n.replace('_', '').replace('CPhos', 'cPhos') if n == 'CPhos' else n.replace('_', ''), (score[i], ref[i]),
                fontsize=7.4, xytext=(dx, dy), textcoords='offset points', ha='left' if dx > 0 else 'right', zorder=5)
ax.set_xlabel('deconvolution score from the pooled plate')
ax.set_ylabel('single-catalyst activity, best of the 6 screened\nconditions (reduced area %)')
rho, p = SI['S2.5 benchmark']['spearman_best_of_6']
ax.set_title(f'{len(names)} pooled precatalysts with individual reference measurements\nSpearman $\\rho$ = {rho:.2f} ($P$ = {p:.2f})', fontsize=9)
ax.legend(frameon=False, fontsize=7.8, loc='center right')
fig.tight_layout(); fig.savefig(os.path.join(OUT, 'SI_benchmark.png'), bbox_inches='tight'); plt.close(fig)

# ---------- rank stability (SI S2.6) ----------
st = SI['S2.6 rank stability (400 resamples)']
names_st = [n for n in st if st[n]['top3'] >= 0.05]
fig, ax = plt.subplots(figsize=(5.4, 3.0))
ypos = np.arange(len(names_st))[::-1]
ax.barh(ypos, [st[n]['top3'] for n in names_st], color=C['obs'], height=.62, label='in top three')
ax.barh(ypos, [st[n]['first'] for n in names_st], color=C['fit'], height=.30, label='ranked first')
ax.set_yticks(ypos); ax.set_yticklabels([n.replace('_', '') for n in names_st], fontsize=8); ax.set_xlim(0, 1)
ax.set_xlabel('frequency over 400 resamples'); ax.legend(frameon=False, fontsize=7.5, loc='lower right')
fig.tight_layout(); fig.savefig(os.path.join(OUT, 'SI_rank_stability.png'), bbox_inches='tight'); plt.close(fig)

# ---------- recovery simulations (SI S2.6) ----------
S = pd.read_csv(out('simulation.csv'))
fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.3))
cols = {0.0: '#1f4e79', 0.15: '#4f81bd', 0.30: '#e0803c', 0.50: '#c0504d'}
for ax, metric, ttl in zip(axes, ['top1', 'recall_k'], ['probability that the best catalyst is ranked first',
                                                       'fraction of the $k$ active catalysts in the top $k$']):
    for cv, c in cols.items():
        s = S[(S.gamma == 0.0) & (S.cv == cv)]
        ax.plot(s.k, s[metric], '-o', ms=3.5, c=c, label=f'CV = {int(cv * 100)} %')
    s = S[(S.gamma == 1.0) & (S.cv == 0.15)]
    ax.plot(s.k, s[metric], '--s', ms=3.5, c='#6a3d9a', label='CV = 15 %, max model')
    ax.set_xlabel('number of active catalysts $k$'); ax.set_ylim(0, 1.03); ax.set_title(ttl, fontsize=8.5)
axes[0].set_ylabel('probability'); axes[1].set_ylabel('recall')
axes[0].legend(frameon=False, fontsize=7.5, loc='lower left')
fig.tight_layout(); fig.savefig(os.path.join(OUT, 'SI_simulation.png'), bbox_inches='tight'); plt.close(fig)

# ---------- XantPhos titration (SI S2.2) ----------
T = pd.read_csv(out('titration.csv'))
rows = list(dict.fromkeys(T['row']))
cols = plt.cm.viridis(np.linspace(0, 0.9, len(rows)))
fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.3))
for c, r in zip(cols, rows):
    s = T[T['row'] == r]
    axes[0].plot(s.xantphos_mol_percent, s.product_area_percent, '-o', ms=3.5, c=c, label=f'row {r}')
    axes[1].plot(s.xantphos_mol_percent, s.relative_to_zero, '-o', ms=3.5, c=c)
v = np.linspace(0, 20, 50)
axes[1].plot(v, 8 / (8 + v), '--', c=C['fit'], lw=2, label='8/(8+$v$)')
axes[0].set_ylabel('product, reduced area %')
axes[1].set_ylabel('relative to the same row at 0 mol%')
axes[1].set_ylim(0, max(1.6, T.relative_to_zero.max() * 1.05))
for ax in axes:
    ax.set_xlabel('XantPhos Pd(allyl)Cl added (mol %)'); ax.set_xticks(sorted(T.xantphos_mol_percent.unique()))
axes[0].legend(frameon=False, fontsize=7.5, ncol=2)
axes[1].legend(frameon=False, fontsize=7.5, loc='upper right')
fig.tight_layout(); fig.savefig(os.path.join(OUT, 'SI_dilution.png'), bbox_inches='tight'); plt.close(fig)
print('figures written to', OUT)
