"""Search the historical pipelines for the settings that reproduce reported scores.

Example (campaign 821, late-stage Suzuki):

    python reproduce.py --export 821_plate1.csv --plate 1 \
        --cats ../Data/cats_names_35.csv --layout horizontal --carbene first \
        --embedding-sm M_E_embeddings_sm.csv --alias rac-Antphos=dppdtbpf \
        --targets "dppp=99.2,BINAP=89.3,Jackie CPhos=75.0,CPhos=69.8,dtbpf=49.3"

Each configuration is scored by the largest absolute difference from the
targets. The aggregate score (average over retained conditions) and each
single-condition reconstruction are both tested, the latter raw (7x) and
normalised so that the plate maximum is 100 (df_complete_data).
"""
import argparse
import itertools
import warnings

import numpy as np
import pandas as pd

from historical_pipelines import (PIPELINES, averaged_signal, mixtures_by_condition, plate_signal,
                                  prior_matrix, read_export, reconstruct)


def parse_targets(text):
    pairs = (item.split('=') for item in text.split(','))
    return {name.strip(): float(value) for name, value in pairs}


def load_embedding(path, alias):
    if path is None:
        return None
    emb = pd.read_csv(path, encoding='utf-16' if open(path, 'rb').read(2) in (b'\xff\xfe', b'\xfe\xff') else 'utf-8',
                      sep=None, engine='python')
    if 'catalyst_name' not in emb.columns:  # hand-made 2D coordinates: first column holds the names
        emb = emb.set_index(emb.columns[0])
        for new, old in alias.items():
            if new not in emb.index and old in emb.index:
                emb.loc[new] = emb.loc[old]
        return emb
    for new, old in alias.items():
        if new not in set(emb['catalyst_name']) and old in set(emb['catalyst_name']):
            row = emb[emb['catalyst_name'] == old].iloc[0].copy()
            row['catalyst_name'] = new
            emb = pd.concat([emb, row.to_frame().T], ignore_index=True)
    return emb


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--export', required=True)
    ap.add_argument('--plate', type=int)
    ap.add_argument('--cats', required=True, help='CSV with a CATALYST column, in pooling order')
    ap.add_argument('--layout', choices=['vertical', 'horizontal'], required=True)
    ap.add_argument('--carbene', choices=['first', 'last'], required=True)
    ap.add_argument('--targets', required=True, help='"name=score,..." as reported')
    ap.add_argument('--pipelines', default=','.join(PIPELINES))
    ap.add_argument('--embedding-sm')
    ap.add_argument('--embedding-bh')
    ap.add_argument('--embedding-neutral', help='2D coordinates used by spotfire_2024 and rxn_category=neutral')
    ap.add_argument('--lambdas', default='0.5,1,5,10')
    ap.add_argument('--sigmas-2024', default='50,100,300', help='sigma grid for the spotfire_2024 prior')
    ap.add_argument('--epsilons', default='',
                    help='extra tolerances to try: fixed values, or fractions of ||y|| ending in y, e.g. "1,10,0.5y"')
    ap.add_argument('--extra-layouts', default='', help='further layouts to try, e.g. "horizontal/first,vertical/last"')
    ap.add_argument('--alias', default='', help='use another catalyst\'s coordinates, e.g. rac-Antphos=dppdtbpf')
    ap.add_argument('--also-as-run-layout', action='store_true',
                    help='also try vertical/last, the only layout supported before November 2025')
    ap.add_argument('--top', type=int, default=15)
    args = ap.parse_args()
    warnings.filterwarnings('ignore')

    cats = list(pd.read_csv(args.cats)['CATALYST'])
    targets = parse_targets(args.targets)
    unknown = set(targets) - set(cats)
    if unknown:
        raise SystemExit(f"Target names not in --cats: {sorted(unknown)}")
    target_idx = np.array([cats.index(n) for n in targets])
    target_val = np.array(list(targets.values()))
    alias = dict(a.split('=') for a in args.alias.split(',') if a)
    embeddings = {k: load_embedding(p, alias) for k, p in
                  (('sm', args.embedding_sm), ('bh', args.embedding_bh), ('neutral', args.embedding_neutral)) if p}
    export = read_export(args.export)

    layouts = [(args.layout, args.carbene)]
    if args.also_as_run_layout and layouts[0] != ('vertical', 'last'):
        layouts.append(('vertical', 'last'))
    for item in filter(None, args.extra_layouts.split(',')):
        if tuple(item.split('/')) not in layouts:
            layouts.append(tuple(item.split('/')))
    lambdas = [float(v) for v in args.lambdas.split(',')]
    extra_eps = [v for v in args.epsilons.split(',') if v]

    def tolerance(eps, signal):
        if eps is None or not eps.endswith('y'):
            return None if eps is None else float(eps)
        return float(eps[:-1]) * np.linalg.norm(signal)

    results = []
    for pipeline in args.pipelines.split(','):
        priors = [None]
        if pipeline == 'spotfire_2024' and 'neutral' in embeddings:
            priors += [('neutral', s, lam) for s in map(float, args.sigmas_2024.split(',')) for lam in lambdas]
        elif pipeline in ('spotfire_aug2025', 'spotfire_current'):
            priors += [(k, None, lam) for k in embeddings for lam in lambdas]
        for agg, smw, (layout, carbene), prior, eps in itertools.product(
                ['last', 'sum'], [True, False], layouts, priors, [None] + extra_eps):
            plate = plate_signal(export, args.plate, agg, True, smw)
            readings = mixtures_by_condition(plate, layout, carbene)
            D, lam = None, 0.0
            if prior:
                kind, sigma, lam = prior
                D = prior_matrix(pipeline, cats, embeddings[kind], kind, sigma)
            y, keep = averaged_signal(readings)
            candidates = {}
            scores, used = reconstruct(y, pipeline, D, lam, tolerance(eps, y))
            if scores is not None:
                candidates['aggregate'] = (scores, used)
            per_cond = []
            for c in range(readings.shape[0]):
                s, e = reconstruct(readings[c], pipeline, D, lam, tolerance(eps, readings[c]))
                per_cond.append(s)
                if s is not None:
                    candidates[f'condition {c + 1}'] = (s, e)
            if all(s is not None for s in per_cond):
                plate_max = max(np.max(s) for s in per_cond) / 7  # df_complete_data uses x, not 7x
                for c, s in enumerate(per_cond):
                    candidates[f'condition {c + 1} normalised'] = (s / 7 / plate_max * 100, None)
            for where, (s, e) in candidates.items():
                err = float(np.max(np.abs(s[target_idx] - target_val)))
                results.append(dict(error=round(err, 3), pipeline=pipeline, signal=where,
                                    product=agg, sm_weights=smw, layout=f'{layout}/{carbene}',
                                    prior='off' if not prior else f'{prior[0]} lambda={prior[2]}' +
                                    (f' sigma={prior[1]}' if prior[1] else ''),
                                    epsilon=None if e is None else round(e, 4),
                                    scores=' '.join(f'{n}={s[cats.index(n)]:.1f}' for n in targets)))
    table = pd.DataFrame(results).sort_values('error')
    with pd.option_context('display.max_colwidth', 120, 'display.width', 250):
        print(table.head(args.top).to_string(index=False))


if __name__ == '__main__':
    main()
