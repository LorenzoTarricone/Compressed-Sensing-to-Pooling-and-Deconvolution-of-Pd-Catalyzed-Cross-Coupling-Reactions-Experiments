"""Paths, catalyst-name mapping and the solver shared by the scripts that reproduce the paper."""
import os

import cvxpy as cp
import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
DATA = os.path.join(REPO, 'Data', 'paper')
POOLED = os.path.join(DATA, 'pooled_plates')
SFD = os.path.join(REPO, 'Spotfire', 'Lorenzo_M_E_deconvolution.sfd')
CATALYST_LISTS = os.path.join(REPO, 'Data')
OUT = os.environ.get('CATMIX_OUTPUT', os.path.join(HERE, 'output'))
FIGURES = os.path.join(OUT, 'figures')
EXPECTED = os.path.join(HERE, 'expected')

# Tolerance schedule of the Spotfire function: 0.01 ||y||_2 times these factors, first feasible solve accepted.
SCALES = (1, 2, 5, 10, 20, 50, 100, 200)
TIE_TOLERANCE = 1e-4   # scores closer than this share a rank

# Catalyst name in the pooled-plate catalyst lists -> name in the single-catalyst reference plates of campaign 334
# (ELN032036-334, plates 15 to 18). P(oTol)_3 was pooled but never measured on its own, so it has no entry.
REFERENCE_NAME = {
    'meCgPPh': 'meCgPPh Pd G3', 'TrixiePhos': 'TrixiePhos Pd G3', 'RuPhos': 'RuPhos Pd(crotyl)Cl',
    'RockPhos': 'RockPhos Pd G6 Br', 'VPhos': 'VPhos Pd G4', 'PPh_3': 'PPh3)2PdCl2',
    'tBuBrettPhos': 'tBuBrettPhos Pd(allyl)OTf', 'CyJohnPhos': 'CyJohnPhos Pd G3', 'CPhos': 'cPhos Pd G3',
    'DPEPhos': 'DPEPhos Pd(allyl)Cl', 'AlPhos': 'AlPhosPd G6 Br', 'Triisobutyl-Phosphatrane': 'triisobutylphosphatrane Pd G6',
    'SPhos': 'SPhos Pd(tBuInd)Cl', 'XantPhos': 'XantPhos Pd(allyl)Cl', 'Me3(OMe)-tBuXPhos': 'Me3(OMe)tBuXPhos Pd G3',
    'Amphos': 'AmPhos Pd(crotyl)Cl', 'Jackie CPhos': 'Jackie CPhos Pd G4', 'dppp': 'dpppPdCl2',
    'Mordal-Phos': 'MorDalPhos Pd G6 Br', 'Jackie Phos': 'JackiePhos Pd G3', 'dppf': 'dppfPdCl2',
    # pooled as the G2 precatalyst, measured on its own as the crotyl chloride: matched by ligand
    'P(tBu)_3': 'P(tBu)3 Pd(crotyl)Cl',
    'GPhos': 'GPhos Pd G6 Br', 'BINAP': 'rac-BINAP Pd G3', 'BippyPhos': 'BippyPhos Pd(allyl)OTf',
    'cataCXlum A': 'CataCXiumA Pd(allyl)Cl', 'QPhos': 'QPhos Pd(crotyl)Cl', 'cBridP': 'cBridP Pd(allyl)Cl',
    'dtbpf': 'dtbpfPdCl2', 'JosiPhosSL J009-1': 'JosiphosSL-J009-1 Pd G3', 'PCy_3': 'P(Cy)3 Pd(crotyl)Cl',
    'XPhos': 'XPhos Pd(crotyl)Cl', 'tBuXPhos': 'tBuXPhos Pd(allyl)OTf',
    # library position 28 of campaign 334 (the reference export writes the name without its opening bracket)
    'dppdtbpf': 'dppdtbpf)PdCl2',
}


def out(*parts):
    """Path inside the output folder, creating the folder if needed."""
    path = os.path.join(OUT, *parts)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    return path


def solve(y, A):
    """L1 reconstruction exactly as in the Spotfire function: x >= 0, first feasible tolerance. Returns 7x."""
    x = cp.Variable(A.shape[1])
    for s in SCALES:
        prob = cp.Problem(cp.Minimize(cp.norm(x, 1)), [cp.norm2(y - A @ x) <= 0.01 * np.linalg.norm(y) * s, x >= 0])
        try:
            prob.solve()
        except cp.error.SolverError:
            continue
        if prob.status in ('optimal', 'optimal_inaccurate'):
            return x.value * 7
    raise RuntimeError('no feasible tolerance in the schedule')


def spearman(a, b):
    """Spearman correlation; scores within TIE_TOLERANCE count as tied, as in the published rank table."""
    r = stats.spearmanr(np.round(np.asarray(a, float), 4), b)
    return round(float(r.correlation), 3), float(f'{r.pvalue:.2g}')
