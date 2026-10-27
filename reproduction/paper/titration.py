"""Step 5: XantPhos titration of SI S2.2 (ELN032036-334, plate 26).

Rows A to H each carry one active precatalyst at fixed loading; XantPhos Pd(allyl)Cl is added at 0, 4, 8, 12, 16
and 20 mol% in duplicate, in adjacent columns (1-2 = 0 mol%, 3-4 = 4 mol%, ..., 11-12 = 20 mol%). Product peaks
are summed per well at the last sampling point, duplicates are averaged, and each row is divided by its own wells
without added XantPhos.

Writes titration.csv (one line per row and loading) to the output folder.

    python titration.py
"""
import os

import numpy as np
import pandas as pd

from common import DATA, out

LOADINGS = [0, 4, 8, 12, 16, 20]


def main():
    x = pd.read_csv(os.path.join(DATA, 'titration_334_plate_26.csv'))
    last = x.SAMPLENUMBER.max()
    d = x[(x.SAMPLENUMBER == last) & (x.COMPOUND == 'Product')]
    M = np.zeros((8, 12))
    for (r, c), v in d.groupby(['PLATEROW', 'PLATECOLUMN'])['AREA_TOTAL_REDUCED'].sum().items():
        M['ABCDEFGH'.index(r), int(c) - 1] = v
    Y = M.reshape(8, 6, 2).mean(axis=2)
    rel = Y / Y[:, :1]
    rows = [dict(row='ABCDEFGH'[i], xantphos_mol_percent=v, product_area_percent=round(float(Y[i, j]), 4),
                 relative_to_zero=round(float(rel[i, j]), 4)) for i in range(8) for j, v in enumerate(LOADINGS)]
    pd.DataFrame(rows).to_csv(out('titration.csv'), index=False)
    supp = 100 * (1 - rel[:, -1])
    print('suppression at 20 mol%% by row: %s; range %.0f-%.0f%%, median %.0f%%; rows rising first: %s'
          % (np.round(supp).astype(int).tolist(), supp.min(), supp.max(), np.median(supp),
             ''.join('ABCDEFGH'[i] for i in range(8) if rel[i, 1] > 1)))


if __name__ == '__main__':
    main()
