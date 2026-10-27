"""Step 2: single-catalyst reference activities for the Buchwald--Hartwig amination (campaign 334).

43 precatalysts were run one per well at 8 mol% across eight solvent/base combinations on plates 15 to 18 of
ELN032036-334, before the pooled plate. For each plate the last sampling point is used (plates 15 and 18 have
one, plates 16 and 17 three); duplicates of a catalyst/condition pair are combined by their maximum. The result
is a 43 x 8 table of reduced product area percentages, used by the additivity test (SI S2.2) and the benchmark
(SI S2.5).

Writes reference_334.csv to the output folder.

    python reference_set.py
"""
import os

import pandas as pd

from common import DATA, out

CONDITIONS = ['MeCN (10vol)/ K3PO4 (3.0 eq, 3.25M, aq)/', 'CPME (10vol)/ KOH (3.0 eq, 3M, aq)/',
              'tAmOH (10vol)/ KOH (3.0 eq, 3M, aq)/', 'PhMe (10vol)/ NaOtBu (3.0 eq)/', 'CPME (10vol)/ NaOtBu (3.0 eq)/',
              'tAmOH (10vol)/ Cs2CO3 (3.0 eq)/', 'PhMe (10vol)/ TMSOK (3.0 eq)/', 'CPME (10vol)/ KOtBu (3.0 eq)/']


def main():
    df = pd.read_csv(os.path.join(DATA, 'reference_334_single_catalyst_plates_15-18.csv'))
    df['catalyst'] = df['ColumnLabel'].str.extract(r'^\d+:\s*(.*?)\s*\(8mol%\)')[0]
    df['condition'] = df['RowLabel'].str.extract(r'^[A-H]:\s*(.*?)\s*$')[0]
    df = df[df.catalyst.notna() & df.SAMPLENUMBER.notna()]
    last = df.groupby('PLATENUMBER').SAMPLENUMBER.transform('max')
    df = df[df.SAMPLENUMBER == last]
    table = df.groupby(['catalyst', 'condition'])['AREA_TOTAL_REDUCED'].max().unstack()
    table = table.reindex(columns=[c for c in CONDITIONS if c in table.columns]).fillna(0.0)
    table.index.name = 'cat'
    assert table.shape == (43, 8), table.shape
    table.to_csv(out('reference_334.csv'))
    print('reference set:', table.shape[0], 'precatalysts x', table.shape[1], 'conditions')


if __name__ == '__main__':
    main()
