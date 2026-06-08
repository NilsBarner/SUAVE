__all__ = ["build_jvl_grids"]

import os
import re
import glob
import pandas as pd


def build_jvl_grids(folder, load=False):
    pat = re.compile(r'_(\d+)_([0-9]+\.[0-9]+)\.csv$')

    def parse(files):
        rows = []
        for f in files:
            m = pat.search(os.path.basename(f))
            if m:
                rows.append((int(m.group(1)), m.group(2), f))
        return rows

    geom = parse(glob.glob(os.path.join(folder, 'geometry_*.csv')))
    mass = parse(glob.glob(os.path.join(folder, 'mass_distr_*.csv')))

    N_eng = sorted({n for n, _, _ in geom})
    PRs = sorted({p for _, p, _ in geom}, key=float)

    def make_df(entries):
        df = pd.DataFrame(index=N_eng, columns=PRs, dtype=object)
        for n, p, f in entries:
            df.at[n, p] = pd.read_csv(f) if load else f
        return df

    return make_df(geom), make_df(mass)


