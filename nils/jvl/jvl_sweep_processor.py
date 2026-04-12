import matplotlib.pyplot as plt
import pandas as pd
import re
import os
import sys
import glob
import numpy as np

###

folder = os.path.join(r'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\nils\jvl\data')  # or any other folder path

# Find all .csv files in above folder containing geometry and mass distribution information

# geometry_files = glob.glob(os.path.join(folder, 'geometry_*.csv'))
# mass_distr_files = glob.glob(os.path.join(folder, 'mass_distr_*.csv'))

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
    PRs   = sorted({p for _, p, _ in geom}, key=float)

    def make_df(entries):
        df = pd.DataFrame(index=N_eng, columns=PRs, dtype=object)
        for n, p, f in entries:
            df.at[n, p] = pd.read_csv(f) if load else f
        return df

    return make_df(geom), make_df(mass)


geom_df, mass_df = build_jvl_grids(folder)

N_eng_range = list(geom_df.head().index)  # list of ints
Prop_PR_des_range = [float(element) for element in list(geom_df.columns)]
# Prop_PR_des_range = list(geom_df.columns)  # list of strings

N_eng_grid, Prop_PR_des_grid = np.meshgrid(N_eng_range, Prop_PR_des_range, indexing='ij')

###

CLtot_grid_df = pd.read_csv(r'C:\Users\nmb48\Documents\GitHub\SUAVE\Tutorials-2.5.0.2\B737_AVL_Tutorial\jvl\CLtot_grid.csv')
CTtot_grid_df = pd.read_csv(r'C:\Users\nmb48\Documents\GitHub\SUAVE\Tutorials-2.5.0.2\B737_AVL_Tutorial\jvl\CTtot_grid.csv')
dy_lower_grid_df = pd.read_csv(r'C:\Users\nmb48\Documents\GitHub\SUAVE\Tutorials-2.5.0.2\B737_AVL_Tutorial\jvl\dy_lower_grid.csv')
dY_upper_grid_df = pd.read_csv(r'C:\Users\nmb48\Documents\GitHub\SUAVE\Tutorials-2.5.0.2\B737_AVL_Tutorial\jvl\dY_upper_grid.csv')

CLtot_grid_array = CLtot_grid_df.to_numpy().T
CTtot_grid_array = CTtot_grid_df.to_numpy().T
dy_lower_grid_array = dy_lower_grid_df.to_numpy().T
dY_upper_grid_array = dY_upper_grid_df.to_numpy().T

###

fig, ax = plt.subplots()

ax.scatter(dy_lower_grid_array, CLtot_grid_array)

plt.show()


