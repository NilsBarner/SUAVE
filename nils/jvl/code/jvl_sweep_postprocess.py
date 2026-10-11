import matplotlib.pyplot as plt
import pandas as pd
import re
import os
import sys
import glob
import numpy as np

from nils.jvl.code.shared_methods import build_jvl_grids

###

folder = os.path.join(r'C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\nils\jvl\data')  # or any other folder path

# Find all .csv files in above folder containing geometry and mass distribution information

# geometry_files = glob.glob(os.path.join(folder, 'geometry_*.csv'))
# mass_distr_files = glob.glob(os.path.join(folder, 'mass_distr_*.csv'))

geom_df, mass_df = build_jvl_grids(folder)
N_eng_range = list(geom_df.head().index)  # list of ints
Prop_PR_des_range = [float(element) for element in list(geom_df.columns)]
# Prop_PR_des_range = list(geom_df.columns)  # list of strings

N_eng_grid, Prop_PR_des_grid = np.meshgrid(N_eng_range, Prop_PR_des_range, indexing='ij')

###

CLtot_grid_df = pd.read_csv(r'C:\Users\nmb48\Documents\GitHub\SUAVE\nils\jvl\data\CLtot_grid.csv')
CTtot_grid_df = pd.read_csv(r'C:\Users\nmb48\Documents\GitHub\SUAVE\nils\jvl\data\CTtot_grid.csv')
dy_lower_grid_df = pd.read_csv(r'C:\Users\nmb48\Documents\GitHub\SUAVE\nils\jvl\data\dy_lower_grid.csv')
dY_upper_grid_df = pd.read_csv(r'C:\Users\nmb48\Documents\GitHub\SUAVE\nils\jvl\data\dY_upper_grid.csv')

CLtot_grid_array = CLtot_grid_df.to_numpy().T
CTtot_grid_array = CTtot_grid_df.to_numpy().T
dy_lower_grid_array = dy_lower_grid_df.to_numpy().T
dY_upper_grid_array = dY_upper_grid_df.to_numpy().T

###

fig, ax = plt.subplots()

ax.scatter(dy_lower_grid_array, CLtot_grid_array)

plt.show()


