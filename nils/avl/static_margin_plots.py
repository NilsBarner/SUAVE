import os
import io
import re
import sys
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from ambiance import Atmosphere
from collections import defaultdict
from matplotlib import gridspec
from matplotlib.cm import ScalarMappable

# from nils.matplotlib_custom_settings import *

ac_type = "regional"  # "regional" or "narrowbody"

#%% Static stability

folder = r"C:\Users\nmb48\Documents\GitHub\SUAVE"
files = glob.glob(os.path.join(
    folder,
    # 'suave_static_stability_outputs_7_*_ATR_72-600.txt'
    'suave_static_stability_outputs_7_*_Airbus_A220-100.txt'
))

# Nested dict: data[i][j] = DataFrame
df_dict = defaultdict(dict)

pattern = re.compile(r"suave_static_stability_outputs_(\d+)_(\d+)\.txt$")

# for f in files:
#     match = pattern.search(os.path.basename(f))
#     _i, _j = map(int, match.groups())
#     # df_dict[_i][_j] = pd.read_csv(f, sep=r"\s+", header=0)
    
#     with open(f, "r") as fh:
#         lines = fh.readlines()

#     # Remove leading '#' from header if present
#     if lines[0].lstrip().startswith("#"):
#         lines[0] = lines[0].lstrip()[1:].lstrip()

#     df_dict[_i][_j] = pd.read_csv(
#         io.StringIO("".join(lines)),
#         sep=r"\s+",
#         header=0
#     )

dfs = []  # <-- add before the loop

for f in files:
    with open(f, "r") as fh:
        lines = fh.readlines()

    # Clean header
    if lines[0].lstrip().startswith("#"):
        lines[0] = lines[0].lstrip()[1:].lstrip()

    df = pd.read_csv(
        io.StringIO("".join(lines)),
        sep=r"\s+",
        header=0
    )

    dfs.append(df)  # <-- collect instead of nested dict

# After loop: concatenate into single dataframe
df_all = pd.concat(dfs, ignore_index=True)
    
titles = ['1. Take-off', '2. Climb', '3. Beginning of cruise', '4. End of cruise', '5. Descent', '6. Landing']

#%%

fig, ax = plt.subplots()
ax.hist(df_all['NP'].to_numpy(), bins=20)
plt.show()