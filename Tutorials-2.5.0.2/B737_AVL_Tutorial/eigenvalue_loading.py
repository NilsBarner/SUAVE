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

#%% Static stability

folder = r"C:/Users/nmb48/"
files = glob.glob(os.path.join(
    folder,
    'suave_dynamic_stability_outputs_[1-7]_*.txt'
))
txt_files = [
    f for f in files
    if 0 <= int(f.rsplit('_', 1)[1].split('.')[0]) <= 99
]

# Nested dict: data[i][j] = DataFrame
df_dict = defaultdict(dict)

pattern = re.compile(r"suave_dynamic_stability_outputs_(\d+)_(\d+)\.txt$")

#%%

# for f in txt_files:
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

N_COLS = 27  # expected number of columns per row

for f in txt_files:
    match = pattern.search(os.path.basename(f))
    _i, _j = map(int, match.groups())

    with open(f, "r") as fh:
        raw_lines = fh.readlines()

    # Clean header
    header = raw_lines[0].lstrip("#").strip()

    data_tokens = []
    rows = []

    for line in raw_lines[1:]:
        if not line.strip():
            continue

        data_tokens.extend(line.split())

        while len(data_tokens) >= N_COLS:
            rows.append(data_tokens[:N_COLS])
            data_tokens = data_tokens[N_COLS:]

    df_dict[_i][_j] = pd.DataFrame(
        rows,
        columns=header.split(),
        dtype=float
    )

    
# %% Write combined output file

out_file = os.path.join(folder, "suave_dynamic_stability_outputs_combined.txt")

with open(out_file, "w") as fh:
    first = True
    for i in sorted(df_dict):
        for j in sorted(df_dict[i]):
            if not first:
                fh.write("\n")  # empty line between blocks
            first = False

            df_dict[i][j].to_csv(
                fh,
                sep=" ",
                index=False,
                float_format="%.6f"
            )

