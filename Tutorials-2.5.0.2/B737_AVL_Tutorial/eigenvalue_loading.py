"""
This script reads all files obeying the pattern
suave_dynamic_stability_outputs_[1-7]_*.txt and
collates them into a single
suave_dynamic_stability_outputs_combined.txt
file that can then be read by
process_dynamic_stability.py.
"""

__all__ = [

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

#%% Search for files following pattern

folder = r"C:/Users/nmb48/"
files = glob.glob(os.path.join(
    folder,
    'suave_dynamic_stability_outputs_[1-7]_*.txt'
))
txt_files = [
    f for f in files
    if 0 <= int(f.rsplit('_', 1)[1].split('.')[0]) <= 99
]
pattern = re.compile(r"suave_dynamic_stability_outputs_(\d+)_(\d+)\.txt$")

#%% Create nested dict where data[i][j] = DataFrame

N_COLS = 27  # expected number of columns per row
df_dict = defaultdict(dict)

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

#%% Write combined output file

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

