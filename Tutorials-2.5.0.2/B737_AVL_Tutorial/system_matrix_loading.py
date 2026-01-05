"""
This script reads all files obeying the pattern
suave_dynamic_stability_matrix_[1-7]_*.txt and
collates them into a single
suave_dynamic_stability_matrix_combined.txt
file that can then be read by
process_dynamic_stability.py.
"""

__all__ = []

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
    'suave_dynamic_stability_matrix_[1-7]_*.txt'
))
txt_files = [
    f for f in files
    if 0 <= int(f.rsplit('_', 1)[1].split('.')[0]) <= 99
]
pattern = re.compile(r"suave_dynamic_stability_matrix_(\d+)_(\d+)\.txt$")

#%% Create nested dict where data[i][j] = DataFrame

df_dict = defaultdict(dict)

for f in txt_files:
    match = pattern.search(os.path.basename(f))
    _i, _j = map(int, match.groups())

    with open(f, "r") as fh:
        lines = fh.readlines()

    # Remove leading '#' from header if present
    if lines[0].lstrip().startswith("#"):
        lines[0] = lines[0].lstrip()[1:].lstrip()

    df = pd.read_csv(
        io.StringIO("".join(lines)),
        sep=r"\s+",
        header=0
    )

    # Drop the '|' column
    if '|' in df.columns:
        df = df.drop(columns='|')

    df_dict[_i][_j] = df

#%% Write combined output file

out_file = os.path.join(folder, "suave_dynamic_stability_matrix_combined.txt")

with open(out_file, "w") as fh:
    first = True
    for i in sorted(df_dict):
        for j in sorted(df_dict[i]):
            df = df_dict[i][j]

            if not first:
                fh.write("\n")
            first = False

            nrows = len(df)

            # Write in blocks
            for row_start in range(0, nrows, 12):
                block = df.iloc[row_start:row_start + 12]

                block.to_csv(
                    fh,
                    sep=" ",
                    index=False,
                    header=(row_start == 0 or row_start % 72 == 0),
                    float_format="%.6f"
                )

                # One empty line after every 12 rows
                fh.write("\n")

                # Extra empty line after every 72 rows
                if (row_start + 12) % 72 == 0:
                    fh.write("\n")
