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

# folder = r"C:/Users/nmb48/"
# folder = r"C:\Users\nmb48\Documents\GitHub\SUAVE"
#####
# NILS: used to generate figures in Chapter 7 of thesis
folder = r"C:\Users\nmb48\Documents\GitHub\SUAVE\nils\avl\data\ch7"
#####
files = glob.glob(os.path.join(
    folder,
    # 'suave_dynamic_stability_matrix_[1-7]_*.txt'
    # 'suave_dynamic_stability_matrix_7_*_ATR_72-600.txt'
    'suave_dynamic_stability_matrix_6_*_ATR_72-600.txt'
    # 'suave_dynamic_stability_matrix_6_*_Airbus_A220-100.txt'
    # 'suave_dynamic_stability_matrix_7_*_Airbus_A220-100.txt'
))
# txt_files = [
#     f for f in files
#     if 0 <= int(f.rsplit('_', 1)[1].split('.')[0]) <= 99
# ]
# pattern = re.compile(r"suave_dynamic_stability_matrix_(\d+)_(\d+)\.txt$")
pattern = re.compile(r"suave_dynamic_stability_matrix_(\d+)_(\d+)_ATR_72-600\.txt$")
# pattern = re.compile(r"suave_dynamic_stability_matrix_(\d+)_(\d+)_Airbus_A220-100\.txt$")

txt_files = [
    f for f in files
    if (m := pattern.search(os.path.basename(f))) and 0 <= int(m.group(2)) <= 99
]

#%% Create nested dict where data[i][j] = DataFrame

# NILS (04.10.2026): the two cells below are wrapped in a function so that process_dynamic_stability.py and
# stability_zones.py can combine any set of studies of either aircraft; running this script as before still
# writes suave_dynamic_stability_matrix_combined.txt from the files selected above
pattern_any = re.compile(r"suave_dynamic_stability_matrix_(\d+)_(\d+)(?:_[^\\/]*)?\.txt$")  # NILS (04.10.2026): any aircraft tag

def combine_system_matrices(txt_files, out_file, pattern=pattern_any):

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

    # Write combined output file

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

    return out_file


if __name__ == "__main__":
    out_file = os.path.join(folder, "suave_dynamic_stability_matrix_combined.txt")
    combine_system_matrices(txt_files, out_file, pattern=pattern)
