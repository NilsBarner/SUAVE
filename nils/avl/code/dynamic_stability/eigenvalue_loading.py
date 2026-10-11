"""
This script reads all files obeying the pattern
suave_dynamic_stability_outputs_[1-7]_*.txt and
collates them into a single
suave_dynamic_stability_outputs_combined.txt
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
    # 'suave_dynamic_stability_outputs_[1-7]_*.txt'
    # "suave_dynamic_stability_outputs_7_*_ATR_72-600.txt",
    "suave_dynamic_stability_outputs_6_*_ATR_72-600.txt",
    # "suave_dynamic_stability_outputs_6_*_Airbus_A220-100.txt",
    # "suave_dynamic_stability_outputs_7_*_Airbus_A220-100.txt",
))
# txt_files = [
#     f for f in files
#     if 0 <= int(f.rsplit('_', 1)[1].split('.')[0]) <= 99
# ]
# pattern = re.compile(r"suave_dynamic_stability_outputs_(\d+)_(\d+)\.txt$")
pattern = re.compile(r"suave_dynamic_stability_outputs_(\d+)_(\d+)_ATR_72-600\.txt$")
# pattern = re.compile(r"suave_dynamic_stability_outputs_(\d+)_(\d+)_Airbus_A220-100\.txt$")

txt_files = [
    f for f in files
    if (m := pattern.search(os.path.basename(f))) and 0 <= int(m.group(2)) <= 99
]

#%% Create nested dict where data[i][j] = DataFrame

# NILS (04.10.2026): the two cells below are wrapped in a function so that process_dynamic_stability.py and
# stability_zones.py can combine any set of studies of either aircraft; running this script as before still
# writes suave_dynamic_stability_outputs_combined.txt from the files selected above
pattern_any = re.compile(r"suave_dynamic_stability_outputs_(\d+)_(\d+)(?:_[^\\/]*)?\.txt$")  # NILS (04.10.2026): any aircraft tag

def combine_eigen_outputs(txt_files, out_file, pattern=pattern_any):

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

        # df_dict[_i][_j] = pd.DataFrame(
        #     rows,
        #     columns=header.split(),
        #     dtype=float
        # )
        # NILS (04.10.2026): fcs_loc/span_loc are the string "eng" in studies 1/2, which dtype=float cannot parse -> nan
        df_dict[_i][_j] = pd.DataFrame(
            rows,
            columns=header.split(),
        ).apply(pd.to_numeric, errors="coerce")

    # Write combined output file

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
                    float_format="%.6f",
                    na_rep="nan",  # NILS (04.10.2026): keep one token per column (default na_rep='' would shift the columns)
                )

    return out_file


if __name__ == "__main__":
    out_file = os.path.join(folder, "suave_dynamic_stability_outputs_combined.txt")
    combine_eigen_outputs(txt_files, out_file, pattern=pattern)
