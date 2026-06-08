#!/usr/bin/env python3
"""
Build a 3-level dictionary from:
  1) one tab-delimited mass distribution .csv file
  2) 25 whitespace-delimited static stability .txt files

Dictionary structure:
    data[sigma_fcs][fcs_loc][row_number] -> one-row pandas DataFrame

Each stored DataFrame contains:
  - one matching row from the .csv file
  - one row from the .txt file
stacked horizontally.
"""

from pathlib import Path
import io
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines


# ---------------------------------------------------------------------
# User inputs
# ---------------------------------------------------------------------
csv_file = Path(r"C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_results_regional_LH2_250326_7.csv")
# csv_file = Path(r"C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_results_narrowbody_LH2_250326_7.csv")
txt_folder = Path(r"C:\Users\nmb48\Documents\GitHub\SUAVE")
txt_pattern = "suave_static_stability_outputs_7_*_ATR_72-600.txt"
# txt_pattern = "suave_static_stability_outputs_7_*_Airbus_A220-100.txt"


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------
def read_mass_distribution_csv(path: Path) -> pd.DataFrame:
    # Your sample is tab-separated, even though the file extension is .csv.
    return pd.read_csv(path)  #, sep="\t")


def read_static_stability_txt(path: Path) -> pd.DataFrame:
    with open(path, "r") as fh:
        lines = fh.readlines()

    # Remove leading "#" from the header line if present
    if lines and lines[0].lstrip().startswith("#"):
        lines[0] = lines[0].lstrip()[1:].lstrip()

    return pd.read_csv(
        io.StringIO("".join(lines)),
        sep=r"\s+",
        header=0,
    )


def find_matching_csv_row(df_csv: pd.DataFrame, sigma_fcs, fcs_loc) -> pd.DataFrame:
    mask = (
        np.isclose(pd.to_numeric(df_csv["sigma_fcs"], errors="coerce"), float(sigma_fcs))
        & np.isclose(pd.to_numeric(df_csv["fcs_loc"], errors="coerce"), float(fcs_loc))
        # & np.isclose(pd.to_numeric(df_csv["span_loc"], errors="coerce"), float(fcs_loc))
    )

    matches = df_csv.loc[mask]
    if len(matches) == 0:
        raise KeyError(f"No CSV row found for sigma_fcs={sigma_fcs}, fcs_loc={fcs_loc}")
    if len(matches) > 1:
        raise ValueError(f"More than one CSV row found for sigma_fcs={sigma_fcs}, fcs_loc={fcs_loc}")

    return matches.iloc[[0]].reset_index(drop=True)


def merge_one_txt_file(df_csv: pd.DataFrame, txt_path: Path) -> tuple[float, float, dict[int, pd.DataFrame]]:
    df_txt = read_static_stability_txt(txt_path)

    sigma_fcs = df_txt.loc[0, "sigma_fcs"]
    fcs_loc = df_txt.loc[0, "fcs_loc"]
    # fcs_loc = df_txt.loc[0, "span_loc"]

    df_csv_row = find_matching_csv_row(df_csv, sigma_fcs, fcs_loc)

    row_dict = {}
    for row_number, (_, txt_row) in enumerate(df_txt.iterrows(), start=1):
        df_txt_row = txt_row.to_frame().T.reset_index(drop=True)
        merged_row = pd.concat([df_csv_row, df_txt_row], axis=1)
        row_dict[row_number] = merged_row

    return float(sigma_fcs), float(fcs_loc), row_dict


# ---------------------------------------------------------------------
# Main build
# ---------------------------------------------------------------------
df_csv = read_mass_distribution_csv(csv_file)

txt_files = sorted(
    txt_folder.glob(txt_pattern),
    key=lambda p: int(re.search(r"outputs_7_(\d+)_", p.name).group(1)),
    # key=lambda p: int(re.search(r"outputs_6_(\d+)_", p.name).group(1)),
)

data = {}

for txt_path in txt_files:
    sigma_fcs, fcs_loc, row_dict = merge_one_txt_file(df_csv, txt_path)

    if sigma_fcs not in data:
        data[sigma_fcs] = {}
    data[sigma_fcs][fcs_loc] = row_dict


# ---------------------------------------------------------------------
# Example access
# ---------------------------------------------------------------------
# one-row dataframe for sigma_fcs=3375, fcs_loc=0.255, txt row 1
# df_example = data[3375.0][0.255][1]
# print(df_example)

#%%

sigma_fcs_range = list(data.keys())
fcs_loc_range = list(data[1500].keys())

sigma_fcs_grid, fcs_loc_grid = np.meshgrid(sigma_fcs_range, fcs_loc_range, indexing='ij')

x_cg_array = np.zeros((len(sigma_fcs_range), len(fcs_loc_range)))
x_np_array = np.zeros_like(x_cg_array)

for i, sigma_fcs in enumerate(sigma_fcs_range):
    for j, fcs_loc in enumerate(fcs_loc_range):
        
        sub_data = data[sigma_fcs][fcs_loc]
        
        x_cg_array[i, j] = sub_data[2]['x_cg']
        x_np_array[i, j] = sub_data[2]['NP']
        

fig, ax = plt.subplots()

ctf = ax.contourf(sigma_fcs_grid, fcs_loc_grid, x_cg_array - x_np_array)
plt.colorbar(ctf)

plt.show()

#%%

sigma_fcs_range = np.array(sorted(data.keys()), dtype=float)
fcs_loc_range = np.array(sorted(next(iter(data.values())).keys()), dtype=float)
sigma_fcs_grid, fcs_loc_grid = np.meshgrid(sigma_fcs_range, fcs_loc_range, indexing='ij')

# Build one diff surface per txt row number (1..6)
diff_by_row = {
    row_number: np.full((len(sigma_fcs_range), len(fcs_loc_range)), np.nan)
    for row_number in range(1, 7)
}

for i, sigma_fcs in enumerate(sigma_fcs_range):
    for j, fcs_loc in enumerate(fcs_loc_range):
        sub_data = data[sigma_fcs][fcs_loc]
        for row_number in range(1, 7):
            diff_by_row[row_number][i, j] = (
                sub_data[row_number]["x_cg"].iloc[0] - sub_data[row_number]["NP"].iloc[0]
            )

fig, ax = plt.subplots(figsize=(6.4, 4.8))

all_contours = []

for row_number in range(1, 7):
    cs = ax.contour(
        sigma_fcs_grid/1e3,
        fcs_loc_grid,
        diff_by_row[row_number],
        levels=[-2.0],
        colors='0.7',
        linewidths=0.8
    )

    # Extract all paths for this contour level
    for path in cs.collections[0].get_paths():
        verts = path.vertices  # shape (N, 2) → [x, y]
        all_contours.append(verts)


# --- Build common x-grid ---
x_common = sigma_fcs_range/1e3
y_min = np.nanmin(fcs_loc_range)
y_max = np.nanmax(fcs_loc_range)

# Initialize envelope at top
y_envelope = np.full_like(x_common, y_max, dtype=float)


# --- Interpolate each contour onto x_common and take minimum ---
for verts in all_contours:
    x_vals = verts[:, 0]
    y_vals = verts[:, 1]

    # Sort for interpolation (important!)
    order = np.argsort(x_vals)
    x_sorted = x_vals[order]
    y_sorted = y_vals[order]

    # Remove duplicate x values (interp requires strictly increasing)
    x_unique, idx = np.unique(x_sorted, return_index=True)
    y_unique = y_sorted[idx]

    # Interpolate onto common grid
    y_interp = np.interp(
        x_common,
        x_unique,
        y_unique,
        left=np.nan,
        right=np.nan
    )

    # Update envelope (take minimum y across contours)
    mask = ~np.isnan(y_interp)
    y_envelope[mask] = np.minimum(y_envelope[mask], y_interp[mask])


# --- Plot envelope and shading ---
ax.plot(x_common, y_envelope, 'k-')

ax.fill_between(x_common, y_min, y_envelope, color='green', alpha=0.3)
ax.fill_between(x_common, y_envelope, y_max, color='red', alpha=0.3)

ax.set_xlabel('Fuel cell system specific power (kW/kg)')
ax.set_ylabel('Fuel cell system position along fuselage')

ax.tick_params(axis='x', top=False)
ax.tick_params(axis='y', right=False)
ax.spines[['right', 'top']].set_visible(False)
ax.tick_params(axis='y', which='minor', right=False)
ax.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False)

ax.set_ylim(y_min, y_max)

plt.show()

#%%

# SM_min = -1.5
SM_min = -1.0

# -----------------------------------------------------------------------------
# Build grids with axes swapped:
#   x-axis -> fcs_loc
#   y-axis -> sigma_fcs / 1e3
# -----------------------------------------------------------------------------
fcs_loc_range = np.array(sorted(next(iter(data.values())).keys()), dtype=float)
sigma_fcs_range = np.array(sorted(data.keys()), dtype=float)

fcs_loc_grid, sigma_fcs_grid = np.meshgrid(
    fcs_loc_range,
    sigma_fcs_range / 1e3,
    indexing='xy'
)

# -----------------------------------------------------------------------------
# Build one diff surface per txt row number (1..6)
# Shape is (len(sigma_fcs_range), len(fcs_loc_range))
# -----------------------------------------------------------------------------
diff_by_row = {
    row_number: np.full((len(sigma_fcs_range), len(fcs_loc_range)), np.nan)
    for row_number in range(1, 7)
}

for i, sigma_fcs in enumerate(sigma_fcs_range):
    for j, fcs_loc in enumerate(fcs_loc_range):
        sub_data = data[sigma_fcs][fcs_loc]
        for row_number in range(1, 7):
            diff_by_row[row_number][i, j] = (
                sub_data[row_number]["x_cg"].iloc[0] - sub_data[row_number]["NP"].iloc[0]
            )

fig, ax = plt.subplots(figsize=(4.8 * 1.1, 4.8))

all_contours = []

# -----------------------------------------------------------------------------
# Plot contours with swapped axes:
#   x = fcs_loc
#   y = sigma_fcs / 1e3
# -----------------------------------------------------------------------------

# =============================================================================
line_styles = [
    # linestyles_dict['dotted'],
    # linestyles_dict['dashed'],
    # linestyles_dict['dashdotted'],
    # linestyles_dict['dotted'],
    # linestyles_dict['dashed'],
    # linestyles_dict['dashdotted'],
    '-',
    ':',
    '--',
    '-.',
    '-',
    ':',
]
contour_handles = []


# Define linestyles and store dummy handles for legend
line_styles = ['-', ':', '--', '-.', '-', ':']
# colors = ['0.75'] * 4 + ['0.1'] * 2 # same as used in contour
colors = np.linspace(0.0, 0.8, 6)
colors = [str(color) for color in colors]
labels = [str(i) for i in range(1, 7)]

legend_handles = [
    mlines.Line2D([], [], color=colors[i], linestyle=line_styles[i], linewidth=1.2)
    for i in range(6)
]
# =============================================================================

for row_number in range(1, 7):
    cs = ax.contour(
        fcs_loc_grid,
        sigma_fcs_grid,
        diff_by_row[row_number],
        levels=[SM_min],
        # colors='0.3',
        # linewidths=1.2,
        # linestyles=line_styles[row_number - 1]
        colors=colors[row_number-1],
        linewidths=1.2,
        linestyles=line_styles[row_number-1],
    )

    # Store one handle per contour for legend
    if len(cs.collections) > 0:
        contour_handles.append(cs.collections[0])

    # Extract paths as before
    for path in cs.collections[0].get_paths():
        all_contours.append(path.vertices)

# -----------------------------------------------------------------------------
# Build common y-grid (because envelope is now x = x(y))
# -----------------------------------------------------------------------------
y_common = sigma_fcs_range / 1e3
x_min = np.nanmin(fcs_loc_range)
x_max = np.nanmax(fcs_loc_range)

# Envelope starts at the far right
x_envelope = np.full_like(y_common, x_max, dtype=float)

# -----------------------------------------------------------------------------
# Interpolate each contour onto y_common and take the minimum x at each y
# -----------------------------------------------------------------------------
for verts in all_contours:
    x_vals = verts[:, 0]
    y_vals = verts[:, 1]

    # Sort by y for interpolation
    order = np.argsort(y_vals)
    y_sorted = y_vals[order]
    x_sorted = x_vals[order]

    # Remove duplicate y values
    y_unique, idx = np.unique(y_sorted, return_index=True)
    x_unique = x_sorted[idx]

    x_interp = np.interp(
        y_common,
        y_unique,
        x_unique,
        left=np.nan,
        right=np.nan
    )

    mask = ~np.isnan(x_interp)
    x_envelope[mask] = np.minimum(x_envelope[mask], x_interp[mask])

# -----------------------------------------------------------------------------
# Plot envelope and shading
# -----------------------------------------------------------------------------
# ax.plot(x_envelope, y_common, 'k-')

ax.fill_betweenx(y_common, x_min, x_envelope, color='green', alpha=0.1)
ax.fill_betweenx(y_common, x_envelope, x_max, color='red', alpha=0.1)

# -----------------------------------------------------------------------------
# Axis labels swapped
# -----------------------------------------------------------------------------
ax.set_xlabel('Fuel cell system position along fuselage', labelpad=20)
ax.set_ylabel('Fuel cell system specific power (kW/kg)', labelpad=20)

# -----------------------------------------------------------------------------
# Axis styling swapped accordingly
# -----------------------------------------------------------------------------
ax.tick_params(axis='x', top=False)
ax.tick_params(axis='y', right=False)

ax.spines[['right', 'top']].set_visible(False)

ax.tick_params(axis='x', which='minor', top=False)
ax.tick_params(axis='y', which='minor', right=False)

ax.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False)

ax.set_xlim(x_min, x_max)
ax.set_ylim(np.nanmin(sigma_fcs_range / 1e3), np.nanmax(sigma_fcs_range / 1e3))

# labels = [f'{i}' for i in range(1, 7)]
labels = ['Take-off', 'Climb', 'Start of cruise', 'End of cruise', 'Descent', 'Landing']

# ax.legend(
#     contour_handles,
#     labels,
#     title='Case',
#     frameon=False,
#     loc='best'
# )
ax.legend(
    legend_handles,
    labels,
    title='Flight condition',
    frameon=False,
    loc='upper left',
    bbox_to_anchor=(0,1),
)

plt.show()