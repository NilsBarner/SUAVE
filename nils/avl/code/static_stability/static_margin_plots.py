"""
This script plots the static stability boundaries at six mission points.
The green area comprises all those specific powers and longitudinal centre
of gravity locations where the aircraft is stable throughout the entire mission.
"""

__all__ = []

import os
import io
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from pathlib import Path


def _contour_segments(cs):
    """NILS (04.10.2026): version-independent extraction of contour line segments. Matplotlib >= 3.8 removed
    ContourSet.collections, and ContourSet.get_paths() returns one (possibly EMPTY) path per level, which made
    np.interp fail with "array of sample points is empty" whenever a criterion is not crossed anywhere in the
    sweep (e.g. SM > SM_min everywhere). ContourSet.allsegs exists in all versions and lists only real segments."""
    return [seg for level_segs in cs.allsegs for seg in level_segs if len(seg) >= 2]

# ---------------------------------------------------------------------
# User inputs
# ---------------------------------------------------------------------

# csv_folder = r"C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\nils\airbus\data"
# txt_folder = Path(r"C:\Users\nmb48\Documents\GitHub\SUAVE\nils\avl\data")
#####
csv_folder = r"C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv"
txt_folder = Path(r"C:\Users\nmb48\Documents\GitHub\SUAVE\nils\avl\data\ch7")
#####

# ATR 72-600
# csv_file = os.path.join(csv_folder, "mass_distr_results_regional_LH2_250326_7.csv")
# txt_pattern = "suave_static_stability_outputs_7_*_ATR_72-600.txt"
#####
csv_file = os.path.join(csv_folder, "mass_distr_results_regional_LH2_041026_6.csv")
txt_pattern = "suave_static_stability_outputs_6_*_ATR_72-600.txt"
#####

# # Airbus A220-100
# csv_file = os.path.join(csv_folder, "mass_distr_results_narrowbody_LH2_250326_7.csv")
# txt_pattern = "suave_static_stability_outputs_7_*_Airbus_A220-100.txt"

# NILS (04.10.2026): location variable swept in the selected study (fuselage studies 3/7: fcs_loc, wing studies 4/6: span_loc)
study_idx = int(re.search(r"outputs_(\d+)_", txt_pattern).group(1))
LOC_COL = {3: "fcs_loc", 4: "span_loc", 6: "span_loc", 7: "fcs_loc"}.get(study_idx, "fcs_loc")
LOC_LABEL = {"fcs_loc": "Fuel cell system position along fuselage", "span_loc": "Fuel cell system position along wing semi-span"}[LOC_COL]

# NILS (03.10.2026, fix A3): stability criteria
# Aft CG limit: minimum static margin (fraction of MAC) as in TASOPT's default horizontal tail SM_min.
# The criterion was previously x_cg - x_NP < -1.5 m (-2.0 m in the second cell), i.e. ~34 % MAC
# for the A220-100 and ~62 % MAC for the ATR 72-600.
SM_min = 0.15
# Forward CG limit: CS-25 does not prescribe a single static margin; the forward limit follows from
# elevator authority (e.g. CS 25.143/25.145 controllability, nose-wheel lift-off and landing flare).
# Proxy used here: trimmed elevator deflection (deg, +ve trailing edge down) must not exceed
# DELTA_E_MIN in nose-up (negative) direction, leaving the remaining travel for manoeuvre/rotation.
DELTA_E_MIN = -20.0  # (deg) ~2/3 of a typical -30 deg nose-up elevator travel (assumption, adjust per aircraft)

#%%

# Build a 3-level dictionary from:
#   1) one tab-delimited mass distribution .csv file
#   2) 25 whitespace-delimited static stability .txt files

# Dictionary structure:
#     data[sigma_fcs][fcs_loc][row_number] -> one-row pandas DataFrame

# Each stored DataFrame contains:
#   - one matching row from the .csv file
#   - one row from the .txt file
# stacked horizontally.

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
        # & np.isclose(pd.to_numeric(df_csv["fcs_loc"], errors="coerce"), float(fcs_loc))
        & np.isclose(pd.to_numeric(df_csv[LOC_COL], errors="coerce"), float(fcs_loc))  # NILS (04.10.2026): fcs_loc is constant (-1) in the wing studies
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
    # fcs_loc = df_txt.loc[0, "fcs_loc"]
    # fcs_loc = df_txt.loc[0, "span_loc"]
    fcs_loc = df_txt.loc[0, LOC_COL]  # NILS (04.10.2026): generic location variable of the selected study

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
    # key=lambda p: int(re.search(r"outputs_7_(\d+)_", p.name).group(1)),
    # key=lambda p: int(re.search(r"outputs_6_(\d+)_", p.name).group(1)),
    key=lambda p: int(re.search(r"outputs_\d+_(\d+)_", p.name).group(1)),  # NILS (04.10.2026): any study index
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

#%% Plot distance between centre of gravity and neutral point as a function of
# specific power and longitudinal location of fuel cell system

sigma_fcs_range = list(data.keys())
fcs_loc_range = list(data[1500].keys())

sigma_fcs_grid, fcs_loc_grid = np.meshgrid(sigma_fcs_range, fcs_loc_range, indexing='ij')

x_cg_array = np.zeros((len(sigma_fcs_range), len(fcs_loc_range)))
x_np_array = np.zeros_like(x_cg_array)
sm_array = np.zeros_like(x_cg_array)  # NILS (03.10.2026, fix A3)

for i, sigma_fcs in enumerate(sigma_fcs_range):
    for j, fcs_loc in enumerate(fcs_loc_range):
        
        sub_data = data[sigma_fcs][fcs_loc]
        
        x_cg_array[i, j] = sub_data[2]['x_cg'].iloc[0]
        x_np_array[i, j] = sub_data[2]['NP'].iloc[0]
        sm_array[i, j] = sub_data[2]['static_margin'].iloc[0]  # NILS (03.10.2026, fix A3): (x_NP - x_cg) / MAC from AVL.py, per-case CG  # NILS (04.10.2026): .iloc[0] here and in the two lines above (one-row frames fail with NumPy >= 2)
        

fig, ax = plt.subplots()

# ctf = ax.contourf(sigma_fcs_grid, fcs_loc_grid, x_cg_array - x_np_array)
ctf = ax.contourf(sigma_fcs_grid, fcs_loc_grid, sm_array)  # NILS (03.10.2026, fix A3): static margin instead of metres
plt.colorbar(ctf)

plt.show()

#%% Plot zones of stability and instability

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
            # diff_by_row[row_number][i, j] = (
            #     sub_data[row_number]["x_cg"].iloc[0] - sub_data[row_number]["NP"].iloc[0]
            # )
            diff_by_row[row_number][i, j] = sub_data[row_number]["static_margin"].iloc[0]  # NILS (03.10.2026, fix A3)

fig, ax = plt.subplots(figsize=(6.4, 4.8))

all_contours = []

for row_number in range(1, 7):
    cs = ax.contour(
        sigma_fcs_grid/1e3,
        fcs_loc_grid,
        diff_by_row[row_number],
        levels=[SM_min],  # NILS (03.10.2026, fix A3): was [-2.0] (metres)
        colors='0.7',
        linewidths=0.8
    )

    # Extract all paths for this contour level
    # for path in cs.get_paths():
    #     verts = path.vertices  # shape (N, 2) → [x, y]
    #     all_contours.append(verts)
    for verts in _contour_segments(cs):  # NILS (04.10.2026): see _contour_segments()
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
ax.set_ylabel(LOC_LABEL)  # NILS (04.10.2026): was 'Fuel cell system position along fuselage'

ax.tick_params(axis='x', top=False)
ax.tick_params(axis='y', right=False)
ax.spines[['right', 'top']].set_visible(False)
ax.tick_params(axis='y', which='minor', right=False)
ax.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False)

ax.set_ylim(y_min, y_max)

plt.show()

#%% Transpose of above function

# SM_min = -1.5
# SM_min = -1.0
# NILS (03.10.2026, fix A3): SM_min now defined under user inputs at the top (0.15 MAC)

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
de_by_row = {  # NILS (03.10.2026, fix A3): trimmed elevator deflection for forward limit
    row_number: np.full((len(sigma_fcs_range), len(fcs_loc_range)), np.nan)
    for row_number in range(1, 7)
}

for i, sigma_fcs in enumerate(sigma_fcs_range):
    for j, fcs_loc in enumerate(fcs_loc_range):
        sub_data = data[sigma_fcs][fcs_loc]
        for row_number in range(1, 7):
            # diff_by_row[row_number][i, j] = (
            #     sub_data[row_number]["x_cg"].iloc[0] - sub_data[row_number]["NP"].iloc[0]
            # )
            diff_by_row[row_number][i, j] = sub_data[row_number]["static_margin"].iloc[0]  # NILS (03.10.2026, fix A3)
            if "delta_e" in sub_data[row_number].columns:
                de_by_row[row_number][i, j] = sub_data[row_number]["delta_e"].iloc[0]

fig, ax = plt.subplots(figsize=(4.8 * 1.1, 4.8))

all_contours = []

# -----------------------------------------------------------------------------
# Plot contours with swapped axes:
#   x = fcs_loc
#   y = sigma_fcs / 1e3
# -----------------------------------------------------------------------------

line_styles = [
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
colors = np.linspace(0.0, 0.8, 6)
colors = [str(color) for color in colors]
labels = [str(i) for i in range(1, 7)]

legend_handles = [
    mlines.Line2D([], [], color=colors[i], linestyle=line_styles[i], linewidth=1.2)
    for i in range(6)
]

for row_number in range(1, 7):
    cs = ax.contour(
        fcs_loc_grid,
        sigma_fcs_grid,
        diff_by_row[row_number],
        levels=[SM_min],
        colors=colors[row_number-1],
        linewidths=1.2,
        linestyles=line_styles[row_number-1],
    )

    # Store one handle per contour for legend
    # if len(cs.collections) > 0:  # NILS (04.10.2026): ContourSet.collections removed in Matplotlib >= 3.8 (legend uses legend_handles)
    #     contour_handles.append(cs.collections[0])

    # Extract paths as before
    # for path in cs.collections[0].get_paths():
    #     all_contours.append(path.vertices)
    for verts in _contour_segments(cs):  # NILS (04.10.2026): see _contour_segments()
        all_contours.append(verts)

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

# NILS (03.10.2026, fix A3): forward limit from elevator authority (trimmable where delta_e >= DELTA_E_MIN,
# i.e. aft of the contour); envelope = most restrictive (largest x) contour across all six mission points
x_envelope_fwd = np.full_like(y_common, x_min, dtype=float)
for row_number in range(1, 7):
    if np.all(np.isnan(de_by_row[row_number])):
        continue
    cs_fwd = ax.contour(
        fcs_loc_grid,
        sigma_fcs_grid,
        de_by_row[row_number],
        levels=[DELTA_E_MIN],
        colors=colors[row_number-1],
        linewidths=0.6,
        linestyles=line_styles[row_number-1],
    )
    # for path in cs_fwd.get_paths():
    #     verts = path.vertices
    for verts in _contour_segments(cs_fwd):  # NILS (04.10.2026): see _contour_segments()
        order = np.argsort(verts[:, 1])
        y_unique, idx = np.unique(verts[order, 1], return_index=True)
        x_unique = verts[order, 0][idx]
        x_interp = np.interp(y_common, y_unique, x_unique, left=np.nan, right=np.nan)
        mask = ~np.isnan(x_interp)
        x_envelope_fwd[mask] = np.maximum(x_envelope_fwd[mask], x_interp[mask])

# NILS (04.10.2026): no contour exists where a criterion is violated over the whole range of locations,
# in which case the envelopes above would wrongly remain at "all feasible" for that specific power
for i in range(len(y_common)):
    for row_number in range(1, 7):
        if np.all(diff_by_row[row_number][i, :] < SM_min):
            x_envelope[i] = x_min
        if np.all(de_by_row[row_number][i, :] < DELTA_E_MIN):
            x_envelope_fwd[i] = x_max

# ax.fill_betweenx(y_common, x_min, x_envelope, color='green', alpha=0.1)
# ax.fill_betweenx(y_common, x_envelope, x_max, color='red', alpha=0.1)
ax.fill_betweenx(y_common, x_min, x_envelope_fwd, color='red', alpha=0.1)  # untrimmable (forward limit)
ax.fill_betweenx(y_common, x_envelope_fwd, x_envelope, where=(x_envelope > x_envelope_fwd), color='green', alpha=0.1)  # stable and trimmable
ax.fill_betweenx(y_common, np.maximum(x_envelope, x_envelope_fwd), x_max, color='red', alpha=0.1)  # SM < SM_min (aft limit)

# -----------------------------------------------------------------------------
# Axis labels swapped
# -----------------------------------------------------------------------------
ax.set_xlabel(LOC_LABEL, labelpad=20)  # NILS (04.10.2026): was 'Fuel cell system position along fuselage'
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