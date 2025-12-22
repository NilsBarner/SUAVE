"""
This script post-processes the result of the dynamic stability
analysis conducted using the TASOPT-SUAVE-AVL wrapper.
"""

__all__ = []

import os
import re
import sys
import copy
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from collections import defaultdict
from ambiance import Atmosphere
from matplotlib import gridspec
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

from matplotlib_custom_settings import *

#%%

import os
import re
import glob
import numpy as np
import pandas as pd

# -------------------------
# utils: build the dataset
# -------------------------
def build_dyn_dataset(root_dir):
    """
    Scan root_dir for files named
      suave_dynamic_stability_outputs_<x>_<y>.txt
    Read each file, annotate with x (int), file_y (int) and row (int),
    and return a single DataFrame with a MultiIndex:
    (x, sigma_fcs, span_loc, fcs_loc, wing_frac, nacelle_frac, row)
    """
    pattern = os.path.join(root_dir, "suave_dynamic_stability_outputs_*.txt")
    files = sorted(glob.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No matching files in {root_dir}")

    rows = []
    filename_re = re.compile(r"suave_dynamic_stability_outputs_(\d+)_(\d+)\.txt$")
    for fn in files:
        m = filename_re.search(os.path.basename(fn))
        if not m:
            continue
        x = int(m.group(1))
        y = int(m.group(2))
        # read as whitespace-delimited table
        df = pd.read_csv(fn, sep=r'\s+', engine='python', header=0)
        # annotate meta columns
        df = df.copy()
        df['x'] = x
        df['file_y'] = y
        df['row'] = np.arange(len(df))  # row number within this file
        rows.append(df)

    df_all = pd.concat(rows, ignore_index=True, sort=False)

    # ensure numeric columns have numeric dtype
    numeric_cols = [c for c in df_all.columns if c not in ('x','file_y','row')]
    df_all[numeric_cols] = df_all[numeric_cols].apply(pd.to_numeric, errors='coerce')

    # create the MultiIndex requested:
    key_cols = ['x', 'sigma_fcs', 'span_loc', 'fcs_loc', 'wing_frac', 'nacelle_frac', 'row']
    df_all.set_index(key_cols, inplace=True)

    return df_all


def _expand_key(x, first5, row):
    """
    Replace None with slice(None) for MultiIndex lookup.
    """
    xk = slice(None) if x is None else x
    f5 = tuple(slice(None) if v is None else v for v in first5)
    rk = slice(None) if row is None else row
    return (xk, *f5, rk)


# -------------------------
# small selection helpers
# -------------------------
RE_COLS = [f"Re{i}" for i in range(1,9)]
IM_COLS = [f"Im{i}" for i in range(1,9)]

def get_eigen_parts(df_all, x, first5, row, part='Re', idx=None, exact=True, tol=1e-9):
    assert part in ('Re','Im'), "part must be 'Re' or 'Im'"
    cols = RE_COLS if part == 'Re' else IM_COLS

    key = _expand_key(x, first5, row)

    # ---- wildcard path (fast, native MultiIndex slicing) ----
    if x is None or None in first5 or row is None:
        out = df_all.loc[key, cols]
        vals = out.to_numpy(dtype=float)
        if idx is None:
            return vals
        else:
            assert 1 <= int(idx) <= 8
            return vals[:, int(idx)-1]

    # ---- original exact / fuzzy logic (unchanged) ----
    try:
        ser = df_all.loc[key]
    except KeyError:
        if exact:
            raise

        x_mask = (df_all.index.get_level_values('x') == x)
        sigma = df_all.index.get_level_values('sigma_fcs').astype(float)
        span  = df_all.index.get_level_values('span_loc').astype(float)
        fcs   = df_all.index.get_level_values('fcs_loc').astype(float)
        wing  = df_all.index.get_level_values('wing_frac').astype(float)
        nac   = df_all.index.get_level_values('nacelle_frac').astype(float)
        row_vals = df_all.index.get_level_values('row').astype(int)

        cond = (
            x_mask &
            (np.abs(sigma - first5[0]) <= tol) &
            (np.abs(span  - first5[1]) <= tol) &
            (np.abs(fcs   - first5[2]) <= tol) &
            (np.abs(wing  - first5[3]) <= tol) &
            (np.abs(nac   - first5[4]) <= tol) &
            (row_vals == row)
        )
        matches = df_all[cond]
        if len(matches) == 0:
            raise KeyError("No matching row found (even fuzzy).")
        ser = matches.iloc[0]

    vals = ser[cols].to_numpy(dtype=float)
    return vals if idx is None else float(vals[int(idx)-1])

    
    
def get_flight_condition(df_all, x, first5, row, exact=True, tol=1e-9):
    key = _expand_key(x, first5, row)

    # ---- wildcard path ----
    if x is None or None in first5 or row is None:
        cols = ['AoA','Mach','Beta','h','n','W']
        out = df_all.loc[key, cols]
        return out.to_numpy(dtype=float)

    # ---- original exact / fuzzy logic ----
    try:
        ser = df_all.loc[key]
    except KeyError:
        if exact:
            raise

        x_mask = (df_all.index.get_level_values('x') == x)
        sigma = df_all.index.get_level_values('sigma_fcs').astype(float)
        span  = df_all.index.get_level_values('span_loc').astype(float)
        fcs   = df_all.index.get_level_values('fcs_loc').astype(float)
        wing  = df_all.index.get_level_values('wing_frac').astype(float)
        nac   = df_all.index.get_level_values('nacelle_frac').astype(float)
        row_vals = df_all.index.get_level_values('row').astype(int)

        cond = (
            x_mask &
            (np.abs(sigma - first5[0]) <= tol) &
            (np.abs(span  - first5[1]) <= tol) &
            (np.abs(fcs   - first5[2]) <= tol) &
            (np.abs(wing  - first5[3]) <= tol) &
            (np.abs(nac   - first5[4]) <= tol) &
            (row_vals == row)
        )
        matches = df_all[cond]
        if len(matches) == 0:
            raise KeyError("No matching row found (even fuzzy).")
        ser = matches.iloc[0]

    return {
        'AoA':  float(ser['AoA']),
        'Mach': float(ser['Mach']),
        'Beta': float(ser['Beta']),
        'h':    float(ser['h']),
        'n':    float(ser['n']),
        'W':    float(ser['W']),
    }


def get_mass_parameters(df_all, x, first5, row, exact=True, tol=1e-9):
    """
    Retrieve mass magnitude/location parameters:
    (sigma_fcs, span_loc, fcs_loc, wing_frac, nacelle_frac)

    Supports wildcards (None) for x, first5 entries, and row.
    Returns an (N,5) numpy array.
    """

    # build boolean mask explicitly (robust to partial slicing)
    mask = np.ones(len(df_all), dtype=bool)

    if x is not None:
        mask &= (df_all.index.get_level_values('x') == x)

    names = ['sigma_fcs','span_loc','fcs_loc','wing_frac','nacelle_frac']
    for name, val in zip(names, first5):
        if val is not None:
            vals = df_all.index.get_level_values(name).astype(float)
            mask &= (np.abs(vals - val) <= tol)

    if row is not None:
        mask &= (df_all.index.get_level_values('row') == row)

    sel = df_all[mask]
    if len(sel) == 0:
        raise KeyError("No matching mass-parameter rows found.")

    idx = sel.index

    return np.column_stack([
        idx.get_level_values('sigma_fcs').astype(float),
        idx.get_level_values('span_loc').astype(float),
        idx.get_level_values('fcs_loc').astype(float),
        idx.get_level_values('wing_frac').astype(float),
        idx.get_level_values('nacelle_frac').astype(float),
    ])

#%%

# @NILS: could the complex eigenvalues at low Mach numbers have to do
# witht the wing being stalled there?

# mach_vals = np.array([0.05, 0.15, 0.25, 0.45, 0.65, 0.85])
# aoa_vals = np.array([-2., 0., 2., 5., 7., 10.])

# NOTE: cbar is the MAC (see page 98 in AE3202 Flight Dynamics Lecture Notes)

# Table 3.9 (phugoid)

# cbar = 4.2350  # from "C:\Users\nmb48\avl_files\body_axis_derivatives_case_01_01.txt"
# V = 0.05 * Atmosphere(11e3).speed_of_sound[0]
# n = 1

def zeta(lamda):  # (-)
    return -lamda.real / np.sqrt(lamda.real**2 + lamda.imag**2)  # page 129 in AE3202 Flight Dynamics Lecture Notes

def T0p5(lamda):  # {s}
    return -np.log(0.5) / lamda.real  # pages 125 and 127 in AE3202 Flight Dynamics Lecture Notes

def T2(lamda):  # (s)
    return -T0p5(lamda)  # page 127 in AE3202 Flight Dynamics Lecture Notes

def omega0(lamda):  # (1/s)
    return np.sqrt(lamda.real**2 + lamda.imag**2)  # page 129 in AE3202 Flight Dynamics Lecture Notes

def CAP(lamda, n, alpha):  # (1/s^2)
    return omega0(lamda)**2 / (n / alpha)  # (3.111) in medium_PhDthesis_2013_flightdynconstrconcacmdiscanalsopt_morris

def tau(lamda):  # (s)
    return -np.log(0.5) * T0p5(lamda)  # page 125 in AE3202 Flight Dynamics Lecture Notes

def compute_lambdas_from_zeta_omega(zeta0, omega0):
    
    x = -zeta0 * omega0
    y_mag = omega0 * np.sqrt(max(0.0, 1.0 - zeta0**2))

    cplx = x + 1j * y_mag
    cplx_conj = x - 1j * y_mag
    return cplx, cplx_conj

folder = r"C:/Users/nmb48/"
txt_files = glob.glob(os.path.join(folder, 'suave_dynamic_stability_outputs_[1-7]_[0-99].txt'))  # material_distr_eng_pos_1.mat gives nonsensical results!

# Nested dict: data[i][j] = DataFrame
df_dict = defaultdict(dict)

pattern = re.compile(r"suave_dynamic_stability_outputs_(\d+)_(\d+)\.txt$")

for f in txt_files:
    match = pattern.search(os.path.basename(f))
    if not match:
        continue
    _i, _j = map(int, match.groups())
    df_dict[_i][_j] = pd.read_csv(f, sep=r"\s+", engine="python")
    
#%% Calculate MIL-STD-1797 requirements (section 3.4.2 in medium_PhDthesis_2013_flightdynconstrconcacmdiscanalsopt_morris)

def plot_mil_limits(aoa, mach, beta, h, n, W):
    
    amb = Atmosphere(h)
    V = mach * amb.speed_of_sound[0]
    aoa = np.deg2rad(aoa)
    beta = np.deg2rad(beta)

    N_points = 50
    
    # Lateral stability limits
    zeta_dr_min = 0.08  # 3rd row of Table 3.14 (Class II, Category B, Level I)
    omega_n_dr_min = 0.4  # 3rd row of Table 3.14 (Class II, Category B, Level I)
    T_2_sl_min = 20  # 2nd row of Table 3.13 (Category B, Level I)
    tau_roll_max = 1.4  # 3rd row of Table 3.12 (Class II, Category B, Level I)
    
    # Dutch roll
    
    r_circle_dr_min = omega_n_dr_min
    theta = np.linspace(0, 2 * np.pi, N_points)
    x_circle_dr_min = r_circle_dr_min * np.cos(theta)
    y_circle_dr_min = r_circle_dr_min * np.sin(theta)
    
    # x_vert_dr_max = -zeta_dr_min * r_circle_dr_min * cbar / V
    x_vert_dr_max = -zeta_dr_min * r_circle_dr_min
    
    a_dr = zeta_dr_min
    t_dr = np.linspace(0, r_circle_dr_min / a_dr, N_points)
    x_ray_dr_max = -t_dr * a_dr
    y_ray_dr_max = t_dr * np.sqrt(1 - a_dr**2)
    
    # Spiral
    # x_vert_sl_max = -np.log(0.5) / T_2_sl_min * cbar / V
    x_vert_sl_max = -np.log(0.5) / T_2_sl_min
    
    # Roll
    # x_vert_roll_max = -np.log(0.5)**2 / tau_roll_max * cbar / V
    x_vert_roll_max = -np.log(0.5)**2 / tau_roll_max
    
    ###
    
    # Longitudinal stability limits
    zeta_ph_min = 0.04  # 1st row in Table 3.9 (Class II, Category B, Level I)
    zeta_sp_min = 0.3  # 1st row, 4th column in Table 3.10 (Class II, Category B, Level I)
    zeta_sp_max = 2.0  # 1st row, 5th column in Table 3.10 (Class II, Category B, Level I)
    CAP_sp_min = 0.085  # 1st row, 4th column in Table 3.11 (Class II, Category B, Level I)
    CAP_sp_max = 3.6  # 1st row, 5th column in Table 3.11 (Class II, Category B, Level I)
    
    # Short period
    
    # NOTE!!!: abs() below not in original (3.111) in medium_PhDthesis_2013_flightdynconstrconcacmdiscanalsopt_morris but required to prevent nan when aoa < 0
    omega_n_sp_min = np.sqrt(CAP_sp_min * n / abs(aoa))  # (3.111) in medium_PhDthesis_2013_flightdynconstrconcacmdiscanalsopt_morris rearranged
    omega_n_sp_max = np.sqrt(CAP_sp_max * n / abs(aoa))  # (3.111) in medium_PhDthesis_2013_flightdynconstrconcacmdiscanalsopt_morris rearranged
    r_circle_sp_min = omega_n_sp_min
    r_circle_sp_max = omega_n_sp_max
    theta = np.linspace(0, 2 * np.pi, N_points)
    x_circle_sp_min = r_circle_sp_min * np.cos(theta)
    y_circle_sp_min = r_circle_sp_min * np.sin(theta)
    x_circle_sp_max = r_circle_sp_max * np.cos(theta)
    y_circle_sp_max = r_circle_sp_max * np.sin(theta)
    
    a_sp = zeta_sp_min
    t_sp = np.linspace(0, r_circle_sp_max / a_sp, N_points)
    x_ray_sp_max = -t_sp * a_sp
    y_ray_sp_max = t_sp * np.sqrt(1 - a_sp**2)
    
    # Phugoid
    a_ph = zeta_ph_min
    t_ph = np.linspace(0, 10 / np.sqrt(1 - a_ph**2), N_points)  # ax.get_ylim()[1]
    x_ray_ph_max = -t_ph * a_ph
    y_ray_ph_max = t_ph * np.sqrt(1 - a_ph**2)
    
    return (
        x_circle_dr_min, y_circle_dr_min, x_vert_dr_max, x_ray_dr_max, y_ray_dr_max, r_circle_dr_min, zeta_dr_min,
        x_vert_sl_max,
        x_vert_roll_max,
        
        x_circle_sp_min, y_circle_sp_min, x_circle_sp_max, y_circle_sp_max, x_ray_sp_max, y_ray_sp_max, r_circle_sp_min, r_circle_sp_max, zeta_sp_min,
        x_ray_ph_max, y_ray_ph_max,
    )
    
#%% Calculate MIL-STD-1797 requirements (section 3.4.2 in medium_PhDthesis_2013_flightdynconstrconcacmdiscanalsopt_morris)

def mask_mil_limits(
    ev_roll_list, ev_sp_list, ev_dr_list, ev_sl_list, ev_ph_list,
    fc,
):
    
    aoa, mach, beta, h, n, W = fc.T
    aoa = np.deg2rad(aoa)
    beta = np.deg2rad(beta)
    
    # Lateral stability limits
    zeta_dr_min = 0.08  # 3rd row of Table 3.14 (Class II, Category B, Level I)
    omega_n_dr_min = 0.4  # 3rd row of Table 3.14 (Class II, Category B, Level I)
    T_2_sl_min = 20  # 2nd row of Table 3.13 (Category B, Level I)
    tau_roll_max = 1.4  # 3rd row of Table 3.12 (Class II, Category B, Level I)
    
    zeta_dr = zeta(ev_dr_list)
    omega_n_dr = omega0(ev_dr_list)
    T_2_sl = T2(ev_sl_list)
    tau_roll = tau(ev_roll_list)
    
    # Longitudinal stability limits
    zeta_ph_min = 0.04  # 1st row in Table 3.9 (Class II, Category B, Level I)
    zeta_sp_min = 0.3  # 1st row, 4th column in Table 3.10 (Class II, Category B, Level I)
    zeta_sp_max = 2.0  # 1st row, 5th column in Table 3.10 (Class II, Category B, Level I)
    CAP_sp_min = 0.085  # 1st row, 4th column in Table 3.11 (Class II, Category B, Level I)
    CAP_sp_max = 3.6  # 1st row, 5th column in Table 3.11 (Class II, Category B, Level I)
    
    zeta_ph = zeta(ev_ph_list)
    zeta_sp = zeta(ev_sp_list)
    CAP_sp = CAP(ev_sp_list, n, aoa)
    
    # Dutch roll
    # valid_mask_dr = (zeta_dr > zeta_dr_min) & (omega_n_dr > omega_n_dr_min)
    valid_mask_dr = (zeta_dr > zeta_sp_min)
    
    # Spiral
    valid_mask_sl = (T_2_sl > T_2_sl_min)
    
    # Roll
    valid_mask_roll = (tau_roll < tau_roll_max)
    
    # Short period
    valid_mask_sp = (
        zeta_sp > zeta_sp_min
        # (zeta_sp > zeta_sp_min) & (zeta_sp < zeta_sp_max) &
        # (CAP_sp > CAP_sp_min) & (CAP_sp < CAP_sp_max)
    )
    
    # Phugoid
    valid_mask_ph = (zeta_ph > zeta_ph_min)
    
    ev_dr_valid = ev_dr_list[valid_mask_dr]
    ev_dr_invalid = ev_dr_list[~valid_mask_dr]
    # ev_dr_valid = ev_dr_list[valid_mask_dr & valid_mask_sp]
    # ev_dr_invalid = ev_dr_list[~(valid_mask_dr & valid_mask_sp)]
    # ev_dr_valid = ev_dr_list[valid_mask_sp]
    # ev_dr_invalid = ev_dr_list[~valid_mask_sp]
    
    ev_sl_valid = ev_sl_list[valid_mask_sl]
    ev_sl_invalid = ev_sl_list[~valid_mask_sl]
    
    ev_roll_valid = ev_roll_list[valid_mask_roll]
    ev_roll_invalid = ev_roll_list[~valid_mask_roll]
    
    print('ev_roll_list =', ev_roll_list)
    print()
    
    if np.any(ev_roll_list.imag > 0):
        raise Exception
    # sys.exit()
    
    ev_sp_valid = ev_sp_list[valid_mask_sp]
    ev_sp_invalid = ev_sp_list[~valid_mask_sp]
    # ev_sp_valid = ev_sp_list[valid_mask_sp & valid_mask_dr]
    # ev_sp_invalid = ev_sp_list[~(valid_mask_sp & valid_mask_dr)]
    
    ev_ph_valid = ev_ph_list[valid_mask_ph]
    ev_ph_invalid = ev_ph_list[~valid_mask_ph]
    
    return (
        ev_dr_valid, ev_dr_invalid,
        ev_sl_valid, ev_sl_invalid,
        ev_roll_valid, ev_roll_invalid,
        ev_sp_valid, ev_sp_invalid,
        ev_ph_valid, ev_ph_invalid,
    )

#%%

root_dir = r"C:/Users/nmb48/"   # change to folder containing your .txt files
df_all = build_dyn_dataset(root_dir)

# first5 = (4000.0, None, None, None, None)
# x = None
# row = None

# Re = get_eigen_parts(df_all, x, first5, row, part='Re')
# Im = get_eigen_parts(df_all, x, first5, row, part='Im')
# fc = get_flight_condition(df_all, x, first5, row)

# evs = get_eigen_parts(df_all, x=3, first5=(4000.0, None, None, None, None), row=None)

# evs = get_eigen_parts(df_all, 4, (4000,0.99,-1,1,0.2), row=2)

# print("Flight condition:", fc)
# print("Re eigenvalues:", Re)
# print("Im eigenvalues:", Im)

# sys.exit('Stop.')

#%%

modes_type = 'longitudinal'  # 'longitudinal', 'lateral'

edgecolor = 'none'
linewidth = 0.5

fig = plt.figure(figsize=(15,9), constrained_layout=True)
gs = gridspec.GridSpec(2, 3, width_ratios=[1,1,1], hspace=0.2, height_ratios=[1,1], wspace=0.3)
axes = [fig.add_subplot(gs[_i,_j]) for _i in range(2) for _j in range(3)]

# titles = ['Short period', 'Phugoid', 'Roll subsidence', 'Spiral', 'Dutch roll']
titles = ['1. Take-off', '2. Climb', '3. Beginning of cruise', '4. End of cruise', '5. Descent', '6. Landing']

for row_idx, ax in enumerate(axes):
    
    if modes_type == 'longitudinal':
        
        # Axis settings
        ax.set_xlim((-1.5, 1))
        ax.set_ylim((-0.1, 5))
        ax.set_xticks([-1, -0.5])
        # ax.set_yticks([-5, -4, -3, -2, -1, 1, 2, 3, 4, 5])
        ax.set_yticks([1, 2, 3, 4, 5])
        
        # Figure legend
        markers=['o', 'd']
        legend_marker_colours = [colors[0], colors[1], colors[2], colors[0], colors[1], colors[2]]
        legend_markers = ['o', 'o', 'o', 'd', 'd', 'd']
        legend_marker_labels = [
            'Short-period roots nacelle', 'Short-period roots wing', 'Short-period roots fuselage',
            'Phugoid roots nacelle', 'Phugoid roots wing', 'Phugoid roots fuselage',
        ]
        
        short_period_patch = mpatches.Patch(
            facecolor='none',
            edgecolor='black',
            hatch='//',
            label='Short-period space',
            alpha=0.1,
        )
        phugoid_patch = mpatches.Patch(
            facecolor='lightgrey',
            edgecolor='none',
            label='Phugoid space',
            alpha=0.1,
        )
        marker_entries = [
            Line2D(
                [], [], marker=legend_markers[_], color='none', markerfacecolor=legend_marker_colours[_],
                markeredgecolor=legend_marker_colours[_], markersize=8, linestyle='None',
                label=legend_marker_labels[_],
            ) for _ in range(len(legend_marker_colours))
        ]
        handles = [short_period_patch, phugoid_patch] + marker_entries
        fig.legend(
            handles=handles,
            loc='center left',
            bbox_to_anchor=(0.9, 0.5),
            frameon=False
        )
        
        # Inset axis and settings
        axins = inset_axes(
            ax,
            width="100%",
            height="100%",
            bbox_to_anchor=(0.65, 0.65, 0.35, 0.35),  # (x0, y0, w, h) in ax coords
            bbox_transform=ax.transAxes,
            borderpad=0
        )
        axins.set_xlim(-0.045, 0.0)
        axins.set_ylim(0, 0.75)
        axins.axvline(0.0, color='k')
        axins.axhline(0.0, color='k')
        axins.tick_params(labelleft=False, labelbottom=False)
        axins.tick_params(axis='y', which='both', length=0)
        axins.tick_params(axis='x', which='both', length=0)
        mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")
        
    elif modes_type == 'lateral':    
        markers=['o', 'd', '^']
        
        y_lim_lat = 5
        
        ax.set_xlim((-1.35, 0.1))
        ax.set_ylim((-0.1, y_lim_lat))
        ax.set_xticks([-1, -0.5])
        ax.set_yticks([1, 2])
        
        # Add second x-axis left of break to plot roll subsidence roots    
        gap = 0.05  # fraction of ax width
        ax_left = ax.inset_axes(
            [-0.25 - gap, 0.0, 0.25, 1.0],  # <-- shifted left
            sharey=ax
        )
        ax_left.set_xlim(-65, -4)
        ax_left.set_xticks([-65, -5])
        ax_left.set_ylim(-0.1, y_lim_lat)
        ax_left.spines[['left', 'right','top']].set_visible(False)
        ax_left.tick_params(axis='y', which='both', right=False, length=0)
        ax_left.tick_params(axis='x', which='both', length=0)
        ax_left.spines['bottom'].set_position(('data', 0.0))
        ax_left.tick_params(labelright=False)
        ax_left.tick_params(labelleft=False)
    
    # Extract eigenvalue and flight condition data
    config_tuple = (None, None, None, None, None)
    study_idx = None
    ev_real = get_eigen_parts(df_all, study_idx, config_tuple, row_idx, part='Re')
    ev_imag = get_eigen_parts(df_all, study_idx, config_tuple, row_idx, part='Im')
    fc = get_flight_condition(df_all, study_idx, config_tuple, row_idx)
    mp = get_mass_parameters(df_all, study_idx, config_tuple, row_idx)
    
    # ax.scatter(ev_real, ev_imag, marker='.')#, c=np.tile(mp[:,1], (1,8)))
    
    sigma_fcs_unique = np.unique(mp[:, 0])
    # sigma_fcs_unique = [1500, 4000]
    study_idx_list = [3, 4, 6]
    
    for j, study_idx in enumerate(study_idx_list):
    
        for sigma_fcs in sigma_fcs_unique:
            
            # Extract eigenvalue and flight condition data
            config_tuple = (sigma_fcs, None, None, None, None)
            study_idx = study_idx
            ev_real_filtered = get_eigen_parts(df_all, study_idx, config_tuple, row_idx, part='Re')
            ev_imag_filtered = get_eigen_parts(df_all, study_idx, config_tuple, row_idx, part='Im')
            fc_filtered = get_flight_condition(df_all, study_idx, config_tuple, row_idx)
            # print('fc_filtered =', fc_filtered)
            mp_filtered = get_mass_parameters(df_all, study_idx, config_tuple, row_idx)
            
            ax.scatter(ev_real_filtered, ev_imag_filtered, marker='.', color=colors[j])#, c=np.tile(mp[:,1], (1,8)))
            
            if modes_type == 'lateral':
                ax_left.scatter(ev_real_filtered, ev_imag_filtered, marker='.', color=colors[j])#, c=np.tile(mp[:,1], (1,8)))
            
            # =============================================================================
            ev_roll_list = ev_real_filtered[:,0] + 1j * ev_imag_filtered[:,0]
            ev_sp_list = ev_real_filtered[:,1] + 1j * ev_imag_filtered[:,1]
            ev_dr_list = ev_real_filtered[:,3] + 1j * ev_imag_filtered[:,3]
            # ev_sp_list = ev_real_filtered[:,3] + 1j * ev_imag_filtered[:,3]
            # ev_dr_list = ev_real_filtered[:,1] + 1j * ev_imag_filtered[:,1]
            ev_sl_list = ev_real_filtered[:,5] + 1j * ev_imag_filtered[:,5]
            ev_ph_list = ev_real_filtered[:,6:8] + 1j * ev_imag_filtered[:,6:8]
            
            ev_1_list = ev_roll_list
            ev_2_list = ev_dr_list
            ev_3_list = ev_sp_list
            ev_4_list = ev_sl_list
            ev_5_list = ev_ph_list
            
            (ev_dr_valid, ev_dr_invalid,
            ev_sl_valid, ev_sl_invalid,
            ev_roll_valid, ev_roll_invalid,
            ev_sp_valid, ev_sp_invalid,
            ev_ph_valid, ev_ph_invalid) = \
            mask_mil_limits(
                ev_roll_list, ev_sp_list, ev_dr_list, ev_sl_list, ev_ph_list,
                fc_filtered,
            )
            
            alpha_valid = 0.2
            
            ax.scatter(ev_dr_valid.real, ev_dr_valid.imag, marker='o', facecolor=opaque_color_from_hex(colors[j], alpha=alpha_valid), edgecolor=edgecolor, linewidth=linewidth)
            ax.scatter(ev_dr_invalid.real, ev_dr_invalid.imag, marker='o', facecolor=colors[j], edgecolor=edgecolor, linewidth=linewidth)
            
            ax.scatter(ev_sl_valid.real, ev_sl_valid.imag, marker='s', facecolor=opaque_color_from_hex(colors[j], alpha=alpha_valid), edgecolor=edgecolor, linewidth=linewidth)
            ax.scatter(ev_sl_invalid.real, ev_sl_invalid.imag, marker='s', facecolor=colors[j], edgecolor=edgecolor, linewidth=linewidth)
            
            # ax.scatter(ev_roll_valid.real, ev_roll_valid.imag, marker='d', facecolor=opaque_color_from_hex(colors[j], alpha=alpha_valid), edgecolor=edgecolor, linewidth=linewidth)
            # ax.scatter(ev_roll_invalid.real, ev_roll_invalid.imag, marker='d', facecolor=colors[j], edgecolor=edgecolor, linewidth=linewidth)
            
            ax.scatter(ev_sp_valid.real, ev_sp_valid.imag, marker='^', facecolor=opaque_color_from_hex(colors[j], alpha=alpha_valid), edgecolor=edgecolor, linewidth=linewidth)
            ax.scatter(ev_sp_invalid.real, ev_sp_invalid.imag, marker='^', facecolor=colors[j], edgecolor=edgecolor, linewidth=linewidth)
            
            ax.scatter(ev_ph_valid.real, ev_ph_valid.imag, marker='*', facecolor=opaque_color_from_hex(colors[j], alpha=alpha_valid), edgecolor=edgecolor, linewidth=linewidth)
            ax.scatter(ev_ph_invalid.real, ev_ph_invalid.imag, marker='*', facecolor=colors[j], edgecolor=edgecolor, linewidth=linewidth)
            
            if modes_type == 'lateral':
                
                ax_left.scatter(ev_dr_valid.real, ev_dr_valid.imag, marker='o', facecolor=opaque_color_from_hex(colors[j], alpha=alpha_valid), edgecolor=edgecolor, linewidth=linewidth)
                ax_left.scatter(ev_dr_invalid.real, ev_dr_invalid.imag, marker='o', facecolor=colors[j], edgecolor=edgecolor, linewidth=linewidth)
                
                ax_left.scatter(ev_sl_valid.real, ev_sl_valid.imag, marker='s', facecolor=opaque_color_from_hex(colors[j], alpha=alpha_valid), edgecolor=edgecolor, linewidth=linewidth)
                ax_left.scatter(ev_sl_invalid.real, ev_sl_invalid.imag, marker='s', facecolor=colors[j], edgecolor=edgecolor, linewidth=linewidth)
                
                ax_left.scatter(ev_roll_valid.real, ev_roll_valid.imag, marker='s', facecolor=opaque_color_from_hex(colors[j], alpha=alpha_valid), edgecolor=edgecolor, linewidth=linewidth)
                ax_left.scatter(ev_roll_invalid.real, ev_roll_invalid.imag, marker='s', facecolor=colors[j], edgecolor=edgecolor, linewidth=linewidth)
                
                ax_left.scatter(ev_sp_valid.real, ev_sp_valid.imag, marker='^', facecolor=opaque_color_from_hex(colors[j], alpha=alpha_valid), edgecolor=edgecolor, linewidth=linewidth)
                ax_left.scatter(ev_sp_invalid.real, ev_sp_invalid.imag, marker='^', facecolor=colors[j], edgecolor=edgecolor, linewidth=linewidth)
                
                ax_left.scatter(ev_ph_valid.real, ev_ph_valid.imag, marker='*', facecolor=opaque_color_from_hex(colors[j], alpha=alpha_valid), edgecolor=edgecolor, linewidth=linewidth)
                ax_left.scatter(ev_ph_invalid.real, ev_ph_invalid.imag, marker='*', facecolor=colors[j], edgecolor=edgecolor, linewidth=linewidth)
            
            # =============================================================================
            
            # Calculate MIL-constraints from flight conditions
            aoa = fc_filtered[0,0]
            mach = fc_filtered[0,1]
            beta = fc_filtered[0,2]
            h = fc_filtered[0,3]
            n = fc_filtered[0,4]
            W = fc_filtered[0,5]
            (
                x_circle_dr_min, y_circle_dr_min, x_vert_dr_max, x_ray_dr_max, y_ray_dr_max, r_circle_dr_min, zeta_dr_min,
                x_vert_sl_max,
                x_vert_roll_max,
                
                x_circle_sp_min, y_circle_sp_min, x_circle_sp_max, y_circle_sp_max, x_ray_sp_max, y_ray_sp_max, r_circle_sp_min, r_circle_sp_max, zeta_sp_min,
                x_ray_ph_max, y_ray_ph_max,
            ) = \
            plot_mil_limits(aoa, mach, beta, h, n, W)
            
            # Plot MIL-constraints
            
            if modes_type == 'longitudinal':
                
                # Short period
                ax.plot(x_circle_sp_min, y_circle_sp_min, color=opaque_color_from_hex(mcolors.to_hex("black"), alpha=0.2))
                ax.plot(x_circle_sp_max, y_circle_sp_max, color=opaque_color_from_hex(mcolors.to_hex("black"), alpha=0.2))
                ax.plot(x_ray_sp_max, y_ray_sp_max, '-', color=opaque_color_from_hex(mcolors.to_hex("black"), alpha=0.2))
                ax.plot(x_ray_sp_max, -y_ray_sp_max, '-', color=opaque_color_from_hex(mcolors.to_hex("black"), alpha=0.2))
                
                # Angular range for the wedge
                theta_sp_min = np.arccos(zeta_sp_min)
                angles_sp = np.linspace(np.pi - theta_sp_min, np.pi + theta_sp_min, 200)
        
                # Outer arc
                x_wedge_sp_max = r_circle_sp_max * np.cos(angles_sp)
                y_wedge_sp_max = r_circle_sp_max * np.sin(angles_sp)
        
                # Inner arc (reverse direction so polygon closes properly)
                x_wedge_sp_min = r_circle_sp_min * np.cos(angles_sp[::-1])
                y_wedge_sp_min = r_circle_sp_min * np.sin(angles_sp[::-1])
        
                # Combine into a single polygon
                x_wedge_sp = np.concatenate([x_wedge_sp_min, x_wedge_sp_max])
                y_wedge_sp = np.concatenate([y_wedge_sp_min, y_wedge_sp_max])
                ax.fill(x_wedge_sp, y_wedge_sp, facecolor='None', edgecolor=opaque_color_from_hex(mcolors.to_hex("black"), alpha=0.2), hatch='//')
                
                # Phugoid
                ax.plot(x_ray_ph_max, y_ray_ph_max, '-', color='black')
                ax.plot(x_ray_ph_max, -y_ray_ph_max, '-', color='black')
                
                x_stacked_ray_ph_max = np.hstack((np.flip(x_ray_ph_max), x_ray_ph_max))
                y_stacked_ray_ph_max = np.hstack((np.flip(y_ray_ph_max), -y_ray_ph_max))
                ax.fill_betweenx(
                    y=y_stacked_ray_ph_max, x1=np.ones_like(x_stacked_ray_ph_max) * -3, x2=x_stacked_ray_ph_max,
                    facecolor=opaque_color_from_hex(mcolors.to_hex("black"), alpha=0.1), edgecolor='None', zorder=-1,
                )
                axins.fill_betweenx(
                    y=y_stacked_ray_ph_max, x1=np.ones_like(x_stacked_ray_ph_max) * -3, x2=x_stacked_ray_ph_max,
                    facecolor=opaque_color_from_hex(mcolors.to_hex("black"), alpha=0.1), edgecolor='None', zorder=-1,
                )
         
            elif modes_type == 'lateral':
                
                alpha_constr = 0.1
            
                # Roll subsidence
                # ax.axvline(x_vert_roll_max, color='black', alpha=alpha_constr)
                y_stacked_ray_roll_max = [-y_lim_lat, y_lim_lat]
                ax.fill_betweenx(
                    # y=y_stacked_ray_roll_max, x1=np.ones_like(x_stacked_ray_roll_max) * ax.get_xlim()[0], x2=x_stacked_ray_roll_max,
                    y=y_stacked_ray_roll_max, x1=-y_lim_lat, x2=x_vert_roll_max,
                    # facecolor='lightgrey', edgecolor='None', alpha=0.5, zorder=-1,
                    facecolor='None', edgecolor='black', hatch='..', alpha=0.1, zorder=-1,
                )
                
                # Spiral
                # ax.axvline(x_vert_sl_max, color='black', alpha=alpha_constr)
                y_stacked_ray_sl_max = [-y_lim_lat, y_lim_lat]
                ax.fill_betweenx(
                    # y=y_stacked_ray_sl_max, x1=np.ones_like(x_stacked_ray_sl_max) * ax.get_xlim()[0], x2=x_stacked_ray_sl_max,
                    y=y_stacked_ray_sl_max, x1=-y_lim_lat, x2=x_vert_sl_max,
                    facecolor='lightgrey', edgecolor='None', alpha=0.5, zorder=-1,
                )
                
                # Dutch roll
                # ax.plot(x_circle_dr_min, y_circle_dr_min, color='black', alpha=alpha_constr)
                # ax.axvline(x_vert_dr_max, color='black', alpha=alpha_constr)
                # ax.plot(x_ray_dr_max, y_ray_dr_max, '-', color='black', alpha=alpha_constr)
                # ax.plot(x_ray_dr_max, -y_ray_dr_max, '-', color='black', alpha=alpha_constr)
                
                # Angular range for the wedge
                theta_dr_min = np.arccos(zeta_dr_min)
                angles_dr = np.linspace(np.pi - theta_dr_min, np.pi + theta_dr_min, 200)
        
                # Outer arc
                x_wedge_dr_max = r_circle_dr_min * np.cos(angles_dr)
                y_wedge_dr_max = r_circle_dr_min * np.sin(angles_dr)
        
                x_dr_ray_vert_intersect = -zeta_dr_min * r_circle_dr_min
                y_dr_ray_vert_intersect = r_circle_dr_min * np.sqrt(1 - zeta_dr_min**2)
        
                # Combine into a single polygon
                # x_wedge_dr = np.concatenate([x_wedge_dr_max, [x_dr_ray_vert_intersect, x_dr_ray_vert_intersect]])
                # y_wedge_dr = np.concatenate([y_wedge_dr_max, [-y_dr_ray_vert_intersect, y_dr_ray_vert_intersect]])
                # ax.fill(x_wedge_dr, y_wedge_dr, facecolor='None', edgecolor='black', hatch='////', alpha=0.1, zorder=-1)
                x_wedge_dr = np.concatenate([
                    x_wedge_dr_max, [x_ray_dr_max[-1], -y_lim_lat, -y_lim_lat, x_ray_dr_max[-1]],
                ])
                y_wedge_dr = np.concatenate([
                    y_wedge_dr_max, [-y_ray_dr_max[-1], -y_ray_dr_max[-1], y_ray_dr_max[-1], y_ray_dr_max[-1]],
                ])
                ax.fill(x_wedge_dr, y_wedge_dr, facecolor='None', edgecolor='black', hatch='//', alpha=0.1, zorder=-1)
            
            """
            # # Plot lateral roots
            # if modes_type == 'lateral':
            #     ax.scatter(ev_1_list.real, ev_1_list.imag, marker='o', color=colors[j], edgecolor=edgecolor, linewidth=linewidth, zorder=100)
            #     ax_left.scatter(ev_1_list.real, ev_1_list.imag, marker='o', color=colors[j], edgecolor=edgecolor, linewidth=linewidth, zorder=100)
            
            #     ax.scatter(ev_4_list.real, ev_4_list.imag, marker='^', color=colors[j], edgecolor=edgecolor, linewidth=linewidth, zorder=100)
            #     ax_left.scatter(ev_4_list.real, ev_4_list.imag, marker='^', color=colors[j], edgecolor=edgecolor, linewidth=linewidth, zorder=100)
            
            # # Plot longitudinal roots
            # elif modes_type == 'longitudinal':
            #     ax.scatter(ev_5_list.real, ev_5_list.imag, marker='*', color=colors[j], edgecolor=edgecolor, linewidth=linewidth, zorder=100)
            #     axins.scatter(ev_5_list.real, ev_5_list.imag, marker='*', color=colors[j], edgecolor=edgecolor, linewidth=linewidth, zorder=100)
                
            # Disentangle short-period and dutch rolle modes (substantial cross-coupling)
                
            def classify_modes(eigs_pos, labels_guess):
                
                z = eigs_pos
                lab = labels_guess
            
                # Compute "cloud centres"
                c0 = z[lab == 0].mean()
                c1 = z[lab == 1].mean()
            
                # Assign by nearest centre
                d0 = np.abs(z - c0)
                d1 = np.abs(z - c1)
            
                return (d1 < d0).astype(int)
    
            # Start by assuming that ev_2_list and ev_3_list are correctly separated
            eigs = np.concatenate([np.array(ev_2_list), np.array(ev_3_list)])
            labels_guess = np.concatenate([np.ones(len(ev_2_list)), np.zeros(len(ev_3_list))])
            labels_fixed = classify_modes(eigs, labels_guess)
            summed_abs_diff = sum(np.abs(labels_fixed - labels_guess))
            
            # sys.exit('Stop here.')
            
            # Only reorder if above assumption turns out to be incorrect
            if summed_abs_diff > 0:
                
                # Original group sizes
                n2 = len(ev_2_list)
                
                # How many from each original list went into each label?
                count_2_label1 = np.sum(labels_fixed[:n2] == 1)
                count_2_label0 = np.sum(labels_fixed[:n2] == 0)
                count_3_label0 = np.sum(labels_fixed[n2:] == 0)
                count_3_label1 = np.sum(labels_fixed[n2:] == 1)
                
                # Decide mapping by majority vote
                if count_2_label0 + count_3_label1 <= count_2_label1 + count_3_label0:
                    ev_2_list_sorted = eigs[labels_fixed == 1]
                    ev_3_list_sorted = eigs[labels_fixed == 0]
                else:
                    ev_2_list_sorted = eigs[labels_fixed == 0]
                    ev_3_list_sorted = eigs[labels_fixed == 1]
                
            else:
                ev_2_list_sorted = np.array(ev_2_list)
                ev_3_list_sorted = np.array(ev_3_list)
                
            # Once disentangled, the two clouds may still be swapped. Note that the imaginary
            # part corresponding to the short pediod should always be greater than that
            # corresponding to the dutch roll (flip below if not the case)
            if np.average(np.abs(ev_2_list_sorted.imag)) > np.average(np.abs(ev_3_list_sorted.imag)):
                _ev_2_list_sorted = copy.deepcopy(ev_2_list_sorted)
                _ev_3_list_sorted = copy.deepcopy(ev_3_list_sorted)
                ev_2_list_sorted = copy.deepcopy(_ev_3_list_sorted)
                ev_3_list_sorted = copy.deepcopy(_ev_2_list_sorted)
            
            # # # Blue markers (nacelle) are never misordered
            # # if j == 0:
            # #     if modes_type == 'lateral':
            # #         ax.scatter(np.array(ev_2_list).real, np.array(ev_2_list).imag, marker='s', color=colors[j], edgecolor=edgecolor, linewidth=linewidth)
            # #         ax.scatter(np.array(ev_2_list).real, -np.array(ev_2_list).imag, marker='s', color=colors[j], edgecolor=edgecolor, linewidth=linewidth)
            # #     elif modes_type == 'longitudinal':
            # #         ax.scatter(np.array(ev_3_list).real, np.array(ev_3_list).imag, marker='d', color=colors[j], edgecolor=edgecolor, linewidth=linewidth)
            # #         ax.scatter(np.array(ev_3_list).real, -np.array(ev_3_list).imag, marker='d', color=colors[j], edgecolor=edgecolor, linewidth=linewidth)
            # # elif (j == 1 or j == 2):
            # if modes_type == 'lateral':
            #     ax.scatter(ev_2_list_sorted.real, ev_2_list_sorted.imag, marker='s', color=colors[j], edgecolor=edgecolor, linewidth=linewidth)
            #     ax.scatter(ev_2_list_sorted.real, -ev_2_list_sorted.imag, marker='s', color=colors[j], edgecolor=edgecolor, linewidth=linewidth)
            # elif modes_type == 'longitudinal':
            #     ax.scatter(ev_3_list_sorted.real, ev_3_list_sorted.imag, marker='d', color=colors[j], edgecolor=edgecolor, linewidth=linewidth)
            #     ax.scatter(ev_3_list_sorted.real, -ev_3_list_sorted.imag, marker='d', color=colors[j], edgecolor=edgecolor, linewidth=linewidth)
            """
            
    ax.set_title(titles[row_idx], pad=10)
    ax.spines[['right','top']].set_visible(False)
    ax.tick_params(axis='y', which='both', right=False, length=0)
    ax.tick_params(axis='x', which='both', length=0)
    ax.spines['bottom'].set_position(('data', 0.0))
    ax.spines['left'].set_position(('data', 0.0))
    

plt.show()


