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

# import os
# import re
# import glob
# import numpy as np
# import pandas as pd

# # -------------------------
# # utils: build the dataset
# # -------------------------
# def build_dyn_dataset(root_dir):
#     """
#     Scan root_dir for files named
#       suave_dynamic_stability_outputs_<x>_<y>.txt
#     Read each file, annotate with x (int), file_y (int) and row (int),
#     and return a single DataFrame with a MultiIndex:
#     (x, sigma_fcs, span_loc, fcs_loc, wing_frac, nacelle_frac, row)
#     """
#     pattern = os.path.join(root_dir, "suave_dynamic_stability_outputs_*.txt")
#     files = sorted(glob.glob(pattern))
#     if not files:
#         raise FileNotFoundError(f"No matching files in {root_dir}")

#     rows = []
#     filename_re = re.compile(r"suave_dynamic_stability_outputs_(\d+)_(\d+)\.txt$")
#     for fn in files:
#         m = filename_re.search(os.path.basename(fn))
#         if not m:
#             continue
#         x = int(m.group(1))
#         y = int(m.group(2))
#         # read as whitespace-delimited table
#         df = pd.read_csv(fn, sep=r'\s+', engine='python', header=0)
#         # annotate meta columns
#         df = df.copy()
#         df['x'] = x
#         df['file_y'] = y
#         df['row'] = np.arange(len(df))  # row number within this file
#         rows.append(df)

#     df_all = pd.concat(rows, ignore_index=True, sort=False)

#     # ensure numeric columns have numeric dtype
#     numeric_cols = [c for c in df_all.columns if c not in ('x','file_y','row')]
#     df_all[numeric_cols] = df_all[numeric_cols].apply(pd.to_numeric, errors='coerce')

#     # create the MultiIndex requested:
#     key_cols = ['x', 'sigma_fcs', 'span_loc', 'fcs_loc', 'wing_frac', 'nacelle_frac', 'row']
#     df_all.set_index(key_cols, inplace=True)

#     return df_all


# def _expand_key(x, first5, row):
#     """
#     Replace None with slice(None) for MultiIndex lookup.
#     """
#     xk = slice(None) if x is None else x
#     f5 = tuple(slice(None) if v is None else v for v in first5)
#     rk = slice(None) if row is None else row
#     return (xk, *f5, rk)


# # -------------------------
# # small selection helpers
# # -------------------------
# RE_COLS = [f"Re{i}" for i in range(1,9)]
# IM_COLS = [f"Im{i}" for i in range(1,9)]

# def get_eigen_parts(df_all, x, first5, row, part='Re', idx=None, exact=True, tol=1e-9):
#     assert part in ('Re','Im'), "part must be 'Re' or 'Im'"
#     cols = RE_COLS if part == 'Re' else IM_COLS

#     key = _expand_key(x, first5, row)

#     # ---- wildcard path (fast, native MultiIndex slicing) ----
#     if x is None or None in first5 or row is None:
#         out = df_all.loc[key, cols]
#         vals = out.to_numpy(dtype=float)
#         if idx is None:
#             return vals
#         else:
#             assert 1 <= int(idx) <= 8
#             return vals[:, int(idx)-1]

#     # ---- original exact / fuzzy logic (unchanged) ----
#     try:
#         ser = df_all.loc[key]
#     except KeyError:
#         if exact:
#             raise

#         x_mask = (df_all.index.get_level_values('x') == x)
#         sigma = df_all.index.get_level_values('sigma_fcs').astype(float)
#         span  = df_all.index.get_level_values('span_loc').astype(float)
#         fcs   = df_all.index.get_level_values('fcs_loc').astype(float)
#         wing  = df_all.index.get_level_values('wing_frac').astype(float)
#         nac   = df_all.index.get_level_values('nacelle_frac').astype(float)
#         row_vals = df_all.index.get_level_values('row').astype(int)

#         cond = (
#             x_mask &
#             (np.abs(sigma - first5[0]) <= tol) &
#             (np.abs(span  - first5[1]) <= tol) &
#             (np.abs(fcs   - first5[2]) <= tol) &
#             (np.abs(wing  - first5[3]) <= tol) &
#             (np.abs(nac   - first5[4]) <= tol) &
#             (row_vals == row)
#         )
#         matches = df_all[cond]
#         if len(matches) == 0:
#             raise KeyError("No matching row found (even fuzzy).")
#         ser = matches.iloc[0]

#     vals = ser[cols].to_numpy(dtype=float)
#     return vals if idx is None else float(vals[int(idx)-1])

    
    
# def get_flight_condition(df_all, x, first5, row, exact=True, tol=1e-9):
#     key = _expand_key(x, first5, row)

#     # ---- wildcard path ----
#     if x is None or None in first5 or row is None:
#         cols = ['AoA','Mach','Beta','h','n','W']
#         out = df_all.loc[key, cols]
#         return out.to_numpy(dtype=float)

#     # ---- original exact / fuzzy logic ----
#     try:
#         ser = df_all.loc[key]
#     except KeyError:
#         if exact:
#             raise

#         x_mask = (df_all.index.get_level_values('x') == x)
#         sigma = df_all.index.get_level_values('sigma_fcs').astype(float)
#         span  = df_all.index.get_level_values('span_loc').astype(float)
#         fcs   = df_all.index.get_level_values('fcs_loc').astype(float)
#         wing  = df_all.index.get_level_values('wing_frac').astype(float)
#         nac   = df_all.index.get_level_values('nacelle_frac').astype(float)
#         row_vals = df_all.index.get_level_values('row').astype(int)

#         cond = (
#             x_mask &
#             (np.abs(sigma - first5[0]) <= tol) &
#             (np.abs(span  - first5[1]) <= tol) &
#             (np.abs(fcs   - first5[2]) <= tol) &
#             (np.abs(wing  - first5[3]) <= tol) &
#             (np.abs(nac   - first5[4]) <= tol) &
#             (row_vals == row)
#         )
#         matches = df_all[cond]
#         if len(matches) == 0:
#             raise KeyError("No matching row found (even fuzzy).")
#         ser = matches.iloc[0]

#     return {
#         'AoA':  float(ser['AoA']),
#         'Mach': float(ser['Mach']),
#         'Beta': float(ser['Beta']),
#         'h':    float(ser['h']),
#         'n':    float(ser['n']),
#         'W':    float(ser['W']),
#     }


# def get_mass_parameters(df_all, x, first5, row, exact=True, tol=1e-9):
#     """
#     Retrieve mass magnitude/location parameters:
#     (sigma_fcs, span_loc, fcs_loc, wing_frac, nacelle_frac)

#     Supports wildcards (None) for x, first5 entries, and row.
#     Returns an (N,5) numpy array.
#     """

#     # build boolean mask explicitly (robust to partial slicing)
#     mask = np.ones(len(df_all), dtype=bool)

#     if x is not None:
#         mask &= (df_all.index.get_level_values('x') == x)

#     names = ['sigma_fcs','span_loc','fcs_loc','wing_frac','nacelle_frac']
#     for name, val in zip(names, first5):
#         if val is not None:
#             vals = df_all.index.get_level_values(name).astype(float)
#             mask &= (np.abs(vals - val) <= tol)

#     if row is not None:
#         mask &= (df_all.index.get_level_values('row') == row)

#     sel = df_all[mask]
#     if len(sel) == 0:
#         raise KeyError("No matching mass-parameter rows found.")

#     idx = sel.index

#     return np.column_stack([
#         idx.get_level_values('sigma_fcs').astype(float),
#         idx.get_level_values('span_loc').astype(float),
#         idx.get_level_values('fcs_loc').astype(float),
#         idx.get_level_values('wing_frac').astype(float),
#         idx.get_level_values('nacelle_frac').astype(float),
#     ])

# #%%

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
    return -np.log(2) / lamda.real  # page 127 in AE3202 Flight Dynamics Lecture Notes

def T2(lamda):  # (s)
    return -T0p5(lamda)  # page 127 in AE3202 Flight Dynamics Lecture Notes

def omega0(lamda):  # (1/s)
    return np.sqrt(lamda.real**2 + lamda.imag**2)  # page 129 in AE3202 Flight Dynamics Lecture Notes

def CAP(lamda, n, alpha):  # (1/s^2)
    return omega0(lamda)**2 / (n / alpha)  # (3.111) in medium_PhDthesis_2013_flightdynconstrconcacmdiscanalsopt_morris

def tau(lamda):  # (s)
    return T0p5(lamda) / np.log(2)  # page 125 in AE3202 Flight Dynamics Lecture Notes

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

    N_points = 20
    
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
    
    x_vert_dr_max = -zeta_dr_min * r_circle_dr_min
    
    a_dr = zeta_dr_min
    t_dr = np.linspace(0, r_circle_dr_min / a_dr, N_points)
    x_ray_dr_max = -t_dr * a_dr
    y_ray_dr_max = t_dr * np.sqrt(1 - a_dr**2)
    
    # Spiral
    x_vert_sl_max = -np.log(0.5) / T_2_sl_min
    
    # Roll
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
    t_ph = np.linspace(0, 10 / np.sqrt(1 - a_ph**2), N_points)
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
    mode,
    df_mode,
):
    
    ev_array = df_mode['Re'].to_numpy() + 1j * df_mode['Im'].to_numpy()
    sigma_fcs_array = df_mode['sigma_fcs'].to_numpy()
    fcs_loc_array = df_mode['fcs_loc'].to_numpy()
    span_loc_array = df_mode['span_loc'].to_numpy()
    aoa = np.deg2rad(df_mode['AoA'].to_numpy())
    mach = df_mode['Mach'].to_numpy()
    beta = np.deg2rad(df_mode['Beta'].to_numpy())
    h = df_mode['h'].to_numpy()
    n = df_mode['n'].to_numpy()
    W = df_mode['W'].to_numpy()
    
    # Longitudinal stability limits
    
    # Short period
    if mode == 'short-period':
        zeta_sp_min = 0.3  # 1st row, 4th column in Table 3.10 (Class II, Category B, Level I)
        zeta_sp_max = 2.0  # 1st row, 5th column in Table 3.10 (Class II, Category B, Level I)
        CAP_sp_min = 0.085  # 1st row, 4th column in Table 3.11 (Class II, Category B, Level I)
        CAP_sp_max = 3.6  # 1st row, 5th column in Table 3.11 (Class II, Category B, Level I)
        zeta_sp = zeta(ev_array)
        CAP_sp = CAP(ev_array, n, aoa)
        valid_mask = (
            (zeta_sp > zeta_sp_min) & (zeta_sp < zeta_sp_max) &
            (CAP_sp > CAP_sp_min) & (CAP_sp < CAP_sp_max)
        )
        
    # Phugoid
    elif mode == 'phugoid':
        zeta_ph_min = 0.04  # 1st row in Table 3.9 (Class II, Category B, Level I)
        zeta_ph = zeta(ev_array)
        valid_mask = (zeta_ph > zeta_ph_min)
        
    # Lateral stability limits
    
    # Roll
    elif mode == 'roll':
        tau_roll_max = 1.4  # 3rd row of Table 3.12 (Class II, Category B, Level I)
        tau_roll = tau(ev_array)
        valid_mask = (tau_roll < tau_roll_max)
        
    # Spiral
    elif mode == 'spiral':
        T_2_sl_min = 20  # 2nd row of Table 3.13 (Category B, Level I)
        T_2_sl = T2(ev_array)
        valid_mask = (T_2_sl > T_2_sl_min) | (T_2_sl < 0)  # only applies to positive real parts
        
    # Dutch roll
    elif mode == 'dutch-roll':
        zeta_dr_min = 0.08  # 3rd row of Table 3.14 (Class II, Category B, Level I)
        omega_n_dr_min = 0.4  # 3rd row of Table 3.14 (Class II, Category B, Level I)
        zeta_dr = zeta(ev_array)
        omega_n_dr = omega0(ev_array)
        valid_mask = (zeta_dr > zeta_dr_min) & (omega_n_dr > omega_n_dr_min)
        
    ev_valid = ev_array[valid_mask]
    ev_invalid = ev_array[~valid_mask]
    sigma_fcs_valid = sigma_fcs_array[valid_mask]
    sigma_fcs_invalid = sigma_fcs_array[~valid_mask]
    fcs_loc_valid = fcs_loc_array[valid_mask]
    fcs_loc_invalid = fcs_loc_array[~valid_mask]
    span_loc_valid = span_loc_array[valid_mask]
    span_loc_invalid = span_loc_array[~valid_mask]
        
    return (
        ev_valid, ev_invalid,
        sigma_fcs_valid, sigma_fcs_invalid,
        fcs_loc_valid, fcs_loc_invalid,
        span_loc_valid, span_loc_invalid,
    )

#%%

from eigenvalue_classification_physics_based import match_evals_to_emodes

matfile = r"C:/Users/nmb48/suave_dynamic_stability_matrix_combined.txt"
eigfile = r"C:/Users/nmb48/suave_dynamic_stability_outputs_combined.txt"
res = match_evals_to_emodes(matfile, eigfile)

#%%

modes_type = 'longitudinal'  # 'longitudinal', 'lateral'

edgecolor = 'k'
linewidth = 0.5
alpha_invalid = 0.2
alpha_valid = 1.0

fig = plt.figure(figsize=(15,9))
gs = gridspec.GridSpec(
    2, 3,
    figure=fig,
    left=0.05,  # for longitudinal
    # left=0.1125,  # for lateral
    # right=0.8,  # for longitudinal
    right=0.775,  # for lateral
    bottom=0.1125,
    top=0.9,
    hspace=0.3,
    wspace=0.2,  # for longitudinal
    # wspace=0.4,  # for lateral
)
axes = [fig.add_subplot(gs[_i,_j]) for _i in range(2) for _j in range(3)]

titles = ['1. Take-off', '2. Climb', '3. Beginning of cruise', '4. End of cruise', '5. Descent', '6. Landing']

# =============================================================================
# init_dict = {
#     cat: {
#         mode: {rid: None for rid in [1, 2, 3, 4, 5, 6]}
#         for mode in ['short-period', 'phugoid', 'roll', 'dutch-roll', 'spiral']
#     }
#     for cat in ['fuselage','wing','nacelle']
# }
# =============================================================================

init_dict = {
    cat: {
        mode: {
            rid: {
                var: None for var in ['ev', 'sigma_fcs', 'fcs_loc', 'span_loc']
            } 
            for rid in [1, 2, 3, 4, 5, 6]
        }
        for mode in ['short-period', 'phugoid', 'roll', 'dutch-roll', 'spiral']
    }
    for cat in ['fuselage','wing','nacelle']
}
ev_valid_dict = copy.deepcopy(init_dict)
ev_invalid_dict = copy.deepcopy(init_dict)
ev_valid_sp_dict = copy.deepcopy(init_dict)
ev_invalid_sp_dict = copy.deepcopy(init_dict)
ev_valid_ph_dict = copy.deepcopy(init_dict)
ev_invalid_ph_dict = copy.deepcopy(init_dict)
ev_valid_roll_dict = copy.deepcopy(init_dict)
ev_invalid_roll_dict = copy.deepcopy(init_dict)
ev_valid_sl_dict = copy.deepcopy(init_dict)
ev_invalid_sl_dict = copy.deepcopy(init_dict)
ev_valid_dr_dict = copy.deepcopy(init_dict)
ev_invalid_dr_dict = copy.deepcopy(init_dict)

for row_idx, ax in enumerate(axes):
    rid = row_idx + 1
    
    if modes_type == 'longitudinal':
        modes_sublist = ['short-period', 'phugoid']
        markers = ['o', 'd']
        
        # Axis settings
        ax.set_xlim((-1.5, 1))
        ax.set_ylim((-0.1, 5))
        ax.set_xticks([-1, -0.5])
        # ax.set_yticks([-5, -4, -3, -2, -1, 1, 2, 3, 4, 5])
        ax.set_yticks([1, 2, 3, 4, 5])
        
        # Inset axis and settings
        axins = inset_axes(
            ax,
            width="100%",
            height="100%",
            bbox_to_anchor=(0.65, 0.65, 0.35, 0.35),  # (x0, y0, w, h) in ax coords
            bbox_transform=ax.transAxes,
            borderpad=0
        )
        axins.axvline(0.0, color='k')
        axins.axhline(0.0, color='k')
        axins.tick_params(labelleft=False, labelbottom=False)
        axins.tick_params(axis='y', which='both', length=0)
        axins.tick_params(axis='x', which='both', length=0)
        mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")
        
    elif modes_type == 'lateral':
        modes_sublist = ['roll', 'dutch-roll', 'spiral']
        markers = ['o', 'd', '^']
        
        y_lim_lat = 3
        
        ax.set_xlim((-1.35, 0.1))
        ax.set_ylim((-0.1, y_lim_lat))
        ax.set_xticks([-1, -0.5])
        ax.set_yticks([1, 2])
        
        # Add second x-axis left of break to plot roll subsidence roots    
        gap = 0.05  # fraction of ax width
        axleft = ax.inset_axes(
            [-0.25 - gap, 0.0, 0.2, 1.0],  # <-- shifted left
            sharey=ax
        )
        axleft.set_xlim(-65, -4)
        axleft.set_xticks([-65, -5])
        axleft.set_ylim(-0.1, y_lim_lat)
        axleft.spines[['left', 'right','top']].set_visible(False)
        axleft.tick_params(axis='y', which='both', right=False, length=0)
        axleft.tick_params(axis='x', which='both', length=0)
        axleft.spines['bottom'].set_position(('data', 0.0))
        axleft.tick_params(labelright=False)
        axleft.tick_params(labelleft=False)
    
    x_circle_dr_min_list = []
    y_circle_dr_min_list = []
    x_vert_dr_max_list = []
    x_ray_dr_max_list = []
    y_ray_dr_max_list = []
    r_circle_dr_min_list = []
    zeta_dr_min_list = []
    x_vert_sl_max_list = []
    x_vert_roll_max_list = []
    
    x_circle_sp_min_list = []
    y_circle_sp_min_list = []
    x_circle_sp_max_list = []
    y_circle_sp_max_list = []
    x_ray_sp_max_list = []
    y_ray_sp_max_list = []
    r_circle_sp_min_list = []
    r_circle_sp_max_list = []
    zeta_sp_min_list = []
    x_ray_ph_max_list = []
    y_ray_ph_max_list = []
    
    for j, cat in enumerate(['fuselage','wing','nacelle']):
        
        for k, mode in enumerate(modes_sublist):
            
            df_cat = res[cat][mode][rid]
            (
            ev_valid, ev_invalid,
            sigma_fcs_valid, sigma_fcs_invalid,
            fcs_loc_valid, fcs_loc_invalid,
            span_loc_valid, span_loc_invalid,
            ) = \
            mask_mil_limits(
                mode,
                df_cat,
            )
            
            ev_valid_dict[cat][mode][rid]['ev'] = ev_valid
            ev_invalid_dict[cat][mode][rid]['ev'] = ev_invalid
            ev_valid_dict[cat][mode][rid]['sigma_fcs'] = sigma_fcs_valid
            ev_invalid_dict[cat][mode][rid]['sigma_fcs'] = sigma_fcs_invalid
            ev_valid_dict[cat][mode][rid]['fcs_loc'] = fcs_loc_valid
            ev_invalid_dict[cat][mode][rid]['fcs_loc'] = fcs_loc_invalid
            ev_valid_dict[cat][mode][rid]['span_loc'] = span_loc_valid
            ev_invalid_dict[cat][mode][rid]['span_loc'] = span_loc_invalid
            
            ax.scatter(ev_valid.real, ev_valid.imag, marker=markers[k], facecolor=opaque_color_from_hex(colors[j], alpha=alpha_valid), edgecolor=opaque_color_from_hex(mcolors.to_hex(edgecolor), alpha=alpha_valid), linewidth=linewidth, zorder=100, clip_on=True)
            ax.scatter(ev_invalid.real, ev_invalid.imag, marker=markers[k], facecolor=opaque_color_from_hex(colors[j], alpha=alpha_invalid), edgecolor=opaque_color_from_hex(mcolors.to_hex(edgecolor), alpha=alpha_invalid), linewidth=linewidth, zorder=100, clip_on=True)
        
            if mode == 'short-period':
                ev_valid_sp_dict[cat][mode][rid]['ev'] = ev_valid
                ev_invalid_sp_dict[cat][mode][rid]['ev'] = ev_invalid
                # ev_valid_sp_dict[cat][mode][rid]['sigma_fcs'] = sigma_fcs_valid
                # ev_invalid_sp_dict[cat][mode][rid]['sigma_fcs'] = sigma_fcs_invalid
                # ev_valid_sp_dict[cat][mode][rid]['fcs_loc'] = fcs_loc_valid
                # ev_invalid_sp_dict[cat][mode][rid]['fcs_loc'] = fcs_loc_invalid
                # ev_valid_sp_dict[cat][mode][rid]['span_loc'] = span_loc_valid
                # ev_invalid_sp_dict[cat][mode][rid]['span_loc'] = span_loc_invalid
                # pass
            elif mode == 'phugoid':
                axins.scatter(ev_valid.real, ev_valid.imag, marker=markers[k], facecolor=opaque_color_from_hex(colors[j], alpha=alpha_valid), edgecolor=opaque_color_from_hex(mcolors.to_hex(edgecolor), alpha=alpha_valid), linewidth=linewidth, zorder=100, clip_on=False)
                axins.scatter(ev_invalid.real, ev_invalid.imag, marker=markers[k], facecolor=opaque_color_from_hex(colors[j], alpha=alpha_invalid), edgecolor=opaque_color_from_hex(mcolors.to_hex(edgecolor), alpha=alpha_invalid), linewidth=linewidth, zorder=100, clip_on=False)
                ev_valid_ph_dict[cat][mode][rid]['ev'] = ev_valid
                ev_invalid_ph_dict[cat][mode][rid]['ev'] = ev_invalid
                pass
            elif mode == 'roll':
                axleft.scatter(ev_valid.real, ev_valid.imag, marker=markers[k], facecolor=opaque_color_from_hex(colors[j], alpha=alpha_valid), edgecolor=opaque_color_from_hex(mcolors.to_hex(edgecolor), alpha=alpha_valid), linewidth=linewidth, zorder=100, clip_on=False)
                axleft.scatter(ev_invalid.real, ev_invalid.imag, marker=markers[k], facecolor=opaque_color_from_hex(colors[j], alpha=alpha_invalid), edgecolor=opaque_color_from_hex(mcolors.to_hex(edgecolor), alpha=alpha_invalid), linewidth=linewidth, zorder=100, clip_on=False)
                # ev_valid_roll_dict[cat][mode][rid]['ev'] = ev_valid
                # ev_invalid_roll_dict[cat][mode][rid]['ev'] = ev_invalid
                pass
            elif mode == 'spiral':
                # ev_valid_sl_dict[cat][mode][rid]['ev'] = ev_valid
                # ev_invalid_sl_dict[cat][mode][rid]['ev'] = ev_invalid
                pass
            elif mode == 'dutch-roll':
                # ev_valid_dr_dict[cat][mode][rid]['ev'] = ev_valid
                # ev_invalid_dr_dict[cat][mode][rid]['ev'] = ev_invalid
                pass
            
            aoa_array = df_cat['AoA'].to_numpy()
            mach_array = df_cat['Mach'].to_numpy()
            beta_array = df_cat['Beta'].to_numpy()
            h_array = df_cat['h'].to_numpy()
            n_array = df_cat['n'].to_numpy()
            W_array = df_cat['W'].to_numpy()
            
            for l in range(len(df_cat)):
            
                # Calculate MIL-constraints from flight conditions
                aoa = aoa_array[l]
                mach = mach_array[l]
                beta = beta_array[l]
                h = h_array[l]
                n = n_array[l]
                W = W_array[l]
                (
                    x_circle_dr_min, y_circle_dr_min, x_vert_dr_max, x_ray_dr_max, y_ray_dr_max, r_circle_dr_min, zeta_dr_min,
                    x_vert_sl_max,
                    x_vert_roll_max,
                    
                    x_circle_sp_min, y_circle_sp_min, x_circle_sp_max, y_circle_sp_max, x_ray_sp_max, y_ray_sp_max, r_circle_sp_min, r_circle_sp_max, zeta_sp_min,
                    x_ray_ph_max, y_ray_ph_max,
                ) = \
                plot_mil_limits(aoa, mach, beta, h, n, W)
                
                x_circle_dr_min_list.append(x_circle_dr_min)
                y_circle_dr_min_list.append(y_circle_dr_min)
                x_vert_dr_max_list.append(x_vert_dr_max)
                x_ray_dr_max_list.append(x_ray_dr_max)
                y_ray_dr_max_list.append(y_ray_dr_max)
                r_circle_dr_min_list.append(r_circle_dr_min)
                zeta_dr_min_list.append(zeta_dr_min)
                x_vert_sl_max_list.append(x_vert_sl_max)
                x_vert_roll_max_list.append(x_vert_roll_max)
                
                x_circle_sp_min_list.append(x_circle_sp_min)
                y_circle_sp_min_list.append(y_circle_sp_min)
                x_circle_sp_max_list.append(x_circle_sp_max)
                y_circle_sp_max_list.append(y_circle_sp_max)
                x_ray_sp_max_list.append(x_ray_sp_max)
                y_ray_sp_max_list.append(y_ray_sp_max)
                r_circle_sp_min_list.append(r_circle_sp_min)
                r_circle_sp_max_list.append(r_circle_sp_max)
                zeta_sp_min_list.append(zeta_sp_min)
                x_ray_ph_max_list.append(x_ray_ph_max)
                y_ray_ph_max_list.append(y_ray_ph_max)
                
    # Plot MIL-constraints
    
    if modes_type == 'longitudinal':
        
        x_circle_dr_min_array = np.hstack(x_circle_dr_min_list)
        y_circle_dr_min_array = np.hstack(y_circle_dr_min_list)
        x_vert_dr_max_array = np.hstack(x_vert_dr_max_list)
        x_ray_dr_max_array = np.hstack(x_ray_dr_max_list)
        y_ray_dr_max_array = np.hstack(y_ray_dr_max_list)
        r_circle_dr_min_array = np.hstack(r_circle_dr_min_list)
        zeta_dr_min_array = np.hstack(zeta_dr_min_list)
        x_vert_sl_max_array = np.hstack(x_vert_sl_max_list)
        x_vert_roll_max_array = np.hstack(x_vert_roll_max_list)
        
        x_circle_sp_min_array = np.hstack(x_circle_sp_min_list)
        y_circle_sp_min_array = np.hstack(y_circle_sp_min_list)
        x_circle_sp_max_array = np.hstack(x_circle_sp_max_list)
        y_circle_sp_max_array = np.hstack(y_circle_sp_max_list)
        x_ray_sp_max_array = np.hstack(x_ray_sp_max_list)
        y_ray_sp_max_array = np.hstack(y_ray_sp_max_list)
        r_circle_sp_min_array = np.hstack(r_circle_sp_min_list)
        r_circle_sp_max_array = np.hstack(r_circle_sp_max_list)
        zeta_sp_min_array = np.hstack(zeta_sp_min_list)
        x_ray_ph_max_array = np.hstack(x_ray_ph_max_list)
        y_ray_ph_max_array = np.hstack(y_ray_ph_max_list)
        
        # Short period
        # ax.plot(x_circle_sp_min, y_circle_sp_min, color=opaque_color_from_hex(mcolors.to_hex("black"), alpha=0.2))
        # ax.plot(x_circle_sp_max, y_circle_sp_max, color=opaque_color_from_hex(mcolors.to_hex("black"), alpha=0.2))
        # ax.plot(x_ray_sp_max, y_ray_sp_max, '-', color=opaque_color_from_hex(mcolors.to_hex("black"), alpha=0.2))
        # ax.plot(x_ray_sp_max, -y_ray_sp_max, '-', color=opaque_color_from_hex(mcolors.to_hex("black"), alpha=0.2))
        
        # Angular range for the wedge
        theta_sp_min = np.arccos(zeta_sp_min)
        angles_sp = np.linspace(np.pi - theta_sp_min, np.pi + theta_sp_min, 20)

        # Outer arc
        x_wedge_sp_min_max = np.min(r_circle_sp_max_array) * np.cos(angles_sp)
        y_wedge_sp_min_max = np.min(r_circle_sp_max_array) * np.sin(angles_sp)
        x_wedge_sp_max_max = np.max(r_circle_sp_max_array) * np.cos(angles_sp)
        y_wedge_sp_max_max = np.max(r_circle_sp_max_array) * np.sin(angles_sp)

        # Inner arc (reverse direction so polygon closes properly)
        x_wedge_sp_min_min = np.min(r_circle_sp_min_array) * np.cos(angles_sp[::-1])
        y_wedge_sp_min_min = np.min(r_circle_sp_min_array) * np.sin(angles_sp[::-1])
        x_wedge_sp_max_min = np.max(r_circle_sp_min_array) * np.cos(angles_sp[::-1])
        y_wedge_sp_max_min = np.max(r_circle_sp_min_array) * np.sin(angles_sp[::-1])

        # Main wedge
        x_wedge_sp = np.concatenate([x_wedge_sp_min_max, x_wedge_sp_max_min])
        y_wedge_sp = np.concatenate([y_wedge_sp_min_max, y_wedge_sp_max_min])
        ax.fill(x_wedge_sp, y_wedge_sp, facecolor='None', edgecolor=opaque_color_from_hex(mcolors.to_hex("black"), alpha=0.2), hatch='//')
        
        # Inner wedge
        x_wedge_sp = np.concatenate([x_wedge_sp_max_max, np.flip(x_wedge_sp_min_max)])
        y_wedge_sp = np.concatenate([y_wedge_sp_max_max, np.flip(y_wedge_sp_min_max)])
        ax.fill(x_wedge_sp, y_wedge_sp, facecolor='none', edgecolor=opaque_color_from_hex(mcolors.to_hex("black"), alpha=0.2), hatch='//')
        
        # Outer wedge
        x_wedge_sp = np.concatenate([x_wedge_sp_max_min, np.flip(x_wedge_sp_min_min)])
        y_wedge_sp = np.concatenate([y_wedge_sp_max_min, np.flip(y_wedge_sp_min_min)])
        ax.fill(x_wedge_sp, y_wedge_sp, facecolor='none', edgecolor=opaque_color_from_hex(mcolors.to_hex("black"), alpha=0.2), hatch='//')
        
        # Phugoid
        # ax.plot(x_ray_ph_max_list[0], y_ray_ph_max_list[0], '-', color='black')
        # ax.plot(x_ray_ph_max_list[0], -y_ray_ph_max_list[0], '-', color='black')
        
        x_stacked_ray_ph_max = np.hstack((np.flip(x_ray_ph_max_list[0]), x_ray_ph_max_list[0]))
        y_stacked_ray_ph_max = np.hstack((np.flip(y_ray_ph_max_list[0]), -y_ray_ph_max_list[0]))
        ax.fill_betweenx(
            y=y_stacked_ray_ph_max, x1=np.ones_like(x_stacked_ray_ph_max) * -3, x2=x_stacked_ray_ph_max,
            facecolor=opaque_color_from_hex(mcolors.to_hex("black"), alpha=0.1), edgecolor='None', zorder=-1,
        )
 
    elif modes_type == 'lateral':
        
        # NOTE: the lateral MIL constraints are insensitive to flight condition!
        
        alpha_constr = 0.1
    
        # Roll subsidence
        # ax.axvline(x_vert_roll_max, color='black', alpha=alpha_constr)
        y_stacked_ray_roll_max = [-y_lim_lat, y_lim_lat]
        ax.fill_betweenx(
            y=y_stacked_ray_roll_max, x1=-y_lim_lat, x2=x_vert_roll_max,
            facecolor='None', edgecolor='black', hatch='..', alpha=0.1, zorder=-1,
        )
        axleft.fill_betweenx(
            y=y_stacked_ray_roll_max, x1=-65, x2=5,
            facecolor='None', edgecolor='black', hatch='..', alpha=0.1, zorder=-1,
        )
        
        # Spiral
        # ax.axvline(x_vert_sl_max, color='black', alpha=alpha_constr)
        y_stacked_ray_sl_max = [-y_lim_lat, y_lim_lat]
        ax.fill_betweenx(
            y=y_stacked_ray_sl_max, x1=-y_lim_lat, x2=x_vert_sl_max,
            facecolor='lightgrey', edgecolor='None', alpha=0.5, zorder=-1,
        )
        axleft.fill_betweenx(
            y=y_stacked_ray_sl_max, x1=-65, x2=5,
            facecolor='lightgrey', edgecolor='None', alpha=0.5, zorder=-1,
        )
        
        # Dutch roll
        # ax.plot(x_circle_dr_min, y_circle_dr_min, color='black', alpha=alpha_constr)
        # ax.axvline(x_vert_dr_max, color='black', alpha=alpha_constr)
        # ax.plot(x_ray_dr_max, y_ray_dr_max, '-', color='black', alpha=alpha_constr)
        # ax.plot(x_ray_dr_max, -y_ray_dr_max, '-', color='black', alpha=alpha_constr)
        
        # Angular range for the wedge
        theta_dr_min = np.arccos(zeta_dr_min)
        angles_dr = np.linspace(np.pi - theta_dr_min, np.pi + theta_dr_min, 20)

        # Outer arc
        x_wedge_dr_max = r_circle_dr_min * np.cos(angles_dr)
        y_wedge_dr_max = r_circle_dr_min * np.sin(angles_dr)

        x_dr_ray_vert_intersect = -zeta_dr_min * r_circle_dr_min
        y_dr_ray_vert_intersect = r_circle_dr_min * np.sqrt(1 - zeta_dr_min**2)

        # Combine into a single polygon
        x_wedge_dr = np.concatenate([
            x_wedge_dr_max, [x_ray_dr_max[-1], -y_lim_lat, -y_lim_lat, x_ray_dr_max[-1]],
        ])
        y_wedge_dr = np.concatenate([
            y_wedge_dr_max, [-y_ray_dr_max[-1], -y_ray_dr_max[-1], y_ray_dr_max[-1], y_ray_dr_max[-1]],
        ])
        ax.fill(x_wedge_dr, y_wedge_dr, facecolor='None', edgecolor='black', hatch='//', alpha=0.1, zorder=-1)
        x_wedge_dr = np.concatenate([
            x_wedge_dr_max, [x_ray_dr_max[-1], -65, -65, x_ray_dr_max[-1]],
        ])
        y_wedge_dr = np.concatenate([
            y_wedge_dr_max, [-y_ray_dr_max[-1], -y_ray_dr_max[-1], y_ray_dr_max[-1], y_ray_dr_max[-1]],
        ])
        axleft.fill(x_wedge_dr, y_wedge_dr, facecolor='None', edgecolor='black', hatch='//', alpha=0.1, zorder=-1)
            
    if modes_type == 'longitudinal':
        
        axins.fill_betweenx(
            y=y_stacked_ray_ph_max, x1=np.ones_like(x_stacked_ray_ph_max) * -3, x2=x_stacked_ray_ph_max,
            facecolor=opaque_color_from_hex(mcolors.to_hex("black"), alpha=0.1), edgecolor='None', zorder=-1,
        )
        
        arrays = [
            np.ravel(mode_dict[rid]['ev'])
            for cat_dict in ev_valid_ph_dict.values()
            for mode_dict in cat_dict.values()
            if mode_dict[rid]['ev'] is not None and np.asarray(mode_dict[rid]['ev']).size > 0
        ]
        ev_valid_phugoid_array = (
            np.concatenate(arrays)
            if len(arrays) > 0
            else np.empty(0, dtype=np.complex128)
        )
        
        arrays = [
            np.ravel(mode_dict[rid]['ev'])
            for cat_dict in ev_invalid_ph_dict.values()
            for mode_dict in cat_dict.values()
            if mode_dict[rid]['ev'] is not None and np.asarray(mode_dict[rid]['ev']).size > 0
        ]
        ev_invalid_phugoid_array = (
            np.concatenate(arrays)
            if len(arrays) > 0
            else np.empty(0, dtype=np.complex128)
        )
        
        ev_real_phugoid_min = np.min(np.hstack((ev_valid_phugoid_array.real, ev_invalid_phugoid_array.real)))
        ev_real_phugoid_max = np.max(np.hstack((ev_valid_phugoid_array.real, ev_invalid_phugoid_array.real)))
        ev_imag_phugoid_min = np.min(np.hstack((ev_valid_phugoid_array.imag, ev_invalid_phugoid_array.imag)))
        ev_imag_phugoid_max = np.max(np.hstack((ev_valid_phugoid_array.imag, ev_invalid_phugoid_array.imag)))
        
        axins.set_xlim(ev_real_phugoid_min, ev_real_phugoid_max)
        axins.set_ylim(ev_imag_phugoid_min, ev_imag_phugoid_max)
        add_margin(axins, m=0.5)
            
    ax.set_title(titles[row_idx], pad=20)
    ax.spines[['right','top']].set_visible(False)
    ax.tick_params(axis='y', which='both', right=False, length=0)
    ax.tick_params(axis='x', which='both', length=0)
    ax.spines['bottom'].set_position(('data', 0.0))
    ax.spines['left'].set_position(('data', 0.0))
    
# Figure legend

if modes_type == 'longitudinal':
    
    legend_ax = fig.add_axes([0.8, 0.0, 0.225, 1.0])
    legend_marker_colours = [colors[0], colors[1], colors[2], colors[0], colors[1], colors[2]]
    legend_markers = ['o', 'o', 'o', 's', 's', 's']
    legend_marker_labels = [
        'Short period nacelle', 'Short period wing', 'Short period fuselage',
        'Phugoid nacelle', 'Phugoid wing', 'Phugoid fuselage',
    ]
    legend_elements = [
        mpatches.Patch(
            facecolor='none',
            edgecolor='black',
            hatch='//',
            label='Short period space',
            alpha=0.1,
        ),
        mpatches.Patch(
            facecolor='black',
            edgecolor='none',
            label='Phugoid space',
            alpha=0.1,
        )
    ] + [
        Line2D(
            [], [], marker=legend_markers[_], color='none', markerfacecolor=legend_marker_colours[_],
            markeredgecolor=legend_marker_colours[_], markersize=8, linestyle='None',
            label=legend_marker_labels[_],
        ) for _ in range(len(legend_marker_colours))
    ]
    legend_ax.legend(
        handles=legend_elements, loc='center right', bbox_to_anchor=(1.0, 0.5), ncol=1, labelspacing=1.0, frameon=False,
    )
elif modes_type == 'lateral':
    
    legend_ax = fig.add_axes([0.8, 0.0, 0.225, 1.0])
    legend_marker_colours = [colors[0], colors[1], colors[2], colors[0], colors[1], colors[2], colors[0], colors[1], colors[2]]
    legend_markers = ['o', 'o', 'o', 's', 's', 's', 'd', 'd', 'd']
    legend_marker_labels = [
        'Roll subsidence nacelle', 'Roll subsidence wing', 'Roll subsidence fuselage',
        'Spiral nacelle', 'Spiral wing', 'Spiral fuselage',
        'Dutch roll nacelle', 'Dutch roll wing', 'Dutch roll fuselage',
    ]
    legend_elements = [
        mpatches.Patch(
            facecolor='black',
            edgecolor='none',
            label='Spiral space',
            alpha=0.1,
        ),
        mpatches.Patch(
            facecolor='none',
            edgecolor='black',
            hatch='..',
            label='Roll subsidence space',
            alpha=0.1,
        ),
        mpatches.Patch(
            facecolor='none',
            edgecolor='black',
            hatch='//',
            label='Dutch roll space',
            alpha=0.1,
        )
    ] + [
        Line2D(
            [], [], marker=legend_markers[_], color='none', markerfacecolor=legend_marker_colours[_],
            markeredgecolor=legend_marker_colours[_], markersize=8, linestyle='None',
            label=legend_marker_labels[_],
        ) for _ in range(len(legend_marker_colours))
    ]
    legend_ax.legend(
        handles=legend_elements, loc='center right', bbox_to_anchor=(1.0, 0.5), ncol=1, labelspacing=1.0, frameon=False,
    )

legend_ax.spines[['left', 'right', 'top', 'bottom']].set_visible(False)
legend_ax.tick_params(
    axis='both',
    which='both',
    right=False,
    labelright=False,
    left=False,
    labelleft=False,
    length=0
)
legend_ax.set_xticks([])
legend_ax.set_yticks([])

    
fig.text(0.45, 0.025, 'Real component (1/s)', ha='center', va='bottom')
fig.text(0.025, 0.5, 'Imaginary component (1/s)', ha='right', va='center', rotation='vertical')
# plt.savefig('longitudinal_roots.png', format='png', dpi=600)
# plt.savefig('lateral_roots.png', format='png', dpi=600)

plt.show()

#%%

for j, cat in enumerate(['fuselage','wing','nacelle']):
    
    for k, mode in enumerate(modes_sublist):
        
        fig, ax = plt.subplots()
    
        ev_invalid_sigma_fcs_list = []
        ev_invalid_y_list = []
    
        for row_idx in range(6):
            rid = row_idx + 1
            
            ev_invalid_sigma_fcs_list.append(ev_invalid_dict[cat][mode][rid]['sigma_fcs'])
            if cat == 'fuselage':
                ev_invalid_y_list.append(ev_invalid_dict[cat][mode][rid]['fcs_loc'])
            elif cat == 'wing':
                ev_invalid_y_list.append(ev_invalid_dict[cat][mode][rid]['span_loc'])
            elif cat == 'nacelle':
                ev_invalid_y_list.append(np.zeros_like(ev_invalid_dict[cat][mode][rid]['span_loc']))
    
        ax.scatter(
            np.hstack(ev_invalid_sigma_fcs_list),
            np.hstack(ev_invalid_y_list),
        )
        
        x_min_list = []
        for y_var in np.unique(np.hstack(ev_invalid_y_list)):
            x_min = np.min(
                np.hstack(ev_invalid_sigma_fcs_list)[np.hstack(ev_invalid_y_list) == y_var]
            )
            x_min_list.append(x_min)
        print(x_min_list)
        
        y_min_list = []
        for sigma_fcs in np.unique(np.hstack(ev_invalid_sigma_fcs_list)):
            y_min = np.min(
                np.hstack(ev_invalid_y_list)[np.hstack(ev_invalid_sigma_fcs_list) == sigma_fcs]
            )
            y_min_list.append(y_min)
        ax.plot(np.unique(np.hstack(ev_invalid_sigma_fcs_list)), y_min_list)
        print()
        
        ax.set_xlim(1.5e3, 4e3)
        
        if (cat == 'fuselage' and mode == 'phugoid'):
            plt.savefig('TEMP.svg', format='svg')

plt.show()











