"""
This script post-processes the results of the dynamic stability
analysis conducted using the TASOPT-SUAVE-AVL wrapper.

The read 
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
from matplotlib.ticker import FuncFormatter

from nils.matplotlib_custom_settings import *

#%% Import eigenvalue and A and B system matrix data and associate with 5 eigenmodes

from nils.avl.code.dynamic_stability.eigenvalue_classification import match_evals_to_emodes

# matfile = r"C:/Users/nmb48/suave_dynamic_stability_matrix_combined.txt"
# eigfile = r"C:/Users/nmb48/suave_dynamic_stability_outputs_combined.txt"
# matfile = r"C:/Users/nmb48/Documents/GitHub/SUAVE/nils/avl/data/suave_dynamic_stability_matrix_combined.txt"
# eigfile = r"C:/Users/nmb48/Documents/GitHub/SUAVE/nils/avl/data/suave_dynamic_stability_outputs_combined.txt"
#####
# matfile = r"C:/Users/nmb48/Documents/GitHub/SUAVE/suave_dynamic_stability_matrix_combined.txt"
# eigfile = r"C:/Users/nmb48/Documents/GitHub/SUAVE/suave_dynamic_stability_outputs_combined.txt"
#####
# NILS (04.10.2026): Figs. 4.9 and 4.10 for either aircraft with the loading scenarios of all mass distribution
# studies in one figure. The combined files are written here from the per-case SUAVE output files (as
# eigenvalue_loading.py and system_matrix_loading.py do), so the glob patterns no longer need editing by hand.
# Usage: python -m nils.avl.code.dynamic_stability.process_dynamic_stability [aircraft] [modes_type]
from nils.avl.code.dynamic_stability.eigenvalue_loading import combine_eigen_outputs
from nils.avl.code.dynamic_stability.system_matrix_loading import combine_system_matrices
from nils.avl.code.stability_zones import MIL_LIMITS, flight_phase_category, COLOR_OK, COLOR_FAIL  # NILS (04.10.2026): single source of the criteria and colours (Figs. 4.8-4.11)

#####
# NILS: used to generate figures in Chapter 7 of thesis
SUAVE_ROOT = os.environ.get("SUAVE_ROOT", r"C:\Users\nmb48\Documents\GitHub\SUAVE\nils\avl\data\ch7")
#####
aircraft = sys.argv[1] if len(sys.argv) > 1 and sys.argv[1] else 'ATR_72-600'  # 'ATR_72-600' or 'Airbus_A220-100'
#####
# aircraft = 'ATR_72-600'
aircraft = 'Airbus_A220-100'  # 'ATR_72-600' or 'Airbus_A220-100'
#####
study_list = [1, 2, 5, 6, 7]  # studies 3 and 4 are identical to 7 and 6
save_fig = False  # False

eig_files = [f for s in study_list for f in glob.glob(os.path.join(SUAVE_ROOT, f"suave_dynamic_stability_outputs_{s}_*_{aircraft}.txt"))]
mat_files = [f for s in study_list for f in glob.glob(os.path.join(SUAVE_ROOT, f"suave_dynamic_stability_matrix_{s}_*_{aircraft}.txt"))]
eigfile = combine_eigen_outputs(eig_files, os.path.join(SUAVE_ROOT, f"suave_dynamic_stability_outputs_combined_{aircraft}.txt"))
matfile = combine_system_matrices(mat_files, os.path.join(SUAVE_ROOT, f"suave_dynamic_stability_matrix_combined_{aircraft}.txt"))
res = match_evals_to_emodes(matfile, eigfile, plot_matching=True)

#%% Time response properties derived from eigenvalues

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

# NILS (03.10.2026, fix A9): in (3.111), n/alpha denotes n_alpha = dn/dalpha (load factor per radian of AoA),
# not the ratio of trim load factor to trim AoA (which is wrong whenever alpha_0 != 0 and needed abs(aoa)).
# n_alpha = |Z_w| U0 / g is taken from the AVL system matrix in eigenvalue_classification.py.
def CAP_nalpha(omega_n, n_alpha):  # (1/s^2)
    return omega_n**2 / n_alpha

# NILS (03.10.2026, fixes A8 and A10): MIL-F-8785C / MIL-STD-1797A Level 1 limits for Class II (land-based)
# aircraft. Take-off and landing (mission points 1 and 6) are terminal flight phases (Category C), climb,
# cruise and descent are Category B. The Dutch roll limit was previously zeta >= 0.15 (Category B requires
# zeta >= 0.08 AND zeta*omega_n >= 0.15 rad/s AND omega_n >= 0.4 rad/s).
# NILS (04.10.2026): MIL_LIMITS and flight_phase_category() moved to stability_zones.py (imported above)
# MIL_LIMITS = {
#     'B': dict(zeta_sp_min=0.30, zeta_sp_max=2.00, CAP_sp_min=0.085, CAP_sp_max=3.6,
#               zeta_ph_min=0.04, tau_roll_max=1.4, T_2_sl_min=20.0,
#               zeta_dr_min=0.08, zeta_omega_dr_min=0.15, omega_n_dr_min=0.4),
#     'C': dict(zeta_sp_min=0.35, zeta_sp_max=1.30, CAP_sp_min=0.16, CAP_sp_max=3.6,
#               zeta_ph_min=0.04, tau_roll_max=1.4, T_2_sl_min=12.0,
#               zeta_dr_min=0.08, zeta_omega_dr_min=0.10, omega_n_dr_min=0.4),
# }
#
# def flight_phase_category(rid):
#     """NILS (03.10.2026, fix A10): rid = 1..6 (take-off, climb, start of cruise, end of cruise, descent, landing)"""
#     return 'C' if rid in (1, 6) else 'B'

def tau(lamda):  # (s)
    return T0p5(lamda) / np.log(2)  # page 125 in AE3202 Flight Dynamics Lecture Notes

def compute_lambdas_from_zeta_omega(zeta0, omega0):
    
    x = -zeta0 * omega0
    y_mag = omega0 * np.sqrt(max(0.0, 1.0 - zeta0**2))

    cplx = x + 1j * y_mag
    cplx_conj = x - 1j * y_mag
    return cplx, cplx_conj

#%% Plot MIL-STD-1797 requirements (section 3.4.2 in medium_PhDthesis_2013_flightdynconstrconcacmdiscanalsopt_morris)

def plot_mil_limits(aoa, mach, beta, h, n, W, n_alpha=None, category='B'):  # NILS (03.10.2026, fixes A9 and A10): added n_alpha and category
    
    lim = MIL_LIMITS[category]  # NILS (03.10.2026, fix A10)
    amb = Atmosphere(h)
    V = mach * amb.speed_of_sound[0]
    aoa = np.deg2rad(aoa)
    beta = np.deg2rad(beta)

    N_points = 20
    
    # Lateral stability limits
    # zeta_dr_min = 0.08  # 3rd row of Table 3.14 (Class II, Category B, Level I)
    # =============================================================================
    # 26.03.2026
    # zeta_dr_min = 0.15
    # =============================================================================
    # omega_n_dr_min = 0.4  # 3rd row of Table 3.14 (Class II, Category B, Level I)
    # T_2_sl_min = 20  # 2nd row of Table 3.13 (Category B, Level I)
    # tau_roll_max = 1.4  # 3rd row of Table 3.12 (Class II, Category B, Level I)
    # NILS (03.10.2026, fixes A8 and A10): limits by flight phase category (shading omits zeta*omega_n for now)
    zeta_dr_min = lim['zeta_dr_min']
    omega_n_dr_min = lim['omega_n_dr_min']
    T_2_sl_min = lim['T_2_sl_min']
    tau_roll_max = lim['tau_roll_max']
    
    # Dutch roll
    
    r_circle_dr_min = omega_n_dr_min
    theta = np.linspace(0, 2 * np.pi, N_points)
    x_circle_dr_min = r_circle_dr_min * np.cos(theta)
    y_circle_dr_min = r_circle_dr_min * np.sin(theta)
    
    # print('x_circle_dr_min, y_circle_dr_min =', x_circle_dr_min, y_circle_dr_min)
    
    x_vert_dr_max = -zeta_dr_min * r_circle_dr_min
    
    a_dr = zeta_dr_min
    t_dr = np.linspace(0, r_circle_dr_min / a_dr, N_points)
    x_ray_dr_max = -t_dr * a_dr
    y_ray_dr_max = t_dr * np.sqrt(1 - a_dr**2)
    
    # Spiral
    x_vert_sl_max = -np.log(0.5) / T_2_sl_min
    
    # Roll
    # x_vert_roll_max = -np.log(0.5)**2 / tau_roll_max
    x_vert_roll_max = -1.0 / tau_roll_max  # NILS (04.10.2026, fix B2): tau = -1/Re (see tau() above), i.e. Re <= -1/tau_max
    
    ###
    
    # Longitudinal stability limits
    zeta_ph_min = lim['zeta_ph_min']  # 1st row in Table 3.9 (Class II, Category B, Level I)  # NILS (03.10.2026, fix A10): by category
    zeta_sp_min = lim['zeta_sp_min']  # 1st row, 4th column in Table 3.10 (Class II, Category B, Level I)  # NILS (03.10.2026, fix A10): by category
    zeta_sp_max = lim['zeta_sp_max']  # 1st row, 5th column in Table 3.10 (Class II, Category B, Level I)  # NILS (03.10.2026, fix A10): by category
    CAP_sp_min = lim['CAP_sp_min']  # 1st row, 4th column in Table 3.11 (Class II, Category B, Level I)  # NILS (03.10.2026, fix A10): by category
    CAP_sp_max = lim['CAP_sp_max']  # 1st row, 5th column in Table 3.11 (Class II, Category B, Level I)  # NILS (03.10.2026, fix A10): by category
    
    # Short period
    
    # NOTE: abs() below not in original (3.111) in medium_PhDthesis_2013_flightdynconstrconcacmdiscanalsopt_morris but required to prevent nan when aoa < 0
    # omega_n_sp_min = np.sqrt(CAP_sp_min * n / abs(aoa))  # (3.111) in medium_PhDthesis_2013_flightdynconstrconcacmdiscanalsopt_morris rearranged
    # omega_n_sp_max = np.sqrt(CAP_sp_max * n / abs(aoa))  # (3.111) in medium_PhDthesis_2013_flightdynconstrconcacmdiscanalsopt_morris rearranged
    omega_n_sp_min = np.sqrt(CAP_sp_min * n_alpha)  # NILS (03.10.2026, fix A9): CAP = omega_n^2 / n_alpha rearranged
    omega_n_sp_max = np.sqrt(CAP_sp_max * n_alpha)  # NILS (03.10.2026, fix A9): CAP = omega_n^2 / n_alpha rearranged
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
    
#%% Mask MIL-STD-1797 requirements (section 3.4.2 in medium_PhDthesis_2013_flightdynconstrconcacmdiscanalsopt_morris)

def mask_mil_limits(
    mode,
    df_mode,
    rid=None,  # NILS (03.10.2026, fix A10): mission point 1..6 determines flight phase category
    # =============================================================================
    limit_subdict=None,
    # =============================================================================
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
    # NILS (03.10.2026, fixes A9 and A11): damping ratio and natural frequency from eigenvalue_classification.py
    # (correct for overdamped real pairs, for which zeta(ev) = 1 per root), and n_alpha from the A matrix
    zeta_mode = df_mode['zeta'].to_numpy()
    wn_mode = df_mode['wn'].to_numpy()
    n_alpha = df_mode['n_alpha'].to_numpy()
    lim = MIL_LIMITS[flight_phase_category(rid)]  # NILS (03.10.2026, fix A10)
    
    # # =============================================================================
    # limit_dict = {}
    # # =============================================================================
    
    # Longitudinal stability limits
    
    # Short period
    if mode == 'short-period':
        zeta_sp_min = lim['zeta_sp_min']  # NILS (03.10.2026, fix A10): was 0.3 (Category B) for all points
        zeta_sp_max = lim['zeta_sp_max']  # NILS (03.10.2026, fix A10): was 2.0 (Category B) for all points
        CAP_sp_min = lim['CAP_sp_min']  # NILS (03.10.2026, fix A10): was 0.085 (Category B) for all points
        CAP_sp_max = lim['CAP_sp_max']  # NILS (03.10.2026, fix A10): was 3.6 for all points
        # zeta_sp = zeta(ev_array)
        # CAP_sp = CAP(ev_array, n, aoa)
        zeta_sp = zeta_mode  # NILS (03.10.2026, fix A11)
        CAP_sp = CAP_nalpha(wn_mode, n_alpha)  # NILS (03.10.2026, fix A9)
        valid_mask = (
            (zeta_sp > zeta_sp_min) & (zeta_sp < zeta_sp_max) &
            (CAP_sp > CAP_sp_min) & (CAP_sp < CAP_sp_max)
        )
        # =============================================================================
        limit_subdict['zeta_sp_min'] = (zeta_sp - zeta_sp_min) / zeta_sp_min
        limit_subdict['zeta_sp_max'] = (zeta_sp_max - zeta_sp) / zeta_sp_max
        limit_subdict['CAP_sp_min'] = (CAP_sp - CAP_sp_min) / CAP_sp_min
        limit_subdict['CAP_sp_max'] = (CAP_sp_max - CAP_sp) / CAP_sp_max
        # =============================================================================
        
    # Phugoid
    elif mode == 'phugoid':
        zeta_ph_min = lim['zeta_ph_min']  # 1st row in Table 3.9 (Class II, Category B, Level I)
        # zeta_ph = zeta(ev_array)
        zeta_ph = zeta_mode  # NILS (03.10.2026, fix A11)
        valid_mask = (zeta_ph > zeta_ph_min)
        # =============================================================================
        limit_subdict['zeta_ph_min'] = (zeta_ph - zeta_ph_min) / zeta_ph_min
        # limit_subdict['zeta_ph_min'] = zeta_ph / zeta_ph_min
        # =============================================================================
        
    # Lateral stability limits
    
    # Roll
    elif mode == 'roll':
        tau_roll_max = lim['tau_roll_max']  # 3rd row of Table 3.12 (Class II, Category B, Level I)
        tau_roll = tau(ev_array)
        valid_mask = (tau_roll < tau_roll_max)
        # =============================================================================
        limit_subdict['tau_roll_max'] = (tau_roll_max - tau_roll) / tau_roll_max
        # =============================================================================
        
    # Spiral
    elif mode == 'spiral':
        T_2_sl_min = lim['T_2_sl_min']  # NILS (03.10.2026, fix A10): 20 s (Category B), 12 s (Category C)
        T_2_sl = T2(ev_array)
        valid_mask = (T_2_sl > T_2_sl_min) | (T_2_sl < 0)  # only applies to positive real parts
        # =============================================================================
        limit_subdict['T_2_sl_min'] = (T_2_sl - T_2_sl_min) / T_2_sl_min
        # =============================================================================
        
    # Dutch roll
    elif mode == 'dutch-roll':
        # zeta_dr_min = 0.08  # 3rd row of Table 3.14 (Class II, Category B, Level I)
        # =============================================================================
        # 26.03.2026
        # zeta_dr_min = 0.15
        # =============================================================================
        # omega_n_dr_min = 0.4  # 3rd row of Table 3.14 (Class II, Category B, Level I)
        # NILS (03.10.2026, fix A8): zeta >= 0.08, zeta*omega_n >= 0.15 (B) / 0.10 (C) rad/s, omega_n >= 0.4 rad/s
        zeta_dr_min = lim['zeta_dr_min']
        zeta_omega_dr_min = lim['zeta_omega_dr_min']
        omega_n_dr_min = lim['omega_n_dr_min']
        zeta_dr = zeta(ev_array)
        omega_n_dr = omega0(ev_array)
        # valid_mask = (zeta_dr > zeta_dr_min) & (omega_n_dr > omega_n_dr_min)
        valid_mask = (zeta_dr > zeta_dr_min) & (zeta_dr * omega_n_dr > zeta_omega_dr_min) & (omega_n_dr > omega_n_dr_min)
        # print('zeta_dr, zeta_dr_min =', zeta_dr, zeta_dr_min)
        # print('valid_mask =', valid_mask)
        # raise Exception()
        # =============================================================================
        limit_subdict['zeta_dr_min'] = (zeta_dr - zeta_dr_min) / zeta_dr_min
        limit_subdict['zeta_omega_dr_min'] = (zeta_dr * omega_n_dr - zeta_omega_dr_min) / zeta_omega_dr_min
        limit_subdict['omega_n_dr_min'] = (omega_n_dr - omega_n_dr_min) / omega_n_dr_min
        # =============================================================================
        
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
        # # =============================================================================
        # limit_dict,
        # # =============================================================================
    )

#%% Plot root loci of longitudinal or lateral eigenmotions

# import matplotlib as mpl
# fontsize = 25
# params = {
#     'axes.labelsize': fontsize,
#     'axes.titlesize': fontsize,
#     'font.size': fontsize,
#     'legend.fontsize': fontsize,
#     'legend.title_fontsize': fontsize,
#     'xtick.labelsize': fontsize,
#     'ytick.labelsize': fontsize,
# }
# mpl.rcParams.update(params)

# User input
modes_type = 'longitudinal'  # 'longitudinal', 'lateral'
modes_type = sys.argv[2] if len(sys.argv) > 2 else modes_type  # NILS (04.10.2026)

edgecolor = 'k'
linewidth = 0.5
alpha_invalid = 0.2
alpha_valid = 1.0

# fig = plt.figure(figsize=(15,9))
# gs = gridspec.GridSpec(
#     2, 3,
#     figure=fig,
#     left=0.0675,  # for longitudinal
#     # right=0.8,  # for longitudinal
#     # right=0.775,  # for lateral
#     right=0.95,
#     # right=1.0,
#     # bottom=0.1375,
#     bottom=0.125,
#     top=0.975,
#     hspace=0.2,
#     wspace=0.2,  # for longitudinal
#     # wspace=0.4,  # for lateral
# )
# axes = [fig.add_subplot(gs[_i,_j]) for _i in range(2) for _j in range(3)]
# titles = ['1. Take-off', '2. Climb', '3. Start of\ncruise', '4. End of\ncruise', '5. Descent', '6. Landing']

# =============================================================================
axes = np.ones(6)
fig, ax = plt.subplots(figsize=(7, 5))
# =============================================================================

# Initialise data dictionary
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
    # for cat in ['fuselage', 'wing', 'nacelle']
    for cat in res.keys()  # NILS (04.10.2026, fix B4): categories from eigenvalue_classification.py (incl. 'mixed')
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

# NILS (04.10.2026, fix B3): left x-limit of the lateral plots from the fastest roll root (roll roots between
# -1.4 and -3.2 1/s were clipped by the fixed xlim of -2.0)
roll_re = [res[cat]['roll'][rid]['Re'].to_numpy() for cat in res for rid in range(1, 7) if len(res[cat]['roll'][rid]) > 0]
x_lim_lat_left = min(-2.0, 1.1 * np.nanmin(np.concatenate(roll_re))) if len(roll_re) > 0 else -2.0

# =============================================================================
if modes_type == 'longitudinal':
    modes_sublist = ['short-period', 'phugoid']
    markers = ['o', 'd']
    
    # Axis settings
    ax.set_xlim((-1.5, 1))
    ax.set_ylim((-0.1, 5))
    ax.set_xticks([-1.5, -1, -0.5, 0, 0.5, 1])
    # ax.set_yticks([-5, -4, -3, -2, -1, 1, 2, 3, 4, 5])  #if want to plot both positive and negative imaginary axis
    ax.set_yticks([1, 2, 3, 4, 5])  #if want to plot positive imaginary axis only
    
    # Inset axis and settings
    axins = inset_axes(
        ax,
        width="100%",
        height="100%",
        bbox_to_anchor=(0.65, 0.1, 0.35, 0.35),  # (x0, y0, w, h) in ax coords
        bbox_transform=ax.transAxes,
        borderpad=0
    )
    axins.axvline(0.0, color='k')
    axins.axhline(0.0, color='k')
    # axins.tick_params(labelleft=False, labelbottom=False)
    # axins.tick_params(labelleft=False, labelright=True, labelbottom=True, labelsize=0.5 * fontsize, pad=2)  # NILS (04.10.2026, fix B6): scale of the phugoid inset (right, clear of the main y-axis)
    # axins.xaxis.set_major_locator(plt.MaxNLocator(2))
    # axins.yaxis.set_major_locator(plt.MaxNLocator(2))
    # axins.set_title('Phugoid', fontsize=0.5 * fontsize, pad=2)  # NILS (04.10.2026, fix B6)
    # axins.tick_params(axis='y', which='both', length=0)
    # axins.tick_params(axis='x', which='both', length=0)
    axins.set_xticks([])
    axins.set_yticks([])
    mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")
    
    # =============================================================================
    ev_real_phugoid_min_list = []
    ev_real_phugoid_max_list = []
    ev_imag_phugoid_min_list = []
    ev_imag_phugoid_max_list = []
    # =============================================================================
    
elif modes_type == 'lateral':
    modes_sublist = ['roll', 'dutch-roll', 'spiral']
    markers = ['o', 'd', '^']
    
    # y_lim_lat = 3
    
    # # ax.set_xlim((-2.0, 0.1))
    # ax.set_xlim((x_lim_lat_left, 0.1))  # NILS (04.10.2026, fix B3)
    # ax.set_ylim((-0.1, y_lim_lat))
    # # ax.set_xticks([-1, -0.5])
    # ax.set_xticks([t for t in [-3, -2, -1] if t > x_lim_lat_left])  # NILS (04.10.2026, fix B3)
    # ax.set_yticks([1, 2])
    # =============================================================================
    y_lim_lat = 5
    
    # Axis settings
    ax.set_xlim((-1.5, 1))
    ax.set_ylim((-0.1, 5))
    ax.set_xticks([-1.5, -1, -0.5, 0, 0.5, 1])
    # ax.set_yticks([-5, -4, -3, -2, -1, 1, 2, 3, 4, 5])  #if want to plot both positive and negative imaginary axis
    ax.set_yticks([1, 2, 3, 4, 5])  #if want to plot positive imaginary axis only
    # =============================================================================
    '''
    # Add second x-axis left of break to plot roll subsidence roots    
    gap = 0.2  # fraction of ax width
    axleft = ax.inset_axes(
        # [-0.25 - gap, 0.0, 0.2, 1.0],  # <-- shifted left
        [-0.1 - gap, 0.0, 0.35, 1.0],  # <-- shifted left
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
    '''
    # =============================================================================
    # Inset axis and settings
    axins = inset_axes(
        ax,
        width="100%",
        height="100%",
        bbox_to_anchor=(0.65, 0.1, 0.35, 0.35),  # (x0, y0, w, h) in ax coords
        bbox_transform=ax.transAxes,
        borderpad=0
    )
    axins.axvline(0.0, color='k')
    axins.axhline(0.0, color='k')
    # axins.tick_params(labelleft=False, labelbottom=False)
    # axins.tick_params(labelleft=False, labelright=True, labelbottom=True, labelsize=0.5 * fontsize, pad=2)  # NILS (04.10.2026, fix B6): scale of the phugoid inset (right, clear of the main y-axis)
    # axins.xaxis.set_major_locator(plt.MaxNLocator(2))
    # axins.yaxis.set_major_locator(plt.MaxNLocator(2))
    # axins.set_title('Phugoid', fontsize=0.5 * fontsize, pad=2)  # NILS (04.10.2026, fix B6)
    # axins.tick_params(axis='y', which='both', length=0)
    # axins.tick_params(axis='x', which='both', length=0)
    axins.set_xticks([])
    axins.set_yticks([])
    mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")
    # =============================================================================
    
    # =============================================================================
    ev_real_spiral_min_list = []
    ev_real_spiral_max_list = []
    ev_imag_spiral_min_list = []
    ev_imag_spiral_max_list = []
    # =============================================================================
# =============================================================================

# =============================================================================
def create_empty_dict(d):
    return {
        key: create_empty_dict(value) if isinstance(value, dict) else {}
        for key, value in d.items()
    }

limit_dict = create_empty_dict(res)
print(limit_dict)
# =============================================================================

for row_idx, _ax in enumerate(axes):
    rid = row_idx + 1
    
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
    
    # for j, cat in enumerate(['fuselage','wing','nacelle']):
    for j, cat in enumerate(res.keys()):  # NILS (04.10.2026, fix B4)
    # =============================================================================
    # for j, cat in enumerate(['wing']):
    # for j, cat in enumerate(['fuselage']):
    # =============================================================================
    
        print('cat =', cat)
        
        for k, mode in enumerate(modes_sublist):
        # =============================================================================
        # for k, mode in enumerate(['dutch-roll']):
        # for k, mode in enumerate(['roll']):
        # for k, mode in enumerate(['spiral']):
        # =============================================================================
        
            print('mode =', mode)
            
            if len(res[cat][mode][rid]) == 0:  # NILS (04.10.2026): category not present in the selected studies
                continue
            # =============================================================================
            limit_subdict = limit_dict[cat][mode][rid]
            # =============================================================================
            # print('cat, mode, rid =', cat, mode, rid)
            df_cat = res[cat][mode][rid]
            (
            ev_valid, ev_invalid,
            sigma_fcs_valid, sigma_fcs_invalid,
            fcs_loc_valid, fcs_loc_invalid,
            span_loc_valid, span_loc_invalid,
            # # =============================================================================
            # limit_dict,
            # # =============================================================================
            ) = \
            mask_mil_limits(
                mode,
                df_cat,
                rid,  # NILS (03.10.2026, fix A10)
                # =============================================================================
                limit_subdict,
                # =============================================================================
            )
            
            ev_valid_dict[cat][mode][rid]['ev'] = ev_valid
            ev_invalid_dict[cat][mode][rid]['ev'] = ev_invalid
            ev_valid_dict[cat][mode][rid]['sigma_fcs'] = sigma_fcs_valid
            ev_invalid_dict[cat][mode][rid]['sigma_fcs'] = sigma_fcs_invalid
            ev_valid_dict[cat][mode][rid]['fcs_loc'] = fcs_loc_valid
            ev_invalid_dict[cat][mode][rid]['fcs_loc'] = fcs_loc_invalid
            ev_valid_dict[cat][mode][rid]['span_loc'] = span_loc_valid
            ev_invalid_dict[cat][mode][rid]['span_loc'] = span_loc_invalid
            
            # ax.scatter(ev_valid.real, ev_valid.imag, marker=markers[k], facecolor=opaque_color_from_hex(colors[j], alpha=alpha_valid), edgecolor=opaque_color_from_hex(mcolors.to_hex(edgecolor), alpha=alpha_valid), linewidth=linewidth, zorder=100, clip_on=True)
            # ax.scatter(ev_invalid.real, ev_invalid.imag, marker=markers[k], facecolor=opaque_color_from_hex(colors[j], alpha=alpha_invalid), edgecolor=opaque_color_from_hex(mcolors.to_hex(edgecolor), alpha=alpha_invalid), linewidth=linewidth, zorder=100, clip_on=True)
            ###
            ax.scatter(ev_valid.real, ev_valid.imag, marker=markers[k], facecolor=COLOR_OK, edgecolor=edgecolor, linewidth=linewidth, zorder=100, clip_on=True)  # NILS (04.10.2026): COLOR_OK = 'green' (stability_zones.py)
            ax.scatter(ev_invalid.real, ev_invalid.imag, marker=markers[k], facecolor=COLOR_FAIL, edgecolor=edgecolor, linewidth=linewidth, zorder=100, clip_on=True)  # NILS (04.10.2026): COLOR_FAIL = 'red'
            ###
        
            if mode == 'short-period':
                ev_valid_sp_dict[cat][mode][rid]['ev'] = ev_valid
                ev_invalid_sp_dict[cat][mode][rid]['ev'] = ev_invalid
                # ev_valid_sp_dict[cat][mode][rid]['sigma_fcs'] = sigma_fcs_valid
                # ev_invalid_sp_dict[cat][mode][rid]['sigma_fcs'] = sigma_fcs_invalid
                # ev_valid_sp_dict[cat][mode][rid]['fcs_loc'] = fcs_loc_valid
                # ev_invalid_sp_dict[cat][mode][rid]['fcs_loc'] = fcs_loc_invalid
                # ev_valid_sp_dict[cat][mode][rid]['span_loc'] = span_loc_valid
                # ev_invalid_sp_dict[cat][mode][rid]['span_loc'] = span_loc_invalid
            elif mode == 'phugoid':
                # axins.scatter(ev_valid.real, ev_valid.imag, marker=markers[k], facecolor=opaque_color_from_hex(colors[j], alpha=alpha_valid), edgecolor=opaque_color_from_hex(mcolors.to_hex(edgecolor), alpha=alpha_valid), linewidth=linewidth, zorder=100, clip_on=False)
                # axins.scatter(ev_invalid.real, ev_invalid.imag, marker=markers[k], facecolor=opaque_color_from_hex(colors[j], alpha=alpha_invalid), edgecolor=opaque_color_from_hex(mcolors.to_hex(edgecolor), alpha=alpha_invalid), linewidth=linewidth, zorder=100, clip_on=False)
                ###
                axins.scatter(ev_valid.real, ev_valid.imag, marker=markers[k], facecolor=COLOR_OK, edgecolor=edgecolor, linewidth=linewidth, zorder=100, clip_on=False)
                axins.scatter(ev_invalid.real, ev_invalid.imag, marker=markers[k], facecolor=COLOR_FAIL, edgecolor=edgecolor, linewidth=linewidth, zorder=100, clip_on=False)
                ###
                ev_valid_ph_dict[cat][mode][rid]['ev'] = ev_valid
                ev_invalid_ph_dict[cat][mode][rid]['ev'] = ev_invalid
            elif mode == 'roll':
                '''
                # axleft.scatter(ev_valid.real, ev_valid.imag, marker=markers[k], facecolor=opaque_color_from_hex(colors[j], alpha=alpha_valid), edgecolor=opaque_color_from_hex(mcolors.to_hex(edgecolor), alpha=alpha_valid), linewidth=linewidth, zorder=100, clip_on=False)
                # axleft.scatter(ev_invalid.real, ev_invalid.imag, marker=markers[k], facecolor=opaque_color_from_hex(colors[j], alpha=alpha_invalid), edgecolor=opaque_color_from_hex(mcolors.to_hex(edgecolor), alpha=alpha_invalid), linewidth=linewidth, zorder=100, clip_on=False)
                axleft.scatter(ev_valid.real, ev_valid.imag, marker=markers[k], facecolor='green', edgecolor=edgecolor, linewidth=linewidth, zorder=100, clip_on=False)
                axleft.scatter(ev_invalid.real, ev_invalid.imag, marker=markers[k], facecolor='red', edgecolor=edgecolor, linewidth=linewidth, zorder=100, clip_on=False)
                # ev_valid_roll_dict[cat][mode][rid]['ev'] = ev_valid
                # ev_invalid_roll_dict[cat][mode][rid]['ev'] = ev_invalid
                '''
                pass
            elif mode == 'spiral':
                ev_valid_sl_dict[cat][mode][rid]['ev'] = ev_valid
                ev_invalid_sl_dict[cat][mode][rid]['ev'] = ev_invalid
                # =============================================================================
                axins.scatter(ev_valid.real, ev_valid.imag, marker=markers[k], facecolor=COLOR_OK, edgecolor=edgecolor, linewidth=linewidth, zorder=100, clip_on=False)
                axins.scatter(ev_invalid.real, ev_invalid.imag, marker=markers[k], facecolor=COLOR_FAIL, edgecolor=edgecolor, linewidth=linewidth, zorder=100, clip_on=False)
                # =============================================================================
                # pass
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
            n_alpha_array = df_cat['n_alpha'].to_numpy()  # NILS (03.10.2026, fix A9)
            
            # print('len(df_cat) =', len(df_cat))
            
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
                plot_mil_limits(aoa, mach, beta, h, n, W, n_alpha=n_alpha_array[l], category=flight_phase_category(rid))  # NILS (03.10.2026, fixes A9 and A10)
                
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
        
        # Short period boundaries (commented for clearer plot)
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
        
        # Phugoid boundaries (commented for clearer plot)
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
        # '''
        # Roll subsidence boundaries (commented for clearer plot)
        # ax.axvline(x_vert_roll_max, color='black', alpha=alpha_constr)
        y_stacked_ray_roll_max = [-y_lim_lat, y_lim_lat]
        ax.fill_betweenx(
            # y=y_stacked_ray_roll_max, x1=-y_lim_lat, x2=x_vert_roll_max,
            y=y_stacked_ray_roll_max, x1=x_lim_lat_left, x2=x_vert_roll_max,  # NILS (04.10.2026, fix B3): to the left x-limit
            # facecolor='None', edgecolor='black', hatch='...', alpha=0.1, zorder=-1,
            # =============================================================================
            facecolor=opaque_color_from_hex(mcolors.to_hex("black"), alpha=0.1), edgecolor='None', zorder=-1,
            # =============================================================================
        )
        '''
        axleft.fill_betweenx(
            y=y_stacked_ray_roll_max, x1=-65, x2=5,
            facecolor='None', edgecolor='black', hatch='...', alpha=0.1, zorder=-1,
        )
        '''
        # Spiral boundaries (commented for clearer plot)
        # ax.axvline(x_vert_sl_max, color='black', alpha=alpha_constr)
        y_stacked_ray_sl_max = [-y_lim_lat, y_lim_lat + 1]  # NILS: extend slightly upward to make border invisible
        ax.fill_betweenx(
            # y=y_stacked_ray_sl_max, x1=-y_lim_lat, x2=x_vert_sl_max,
            y=y_stacked_ray_sl_max, x1=x_lim_lat_left, x2=x_vert_sl_max,  # NILS (04.10.2026, fix B3): to the left x-limit
            # facecolor='lightgrey', edgecolor='None', alpha=0.5, zorder=-1,
            # =============================================================================
            facecolor='None', edgecolor=opaque_color_from_hex(mcolors.to_hex("black"), alpha=0.2), hatch='...', zorder=-1,
            # =============================================================================
        )
        '''
        axleft.fill_betweenx(
            y=y_stacked_ray_sl_max, x1=-65, x2=5,
            facecolor='lightgrey', edgecolor='None', alpha=0.5, zorder=-1,
        )
        '''
        # Dutch roll boundaries (commented for clearer plot)
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
        # ax.fill(x_wedge_dr, y_wedge_dr, facecolor='None', edgecolor='black', hatch='//', alpha=0.1, zorder=-1)
        # NILS (04.10.2026): Dutch roll ok region including zeta*omega_n >= 0.15 (B) / 0.10 (C) rad/s, i.e. the
        # intersection of Re <= -(zeta*omega_n)_min, zeta >= zeta_min (wedge) and omega_n >= omega_n_min (circle)
        lim_lat = MIL_LIMITS[flight_phase_category(rid)]
        y_dr_ok = np.linspace(-y_lim_lat, y_lim_lat + 1, 401)  # NILS: extend slightly upward to make border invisible
        x_dr_ok = np.minimum.reduce([
            np.full_like(y_dr_ok, -lim_lat['zeta_omega_dr_min']),
            -np.abs(y_dr_ok) * lim_lat['zeta_dr_min'] / np.sqrt(1.0 - lim_lat['zeta_dr_min']**2),
            -np.sqrt(np.maximum(lim_lat['omega_n_dr_min']**2 - y_dr_ok**2, 0.0)),
        ])
        ax.fill_betweenx(y_dr_ok, x_lim_lat_left, x_dr_ok, facecolor='None', edgecolor=opaque_color_from_hex(mcolors.to_hex("black"), alpha=0.2), hatch='//', zorder=-1)
        x_wedge_dr = np.concatenate([
            x_wedge_dr_max, [x_ray_dr_max[-1], -65, -65, x_ray_dr_max[-1]],
        ])
        y_wedge_dr = np.concatenate([
            y_wedge_dr_max, [-y_ray_dr_max[-1], -y_ray_dr_max[-1], y_ray_dr_max[-1], y_ray_dr_max[-1]],
        ])
        '''
        axleft.fill(x_wedge_dr, y_wedge_dr, facecolor='None', edgecolor='black', hatch='//', alpha=0.1, zorder=-1)
        '''
        
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
        
        # =============================================================================
        ev_real_phugoid_min_list.append(ev_real_phugoid_min)
        ev_real_phugoid_max_list.append(ev_real_phugoid_max)
        ev_imag_phugoid_min_list.append(ev_imag_phugoid_min)
        ev_imag_phugoid_max_list.append(ev_imag_phugoid_max)
        
        axins.set_xlim(np.nanmin(ev_real_phugoid_min_list), np.nanmax(ev_real_phugoid_max_list))
        axins.set_ylim(np.nanmin(ev_imag_phugoid_min_list), np.nanmax(ev_imag_phugoid_max_list))
        
    elif modes_type == 'lateral':
        
        y_stacked_ray_sl_max = [-y_lim_lat, y_lim_lat]
        axins.fill_betweenx(
            # y=y_stacked_ray_sl_max, x1=-y_lim_lat, x2=x_vert_sl_max,
            y=y_stacked_ray_sl_max, x1=x_lim_lat_left, x2=x_vert_sl_max,  # NILS (04.10.2026, fix B3): to the left x-limit
            # facecolor='lightgrey', edgecolor='None', alpha=0.5, zorder=-1,
            # =============================================================================
            facecolor='None', edgecolor=opaque_color_from_hex(mcolors.to_hex("black"), alpha=0.2), hatch='...', zorder=-1,
            # =============================================================================
        )
        
        # <><><><><><><>
        arrays = [
            np.ravel(mode_dict[rid]['ev'])
            for cat_dict in ev_valid_sl_dict.values()
            for mode_dict in cat_dict.values()
            if mode_dict[rid]['ev'] is not None and np.asarray(mode_dict[rid]['ev']).size > 0
        ]
        ev_valid_spiral_array = (
            np.concatenate(arrays)
            if len(arrays) > 0
            else np.empty(0, dtype=np.complex128)
        )
        
        arrays = [
            np.ravel(mode_dict[rid]['ev'])
            for cat_dict in ev_invalid_sl_dict.values()
            for mode_dict in cat_dict.values()
            if mode_dict[rid]['ev'] is not None and np.asarray(mode_dict[rid]['ev']).size > 0
        ]
        ev_invalid_spiral_array = (
            np.concatenate(arrays)
            if len(arrays) > 0
            else np.empty(0, dtype=np.complex128)
        )
        
        ev_real_spiral_min = np.min(np.hstack((ev_valid_spiral_array.real, ev_invalid_spiral_array.real)))
        ev_real_spiral_max = np.max(np.hstack((ev_valid_spiral_array.real, ev_invalid_spiral_array.real)))
        ev_imag_spiral_min = np.min(np.hstack((ev_valid_spiral_array.imag, ev_invalid_spiral_array.imag)))
        ev_imag_spiral_max = np.max(np.hstack((ev_valid_spiral_array.imag, ev_invalid_spiral_array.imag)))
        # <><><><><><><>
        
        ev_real_spiral_min_list.append(ev_real_spiral_min)
        ev_real_spiral_max_list.append(ev_real_spiral_max)
        ev_imag_spiral_min_list.append(ev_imag_spiral_min)
        ev_imag_spiral_max_list.append(ev_imag_spiral_max)
        
        axins.set_xlim(np.nanmin(ev_real_spiral_min_list), np.nanmax(ev_real_spiral_max_list))
        axins.set_ylim(np.nanmin(ev_imag_spiral_min_list), np.nanmax(ev_imag_spiral_max_list))
        # =============================================================================
        
# =============================================================================

add_margin(axins, m=0.5)

# ax.set_title(titles[row_idx], pad=20, loc='left', x=0.25, y=0.85, va='top', ha='center')
ax.spines[['right','top']].set_visible(False)
ax.tick_params(axis='y', which='both', right=False, length=0)
ax.tick_params(axis='x', which='both', length=0)
ax.spines['bottom'].set_position(('data', 0.0))
ax.spines['left'].set_position(('data', 0.0))

# Figure legend

if modes_type == 'longitudinal':
    
    # ax = fig.add_axes([0.8, 0.0, 0.225, 1.0])
    # ax = fig.add_axes([0.0, 1.0, 1.0, 0.1])
    # legend_marker_colours = [colors[0], colors[1], colors[2], colors[0], colors[1], colors[2]]
    legend_marker_colours = ['k', 'k']
    # legend_markers = ['o', 'o', 'o', 's', 's', 's']
    legend_markers = ['o', 'd']
    legend_marker_labels = [
        # 'Short period nacelle', 'Short period wing', 'Short period fuselage',
        # 'Phugoid nacelle', 'Phugoid wing', 'Phugoid fuselage',
        'Short period', 'Phugoid',
    ]
    legend_elements = [
        mpatches.Patch(
            facecolor='none',
            edgecolor='black',
            hatch='//',
            label='Short period ok',
            alpha=0.1,
        ),
        mpatches.Patch(
            facecolor='black',
            edgecolor='none',
            label='Phugoid ok',
            alpha=0.1,
        )
    ] + [
        Line2D(
            [], [], marker=legend_markers[_], color='none', markerfacecolor=legend_marker_colours[_],
            markeredgecolor=legend_marker_colours[_], markersize=8, linestyle='None',
            label=legend_marker_labels[_],
        ) for _ in range(len(legend_marker_colours))
    ]
    # ax.legend(
    #     handles=legend_elements, loc='lower center', bbox_to_anchor=(0.5, 1.0), ncol=4, labelspacing=1.0, frameon=False,
    #     fontsize=14,
    # )
    ax.legend(
        handles=legend_elements,
        loc='upper right',
        bbox_to_anchor=(1.05, 1.05),
        ncol=1,
        handlelength=1.0,
        handleheight=1.0,
        columnspacing=1.0,
        handletextpad=0.5,
        labelspacing=0.5,
        frameon=False,
        fontsize=14,
    )
elif modes_type == 'lateral':
    
    # ax = fig.add_axes([0.8, 0.0, 0.225, 1.0])
    # ax = fig.add_axes([0.0, 1.0, 1.0, 0.1])
    # legend_marker_colours = [colors[0], colors[1], colors[2], colors[0], colors[1], colors[2], colors[0], colors[1], colors[2]]
    legend_marker_colours = ['k', 'k', 'k']
    # legend_markers = ['o', 'o', 'o', 's', 's', 's', 'd', 'd', 'd']
    legend_markers = ['o', 'd', '^']
    legend_marker_labels = [
        # 'Roll subsidence nacelle', 'Roll subsidence wing', 'Roll subsidence fuselage',
        # 'Spiral nacelle', 'Spiral wing', 'Spiral fuselage',
        # 'Dutch roll nacelle', 'Dutch roll wing', 'Dutch roll fuselage',
        # 'Roll', 'Spiral', 'Dutch roll',
        'Roll', 'Dutch roll', 'Spiral',  # NILS (04.10.2026, fix B1): same order as modes_sublist/markers ('o', 'd', '^')
    ]
    legend_elements = [
        mpatches.Patch(
            facecolor='black',
            edgecolor='none',
            # label='Spiral ok',
            label='Roll ok',
            alpha=0.1,
        ),
        mpatches.Patch(
            facecolor='none',
            edgecolor='black',
            hatch='...',
            # label='Roll ok',
            label='Spiral ok',
            alpha=0.1,
        ),
        mpatches.Patch(
            facecolor='none',
            edgecolor='black',
            hatch='//',
            label='Dutch roll ok',
            alpha=0.1,
        )
    ] + [
        Line2D(
            [], [], marker=legend_markers[_], color='none', markerfacecolor=legend_marker_colours[_],
            markeredgecolor=legend_marker_colours[_], markersize=8, linestyle='None',
            label=legend_marker_labels[_],
        ) for _ in range(len(legend_marker_colours))
    ]
    # ax.legend(
    #     handles=legend_elements, loc='lower center', bbox_to_anchor=(0.5, 0.0), ncol=6, labelspacing=1.0, frameon=False,
    # )
    ax.legend(
        handles=legend_elements,
        loc='upper right',
        bbox_to_anchor=(1.025, 1.05),
        ncol=1,
        handlelength=1.0,
        handleheight=1.0,
        columnspacing=1.0,
        handletextpad=0.5,
        labelspacing=0.5,
        frameon=False,
        fontsize=14,
    )

ax.set_xlabel('Real component [1/s]', labelpad=15)
ax.set_ylabel('Imaginary component [1/s]', labelpad=15)
ax.yaxis.set_label_coords(-0.05, 0.5)
# =============================================================================

r'''
# Figure legend

if modes_type == 'longitudinal':
    
    # legend_ax = fig.add_axes([0.8, 0.0, 0.225, 1.0])
    legend_ax = fig.add_axes([0.0, 1.0, 1.0, 0.1])
    # legend_marker_colours = [colors[0], colors[1], colors[2], colors[0], colors[1], colors[2]]
    legend_marker_colours = ['k', 'k']
    # legend_markers = ['o', 'o', 'o', 's', 's', 's']
    legend_markers = ['o', 'd']
    legend_marker_labels = [
        # 'Short period nacelle', 'Short period wing', 'Short period fuselage',
        # 'Phugoid nacelle', 'Phugoid wing', 'Phugoid fuselage',
        'Short period', 'Phugoid',
    ]
    legend_elements = [
        mpatches.Patch(
            facecolor='none',
            edgecolor='black',
            hatch='//',
            label='Short period ok',
            alpha=0.1,
        ),
        mpatches.Patch(
            facecolor='black',
            edgecolor='none',
            label='Phugoid ok',
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
        handles=legend_elements, loc='lower center', bbox_to_anchor=(0.5, 0.0), ncol=4, labelspacing=1.0, frameon=False,
    )
elif modes_type == 'lateral':
    
    # legend_ax = fig.add_axes([0.8, 0.0, 0.225, 1.0])
    legend_ax = fig.add_axes([0.0, 1.0, 1.0, 0.1])
    # legend_marker_colours = [colors[0], colors[1], colors[2], colors[0], colors[1], colors[2], colors[0], colors[1], colors[2]]
    legend_marker_colours = ['k', 'k', 'k']
    # legend_markers = ['o', 'o', 'o', 's', 's', 's', 'd', 'd', 'd']
    legend_markers = ['o', 'd', '^']
    legend_marker_labels = [
        # 'Roll subsidence nacelle', 'Roll subsidence wing', 'Roll subsidence fuselage',
        # 'Spiral nacelle', 'Spiral wing', 'Spiral fuselage',
        # 'Dutch roll nacelle', 'Dutch roll wing', 'Dutch roll fuselage',
        # 'Roll', 'Spiral', 'Dutch roll',
        'Roll', 'Dutch roll', 'Spiral',  # NILS (04.10.2026, fix B1): same order as modes_sublist/markers ('o', 'd', '^')
    ]
    legend_elements = [
        mpatches.Patch(
            facecolor='black',
            edgecolor='none',
            label='Spiral ok',
            alpha=0.1,
        ),
        mpatches.Patch(
            facecolor='none',
            edgecolor='black',
            hatch='...',
            label='Roll ok',
            alpha=0.1,
        ),
        mpatches.Patch(
            facecolor='none',
            edgecolor='black',
            hatch='//',
            label='Dutch roll ok',
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
        handles=legend_elements, loc='lower center', bbox_to_anchor=(0.5, 0.0), ncol=6, labelspacing=1.0, frameon=False,
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
'''
    
# fig.text(0.5, 0.025, 'Real component (1/s)', ha='center', va='bottom')
# fig.text(0.025, 0.5, 'Imaginary component (1/s)', ha='right', va='center', rotation='vertical')
# plt.savefig('longitudinal_roots.png', format='png', dpi=600)
# plt.savefig('lateral_roots.png', format='png', dpi=600)
if save_fig:  # NILS (04.10.2026)
    plt.savefig(f'{modes_type}_roots_{aircraft}.png', format='png', dpi=300)
    plt.savefig(f'{modes_type}_roots_{aircraft}.svg', format='svg')

plt.show()

#%%

# =============================================================================
def create_empty_dict(d):
    return {
        key: create_empty_dict(value) if isinstance(value, dict) else {}
        for key, value in d.items()
    }

limit_dict = create_empty_dict(res)

fig, ax = plt.subplots()
# =============================================================================

for row_idx, _ax in enumerate(axes):
    rid = row_idx + 1
    
    for j, cat in enumerate(res.keys()):  # NILS (04.10.2026, fix B4)
    
        for k, mode in enumerate(modes_sublist):
            
            if len(res[cat][mode][rid]) == 0:  # NILS (04.10.2026): category not present in the selected studies
                continue
            # =============================================================================
            limit_subdict = limit_dict[cat][mode][rid]
            # =============================================================================
            df_cat = res[cat][mode][rid]
            (
            ev_valid, ev_invalid,
            sigma_fcs_valid, sigma_fcs_invalid,
            fcs_loc_valid, fcs_loc_invalid,
            span_loc_valid, span_loc_invalid,
            # # =============================================================================
            # limit_dict,
            # # =============================================================================
            ) = \
            mask_mil_limits(
                mode,
                df_cat,
                rid,  # NILS (03.10.2026, fix A10)
                # =============================================================================
                limit_subdict,
                # =============================================================================
            )
            # print('limit_dict =', limit_dict)
            
            
            if cat == 'fuselage':
                
                # Longitudinal
                for (key, value) in limit_subdict.items():
                    print('key =', key)
                    if 'zeta_ph_min' in key:
                    # # if 'CAP_sp' in key:
                    # # if 'zeta_sp' in key:
                    # if 'CAP_sp_min' in key:
                    # if 'zeta_sp_min' in key:
                        
                        # _mask = df_cat['sigma_fcs'].values == 1500.
                        _mask = df_cat['sigma_fcs'].values == 4000.
                        ax.plot(df_cat['fcs_loc'].values[_mask], value[_mask], color=colors[row_idx])
                        
                # # Lateral
                # for (key, value) in limit_subdict.items():
                #     print('key =', key)
                #     # if 'tau_roll_max' in key:
                #     # if 'T_2_sl_min' in key:
                #     if 'dr' in key:
                        
                #         _mask = df_cat['sigma_fcs'].values == 1500.
                #         # _mask = df_cat['sigma_fcs'].values == 4000.
                #         ax.plot(df_cat['fcs_loc'].values[_mask], value[_mask], color=colors[row_idx])
                        
            # if cat == 'wing':
                
            #     # Longitudinal
            #     for (key, value) in limit_subdict.items():
            #         print('key =', key)
            #         # if 'zeta_ph_min' in key:
            #         # # if 'CAP_sp' in key:
            #         # # if 'zeta_sp' in key:
            #         if 'CAP_sp_min' in key:
            #         # if 'zeta_sp_min' in key:
                        
            #             _mask = df_cat['sigma_fcs'].values == 1500.
            #             # _mask = df_cat['sigma_fcs'].values == 4000.
            #             ax.plot(df_cat['span_loc'].values[_mask], value[_mask], color=colors[row_idx])
                        
                # # Lateral
                # for (key, value) in limit_subdict.items():
                #     print('key =', key)
                #     # if 'tau_roll_max' in key:
                #     # if 'T_2_sl_min' in key:
                #     if 'dr' in key:
                        
                #         _mask = df_cat['sigma_fcs'].values == 1500.
                #         # _mask = df_cat['sigma_fcs'].values == 4000.
                #         ax.plot(df_cat['span_loc'].values[_mask], value[_mask], color=colors[row_idx])
          
# ax.set_ylim(-200, 200)
ax.axhline(0)
plt.show()
            
# sys.exit()

#%%

# =============================================================================
def create_empty_dict(d):
    return {
        key: create_empty_dict(value) if isinstance(value, dict) else {}
        for key, value in d.items()
    }

limit_dict = create_empty_dict(res)

# Store the minimum finite y-value at each x-value across all row_idx.
# Separate minima are maintained for each (cat, branch, mode, limit_name).
min_limit_data = {}
min_row_data = {}

# Configure the four figures.
if modes_type == 'longitudinal':
    plot_specs = {
        'wing': {
            'conservative': {
                'limits': ('zeta_ph_min', 'CAP_sp_min', 'zeta_sp_min'),
                'labels': ('$\zeta_\mathrm{ph,min}$', '$\mathrm{CAP_{sp,min}}$', '$\zeta_\mathrm{sp,min}$'),
                'sigma_fcs': 1500,
                'x_col': 'span_loc',
                'xlabel': 'Position along half-span [-]',
            },
            'advanced': {
                'limits': ('zeta_ph_min', 'CAP_sp_min', 'zeta_sp_min'),
                'sigma_fcs': 4000.,
                'x_col': 'span_loc',
                'xlabel': 'Position along half-span [-]',
            },
        },
        'fuselage': {
            'conservative': {
                'limits': ('zeta_ph_min', 'CAP_sp_min', 'zeta_sp_min'),
                'sigma_fcs': 1500,
                'x_col': 'fcs_loc',
                'xlabel': 'Position along fuselage [-]',
            },
            'advanced': {
                'limits': ('zeta_ph_min', 'CAP_sp_min', 'zeta_sp_min'),
                'sigma_fcs': 4000.,
                'x_col': 'fcs_loc',
                'xlabel': 'Position along fuselage [-]',
            },
        },
    }
elif modes_type == 'lateral':
    plot_specs = {
        'wing': {
            'conservative': {
                # 'limits': ('tau_roll_max', 'T_2_sl_min', 'dr'),
                # 'limits': ('tau_roll_max', 'dr'),
                'limits': ('tau_roll_max', 'zeta_dr_min', 'zeta_omega_dr_min', 'omega_n_dr_min'),
                # 'labels': (r'$\tau_\mathrm{roll,max}$', '$T_{2_\mathrm{sl,min}}$', '$\zeta_\mathrm{dr,min}$', '$\zeta \omega_\mathrm{dr,min}$', '$\omega_{n_\mathrm{dr,min}}$'),
                'labels': (r'$\tau_\mathrm{roll,max}$', '$\zeta_\mathrm{dr,min}$', '$\zeta \omega_\mathrm{dr,min}$', '$\omega_{n_\mathrm{dr,min}}$'),
                'sigma_fcs': 1500,
                'x_col': 'span_loc',
                'xlabel': 'Position along half-span [-]',
            },
            'advanced': {
                # 'limits': ('tau_roll_max', 'T_2_sl_min', 'dr'),
                # 'limits': ('tau_roll_max', 'dr'),
                'limits': ('tau_roll_max', 'zeta_dr_min', 'zeta_omega_dr_min', 'omega_n_dr_min'),
                'sigma_fcs': 4000.,
                'x_col': 'span_loc',
                'xlabel': 'Position along half-span [-]',
            },
        },
        'fuselage': {
            'conservative': {
                # 'limits': ('tau_roll_max', 'T_2_sl_min', 'dr'),
                # 'limits': ('tau_roll_max', 'dr'),
                'limits': ('tau_roll_max', 'zeta_dr_min', 'zeta_omega_dr_min', 'omega_n_dr_min'),
                'sigma_fcs': 1500,
                'x_col': 'fcs_loc',
                'xlabel': 'Position along fuselage [-]',
            },
            'advanced': {
                # 'limits': ('tau_roll_max', 'T_2_sl_min', 'dr'),
                # 'limits': ('tau_roll_max', 'dr'),
                'limits': ('tau_roll_max', 'zeta_dr_min', 'zeta_omega_dr_min', 'omega_n_dr_min'),
                'sigma_fcs': 4000.,
                'x_col': 'fcs_loc',
                'xlabel': 'Position along fuselage [-]',
            },
        },
    }

# Create the four figures.
# =============================================================================
fig, axes = plt.subplots(1, 4, figsize=(14, 4), sharey=True)
# =============================================================================
figures = {
    (cat, branch): axes[idx]
    for idx, (cat, branch) in enumerate(
        (cat, branch)
        for cat, branches in plot_specs.items()
        for branch in branches
    )
}
# =============================================================================

for row_idx, _ax in enumerate(axes):
    rid = row_idx + 1

    for j, cat in enumerate(res.keys()):  # NILS (04.10.2026, fix B4)

        for k, mode in enumerate(modes_sublist):

            if len(res[cat][mode][rid]) == 0:  # NILS (04.10.2026): category not present in the selected studies
                continue

            # =============================================================================
            limit_subdict = limit_dict[cat][mode][rid]
            # =============================================================================

            df_cat = res[cat][mode][rid]

            (
                ev_valid, ev_invalid,
                sigma_fcs_valid, sigma_fcs_invalid,
                fcs_loc_valid, fcs_loc_invalid,
                span_loc_valid, span_loc_invalid,
                # # =============================================================================
                # limit_dict,
                # # =============================================================================
            ) = \
            mask_mil_limits(
                mode,
                df_cat,
                rid,  # NILS (03.10.2026, fix A10)
                # =============================================================================
                limit_subdict,
                # =============================================================================
            )

            # print('limit_dict =', limit_dict)

            if cat not in plot_specs:
                continue

            # Process longitudinal and lateral limits using the same logic.
            for branch, spec in plot_specs[cat].items():

                _mask = df_cat['sigma_fcs'].values == spec['sigma_fcs']
                x_values_all = df_cat[spec['x_col']].values[_mask]

                for key, value in limit_subdict.items():

                    # Identify the limit to which this key belongs.
                    limit_name = next(
                        (
                            limit
                            for limit in spec['limits']
                            if limit in key
                        ),
                        None,
                    )

                    if limit_name is None:
                        continue

                    print('key =', key)

                    y_values = np.asarray(value)[_mask]

                    # Retain only pairs for which both x and y are finite.
                    _finite = (
                        np.isfinite(x_values_all)
                        & np.isfinite(y_values)
                    )

                    x_values = x_values_all[_finite]
                    y_values = y_values[_finite]

                    # Accumulate the minimum y-value at each x-value.
                    # The aggregation spans all row_idx values and all
                    # matching inner/outer branches.
                    group = (cat, branch, mode, limit_name)
                    min_by_x = min_limit_data.setdefault(group, {})
                    row_by_x = min_row_data.setdefault(group, {})

                    for x_value, y_value in zip(x_values, y_values):
                        x_value = float(x_value)
                        y_value = float(y_value)

                        if (
                            x_value not in min_by_x
                            or y_value < min_by_x[x_value]
                        ):
                            min_by_x[x_value] = y_value
                            row_by_x[x_value] = row_idx


# =============================================================================
# fig, axes = plt.subplots(1, 4, figsize=(14, 5))

# Plot the minimum envelopes after processing all row_idx values.
for i, ((cat, branch, mode, limit_name), min_by_x) in enumerate(min_limit_data.items()):

    if not min_by_x:
        continue

    # Sort the x-values to ensure the curve is plotted in ascending order.
    x_plot = np.array(sorted(min_by_x))
    y_plot = np.array([min_by_x[x_value] for x_value in x_plot])

    ax = figures[(cat, branch)]

    row_plot = np.array([
        min_row_data[(cat, branch, mode, limit_name)][x_value]
        for x_value in x_plot
    ])

    # Break the line whenever the critical row_idx changes.
    switch_idx = np.flatnonzero(row_plot[:-1] != row_plot[1:])

    x_line = []
    y_line = []

    for j, (x_value, y_value) in enumerate(zip(x_plot, y_plot)):
        x_line.append(x_value)
        y_line.append(y_value)

        if j < len(x_plot) - 1 and row_plot[j] != row_plot[j + 1]:
            x_line.append(np.nan)
            y_line.append(np.nan)

    line, = ax.plot(
        x_line,
        y_line,
        label=plot_specs[cat][branch].get(
            'labels', plot_specs['wing']['conservative']['labels']
        )[plot_specs[cat][branch]['limits'].index(limit_name)],
    )

    # Mark both endpoints of every line break.
    if switch_idx.size:
        marker_idx = np.unique(np.concatenate((switch_idx, switch_idx + 1)))

        ax.scatter(
            x_plot[marker_idx],
            y_plot[marker_idx],
            marker='o',
            s=36,
            facecolors='white',
            edgecolors=line.get_color(),
            linewidths=1.2,
            zorder=line.get_zorder() + 1,
            clip_on=False,
        )

lower_y_lim_list = []

# Finalise all four figures.
for _i, ((cat, branch), ax) in enumerate(figures.items()):

    ax.set_title(f'{branch.capitalize()} system in {cat}', pad=15)
    ax.set_xlabel(plot_specs[cat][branch]['xlabel'], labelpad=15)
    
    # =============================================================================
    # Shade intervals where every curve is above zero.
    curves = []
    for line in ax.get_lines():
        x = np.asarray(line.get_xdata(), dtype=float)
        y = np.asarray(line.get_ydata(), dtype=float)
    
        _finite = np.isfinite(x) & np.isfinite(y)
        x, y = x[_finite], y[_finite]
    
        _order = np.argsort(x)
        x, y = x[_order], y[_order]
    
        if len(x) >= 2:
            curves.append((x, y))
    
    if curves:
        # Restrict to the x-range shared by every curve.
        x_min = max(x[0] for x, y in curves)
        x_max = min(x[-1] for x, y in curves)
    
        if x_min < x_max:
            x_breaks = [x_min, x_max]
    
            for x, y in curves:
                x_breaks.extend(x[(x > x_min) & (x < x_max)])
    
                # Include zero crossings so that each interval has a
                # consistent sign for every curve.
                _cross = y[:-1] * y[1:] < 0
                x_cross = (
                    x[:-1][_cross]
                    - y[:-1][_cross]
                    * (x[1:][_cross] - x[:-1][_cross])
                    / (y[1:][_cross] - y[:-1][_cross])
                )
                x_breaks.extend(
                    x_cross[(x_cross > x_min) & (x_cross < x_max)]
                )
    
            x_breaks = np.unique(x_breaks)
            interval_start = None
    
            for x0, x1 in zip(x_breaks[:-1], x_breaks[1:]):
                x_mid = (x0 + x1) / 2
                _positive = all(
                    np.interp(x_mid, x, y) > 0
                    for x, y in curves
                )
    
                if _positive and interval_start is None:
                    interval_start = x0
                elif not _positive and interval_start is not None:
                    ax.axvspan(
                        interval_start, x0,
                        # color='green', alpha=0.15, zorder=0,
                        hatch='///', facecolor='None', edgecolor='lightgrey',
                        clip_on=False,
                    )
                    interval_start = None
    
            if interval_start is not None:
                ax.axvspan(
                    interval_start, x_breaks[-1],
                    # color='green', alpha=0.15, zorder=0,
                    hatch='///', facecolor='None', edgecolor='lightgrey',
                    clip_on=False,
                )
    # =============================================================================
    
    ax.axhline(0, color='k')
    ax.set_xlim(0, 1)
    ###
    ax.set_xticks([0.00, 0.25, 0.50, 0.75, 1.00])
    ax.xaxis.set_major_formatter(
        FuncFormatter(
            lambda x, pos:
                '0' if np.isclose(x, 0.0) else
                '1' if np.isclose(x, 1.0) else
                f'{x:.2f}'
        )
    )
    ax.spines[['right', 'top', 'bottom']].set_visible(False)
    ax.tick_params(axis='y', which='both', right=False, length=0)
    ax.tick_params(axis='x', which='both', length=0)
    ###

    if (_i == 0 and ax.get_legend_handles_labels()[0]):
        ax.legend(frameon=False, fontsize=14)
    # if (_i == 1 and ax.get_legend_handles_labels()[0]):
    #     ax.legend(frameon=False, fontsize=14, loc='lower left', bbox_to_anchor=(0, -0.045), ncol=2)
        
    lower_y_lim_list.append(ax.get_ylim()[0])
        
for _j, ax in enumerate(axes):
    # ax.fill_between(x=[0, 1], y1=np.nanmin(lower_y_lim_list), y2=0, hatch='///', facecolor='None', edgecolor='lightgrey')
    ax.set_ylim(bottom=np.nanmin(lower_y_lim_list))
    if _j == 0:
        ax.set_ylabel('$(x-x_\mathrm{min})/x_\mathrm{min}$ or\n$(x_\mathrm{max}-x)/x_\mathrm{max}$', labelpad=15)
        
plt.show()
# =============================================================================
