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

def T0p5(lamda, V, cbar):  # {s}
    return -np.log(0.5) / lamda.real  # pages 125 and 127 in AE3202 Flight Dynamics Lecture Notes

def T2(lamda):  # (s)
    return -T0p5(lamda)  # page 127 in AE3202 Flight Dynamics Lecture Notes

def omega0(lamda, V, cbar):  # (1/s)
    return np.sqrt(lamda.real**2 + lamda.imag**2)  # page 129 in AE3202 Flight Dynamics Lecture Notes

def CAP(lamda, V, cbar, n, alpha):  # (1/s^2)
    return omega0(lamda, V, cbar)**2 / (n / alpha)  # (3.111) in medium_PhDthesis_2013_flightdynconstrconcacmdiscanalsopt_morris

def tau(lamda, V, cbar):  # (s)
    return -np.log(0.5) * T0p5(lamda, V, cbar)  # page 125 in AE3202 Flight Dynamics Lecture Notes

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
    i, j = map(int, match.groups())
    df_dict[i][j] = pd.read_csv(f, sep=r"\s+", engine="python")
    
#%% Calculate MIL-STD-1797 requirements (section 3.4.2 in medium_PhDthesis_2013_flightdynconstrconcacmdiscanalsopt_morris)

def apply_mil_limits(aoa, mach, beta, h, n, W):
    
    amb = Atmosphere(h)
    V = mach * amb.speed_of_sound[0]

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

#%%

# def classify_avl_modes_clean(eigvals):
#     # eigvals = np.linalg.eigvals(A)

#     tol = 1e-6

#     # # Remove pure kinematic integrators
#     # eigvals = [l for l in eigvals if abs(l) > tol]

#     used = set()
#     # modes = []
#     mode_dict = {}

#     for i, lam in enumerate(eigvals):
#         if i in used:
#             continue

#         if abs(np.imag(lam)) > tol:
#             # complex conjugate pair
#             for j, lam2 in enumerate(eigvals):
#                 if j != i and abs(lam2 - np.conj(lam)) < tol:
#                     used.update({i, j})
#                     pair = [lam, lam2]
#                     break

#             wn = np.sqrt(np.real(lam)**2 + np.imag(lam)**2)

#             if wn > 1.0:
#                 mode = 'short_period'
#             # elif wn > 0.3:
#             elif wn > 0.1:
#                 mode = 'dutch_roll'
#             else:
#                 mode = 'phugoid'

#         else:
#             # real root
#             used.add(i)
#             pair = [lam]

#             if abs(np.real(lam)) > 1.0:
#                 mode = 'roll_subsidence'
#             else:
#                 mode = 'spiral'

#         # modes.append((mode, pair))
#         mode_dict[mode] = pair

#     return mode_dict


# A = np.array(case_res.stability.system_matrix)[:, :12]
# modes = classify_avl_modes_clean(A)

# for m in modes:
#     print(m[0], m[1])

#%%

modes_type = 'longitudinal'  # 'longitudinal', 'lateral'

fig = plt.figure(figsize=(15,9), constrained_layout=True)
gs = gridspec.GridSpec(2, 3, width_ratios=[1,1,1], hspace=0.2, height_ratios=[1,1], wspace=0.3)
axes = [fig.add_subplot(gs[i,j]) for i in range(2) for j in range(3)]

# titles = ['Short period', 'Phugoid', 'Roll subsidence', 'Spiral', 'Dutch roll']
titles = ['1. Take-off', '2. Climb', '3. Beginning of cruise', '4. End of cruise', '5. Descent', '6. Landing']

for i, ax in enumerate(axes):
    
    # if i == 0:
    #     ax.set_xlabel('Re')
    #     ax.xaxis.set_label_coords(1.04, 0.515)
    #     ax.set_ylabel('Im', rotation=0)
    #     ax.yaxis.set_label_coords(0.52, 1.02)
    
    if modes_type == 'longitudinal':
        markers=['o', 'd']
        legend_marker_colours = [colors[0], colors[1], colors[2], colors[0], colors[1], colors[2]]
        legend_markers = ['o', 'o', 'o', 'd', 'd', 'd']
        legend_marker_labels = [
            'Short-period roots nacelle', 'Short-period roots wing', 'Short-period roots fuselage',
            'Phugoid roots nacelle', 'Phugoid roots wing', 'Phugoid roots fuselage',
        ]
        
        ax.set_xlim((-1.5, 1))
        ax.set_ylim((-0.1, 5))
        ax.set_xticks([-1, -0.5])
        # ax.set_yticks([-5, -4, -3, -2, -1, 1, 2, 3, 4, 5])
        ax.set_yticks([1, 2, 3, 4, 5])
        
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
        # axins.spines['bottom'].set_position(('data', 0.0))
        # axins.spines['left'].set_position(('data', 0.0))
        mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")
        
    elif modes_type == 'lateral':    
        markers=['o', 'd', '^']
        
        ax.set_xlim((-1.5, 0.1))
        ax.set_ylim((-0.1, 2.5))
        ax.set_xticks([-1, -0.5])
        ax.set_yticks([1, 2])
        
        y_lim_lat = 3
        
    # =============================================================================
    if modes_type == 'lateral':
        # ax_left = ax.inset_axes([0.0, 0.0, 0.25, 1.0], sharey=ax)
        gap = 0.03  # fraction of ax width
        ax_left = ax.inset_axes(
            [-0.2 - gap, 0.0, 0.2, 1.0],  # <-- shifted left
            sharey=ax
        )
        ax_left.set_xlim(-65, -4)
        ax_left.set_ylim(-0.1, 2.5)
        # ax_left.spines['right'].set_visible(False)
        # ax.spines['left'].set_visible(False)
        # ax_left.tick_params(labelright=False)
        
        ax_left.spines[['left', 'right','top']].set_visible(False)
        ax_left.tick_params(axis='y', which='both', right=False, length=0)
        ax_left.tick_params(axis='x', which='both', length=0)
        ax_left.spines['bottom'].set_position(('data', 0.0))
        ax_left.tick_params(labelright=False)
        ax_left.tick_params(labelleft=False)
    
        # d = .015
        # kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
        # ax.plot((-d, +d), (-d, +d), **kwargs)
        # # ax.plot((-d, +d), (1-d, 1+d), **kwargs)
    
        # kwargs = dict(transform=ax_left.transAxes, color='k', clip_on=False)
        # ax_left.plot((1-d, 1+d), (-d, +d), **kwargs)
        # # ax_left.plot((1-d, 1+d), (1-d, 1+d), **kwargs)
        
        # # --- draw x-axis break marks at y = 0 (DATA coordinates) ---
        # d = 0.1  # size of diagonal in data units
        
        # ax.plot([-3-d, -3+d], [-d, d], color='k', clip_on=False)
        # ax_left.plot([-5-d, -5+d], [-d, d], color='k', clip_on=False)
    # =============================================================================
    
    _row = df_dict[3][1].iloc[i]
    
    aoa = np.deg2rad(_row["AoA"])
    mach = _row["Mach"]
    beta = np.deg2rad(_row["Beta"])
    h = _row["h"]
    n = _row["n"]
    W = _row["W"]
        
    (
        x_circle_dr_min, y_circle_dr_min, x_vert_dr_max, x_ray_dr_max, y_ray_dr_max, r_circle_dr_min, zeta_dr_min,
        x_vert_sl_max,
        x_vert_roll_max,
        
        x_circle_sp_min, y_circle_sp_min, x_circle_sp_max, y_circle_sp_max, x_ray_sp_max, y_ray_sp_max, r_circle_sp_min, r_circle_sp_max, zeta_sp_min,
        x_ray_ph_max, y_ray_ph_max,
    ) = \
    apply_mil_limits(aoa, mach, beta, h, n, W)
    
    if modes_type == 'longitudinal':
        
        # Short period
        # ax.plot(x_circle_sp_min, y_circle_sp_min, color='black')
        # ax.plot(x_circle_sp_max, y_circle_sp_max, color='black')
        # ax.plot(x_ray_sp_max, y_ray_sp_max, '-', color='black')
        # ax.plot(x_ray_sp_max, -y_ray_sp_max, '-', color='black')
        
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
        ax.fill(x_wedge_sp, y_wedge_sp, facecolor='None', edgecolor='black', hatch='//', alpha=0.1)
        
        # Phugoid
        # ax.plot(x_ray_ph_max, y_ray_ph_max, '-', color='black')
        # ax.plot(x_ray_ph_max, -y_ray_ph_max, '-', color='black')
        
        x_stacked_ray_ph_max = np.hstack((np.flip(x_ray_ph_max), x_ray_ph_max))
        y_stacked_ray_ph_max = np.hstack((np.flip(y_ray_ph_max), -y_ray_ph_max))
        ax.fill_betweenx(
            y=y_stacked_ray_ph_max, x1=np.ones_like(x_stacked_ray_ph_max) * -3, x2=x_stacked_ray_ph_max,
            facecolor='lightgrey', edgecolor='None', alpha=0.5, zorder=-1,
        )
        axins.fill_betweenx(
            y=y_stacked_ray_ph_max, x1=np.ones_like(x_stacked_ray_ph_max) * -3, x2=x_stacked_ray_ph_max,
            facecolor='lightgrey', edgecolor='None', alpha=0.5, zorder=-1,
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
    
    real_data = []
    imag_data = []
    
    for j, df_sub_dict in enumerate(df_dict.values()):
        
        # =============================================================================
        # N_points = 5
        N_points = 10
        if j == 0:
            sigma_fcs_range = np.linspace(1.5e3, 4e3, N_points)
            # fcs_loc_range = np.linspace(0.01, 0.99, N_points)
            # x_repeated = np.repeat(sigma_fcs_range, N_points)
            # y_tiled = np.tile(fcs_loc_range,  N_points)
            x_repeated = sigma_fcs_range
        elif j == 1:
            sigma_fcs_range = np.linspace(1.5e3, 4e3, N_points)
            # span_loc_range = np.linspace(0.01, 0.99, N_points)
            # x_repeated = np.repeat(sigma_fcs_range, N_points)
            # y_tiled = np.tile(span_loc_range,  N_points)
            x_repeated = sigma_fcs_range
        elif j == 2:
            sigma_fcs_range = np.linspace(1.5e3, 4e3, N_points)
            x_repeated = sigma_fcs_range
        # =============================================================================
        
        ev_1_list = []
        ev_2_list = []
        ev_3_list = []
        ev_4_list = []
        ev_5_list = []
        
        for k, df in enumerate(df_sub_dict.values()):
            
            row = df.iloc[i]
            
            re_cols = row.filter(regex="_Re")
            im_cols = row.filter(regex="_Im")

            ev_re = row[[f"Re{_}" for _ in range(1, 9)]].to_numpy()
            ev_im = row[[f"Im{_}" for _ in range(1, 9)]].to_numpy()
            ev_cplx = ev_re + 1j * ev_im  # shape (n_rows, 4)
            
            # =============================================================================
            # mode_dict = classify_avl_modes_clean(ev_cplx)
            ax.scatter(ev_cplx.real, ev_cplx.imag, marker='.', color=colors[j])
            ax_left.scatter(ev_cplx.real, ev_cplx.imag, marker='.', color=colors[j])
            # =============================================================================
            
            edgecolor = 'k'
            linewidth = 0.5
            
            if modes_type == 'lateral':
                ax.scatter(ev_cplx.real[0], ev_cplx.imag[0], marker='o', color=colors[j], edgecolor=edgecolor, linewidth=linewidth, zorder=100)
                ax_left.scatter(ev_cplx.real[0], ev_cplx.imag[0], marker='o', color=colors[j], edgecolor=edgecolor, linewidth=linewidth, zorder=100)
            
                ax.scatter(ev_cplx.real[5], ev_cplx.imag[5], marker='^', color=colors[j], edgecolor=edgecolor, linewidth=linewidth, zorder=100)
                ax_left.scatter(ev_cplx.real[5], ev_cplx.imag[5], marker='^', color=colors[j], edgecolor=edgecolor, linewidth=linewidth, zorder=100)
            
            elif modes_type == 'longitudinal':
                ax.scatter(ev_cplx.real[6:8], ev_cplx.imag[6:8], marker='*', color=colors[j], edgecolor=edgecolor, linewidth=linewidth, zorder=100)
            
            ev_1 = ev_cplx[0]
            ev_2 = ev_cplx[1]  # positive conjugate ev
            ev_3 = ev_cplx[3]  # positive conjugate ev
            ev_4 = ev_cplx[5]
            ev_5 = ev_cplx[6:8]
            
            ev_1_list.append(ev_1)
            ev_2_list.append(ev_2)
            ev_3_list.append(ev_3)
            ev_4_list.append(ev_4)
            ev_5_list.append(ev_5)
            
            # =============================================================================
            
            # # =============================================================================
            # ax.scatter(
            #     ev_cplx.real, ev_cplx.imag, marker='.', c=x_repeated[k] * np.ones(8),
            #     vmin=x_repeated.min(),
            #     vmax=x_repeated.max(),
            #     cmap='viridis',
            # )
            # # print(ev_cplx.real, ev_cplx.imag, x_repeated[k] * np.ones(8))
            # # print()
            # # =============================================================================
                
            # if modes_type == 'longitudinal':
            #     mode_list = ['short_period', 'phugoid']
                
            #     for l, mode_name in enumerate(mode_list):
            #         try:
            #             mode_evs = mode_dict[mode_name]
                        
            #             for mode_ev in mode_evs:
            #                 ax.scatter(mode_ev.real, mode_ev.imag, color=colors[j], marker=markers[l], s=20)
                    
            #         except Exception as e:
            #             print(e)
            #             continue
            
            # elif modes_type == 'lateral':
            #     mode_list = ['roll_subsidence', 'spiral', 'dutch_roll']
                
            #     for l, mode_name in enumerate(mode_list):
            #         try:
            #             mode_evs = mode_dict[mode_name]
                        
            #             for mode_ev in mode_evs:
            #                 ax.scatter(mode_ev.real, mode_ev.imag, color=colors[j], marker=markers[l], s=20)
                            
            #         except Exception as e:
            #             print(e)
            #             continue
            
            # =============================================================================
            if modes_type == 'longitudinal':
                
                axins.scatter(ev_cplx.real[6:8], ev_cplx.imag[6:8], marker='*', color=colors[j], edgecolor=edgecolor, linewidth=linewidth, zorder=100)
                # print('ev_cplx.real[6:8], ev_cplx.imag[6:8] =', ev_cplx.real[6:8], ev_cplx.imag[6:8])
                
                
            # =============================================================================
            
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
        
        # Blue markers (nacelle) are never misordered
        if j == 0:
            if modes_type == 'lateral':
                ax.scatter(np.array(ev_2_list).real, np.array(ev_2_list).imag, marker='s', color=colors[j], edgecolor=edgecolor, linewidth=linewidth)
                ax.scatter(np.array(ev_2_list).real, -np.array(ev_2_list).imag, marker='s', color=colors[j], edgecolor=edgecolor, linewidth=linewidth)
            elif modes_type == 'longitudinal':
                ax.scatter(np.array(ev_3_list).real, np.array(ev_3_list).imag, marker='d', color=colors[j], edgecolor=edgecolor, linewidth=linewidth)
                ax.scatter(np.array(ev_3_list).real, -np.array(ev_3_list).imag, marker='d', color=colors[j], edgecolor=edgecolor, linewidth=linewidth)
        elif (j == 1 or j == 2):
            if modes_type == 'lateral':
                ax.scatter(ev_2_list_sorted.real, ev_2_list_sorted.imag, marker='s', color=colors[j], edgecolor=edgecolor, linewidth=linewidth)
                ax.scatter(ev_2_list_sorted.real, -ev_2_list_sorted.imag, marker='s', color=colors[j], edgecolor=edgecolor, linewidth=linewidth)
            elif modes_type == 'longitudinal':
                ax.scatter(ev_3_list_sorted.real, ev_3_list_sorted.imag, marker='d', color=colors[j], edgecolor=edgecolor, linewidth=linewidth)
                ax.scatter(ev_3_list_sorted.real, -ev_3_list_sorted.imag, marker='d', color=colors[j], edgecolor=edgecolor, linewidth=linewidth)
        
    ax.set_title(titles[i], pad=10)
    ax.spines[['right','top']].set_visible(False)
    ax.tick_params(axis='y', which='both', right=False, length=0)
    ax.tick_params(axis='x', which='both', length=0)
    ax.spines['bottom'].set_position(('data', 0.0))
    ax.spines['left'].set_position(('data', 0.0))
    

plt.show()


