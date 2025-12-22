"""
This script post-processes the result of the dynamic stability
analysis conducted using the TASOPT-SUAVE-AVL wrapper.
"""

__all__ = []

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from ambiance import Atmosphere
from matplotlib import gridspec

from matplotlib_custom_settings import *

# @NILS: could the complex eigenvalues at low Mach numbers have to do
# witht the wing being stalled there?

# mach_vals = np.array([0.05, 0.15, 0.25, 0.45, 0.65, 0.85])
# aoa_vals = np.array([-2., 0., 2., 5., 7., 10.])

# NOTE: cbar is the MAC (see page 98 in AE3202 Flight Dynamics Lecture Notes)

# Table 3.9 (phugoid)

cbar = 4.2350  # from "C:\Users\nmb48\avl_files\body_axis_derivatives_case_01_01.txt"
V = 0.05 * Atmosphere(11e3).speed_of_sound[0]
n = 1

def zeta(lamda):
    return lamda.real / np.sqrt(lamda.real**2 + lamda.imag**2)  # page 129 in AE3202 Flight Dynamics Lecture Notes

def T0p5(lamda, V, cbar):
    return np.log(0.5) / lamda.real * cbar / V  # pages 125 and 127 in AE3202 Flight Dynamics Lecture Notes

def T2(lamda):
    return -T0p5(lamda)  # page 127 in AE3202 Flight Dynamics Lecture Notes

def omega0(lamda, V, cbar):
    return np.sqrt(lamda.real**2 + lamda.imag**2) * V / cbar  # page 129 in AE3202 Flight Dynamics Lecture Notes

def CAP(lamda, V, cbar, n, alpha):  # (3.111) in medium_PhDthesis_2013_flightdynconstrconcacmdiscanalsopt_morris
    return omega0(lamda, V, cbar)**2 / (n / alpha)

def tau(lamda, V, cbar):
    return -np.log(0.5) * T0p5(lamda, V, cbar)

def compute_lambdas_from_zeta_omega(zeta0, omega0, V, cbar):

    r0 = omega0 * cbar / V
    x = zeta0 * r0
    y_mag = r0 * np.sqrt(max(0.0, 1.0 - zeta0**2))

    lam_plus = x + 1j * y_mag
    lam_minus = x - 1j * y_mag
    return lam_plus, lam_minus, r0

df = pd.read_csv(r"C:\Users\nmb48\base_dynamic_stability_data.txt", sep=r"\s+")
re_cols = df.filter(regex="_Re")
im_cols = df.filter(regex="_Im")

re = df[[f"Long_Re{i}" for i in range(1, 5)]].to_numpy()
im = df[[f"Long_Im{i}" for i in range(1, 5)]].to_numpy()
long_complex = re + 1j * im  # shape (n_rows, 4)

re = df[[f"Lat_Re{i}" for i in range(1, 5)]].to_numpy()
im = df[[f"Lat_Im{i}" for i in range(1, 5)]].to_numpy()
lat_complex = re + 1j * im  # shape (n_rows, 4)

def pick(arr, ind):
    idx = ind.to_numpy().astype(int)  # convert float → int
    return arr[np.arange(len(arr)), idx]

# Longitudinal modes
shortPeriod_data = pick(long_complex, df["shortPeriodInd"])
phugoid_data = pick(long_complex, df["phugoidInd"])

# Lateral modes
rollSubsistence_data = pick(lat_complex, df["rollSubsistenceInd"])
spiral_data = pick(lat_complex, df["spiralInd"])
dutchRoll_data = pick(lat_complex, df["dutchRollInd"])

data_list = [shortPeriod_data, phugoid_data, rollSubsistence_data, spiral_data, dutchRoll_data]

aoa_vals_flat = df['AoA'].to_numpy()
mach_vals_flat = df['Mach'].to_numpy()
beta_vals_flat = df['Beta'].to_numpy()
h_vals_flat = df['h'].to_numpy()

# aoa_vals = np.array([-2.0, 0.0, 5.0])
# mach_vals = np.array([0.05, 0.45, 0.85]) 
# beta_vals = np.array([-5.0, 0.0, 10.0])
# h_vals = np.linspace(0, 11e3, 3)
aoa_vals = np.unique(aoa_vals_flat)
mach_vals = np.unique(mach_vals_flat) 
beta_vals = np.unique(beta_vals_flat)
h_vals = np.unique(h_vals_flat)

aoa_targets = aoa_vals
mach_targets = mach_vals
beta_targets = beta_vals
h_targets = h_vals

MACH, AOA, BETA, H = np.meshgrid(mach_vals, aoa_vals, beta_vals, h_vals)

#%% Calculate MIL-STD-1797 requirements (section 3.4.2 in medium_PhDthesis_2013_flightdynconstrconcacmdiscanalsopt_morris)

def apply_mil_limits(aoa, mach, beta, h):
    
    amb = Atmosphere(h)
    V = mach * amb.speed_of_sound[0]

    N_points = 50
    
    # Lateral stability limits
    zeta_dr_min = 0.08  # 3rd row of Table 3.14 (Class II, Category B, Level I)
    omega_n_dr_min = 0.4  # 3rd row of Table 3.14 (Class II, Category B, Level I)
    T_2_sl_min = 20  # 2nd row of Table 3.13 (Category B, Level I)
    tau_roll_max = 1.4  # 3rd row of Table 3.12 (Class II, Category B, Level I)
    
    # Dutch roll
    
    _, _, r_circle_dr_min = compute_lambdas_from_zeta_omega(zeta_dr_min, omega_n_dr_min, V, cbar)
    theta = np.linspace(0, 2 * np.pi, N_points)
    x_circle_dr_min = r_circle_dr_min * np.cos(theta)
    y_circle_dr_min = r_circle_dr_min * np.sin(theta)
    
    x_vert_dr_max = -zeta_dr_min * r_circle_dr_min * cbar / V
    
    a_dr = zeta_dr_min
    t_dr = np.linspace(0, r_circle_dr_min / a_dr, N_points)
    x_ray_dr_max = -t_dr * a_dr
    y_ray_dr_max = t_dr * np.sqrt(1 - a_dr**2)
    
    # Spiral
    x_vert_sl_max = -np.log(0.5) / T_2_sl_min * cbar / V
    
    # Roll
    x_vert_roll_min = -np.log(0.5)**2 / tau_roll_max * cbar / V
    
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
    print('omega_n_sp_min, omega_n_sp_max =', omega_n_sp_min, omega_n_sp_max)
    _, _, r_circle_sp_min = compute_lambdas_from_zeta_omega(zeta_sp_min, omega_n_sp_min, V, cbar)
    _, _, r_circle_sp_max = compute_lambdas_from_zeta_omega(zeta_sp_max, omega_n_sp_max, V, cbar)
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
        x_circle_dr_min, y_circle_dr_min, x_vert_dr_max, x_ray_dr_max, y_ray_dr_max,
        x_vert_sl_max,
        x_vert_roll_min,
        
        x_circle_sp_min, y_circle_sp_min, x_circle_sp_max, y_circle_sp_max, x_ray_sp_max, y_ray_sp_max,
        x_ray_ph_max, y_ray_ph_max,
    )

#%%

# aoa_ref = 5.0
aoa_ref = aoa_vals[1]
mach_ref = 0.45
# beta_ref = 0.0
beta_ref = beta_vals[1]
h_ref = 5500

(
    x_circle_dr_min, y_circle_dr_min, x_vert_dr_max, x_ray_dr_max, y_ray_dr_max,
    x_vert_sl_max,
    x_vert_roll_min,
    
    x_circle_sp_min, y_circle_sp_min, x_circle_sp_max, y_circle_sp_max, x_ray_sp_max, y_ray_sp_max,
    x_ray_ph_max, y_ray_ph_max,
) = \
apply_mil_limits(aoa_ref, mach_ref, beta_ref, h_ref)

#%%

fig = plt.figure(figsize=(15,9), constrained_layout=True)
gs = gridspec.GridSpec(2, 3, width_ratios=[1,1,1], hspace=0.2, height_ratios=[1,1], wspace=0.2)
axes = [fig.add_subplot(gs[i,j]) for i in range(2) for j in range(3)]

# =============================================================================
# Remove the upper-right subplot (row 0, col 2 → index 2)
axes[2].remove()

# Create an empty Axes in the same grid position
leg_ax = fig.add_subplot(gs[0,2])
leg_ax.axis("off")     # hide frame

# ---- your legend items here ----
# Example legend:
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([], [], marker='o', linestyle='', label='Low h'),
    Line2D([], [], marker='o', linestyle='', label='High h')
]

leg_ax.legend(handles=legend_elements, loc='center')
# =============================================================================

titles = ['Short period', 'Phugoid', 'Roll subsidence', 'Spiral', 'Dutch roll']

for i, ax in enumerate([axes[j] for j in range(len(axes)) if j != 2]):
    
    data = data_list[i]
    
    mask_aoa = np.where((mach_vals_flat == mach_ref) & (beta_vals_flat == beta_ref) & (h_vals_flat == h_ref))[0]
    filtered_data_aoa = data[mask_aoa]
    
    mask_mach = np.where((aoa_vals_flat == aoa_ref) & (beta_vals_flat == beta_ref) & (h_vals_flat == h_ref))[0]
    filtered_data_mach = data[mask_mach]
    
    mask_beta = np.where((aoa_vals_flat == aoa_ref) & (mach_vals_flat == mach_ref) & (h_vals_flat == h_ref))[0]
    filtered_data_beta = data[mask_beta]
    
    mask_h = np.where((aoa_vals_flat == aoa_ref) & (mach_vals_flat == mach_ref) & (beta_vals_flat == beta_ref))[0]
    filtered_data_h = data[mask_h]
    
    ax.plot(filtered_data_aoa.real, filtered_data_aoa.imag, marker='.', color=colors[0])
    ax.plot(filtered_data_aoa.real, -filtered_data_aoa.imag, marker='.', color=colors[0])
    
    ax.plot(filtered_data_mach.real, filtered_data_mach.imag, marker='.', color=colors[1])
    ax.plot(filtered_data_mach.real, -filtered_data_mach.imag, marker='.', color=colors[1])
    
    ax.plot(filtered_data_beta.real, filtered_data_beta.imag, marker='.', color=colors[2])
    ax.plot(filtered_data_beta.real, -filtered_data_beta.imag, marker='.', color=colors[2])
    
    ax.plot(filtered_data_h.real, filtered_data_h.imag, marker='.', color=colors[3])
    ax.plot(filtered_data_h.real, -filtered_data_h.imag, marker='.', color=colors[3])
    
    if i == 0:  # short period
    
        # # for aoa in aoa_vals:
        # for mach in mach_vals:
            
        #     (
        #         x_circle_dr_min, y_circle_dr_min, x_vert_dr_max, x_ray_dr_max, y_ray_dr_max,
        #         x_vert_sl_max,
        #         x_vert_roll_min,
                
        #         x_circle_sp_min, y_circle_sp_min, x_circle_sp_max, y_circle_sp_max, x_ray_sp_max, y_ray_sp_max,
        #         x_ray_ph_max, y_ray_ph_max,
        #     ) = \
        #     apply_mil_limits(aoa_ref, mach, beta_ref, h_ref)
        #     # apply_mil_limits(aoa, mach_ref, beta_ref, h_ref)
        
        ax.plot(x_circle_sp_min, y_circle_sp_min, color='black')
        ax.plot(x_circle_sp_max, y_circle_sp_max, color='black')
    
        ax.plot(x_ray_sp_max, y_ray_sp_max, '-', color='black')
        ax.plot(x_ray_sp_max, -y_ray_sp_max, '-', color='black')
        
        for aoa in aoa_vals:
            mask_aoa = np.where(aoa_vals_flat == aoa)[0]
            
            ax.scatter(data[mask_aoa].real, data[mask_aoa].imag)
        
    elif i == 1:  # phugoid
        
        ax.plot(x_ray_ph_max, y_ray_ph_max, '-', color='black')
        ax.plot(x_ray_ph_max, -y_ray_ph_max, '-', color='black')
        
    elif i == 2:  # roll-subsidence
    
        ax.axvline(x_vert_roll_min, color='black')
        
    elif i == 3:  # spiral
        
        ax.axvline(x_vert_sl_max, color='black')
        
    elif i == 4:  # dutch roll
        
        ax.plot(x_circle_dr_min, y_circle_dr_min, color='black')
        
        ax.axvline(x_vert_dr_max, color='black')
        
        ax.plot(x_ray_dr_max, y_ray_dr_max, '-', color='black')
        ax.plot(x_ray_dr_max, -y_ray_dr_max, '-', color='black')
        
    ax.set_title(titles[i])
    ax.spines[['right','top']].set_visible(False)
    ax.tick_params(axis='y', which='both', right=False, length=0)
    ax.tick_params(axis='x', which='both', length=0)
    ax.spines['bottom'].set_position(('data', 0.0))
    ax.spines['left'].set_position(('data', 0.0))

plt.show()

sys.exit()

#%%

fig, ax = plt.subplots(figsize=(8, 6))

# for aoa_target in aoa_targets:
#     for mach_target in mach_targets:
    
#         print('aoa_target, mach_target =', aoa_target, mach_target)
    
#         aoa_mask = np.where(aoa_vals_flat == aoa_target)[0]
#         mach_mask = np.where(mach_vals_flat == mach_target)[0]
#         common_mask = np.intersect1d(aoa_mask, mach_mask)
        
#         ax.scatter(rollSubsistence_data.real[common_mask], rollSubsistence_data.imag[common_mask], color=colors[0], marker='o', label='Roll')
#         # ax.scatter(rollSubsistence_data.real[common_mask], np.zeros(len(spiral_data.imag))[common_mask], color=colors[0], marker='o', label='Roll')
#         ax.scatter(spiral_data.real[common_mask], spiral_data.imag[common_mask], color=colors[1], marker='s', label='Spiral')
#         # ax.scatter(spiral_data.real, -spiral_data.imag, color=colors[1], marker='s')
#         # ax.scatter(spiral_data.real[common_mask], np.zeros(len(spiral_data.imag))[common_mask], color=colors[1], marker='s', label='Spiral')
#         ax.scatter(dutchRoll_data.real[common_mask], dutchRoll_data.imag[common_mask], color=colors[2], marker='d', label='Dutch roll')
#         ax.scatter(dutchRoll_data.real[common_mask], -dutchRoll_data.imag[common_mask], color=colors[2], marker='d')

# =============================================================================
# Dutch roll
ax.plot(circle_x, circle_y, label=f"omega0 = {omega0} -> r={r0:.3g}", color='red')
# ax.plot(ray_x, ray_y, '--', label=f"zeta = {_zeta0} (ray)", color='red')
ax.axvline(-_zeta0 * r0 * cbar / V, color='red')

# =============================================================================
a = _zeta0
t = np.linspace(0, r0 / a, 100)
x_ray = -t * a
y_ray = t * np.sqrt(1 - a**2)
ax.plot(x_ray, y_ray, '-', label=f"zeta = {_zeta0} (ray)", color='red')
ax.plot(x_ray, -y_ray, '-', label=f"zeta = {_zeta0} (ray)", color='red')
# =============================================================================

# Spiral
ax.axvline(-np.log(0.5) / _T2 * cbar / V, color='blue')  # allowed to be open-loop unstable
ax.fill_betweenx(y=ax.get_ylim(), x1=ax.get_xlim()[0], x2=-np.log(0.5) / _T2 * cbar / V, facecolor='blue', edgecolor='None', alpha=0.5)

# Roll
ax.axvline(-np.log(0.5)**2 / _tau * cbar / V, color='black')
ax.fill_betweenx(y=ax.get_ylim(), x1=ax.get_xlim()[0], x2=-np.log(0.5)**2 / _tau * cbar / V, facecolor='black', edgecolor='None', alpha=0.5)
# =============================================================================

# =============================================================================
theta = np.arccos(_zeta0)
# x_intersect_inner = _omega0 * _zeta0 * cbar / V
# x_intersect_outer = _omega0_max * _zeta0_max * cbar / V
r_outer = _omega0 * cbar / V

# Angular range for the wedge
angles = np.linspace(np.pi - theta, np.pi + theta, 200)

# Outer arc
x_outer = r_outer * np.cos(angles)
y_outer = r_outer * np.sin(angles)

x_ray_vert_intersect = -_zeta0 * r0 * cbar / V
y_ray_vert_intersect = r0 * cbar / V * np.sqrt(1 - _zeta0**2)

# Combine into a single polygon
x = np.concatenate([x_outer, [x_ray_vert_intersect, x_ray_vert_intersect]])
y = np.concatenate([y_outer, [-y_ray_vert_intersect, y_ray_vert_intersect]])

ax.fill(x, y, facecolor='red', edgecolor='None', alpha=0.5)
# =============================================================================

# # Plot either of the two independent inputs as z-value on colour map

# # common_mask = np.where(aoa_vals_flat == aoa_targets)[0]
# common_mask = np.where(mach_vals_flat == mach_targets)[0]
# z_values = aoa_vals_flat[common_mask]

# ax.scatter(rollSubsistence_data.real[common_mask], rollSubsistence_data.imag[common_mask], marker='o', label='Roll', c=z_values)  # , color=colors[0]
# # ax.scatter(rollSubsistence_data.real[common_mask], np.zeros(len(spiral_data.imag))[common_mask], marker='o', label='Roll', c=z_values)  # , color=colors[0]
# ax.scatter(spiral_data.real[common_mask], spiral_data.imag[common_mask], marker='s', label='Spiral', c=z_values)  # , color=colors[1]
# # ax.scatter(spiral_data.real, -spiral_data.imag, marker='s', c=z_values)  # , color=colors[1]
# # ax.scatter(spiral_data.real[common_mask], np.zeros(len(spiral_data.imag))[common_mask], marker='s', label='Spiral', c=z_values)  # , color=colors[1]
# ax.scatter(dutchRoll_data.real[common_mask], dutchRoll_data.imag[common_mask], marker='d', label='Dutch roll', color=colors[2])
# ax.scatter(dutchRoll_data.real[common_mask], -dutchRoll_data.imag[common_mask], marker='d', color=colors[2])

ax.spines[['right', 'top']].set_visible(False)
ax.tick_params(axis='y', which='both', right=False, length=0)
ax.tick_params(axis='x', which='both', length=0)
ax.spines['bottom'].set_position(('data', 0.0))
ax.spines['left'].set_position(('data', 0.0))

ax.axvspan(-np.inf, 0, color='lightgrey', alpha=1, zorder=1, clip_on=False)

xlim = ax.get_xlim()
ylim = ax.get_ylim()
# ax.fill_between(x=[xlim[0], 0], y1=ylim[0], y2=ylim[1], color='lightgrey', alpha=1, zorder=-1)
ax.set_xlim(xlim)
ax.set_ylim(ylim)

# ax.legend(frameon=False, loc='upper left')

# ax.set_aspect('equal')

plt.show()

# sys.exit()

#%% --- LONGITUDINAL modes ---

_zeta0 = 0.04  # 1st row in Table 3.9 (Class II, Category B, Level I)
_zeta0_min = 0.3  # 1st row, 4th column in Table 3.10 (Class II, Category B, Level I)
_zeta0_max = 2.0  # 1st row, 5th column in Table 3.10 (Class II, Category B, Level I)
_CAP_min = 0.085  # 1st row, 4th column in Table 3.11 (Class II, Category B, Level I)
_CAP_max = 3.6  # 1st row, 5th column in Table 3.11 (Class II, Category B, Level I)

_omega0_min = np.sqrt(_CAP_min * n / aoa_vals[2])  # (3.111) in medium_PhDthesis_2013_flightdynconstrconcacmdiscanalsopt_morris rearranged
_omega0_max = np.sqrt(_CAP_max * n / aoa_vals[2])  # (3.111) in medium_PhDthesis_2013_flightdynconstrconcacmdiscanalsopt_morris rearranged

_, _, r0_min = compute_lambdas_from_zeta_omega(_zeta0_min, _omega0_min, V, cbar)
_, _, r0_max = compute_lambdas_from_zeta_omega(_zeta0_max, _omega0_max, V, cbar)

# --- optional plot of the locus (circle for omega0 and rays for zeta)
theta = np.linspace(0, 2*np.pi, 400)
circle_x_min = r0_min * np.cos(theta)
circle_y_min = r0_min * np.sin(theta)
circle_x_max = r0_max * np.cos(theta)
circle_y_max = r0_max * np.sin(theta)

# ray line (plot as a line through origin with angle arccos(zeta0))
tmax = 2 * r0_min
tmax_min = 2 * r0_min
tmax_max = 2 * r0_max

t = np.linspace(0, tmax, 400)
ray_x = -t * _zeta0
ray_y = t * np.sqrt(max(0.0, 1 - _zeta0**2))
# ray_y = np.sqrt((ray_x / _zeta0)**2 - ray_x)
t_min = np.linspace(-tmax_min, tmax_min, 400)
ray_x_min = -t_min * _zeta0_min
ray_y_min = t_min * np.sqrt(max(0.0, 1 - _zeta0_min**2))
# ray_y_min = np.sqrt((ray_x_min / _zeta0_min)**2 - ray_x_min)
t_max = np.linspace(-tmax_max, tmax_max, 400)
ray_x_max = -t_max * _zeta0_max
ray_y_max = t_max * np.sqrt(max(0.0, 1 - _zeta0_max**2))
# ray_y_max = np.sqrt((ray_x_max / _zeta0_max)**2 - ray_x_max)

fig, ax = plt.subplots(figsize=(8, 6))

# for aoa_target in aoa_targets:
#     for mach_target in mach_targets:
    
#         print('aoa_target, mach_target =', aoa_target, mach_target)
    
#         aoa_mask = np.where(aoa_vals_flat == aoa_target)[0]
#         mach_mask = np.where(mach_vals_flat == mach_target)[0]
#         common_mask = np.intersect1d(aoa_mask, mach_mask)
        
#         ax.scatter(phugoid_data.real[common_mask], phugoid_data.imag[common_mask], color=colors[0], marker='d', label='Phugoid')
#         ax.scatter(phugoid_data.real[common_mask], -phugoid_data.imag[common_mask], color=colors[0], marker='d')
        
#         ax.scatter(shortPeriod_data.real[common_mask], shortPeriod_data.imag[common_mask], color=colors[1], marker='d', label='Short period')
#         ax.scatter(shortPeriod_data.real[common_mask], -shortPeriod_data.imag[common_mask], color=colors[1], marker='d')

# =============================================================================
# Short period
ax.plot(circle_x_min, circle_y_min, label=f"omega0 = {omega0} -> r={r0:.3g}", color='red')
# ax.plot(ray_x_min, ray_y_min, '--', label=f"zeta = {_zeta0} (ray)", color='red')
ax.plot(circle_x_max, circle_y_max, label=f"omega0 = {omega0} -> r={r0:.3g}", color='red')
# ax.plot(ray_x_max, ray_y_max, '--', label=f"zeta = {_zeta0} (ray)", color='red')

# =============================================================================
a = _zeta0_min
t = np.linspace(0, r0_max / a, 100)
x_ray_min = -t * a
y_ray_min = t * np.sqrt(1 - a**2)
ax.plot(x_ray_min, y_ray_min, '-', label=f"zeta = {_zeta0} (ray)", color='red')
ax.plot(x_ray_min, -y_ray_min, '-', label=f"zeta = {_zeta0} (ray)", color='red')
# =============================================================================

# circle_min = plt.Circle((0, 0), r0_min, facecolor='lightgrey', edgecolor='black', hatch='////')
# circle_max = plt.Circle((0, 0), r0_max, facecolor='None', edgecolor='black', hatch='\\\\\\')
# ax.add_patch(circle_min)
# ax.add_patch(circle_max)

# # Phugoid
# ax.plot(ray_x, ray_y, '--', label=f"zeta = {_zeta0} (ray)", color='blue')
# ax.fill_between(x=ray_x, y1=ray_y, y2=ax.get_ylim()[1] * np.ones_like(ray_x), color='lightgrey', alpha=1, zorder=-1)
# =============================================================================

# =============================================================================
a = _zeta0
t = np.linspace(0, ax.get_ylim()[1] / np.sqrt(1 - a**2), 100)
x_ray = -t * a
y_ray = t * np.sqrt(1 - a**2)
ax.plot(x_ray, y_ray, '-', label=f"zeta = {_zeta0} (ray)", color='blue')
ax.plot(x_ray, -y_ray, '-', label=f"zeta = {_zeta0} (ray)", color='blue')

x_stacked = np.hstack((np.flip(x_ray), x_ray))
y_stacked = np.hstack((np.flip(y_ray), -y_ray))
ax.fill_betweenx(y=y_stacked, x1=np.ones_like(y_stacked) * ax.get_xlim()[0], x2=x_stacked, facecolor='blue', edgecolor='None', alpha=0.5)

# =============================================================================

# =============================================================================
theta = np.arccos(_zeta0_min)
# x_intersect_inner = _omega0_min * _zeta0_min * cbar / V
# x_intersect_outer = _omega0_max * _zeta0_max * cbar / V
r_inner = _omega0_min * cbar / V
r_outer = _omega0_max * cbar / V

# Angular range for the wedge
angles = np.linspace(np.pi - theta, np.pi + theta, 200)

# Outer arc
x_outer = r_outer * np.cos(angles)
y_outer = r_outer * np.sin(angles)

# Inner arc (reverse direction so polygon closes properly)
x_inner = r_inner * np.cos(angles[::-1])
y_inner = r_inner * np.sin(angles[::-1])

# Combine into a single polygon
x = np.concatenate([x_outer, x_inner])
y = np.concatenate([y_outer, y_inner])

ax.fill(x, y, facecolor='red', edgecolor='None', alpha=0.5)
# =============================================================================

ax.set_xlabel('Re')
ax.xaxis.set_label_coords(1.04, 0.515)
ax.set_ylabel('Im', rotation=0)
ax.yaxis.set_label_coords(0.52, 1.02)
ax.spines[['right', 'top']].set_visible(False)
ax.tick_params(axis='y', which='both', right=False, length=0)
ax.tick_params(axis='x', which='both', length=0)
ax.spines['bottom'].set_position(('data', 0.0))
ax.spines['left'].set_position(('data', 0.0))

ax.axvspan(-np.inf, 0, color='lightgrey', alpha=1, zorder=1, clip_on=False)

xlim = ax.get_xlim()
ylim = ax.get_ylim()
# ax.fill_between(x=[xlim[0], 0], y1=ylim[0], y2=ylim[1], color='lightgrey', alpha=1, zorder=-1)
ax.set_xlim(xlim)
ax.set_ylim(ylim)

# ax.legend(frameon=False, loc='upper left')

plt.show()

sys.exit()

#%%

A_long = df['Long_A'].to_numpy()
B_long = df['Long_B'].to_numpy()
C_long = df['Long_C'].to_numpy()
D_long = df['Long_D'].to_numpy()
E_long = df['Long_E'].to_numpy()

A_lat = df['Lat_A'].to_numpy()
B_lat = df['Lat_B'].to_numpy()
C_lat = df['Lat_C'].to_numpy()
D_lat = df['Lat_D'].to_numpy()
E_lat = df['Lat_E'].to_numpy()

R_long = B_long * C_long * D_long - A_long * D_long**2 - B_long**2 * E_long
R_lat = B_lat * C_lat * D_lat - A_lat * D_lat**2 - B_lat**2 * E_lat

fig, ax = plt.subplots()
ax.scatter(E_long, R_long)
plt.show()

#%%

fig, ax = plt.subplots()
ax.scatter(E_lat, R_lat)
plt.show()
