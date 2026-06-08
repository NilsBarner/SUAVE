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

from matplotlib_custom_settings import *

#%% Static stability

folder = r"C:/Users/nmb48/"
files = glob.glob(os.path.join(
    folder,
    'suave_static_stability_outputs_[1-7]_*.txt'
))
txt_files = [
    f for f in files
    if 0 <= int(f.rsplit('_', 1)[1].split('.')[0]) <= 99
]

# Nested dict: data[i][j] = DataFrame
df_dict = defaultdict(dict)

pattern = re.compile(r"suave_static_stability_outputs_(\d+)_(\d+)\.txt$")

for f in txt_files:
    match = pattern.search(os.path.basename(f))
    _i, _j = map(int, match.groups())
    # df_dict[_i][_j] = pd.read_csv(f, sep=r"\s+", header=0)
    
    with open(f, "r") as fh:
        lines = fh.readlines()

    # Remove leading '#' from header if present
    if lines[0].lstrip().startswith("#"):
        lines[0] = lines[0].lstrip()[1:].lstrip()

    df_dict[_i][_j] = pd.read_csv(
        io.StringIO("".join(lines)),
        sep=r"\s+",
        header=0
    )
    
titles = ['1. Take-off', '2. Climb', '3. Beginning of cruise', '4. End of cruise', '5. Descent', '6. Landing']

#%%

cmap_cm = plt.get_cmap('Greys')
cmap_cn = plt.get_cmap('Blues')
cmap_cl = plt.get_cmap('Reds')

for i in df_dict.keys():

    df_subdict = df_dict[i]

    fig = plt.figure(figsize=(13, 9), constrained_layout=False)
    gs_outer = gridspec.GridSpec(
        2, 3, figure=fig,
        left=0.1,
        right=0.83,   # reserve space for colorbars
        bottom=0.1,
        top=0.9,
        hspace=0.25,
        wspace=0.3,
    )

    if i in (3, 4):
        Cm_matrix = np.array([df_subdict[k]['Cm_alpha'] for k in range(len(df_subdict))])
        Cn_matrix = np.array([df_subdict[k]['Cn_beta']  for k in range(len(df_subdict))])
        Cl_matrix = np.array([df_subdict[k]['Cl_beta']  for k in range(len(df_subdict))])

        Cm_all = Cm_matrix.ravel()
        Cn_all = Cn_matrix.ravel()
        Cl_all = Cl_matrix.ravel()

        Cm_vmin, Cm_vmax = np.nanmin(Cm_all), np.nanmax(Cm_all)
        Cn_vmin, Cn_vmax = np.nanmin(Cn_all), np.nanmax(Cn_all)
        Cl_vmin, Cl_vmax = np.nanmin(Cl_all), np.nanmax(Cl_all)

        Cm_norm = mcolors.Normalize(vmin=Cm_vmin, vmax=Cm_vmax)
        Cn_norm = mcolors.Normalize(vmin=Cn_vmin, vmax=Cn_vmax)
        Cl_norm = mcolors.Normalize(vmin=Cl_vmin, vmax=Cl_vmax)

        n_levels = 50
        Cm_levels = np.linspace(Cm_vmin, Cm_vmax, n_levels)
        Cn_levels = np.linspace(Cn_vmin, Cn_vmax, n_levels)
        Cl_levels = np.linspace(Cl_vmin, Cl_vmax, n_levels)

    cm_mappable = None
    cn_mappable = None
    cl_mappable = None

    for j in range(6):

        row, col = divmod(j, 3)

        if i == 3:
            y_array = np.array([df_subdict[k]['fcs_loc'] for k in range(len(df_subdict))])[:, j]
        elif i == 4:
            y_array = np.array([df_subdict[k]['span_loc'] for k in range(len(df_subdict))])[:, j]

        sigma_fcs_array = np.array([df_subdict[k]['sigma_fcs'] for k in range(len(df_subdict))])[:, j] / 1e3
        Cm_alpha_array  = np.array([df_subdict[k]['Cm_alpha'] for k in range(len(df_subdict))])[:, j]
        Cn_beta_array   = np.array([df_subdict[k]['Cn_beta'] for k in range(len(df_subdict))])[:, j]
        Cl_beta_array   = np.array([df_subdict[k]['Cl_beta'] for k in range(len(df_subdict))])[:, j]

        if i in (3, 4):
            sigma_grid = np.unique(sigma_fcs_array)
            y_grid = np.unique(y_array)

            Sigma, Y = np.meshgrid(sigma_grid, y_grid)

            # NILS: transpose CRUCIAL, otherwise z-values are misarranged
            Cm_grid = Cm_alpha_array.reshape(len(y_grid), len(sigma_grid)).T
            Cn_grid = Cn_beta_array.reshape(len(y_grid), len(sigma_grid)).T
            Cl_grid = Cl_beta_array.reshape(len(y_grid), len(sigma_grid)).T

            gs_inner = gridspec.GridSpecFromSubplotSpec(
                3, 1,
                subplot_spec=gs_outer[row, col],
                hspace=0.1,
            )

            ax_cm = fig.add_subplot(gs_inner[0])
            ax_cn = fig.add_subplot(gs_inner[1], sharex=ax_cm)
            ax_cl = fig.add_subplot(gs_inner[2], sharex=ax_cm)

            cm = ax_cm.contourf(Sigma, Y, Cm_grid, levels=Cm_levels, cmap=cmap_cm)
            cn = ax_cn.contourf(Sigma, Y, Cn_grid, levels=Cn_levels, cmap=cmap_cn)
            cl = ax_cl.contourf(Sigma, Y, Cl_grid, levels=Cl_levels, cmap=cmap_cl)

            cm.set_clim(Cm_vmin, Cm_vmax)
            cn.set_clim(Cn_vmin, Cn_vmax)
            cl.set_clim(Cl_vmin, Cl_vmax)

            if cm_mappable is None:
                cm_mappable, cn_mappable, cl_mappable = cm, cn, cl

            ax_cm.set_title(titles[j], pad=15)
            
            # Remove spurious axes ticks
            ax_cm.tick_params(labelbottom=False, labelleft=True)
            ax_cn.tick_params(labelbottom=False, labelleft=True)
            ax_cl.tick_params(labelbottom=True, labelleft=True)

            for ax in (ax_cm, ax_cn, ax_cl):
                ax.spines[['right', 'top']].set_visible(False)
                ax.tick_params(axis='x', length=0)
                ax.tick_params(axis='y', length=0)
                ax.set_xticks([1.5, 4])
                ax.set_yticks([0.5, 1])

            ax_cl.tick_params(axis='x', labelbottom=True)

        elif i == 6:
            ax = fig.add_subplot(gs_outer[row, col])

            ax.plot(sigma_fcs_array, Cm_alpha_array)
            ax.plot(sigma_fcs_array, Cn_beta_array)
            ax.plot(sigma_fcs_array, Cl_beta_array)

            ax.set_title(titles[j], pad=15)
            ax.spines[['right', 'top']].set_visible(False)
            ax.tick_params(axis='x', length=0)
            ax.tick_params(axis='y', length=0)
            
    if i in (3, 4):
        cm_sm = ScalarMappable(norm=Cm_norm, cmap=cmap_cm)
        cn_sm = ScalarMappable(norm=Cn_norm, cmap=cmap_cn)
        cl_sm = ScalarMappable(norm=Cl_norm, cmap=cmap_cl)
        cm_sm.set_array([])
        cn_sm.set_array([])
        cl_sm.set_array([])

    if i in (3, 4):

        cbar_width = 0.015
        x_cbar = 0.88
    
        cax_cm = fig.add_axes([x_cbar, 2/3 - 0.03, cbar_width, 0.25])
        cax_cn = fig.add_axes([x_cbar, 1/3, cbar_width, 0.25])
        cax_cl = fig.add_axes([x_cbar, 0.0 + 0.03, cbar_width, 0.25])
    
        cbar_cm = fig.colorbar(cm_sm, cax=cax_cm)
        cbar_cn = fig.colorbar(cn_sm, cax=cax_cn)
        cbar_cl = fig.colorbar(cl_sm, cax=cax_cl)
        cbar_cm.set_label(r'Longitudinal stability, $C_{m,\alpha}$', labelpad=20)
        cbar_cn.set_label(r'Directional stability, $C_{n,\beta}$', labelpad=20)
        cbar_cl.set_label(r'Lateral stability, $C_{l,\beta}$', labelpad=20)
    
    fig.text(
        0.5, 0.02,
        'Specific power (kW/kg)',
        ha='center', va='bottom',
    )
    
    plt.show()


sys.exit()

#%%

# =============================================================================

# --- Load the data ---
data = np.loadtxt(r"C:\Users\nmb48\base_stability_data.txt")      # shape (36, 4)

# --- Extract just one quantity (example: CM = col 0) ---
CM_flat = data[:, 0]
Cm_alpha_flat = data[:, 1]
Cn_beta_flat = data[:, 2]
NP_flat = data[:, 3]
static_margin_flat = data[:, 4]
Cl_beta_flat = data[:, 5]
Cn_r_flat = data[:, 6]
Cl_r_flat = data[:, 7]

# --- Reshape to (6 AoA × 6 Mach) ---
CM = CM_flat.reshape(6, 6)
Cm_alpha = Cm_alpha_flat.reshape(6, 6)
Cn_beta = Cn_beta_flat.reshape(6, 6)
NP = NP_flat.reshape(6, 6)
static_margin = static_margin_flat.reshape(6, 6)
Cl_beta = Cl_beta_flat.reshape(6, 6)
Cn_r = Cn_r_flat.reshape(6, 6)
Cl_r = Cl_r_flat.reshape(6, 6)

# --- Define axes ---
mach_vals = np.array([0.05, 0.15, 0.25, 0.45, 0.65, 0.85])
aoa_vals = np.array([-2., 0., 2., 5., 7., 10.])

# --- Build DataFrame ---
df_CM = pd.DataFrame(CM, index=mach_vals, columns=aoa_vals)
df_Cm_alpha = pd.DataFrame(Cm_alpha, index=mach_vals, columns=aoa_vals)
df_Cn_beta = pd.DataFrame(Cn_beta, index=mach_vals, columns=aoa_vals)
df_NP = pd.DataFrame(NP, index=mach_vals, columns=aoa_vals)
df_static_margin = pd.DataFrame(static_margin, index=mach_vals, columns=aoa_vals)

MACH, AOA = np.meshgrid(mach_vals, aoa_vals)
df_list = [df_CM, df_Cm_alpha, df_Cn_beta, df_NP, df_static_margin]

r"""
# for df in df_list:

#     fig, ax = plt.subplots()
#     cf = ax.contour(MACH, AOA, df.to_numpy())
#     fig.colorbar(cf)
#     plt.show()
    
#%%

fig, ax = plt.subplots()

for i, mach in enumerate(mach_vals):
    ax.plot(aoa_vals, df_list[0].to_numpy()[i,:])

plt.show()

#%%

fig, ax = plt.subplots()

for i, mach in enumerate(mach_vals):
    ax.plot(aoa_vals, df_list[1].to_numpy()[i,:])

plt.show()

#%%

fig, ax = plt.subplots()

for i, mach in enumerate(mach_vals):
    ax.plot(aoa_vals, df_list[2].to_numpy()[i,:])

plt.show()

#%%

# fig, ax = plt.subplots()

# for i, mach in enumerate(mach_vals):
#     ax.plot(aoa_vals, df_list[4].to_numpy()[i,:])

# plt.show()

fig, ax = plt.subplots()

for i, mach in enumerate(aoa_vals):
    ax.plot(mach_vals, df_list[4].to_numpy()[:,i])

plt.show()

# sys.exit()
"""

#%%

E = Cl_beta * Cn_r - Cn_beta * Cl_r
Cl_beta_avg = np.average(Cl_beta)
Cn_beta_avg = np.average(Cn_beta)

fig, ax = plt.subplots(figsize=(8,6))
# ax.scatter(-Cl_beta, Cn_beta, c=E)
# ax.scatter(-Cl_beta, Cn_beta, c=MACH)
ax.scatter(-Cl_beta, Cn_beta, c=AOA)

for i in range(np.shape(Cl_beta)[0]):
    ax.plot(-Cl_beta[i,:], Cn_beta[i,:], color=colors[0])

for j in range(np.shape(Cl_beta)[1]):
    ax.plot(-Cl_beta[:,j], Cn_beta[:,j], color=colors[1])

ax.spines[['right', 'top']].set_visible(False)
ax.tick_params(axis='y', which='both', right=False, length=0)
ax.tick_params(axis='x', which='both', length=0)
# ax.spines['bottom'].set_position(('data', 0.0))
# ax.spines['left'].set_position(('data', 0.0))

plt.show()

#%% This is equivalent to the above code around pick()

# dutchRoll_data = []
# rollSubsistence_data = []
# spiral_data = []
                         
# for i in range(len(df)):
#     _dutchRollInd = df["dutchRollInd"][i].astype(int)
#     _rollSubsistenceInd = df["rollSubsistenceInd"][i].astype(int)
#     _spiralInd = df["spiralInd"][i].astype(int)
    
#     dutchRoll_data.append(df[f"Lat_Re{_dutchRollInd + 1}"][i] + 1j * df[f"Lat_Im{_dutchRollInd + 1}"][i])
#     rollSubsistence_data.append(df[f"Lat_Re{_rollSubsistenceInd + 1}"][i] + 1j * df[f"Lat_Im{_rollSubsistenceInd + 1}"][i])
#     spiral_data.append(df[f"Lat_Re{_spiralInd + 1}"][i] + 1j * df[f"Lat_Im{_spiralInd + 1}"][i])
    
# dutchRoll_data = np.array(dutchRoll_data)
# rollSubsistence_data = np.array(rollSubsistence_data)
# spiral_data = np.array(spiral_data)
