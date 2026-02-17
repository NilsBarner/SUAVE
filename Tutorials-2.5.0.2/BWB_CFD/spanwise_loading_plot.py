__all__ = ["plot_spanwise_loading"]

import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import gridspec
from scipy.interpolate import interp1d

from matplotlib_custom_settings import *

_mesh = pv.read(r"C:\Users\nmb48\Documents\GitHub\SUAVE\Tutorials-2.5.0.2\BWB_CFD\tasopt\base_surface_flow.vtu")

def plot_spanwise_loading(
    # mesh,
    y_main, dCL_main,
):
    
    # Ensure triangles
    mesh = _mesh.triangulate()
    
    # Convert point → cell data (as you already did)
    mesh = mesh.point_data_to_cell_data()
    
    # Extract surface: UnstructuredGrid → PolyData
    surf = mesh.extract_surface()
    
    # Now compute normals (this WILL work)
    surf = surf.compute_normals(
        cell_normals=True,
        point_normals=False,
        auto_orient_normals=True,
        inplace=False,
    )
    
    # Use surf from here on
    centers = surf.cell_centers().points
    areas   = surf.compute_cell_sizes(length=False, area=True).cell_data["Area"]
    normals = surf.cell_data["Normals"]
    Cp      = surf.cell_data["Pressure_Coefficient"]
    
    # Lift direction (global Z)
    lift_dir = np.array([0.0, 0.0, 1.0])
    
    # Differential force per cell (non-dimensional)
    dF = -Cp[:, None] * areas[:, None] * normals
    dL = dF @ lift_dir
    
    # Spanwise coordinate (Y)
    y = centers[:, 1]
    
    # Bin spanwise
    nbins = 50
    y_bins = np.linspace(y.min(), y.max(), nbins + 1)
    y_mid  = 0.5 * (y_bins[:-1] + y_bins[1:])
    
    L_span = np.zeros(nbins)
    for i in range(nbins):
        mask = (y >= y_bins[i]) & (y < y_bins[i+1])
        L_span[i] = dL[mask].sum()
    
    # Normalize (AVL-like)
    L_total = L_span.sum()
    cl_span = L_span / L_total
    
    #%% Plot
    
    y_main_norm = y_main / max(y_main)
    y_mid_norm = y_mid / max(y_mid)
    
    y_main_interp = interp1d(y_main_norm, dCL_main)
    dCL_main_norm = dCL_main / np.trapz(dCL_main, y_main_norm)
    dCL_main_interp_SU2 = y_main_interp(y_mid_norm)
    dCL_main_interp_SU2_norm = dCL_main_interp_SU2 / np.trapz(dCL_main_interp_SU2, y_mid_norm)
    
    # =============================================================================
    cl_span_norm = cl_span / np.trapz(cl_span, y_mid_norm)
    # =============================================================================
    
    fig = plt.figure(figsize=(8,5.5), constrained_layout=True)
    gs = gridspec.GridSpec(
        2, 1, width_ratios=[1], hspace=0.5, height_ratios=[1,0.2], wspace=0.0,
    )
    ax = fig.add_subplot(gs[0, 0])
    
    # ax.plot(y_main, dCL_main, color=colors[0], label='TASOPT.jl')
    ax.plot(y_main_norm, dCL_main_norm, marker=".", color=colors[0], label='TASOPT.jl Trefftz plane vortices', clip_on=False)
    # ax.plot(y_mid, dCL_main_interp_SU2_norm, color=colors[0], label='TASOPT.jl')
    # ax.plot(y_mid, cl_span, marker=".", color=colors[1], label='SU2 CFD')
    ax.plot(y_mid_norm, cl_span_norm, marker=".", color=colors[1], label='SU2 CFD wing surface pressures', clip_on=False)
    
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Normalised spanwise coordinate", labelpad=15)
    ax.set_ylabel("Normalised lift distribution", labelpad=15)
    ax.spines[['right', 'top']].set_visible(False)
    ax.tick_params(axis='y', which='both', right=False, length=0)
    ax.tick_params(axis='x', which='both', length=0)
    # plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,-2))
    
    # ax.legend(frameon=False, loc='upper right')
    ax.legend(frameon=False, loc='lower left')
    
    for spine in ax.spines.values():
        spine.set_zorder(0)
    
    # Plot TASOPT.jl error
    
    ax_2 = fig.add_subplot(gs[1, 0])
    
    alpha = 1.0
    x_values = y_mid
    
    # Change in thickness
    cl_rel_error = (dCL_main_interp_SU2_norm - cl_span_norm) / cl_span_norm
    
    # Create a clip patch in axes coordinates to allow plotting outside axis limits
    clip_rect = patches.Rectangle(
        (0, -1),  # x=0, y slightly below axis bottom (axes coords)
        1,  # full width of axes
        2,  # goes from y=-0.05 to y=1.0
        transform=ax_2.transAxes,
        clip_on=True
    )
    
    # Define fill colour
    # base_hex = ax_2._get_lines.get_next_color()
    # color_active = opaque_color_from_hex(base_hex, alpha=0.8, background='white')
    
    # Positive region
    y_pos = np.where(cl_rel_error[x_values >= 0] > 0, cl_rel_error[x_values >= 0], 0)
    poly_positive = ax_2.fill_between(
        x_values[x_values >= 0] / max(x_values),
        0,
        y_pos,
        facecolor=colors[0],
        alpha=alpha,
    )
    poly_positive.set_clip_path(clip_rect)
        
    # Negative region
    y_neg = np.where(cl_rel_error[x_values >= 0] < 0, cl_rel_error[x_values >= 0], 0)
    poly_negative = ax_2.fill_between(
        x_values[x_values >= 0] / max(x_values),
        0,
        y_neg,
        facecolor=colors[0],
        alpha=alpha,
    )
    poly_negative.set_clip_path(clip_rect)
        
    # ax_2.set_xlim(0, 1)
    # ax_2.set_ylim(bottom=0, top=0.5)
    ax_2.spines['bottom'].set_position(('data', 0.0))
    ax_2.set_xticks([])
    ax_2.set_yticks([])
    ax_2.spines[['left', 'right', 'top']].set_visible(False)
    ax_2.tick_params(axis='y', which='both', right=False, length=0)
    ax_2.tick_params(axis='x', which='both', length=0)
    # ax_2.set_ylabel(
    #     "Error",
    #     rotation='horizontal',   # keep label horizontal
    #     ha='right',              # text aligned to the right (so it goes left of the spine)
    #     # va='center',             # vertical center along the spine
    #     y=0,
    #     labelpad=20,
    # )
    # =============================================================================
    # Example: your axis
    ax_2.spines['bottom'].set_position(('data', 0.0))
    ax_2.spines['left'].set_position(('data', 0.0))
    
    # Place the ylabel centered along the spine in data coordinates
    y_bottom, y_top = ax_2.get_ylim()         # y-limits of the axis
    y_center = 0.5 * (0.0 + y_top)            # center along spine: from bottom of spine (0.0) to top of y-axis
    ax_2.set_ylabel(
        "Error",
        rotation='horizontal',
        ha='right',
        va='center'
    )
    ax_2.yaxis.set_label_coords(-0.1, y_center, transform=ax_2.transData)  # use data coords
    # =============================================================================
    ax_2.set_zorder(-1)

    plt.show()

    return

#%%

if __name__ == '__main__':
    
    mesh = pv.read(r"C:\Users\nmb48\Documents\GitHub\SUAVE\Tutorials-2.5.0.2\BWB_CFD\tasopt\base_surface_flow.vtu")
    
    plot_spanwise_loading(mesh)
