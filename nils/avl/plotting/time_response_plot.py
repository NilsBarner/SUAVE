"""
This script plots 3D modal time response about y=0 and z=0 planes.
It is a qualitative/quantitative plot meant to visualisation only.
"""

__all__ = []

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from matplotlib_custom_settings import *


def plot_modal_response_3d(
    eigenvalues,
    A=1.0,
    phase=0.0,
    n_cycles=3,
    n_time=500,
):
    """
    Plot a 3D modal time response about y=0 and z=0 planes.

    Parameters
    ----------
    eigenvalues : complex
        Eigenvalues sigma + i*omega
    A : float
        Initial amplitude
    phase : float
        Initial phase [rad]
    n_cycles : int
        Number of oscillation cycles
    n_time : int
        Number of time samples
    """
    
    x_offset = -1  # shift axes slightly to the left
    
    # --- 3D plot ---
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')
    
    ax.zaxis.pane.set_alpha(0.1)
    ax.zaxis.pane.set_edgecolor('k')
    
    ax.yaxis.pane.set_alpha(0.1)
    ax.yaxis.pane.set_edgecolor('k')
    
    ax.xaxis.pane.set_visible(False)
    
    j = 0
    
    for i, eigenvalue in enumerate(eigenvalues):

        sigma = np.real(eigenvalue)
        omega = np.imag(eigenvalue)
    
        # --- Time span ---
        t_max = 50
    
        t = np.linspace(0, t_max, n_time)
    
        # --- Modal response ---
        if abs(omega) > 0:
            x = A * np.exp(sigma * t) * np.cos(omega * t + phase)
        else:
            x = A * np.exp(sigma * t)
    
        if (i == 0 or i == 3 or i == 5):
            ax.plot(t,  x, np.zeros_like(x), lw=2, zorder=-i, color=colors[j])
            j += 1
    
        elif (i == 1 or i == 6):
            ax.plot(t, np.zeros_like(x),  x, lw=2, color=colors[j])
            j += 1
    
    # --- Make reference planes transparent ---
    ax.xaxis.pane.set_visible(False)
    ax.yaxis.pane.set_visible(False)
    ax.zaxis.pane.set_visible(False)
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.xaxis.line.set_color((0,0,0,0))
    ax.yaxis.line.set_color((0,0,0,0))
    ax.zaxis.line.set_color((0,0,0,0))
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_zticklabels([])
    
    # --- Draw axes through origin ---
    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()
    z_min, z_max = ax.get_zlim()

    # Axis spines
    # ax.plot([x_offset, t.max()], [0, 0], [0, 0], color='k', lw=1, zorder=-100)
    # ax.plot([x_offset, x_offset], [-1, 1], [0, 0], color='k', lw=1)  # y-axis
    # ax.plot([x_offset, x_offset], [0, 0], [-1, 1], color='k', lw=1)  # z-axis

    # x-y plane at z=0
    xy_rect = np.array([
        [x_offset, -1, 0],
        [t_max, -1, 0],
        [t_max, 1, 0],
        [x_offset, 1, 0]
    ])
    ax.add_collection3d(Poly3DCollection([xy_rect], facecolor='lightgrey', edgecolor='grey', alpha=0.1, lw=1))
    
    # x-z plane at y=0
    xz_rect = np.array([
        [x_offset, 0, -1],
        [t_max, 0, -1],
        [t_max, 0, 1],
        [x_offset, 0, 1]
    ])
    ax.add_collection3d(Poly3DCollection([xz_rect], facecolor='lightgrey', edgecolor='grey', alpha=0.1, lw=1))
    
    # Improve aspect perception
    ax.set_box_aspect((2, 1, 1))
    
    plt.tight_layout()
    plt.savefig('modal_response_raw.svg', format='svg')
    plt.show()

#%%

if __name__ == '__main__':
    
    txt_file = r"C:\Users\nmb48\Documents\GitHub\SUAVE\Tutorials-2.5.0.2\B737_AVL_Tutorial\CUsersnmb48\suave_dynamic_stability_outputs_6_0.txt"
    df = pd.read_csv(txt_file, sep=r"\s+", engine="python")
    
    re = df[[f"Re{i}" for i in range(1, 9)]].to_numpy()
    im = df[[f"Im{i}" for i in range(1, 9)]].to_numpy()
    cplx = re + 1j * im

    plot_modal_response_3d(cplx[0])
        
