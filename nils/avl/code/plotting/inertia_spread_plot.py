r"""
This script plots the spread of the non-zero components of the mass moment of
inertia tensor (I_xx, I_yy, I_zz and I_zx) across the mass distribution
studies 1, 2, 5, 6 and 7 for the ATR 72-600 and the Airbus A220-100, as
written to mass_distr_results_<segment>_LH2_<date>_<study>.csv by
TASOPT.jl-priv/nils/avl/code/mass_distr_exec.jl.

To make both aircraft comparable on one axis, the diagonal terms are shown as
non-dimensional radii of gyration (Roskam, Airplane Design Part V):
    R_x = 2 k_x / b,   R_y = 2 k_y / L,   R_z = 2 k_z / e,   e = (b + L) / 2,
with k_i = sqrt(I_ii / m), wing span b and fuselage length L. The product of
inertia is shown as the inclination of the principal x-axis to the body
x-axis, eps = 0.5 atan2(2 I_zx, I_zz - I_xx). In the x-z plane I_zx only
enters the dynamics through this rotation. AVL axes are used (x aft, z up),
and the sign of I_zx is the same in conventional body axes (x fwd, z down).

Each dot is one sampled design; each bar spans the range of one study.

NILS (10.10.2026): new script.

Usage (Spyder, or from Anaconda Prompt):
(suave) C:\Users\nmb48\Documents\GitHub\SUAVE>python nils\avl\code\plotting\inertia_spread_plot.py
"""

import os
import ast
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from nils.matplotlib_custom_settings import *

# ---------------------------------------------------------------------
# User inputs
# ---------------------------------------------------------------------

TASOPT_ROOT = os.environ.get("TASOPT_ROOT", r"C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv")  # results and geometry .csv files
OUT_DIR = os.environ.get("INERTIA_PLOT_OUT", TASOPT_ROOT)

# Aircraft name -> TASOPT file naming (update the dates after re-running mass_distr_exec.jl)
AIRCRAFT = {
    "ATR 72-600": dict(segment="regional", fuel="LH2", date="041026", color=colors[0], marker="o", linestyle="solid"),
    "A220-100": dict(segment="narrowbody", fuel="LH2", date="041026", color=colors[3], marker="s", linestyle='dashed'),
}
# Study index -> row label (top to bottom)
STUDIES = {
    # # 1: "Nacelles + fuselage (1)",
    # # 2: "Nacelles + wing (2)",
    # 6: "Wing",  # (6)",
    # 7: "Fuselage",  # (7)",
    # 5: r"Wing \& fuselage",  # (5)",
    6: dict(xlabel="Wing", color=colors[1]),
    7: dict(xlabel="Fuselage", color=colors[0]),
    5: dict(xlabel=r"Wing \& fuselage", color=colors[2]),
}
MISSION_POINT = 1  # 1 take-off ... 6 landing (inertia changes with the fuel remaining)
SAVE_FIG = True  # True
USE_CUSTOM_STYLE = True  # nils.matplotlib_custom_settings (LaTeX fonts) if importable

# ---------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------

def _mp(series, mp):
    """Value at mission point mp from the "[a, b, c, d, e, f]" strings written by mass_distr_exec.jl"""
    return np.array([ast.literal_eval(v)[mp - 1] for v in series], dtype=float)


def load_inertia(aircraft, study_idx, mp=MISSION_POINT):
    a = AIRCRAFT[aircraft]
    res = os.path.join(TASOPT_ROOT, f"mass_distr_results_{a['segment']}_{a['fuel']}_{a['date']}_{study_idx}.csv")
    geo = os.path.join(TASOPT_ROOT, f"mass_distr_geometry_{a['segment']}_{a['fuel']}_{a['date']}.csv")
    if not (os.path.exists(res) and os.path.exists(geo)):
        return None
    df = pd.read_csv(res)
    g = pd.read_csv(geo).iloc[0]
    b, L = float(g["wing_spans_projected"]), float(g["body_lengths_total"])
    e = 0.5 * (b + L)
    m = _mp(df["mass_mp"], mp)
    Ixx, Iyy, Izz, Izx = (_mp(df[c], mp) for c in ["Ixx_mp", "Iyy_mp", "Izz_mp", "Izx_mp"])
    return pd.DataFrame({
        "R_x": 2.0 / b * np.sqrt(Ixx / m),
        "R_y": 2.0 / L * np.sqrt(Iyy / m),
        "R_z": 2.0 / e * np.sqrt(Izz / m),
        "eps": np.degrees(0.5 * np.arctan2(2.0 * Izx, Izz - Ixx)),
    })

# ---------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------

PANELS = [
    ("R_x", r"$\bar R_x = 2k_x/b$ [-]"),
    ("R_y", r"$\bar R_y = 2k_y/L$ [-]"),
    ("R_z", r"$\bar R_z = 2k_z/e$ [-]"),
    ("eps", r"Principal axis inclination [$^\circ$]"),  # $\varepsilon$
]


def plot_inertia_spread():
    data = {(ac, s): load_inertia(ac, s) for ac in AIRCRAFT for s in STUDIES}
    studies = [s for s in STUDIES if any(data[(ac, s)] is not None for ac in AIRCRAFT)]
    n_ac = len(AIRCRAFT)
    rng = np.random.default_rng(0)  # reproducible jitter

    fig, axes = plt.subplots(2, 2, figsize=(9.0, 7.0), sharey=True)
    axes = axes.ravel()
    for ax, (key, label) in zip(axes, PANELS):
        for r, s in enumerate(studies):
            y_row = len(studies) - 1 - r
            for k, (ac, a) in enumerate(AIRCRAFT.items()):
                df = data[(ac, s)]
                if df is None:
                    continue
                color = STUDIES[s]["color"]
                y = y_row + 0.18 * (1 - 2 * k / max(n_ac - 1, 1))  # ATR above, A220 below the row centre
                v = df[key].to_numpy()
                # ax.plot([v.min(), v.max()], [y, y], color=a["color"], linewidth=1.0, solid_capstyle="butt",
                #     marker="|", markersize=8, markeredgewidth=1.0, zorder=2)
                ax.plot([v.min(), v.max()], [y, y], color=color, linewidth=1.0, linestyle=a["linestyle"], solid_capstyle="butt",
                    marker="|", markersize=8, markeredgewidth=1.0, zorder=2)
                # ax.scatter(v, y + rng.uniform(-0.06, 0.06, len(v)), s=14, marker=a["marker"], facecolor=a["color"],
                #            edgecolor="k", linewidth=0.3, zorder=3, clip_on=False)
                ax.scatter(v, y + rng.uniform(-0.06, 0.06, len(v)), s=14, marker=a["marker"], facecolor=color,
                    edgecolor="k", linewidth=0.3, zorder=3, clip_on=False)
        ax.set_xlabel(label, labelpad=15)
        ax.set_yticks(range(len(studies)))
        # ax.set_yticklabels([STUDIES[s] for s in studies][::-1])
        ax.set_yticklabels([STUDIES[s]["xlabel"] for s in studies][::-1])
        ax.set_ylim(-0.6, len(studies) - 0.4)
        for yl in np.arange(0.5, len(studies) - 1, 1.0):
            ax.axhline(yl, color="0.85", linewidth=0.6, zorder=0)
        ax.spines[["right", "top"]].set_visible(False)
        ax.tick_params(axis='y', which='both', right=False, length=0)
        ax.tick_params(axis='x', which='both', length=0)
    # handles = [Line2D([], [], color=a["color"], marker=a["marker"], markerfacecolor=a["color"], markeredgecolor="k",
    #                   markeredgewidth=0.3, linewidth=1.0, label=ac) for ac, a in AIRCRAFT.items()]
    handles = [Line2D([], [], color="k", marker=a["marker"], markerfacecolor="k", markeredgecolor="k",
        markeredgewidth=0.3, linestyle=a["linestyle"], linewidth=1.0, label=ac) for ac, a in AIRCRAFT.items()]
    # fig.legend(handles=handles, loc="upper center", ncol=n_ac, frameon=False, bbox_to_anchor=(0.5, 1.0))
    ###
    axes[1].legend(handles=handles, loc="upper right", ncol=1, frameon=False, bbox_to_anchor=(1.05, 1.0), fontsize=14)
    ###
    # fig.tight_layout(rect=(0, 0, 1, 0.95))
    if SAVE_FIG:
        fig.savefig(os.path.join(OUT_DIR, "inertia_spread_stability.png"), format='png', dpi=300)
        fig.savefig(os.path.join(OUT_DIR, "inertia_spread_stability.svg"), format='svg')
    return fig


if __name__ == "__main__":
    if USE_CUSTOM_STYLE:
        try:
            from nils.matplotlib_custom_settings import *  # same fonts as the other thesis figures
        except ImportError:
            pass
    plot_inertia_spread()
    plt.show()
