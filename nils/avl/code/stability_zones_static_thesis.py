r"""
This script generalises Fig. 4.8 (static stability zones) to both aircraft and
all mass distribution studies, adds the corresponding longitudinal and lateral
dynamic stability zone maps (i.e. where the eigenvalues in Figs. 4.9 and 4.10
satisfy MIL-F-8785C / MIL-STD-1797A Level 1 at all six mission points), and
exports the feasible fuel cell system locations to
stability_zones_<aircraft>.csv, which is read by
TASOPT.jl-priv/nils/airbus/stability_zone_plotter.py (Fig. 4.11).

The green area of each map comprises all designs that satisfy the criteria
throughout the entire mission; the lines mark the boundary of each mission
point (line styles and grey shades as in Fig. 4.8).

To be run from Spyder or as module from Anaconda Prompt via:
(suave) C:\Users\nmb48\Documents\GitHub\SUAVE>python -m nils.avl.code.stability_zones

NILS (04.10.2026): new script. The stability criteria (MIL_LIMITS, SM_min,
DELTA_E_MIN) are defined here only and imported by process_dynamic_stability.py.
"""

__all__ = [
    "SM_min", "DELTA_E_MIN", "MIL_LIMITS", "flight_phase_category",
    "COLOR_OK", "COLOR_FAIL", "ZONE_ALPHA", "MP_LABELS", "MP_LINESTYLES", "MP_GREYS",
    "STUDY_INFO", "load_static", "load_dynamic", "static_margins", "dynamic_margins",
    "zone_map", "plot_zone_figure", "export_zones",
]

import os
import io
import re
import json
import glob
import tempfile
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from scipy.interpolate import RegularGridInterpolator, interp1d
from matplotlib.ticker import FuncFormatter

from nils.matplotlib_custom_settings import *

# ---------------------------------------------------------------------
# User inputs
# ---------------------------------------------------------------------

SUAVE_ROOT = os.environ.get("SUAVE_ROOT", r"C:\Users\nmb48\Documents\GitHub\SUAVE\nils\avl\data\ch7")  # folder with suave_*_stability_*.txt files
TASOPT_ROOT = os.environ.get("TASOPT_ROOT", r"C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv")  # folder with mass_distr_results_*.csv files
OUT_DIR = os.environ.get("STABILITY_ZONES_OUT", SUAVE_ROOT)  # figures and stability_zones_<aircraft>.csv

# Aircraft tag in the SUAVE output file names -> TASOPT results file mass_distr_results_<segment>_<fuel>_<date>_<study>.csv
AIRCRAFT = {
    "ATR_72-600": dict(segment="regional", fuel="LH2", date="041026"),
    "Airbus_A220-100": dict(segment="narrowbody", fuel="LH2", date="041026"),  # NOTE: set to the date of the A220 mass_distr_exec.jl run
}
STUDIES = [7, 6, 5, 1, 2]  # studies 3 and 4 are identical to 7 and 6 and are plotted too if their files exist
SAVE_FIGS = False  # save .png figures in OUT_DIR
SHOW_SAMPLES = True  # mark the sampled designs (5 x 5 grid) in the zone maps
N_FINE = 60  # interpolation points between two sampled designs along each axis

# ---------------------------------------------------------------------
# Stability criteria (single source for static_margin_plots.py and process_dynamic_stability.py)
# ---------------------------------------------------------------------

# Aft CG limit: minimum static margin (fraction of MAC) as in TASOPT's default horizontal tail SM_min
SM_min = 0.15
# Forward CG limit: CS-25 does not prescribe a single static margin; the forward limit follows from elevator
# authority (CS 25.143/25.145 controllability, nose-wheel lift-off and landing flare). Proxy used here: the trimmed
# elevator deflection (deg, +ve trailing edge down) must not exceed DELTA_E_MIN in nose-up (negative) direction
DELTA_E_MIN = -20.0  # (deg) ~2/3 of a typical -30 deg nose-up elevator travel (assumption, adjust per aircraft)

# MIL-F-8785C / MIL-STD-1797A Level 1 limits for Class II (land-based) aircraft; take-off and landing (mission
# points 1 and 6) are terminal flight phases (Category C), climb, cruise and descent are Category B
MIL_LIMITS = {
    'B': dict(zeta_sp_min=0.30, zeta_sp_max=2.00, CAP_sp_min=0.085, CAP_sp_max=3.6,
              zeta_ph_min=0.04, tau_roll_max=1.4, T_2_sl_min=20.0,
              zeta_dr_min=0.08, zeta_omega_dr_min=0.15, omega_n_dr_min=0.4),
    'C': dict(zeta_sp_min=0.35, zeta_sp_max=1.30, CAP_sp_min=0.16, CAP_sp_max=3.6,
              zeta_ph_min=0.04, tau_roll_max=1.4, T_2_sl_min=12.0,
              zeta_dr_min=0.08, zeta_omega_dr_min=0.10, omega_n_dr_min=0.4),
}

def flight_phase_category(rid):
    """rid = 1..6 (take-off, climb, start of cruise, end of cruise, descent, landing)"""
    return 'C' if rid in (1, 6) else 'B'

# ---------------------------------------------------------------------
# Plot style shared with Figs. 4.8-4.11
# ---------------------------------------------------------------------

COLOR_OK = colors[1]  # 'green'  # zones and eigenvalues that satisfy the criteria (Figs. 4.8-4.11)
COLOR_FAIL = colors[3]  # 'red'  # zones and eigenvalues that violate them
ZONE_ALPHA = 0.1  # transparency of the zone fill (Fig. 4.8)
MP_LABELS = ['Take-off', 'Climb', 'Start of cruise', 'End of cruise', 'Descent', 'Landing']
MP_LINESTYLES = ['-', ':', '--', '-.', '-', ':']
MP_GREYS = [str(c) for c in np.linspace(0.0, 0.8, 6)]

AXIS_LABELS = {
    'sigma_fcs': 'Fuel cell system specific power (kW/kg)',
    'fcs_loc': 'Fuel cell system position along fuselage',
    'span_loc': 'Fuel cell system position along wing semi-span',
    'nacelle_frac': 'Fuel cell system mass fraction in nacelles',
}
# x: abscissa, y: ordinate (None for one-dimensional sweeps), as in mass_distr_exec.jl
STUDY_INFO = {
    1: dict(name='Nacelles and fuselage', x='nacelle_frac', y=None),
    2: dict(name='Nacelles and wing', x='nacelle_frac', y=None),
    3: dict(name='Fuselage', x='fcs_loc', y='sigma_fcs'),
    4: dict(name='Wing', x='span_loc', y='sigma_fcs'),
    5: dict(name='Wing and fuselage', x='fcs_loc', y='span_loc'),
    6: dict(name='Wing', x='span_loc', y='sigma_fcs'),
    7: dict(name='Fuselage', x='fcs_loc', y='sigma_fcs'),
}
GROUP_TITLES = {
    'static': 'Static stability',
    'longitudinal': 'Longitudinal dynamic stability',
    'lateral': 'Lateral dynamic stability',
}

# ---------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------

def _case_index(path):
    return int(re.search(r"_stability_(?:outputs|matrix)_\d+_(\d+)_", os.path.basename(path)).group(1))

def _stability_files(kind, study_idx, aircraft):
    files = glob.glob(os.path.join(SUAVE_ROOT, f"suave_{kind}_{study_idx}_*_{aircraft}.txt"))
    return sorted(files, key=_case_index)

def read_mass_csv(aircraft, study_idx):
    """One row per case (row index = case index in the SUAVE output file names)"""
    a = AIRCRAFT[aircraft]
    path = os.path.join(TASOPT_ROOT, f"mass_distr_results_{a['segment']}_{a['fuel']}_{a['date']}_{study_idx}.csv")
    if not os.path.exists(path):
        return None
    df = pd.read_csv(path)
    for col in ['sigma_fcs', 'span_loc', 'fcs_loc', 'wing_frac', 'nacelle_frac']:
        df[col] = pd.to_numeric(df[col], errors='coerce')  # "eng" -> nan
    return df

def _read_ws_table(path):
    with open(path, "r") as fh:
        lines = fh.readlines()
    if lines and lines[0].lstrip().startswith("#"):
        lines[0] = lines[0].lstrip()[1:].lstrip()
    return pd.read_csv(io.StringIO("".join(lines)), sep=r"\s+", header=0)

def load_static(aircraft, study_idx):
    """Long table: one row per case and mission point (mp = 1..6) with the static outputs of AVL.py"""
    files = _stability_files("static_stability_outputs", study_idx, aircraft)
    if len(files) == 0:
        return None
    frames = []
    for f in files:
        df = _read_ws_table(f)
        df = df.drop(columns=[c for c in ['sigma_fcs', 'span_loc', 'fcs_loc', 'wing_frac', 'nacelle_frac'] if c in df.columns])
        df.insert(0, 'mp', np.arange(1, len(df) + 1))
        df.insert(0, 'case', _case_index(f))
        frames.append(df)
    return pd.concat(frames, ignore_index=True)

def load_dynamic(aircraft, study_idx):
    """Long table: one row per case and mission point with eigenvalue, zeta, wn and n_alpha of each of the five
    modes, classified by eigenvalue_classification.match_evals_to_emodes (one study per call, so that the case
    index of the combined files equals the case index of the SUAVE output files)"""
    from nils.avl.code.dynamic_stability.eigenvalue_classification import match_evals_to_emodes
    from nils.avl.code.dynamic_stability.eigenvalue_loading import combine_eigen_outputs
    from nils.avl.code.dynamic_stability.system_matrix_loading import combine_system_matrices
    eig_files = _stability_files("dynamic_stability_outputs", study_idx, aircraft)
    mat_files = _stability_files("dynamic_stability_matrix", study_idx, aircraft)
    if len(eig_files) == 0 or len(eig_files) != len(mat_files):
        return None
    with tempfile.TemporaryDirectory() as tmp:
        eigfile = combine_eigen_outputs(eig_files, os.path.join(tmp, "eig.txt"))
        matfile = combine_system_matrices(mat_files, os.path.join(tmp, "mat.txt"))
        res = match_evals_to_emodes(matfile, eigfile, plot_matching=False)
    cases = [_case_index(f) for f in eig_files]
    out = None
    for mode in ['short-period', 'phugoid', 'roll', 'spiral', 'dutch-roll']:
        frames = []
        for rid in range(1, 7):
            df = pd.concat([res[cat][mode][rid] for cat in res if len(res[cat][mode][rid]) > 0], ignore_index=True)
            df = df[['group', 'Re', 'Im', 'zeta', 'wn', 'n_alpha']].copy()
            df['mp'] = rid
            frames.append(df)
        df = pd.concat(frames, ignore_index=True)
        df['case'] = [cases[int(g)] for g in df['group']]
        df = df.drop(columns='group').rename(columns={c: f"{c}_{mode}" for c in ['Re', 'Im', 'zeta', 'wn', 'n_alpha']})
        out = df if out is None else out.merge(df, on=['case', 'mp'])
    return out

# ---------------------------------------------------------------------
# Criteria margins (>= 0 satisfied, < 0 violated, nan violated)
# ---------------------------------------------------------------------

def static_margins(df_static):
    m = df_static[['case', 'mp']].copy()
    m['SM'] = df_static['static_margin'] - SM_min  # aft limit
    if 'delta_e' in df_static.columns:
        m['delta_e'] = df_static['delta_e'] - DELTA_E_MIN  # forward limit (elevator authority)
    return m

def dynamic_margins(df_dyn, group):
    """Margins equivalent to mask_mil_limits() in process_dynamic_stability.py"""
    m = df_dyn[['case', 'mp']].copy()
    lim = pd.DataFrame([MIL_LIMITS[flight_phase_category(rid)] for rid in df_dyn['mp']], index=df_dyn.index)
    if group == 'longitudinal':
        zeta_sp, wn_sp = df_dyn['zeta_short-period'], df_dyn['wn_short-period']
        CAP_sp = wn_sp**2 / df_dyn['n_alpha_short-period']
        m['zeta_sp_min'] = zeta_sp - lim['zeta_sp_min']
        m['zeta_sp_max'] = lim['zeta_sp_max'] - zeta_sp
        m['CAP_sp_min'] = CAP_sp - lim['CAP_sp_min']
        m['CAP_sp_max'] = lim['CAP_sp_max'] - CAP_sp
        m['zeta_ph_min'] = df_dyn['zeta_phugoid'] - lim['zeta_ph_min']
    elif group == 'lateral':
        Re_dr, Im_dr = df_dyn['Re_dutch-roll'], df_dyn['Im_dutch-roll']
        wn_dr = np.hypot(Re_dr, Im_dr)
        zeta_dr = -Re_dr / wn_dr
        m['tau_roll_max'] = -df_dyn['Re_roll'] - 1.0 / lim['tau_roll_max']  # tau = -1/Re <= tau_max
        m['T_2_sl_min'] = np.log(2.0) / lim['T_2_sl_min'] - df_dyn['Re_spiral']  # stable or T_2 = ln2/Re >= T_2_min
        m['zeta_dr_min'] = zeta_dr - lim['zeta_dr_min']
        m['zeta_omega_dr_min'] = zeta_dr * wn_dr - lim['zeta_omega_dr_min']
        m['omega_n_dr_min'] = wn_dr - lim['omega_n_dr_min']
    return m

# ---------------------------------------------------------------------
# Zone maps
# ---------------------------------------------------------------------

def _fine_axis(values):
    """Fine axis that contains every sampled value (linear interpolation is then exact at the samples)"""
    values = np.asarray(values, dtype=float)
    pieces = [np.linspace(values[k], values[k + 1], N_FINE + 1)[:-1] for k in range(len(values) - 1)]
    return np.concatenate(pieces + [values[-1:]])

def zone_map(df_mass, margins, study_idx):
    """Interpolates every criterion margin of every mission point linearly between the sampled designs and
    returns the fine axes, the per-mission-point feasibility (ok_mp[rid]), the overall feasibility (ok_all) and
    the feasibility of each criterion over all mission points (ok_crit). A violated or undefined (nan) margin
    at a sampled design is set to the most negative margin of that criterion before interpolating."""
    info = STUDY_INFO[study_idx]
    xv, yv = info['x'], info['y']
    x = np.sort(df_mass[xv].dropna().unique())
    y = np.sort(df_mass[yv].dropna().unique()) if yv is not None else None
    x_f = _fine_axis(x)
    y_f = _fine_axis(y) if yv is not None else None
    crit_cols = [c for c in margins.columns if c not in ('case', 'mp')]
    ok_mp, ok_crit = {}, {c: None for c in crit_cols}
    for rid in range(1, 7):
        mrow = margins[margins['mp'] == rid].set_index('case')
        ok = None
        for c in crit_cols:
            vals = mrow[c].astype(float)
            finite = vals.to_numpy()[np.isfinite(vals.to_numpy())]
            if len(finite) == 0:
                fill = -1.0
            elif finite.min() < 0.0:
                fill = finite.min()
            else:
                fill = -0.1 * max(finite.max(), 1e-6)
            vals = vals.fillna(fill)
            if yv is None:
                arr = np.array([vals.get(case, np.nan) for case in df_mass.index], dtype=float)
                order = np.argsort(df_mass[xv].to_numpy())
                f = interp1d(df_mass[xv].to_numpy()[order], arr[order])
                field = f(x_f)
            else:
                arr = np.full((len(y), len(x)), np.nan)
                for case, row in df_mass.iterrows():
                    i = np.searchsorted(y, row[yv]); j = np.searchsorted(x, row[xv])
                    arr[i, j] = vals.get(case, np.nan)
                interp = RegularGridInterpolator((y, x), arr, method='linear')
                Yf, Xf = np.meshgrid(y_f, x_f, indexing='ij')
                field = interp(np.stack([Yf.ravel(), Xf.ravel()], axis=-1)).reshape(Yf.shape)
            ok_c = field >= 0.0
            ok = ok_c if ok is None else (ok & ok_c)
            ok_crit[c] = ok_c if ok_crit[c] is None else (ok_crit[c] & ok_c)
        ok_mp[rid] = ok
    ok_all = np.logical_and.reduce([ok_mp[rid] for rid in range(1, 7)])
    return dict(x=x_f, y=y_f, x_samples=x, y_samples=y, xv=xv, yv=yv, ok_mp=ok_mp, ok_all=ok_all, ok_crit=ok_crit)

def _scale(var, v):
    return v / 1e3 if var == 'sigma_fcs' else v

def _style_axes(ax):
    ax.spines[['right', 'top']].set_visible(False)
    ax.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False)

def plot_zone_panel(ax, zm, df_mass, hatch_crit=None):
    """2-D: as Fig. 4.8 (abscissa location, ordinate specific power or second location).
    1-D: one strip per mission point and one for the entire mission."""
    xv, yv = zm['xv'], zm['yv']
    if yv is not None:
        X, Y = np.meshgrid(_scale(xv, zm['x']), _scale(yv, zm['y']))
        ax.contourf(X, Y, zm['ok_all'].astype(float), levels=[-0.5, 0.5, 1.5], colors=[COLOR_FAIL, COLOR_OK], alpha=ZONE_ALPHA)
        if hatch_crit is not None and zm['ok_crit'].get(hatch_crit) is not None and np.any(~zm['ok_crit'][hatch_crit]):
            ax.contourf(X, Y, (~zm['ok_crit'][hatch_crit]).astype(float), levels=[0.5, 1.5], colors='none', hatches=['//'])
        for rid in range(1, 7):
            okf = zm['ok_mp'][rid].astype(float)
            if 0.0 < okf.mean() < 1.0:
                ax.contour(X, Y, okf, levels=[0.5], colors=MP_GREYS[rid - 1], linestyles=MP_LINESTYLES[rid - 1], linewidths=1.2)
        if SHOW_SAMPLES:
            ax.scatter(_scale(xv, df_mass[xv]), _scale(yv, df_mass[yv]), s=6, color='0.3', zorder=5, clip_on=False)
        ax.set_xlim(_scale(xv, zm['x'][0]), _scale(xv, zm['x'][-1]))
        ax.set_ylim(_scale(yv, zm['y'][0]), _scale(yv, zm['y'][-1]))
        ax.set_ylabel(AXIS_LABELS[yv])
    else:
        rows = list(range(1, 7)) + ['all']
        for k, rid in enumerate(rows):
            ok = zm['ok_all'] if rid == 'all' else zm['ok_mp'][rid]
            yk = len(rows) - 1 - k
            ax.fill_between(zm['x'], yk - 0.4, yk + 0.4, where=ok, color=COLOR_OK, alpha=ZONE_ALPHA if rid != 'all' else 3 * ZONE_ALPHA, linewidth=0, step='mid')
            ax.fill_between(zm['x'], yk - 0.4, yk + 0.4, where=~ok, color=COLOR_FAIL, alpha=ZONE_ALPHA if rid != 'all' else 3 * ZONE_ALPHA, linewidth=0, step='mid')
            if rid != 'all':
                ax.plot([zm['x'][0], zm['x'][-1]], [yk - 0.4, yk - 0.4], color=MP_GREYS[rid - 1], linestyle=MP_LINESTYLES[rid - 1], linewidth=1.2)
        if SHOW_SAMPLES:
            ax.scatter(df_mass[xv], np.zeros(len(df_mass)), s=6, color='0.3', zorder=5, clip_on=False)
        ax.set_yticks(range(len(rows)))
        ax.set_yticklabels(['Entire mission'] + MP_LABELS[::-1])
        ax.set_xlim(zm['x'][0], zm['x'][-1])
        ax.set_ylim(-0.6, len(rows) - 0.4)
    ax.set_xlabel(AXIS_LABELS[xv])
    _style_axes(ax)

def plot_zone_figure(aircraft, group, zone_maps, masses):
    """One panel per study (in the order of STUDIES) for one aircraft and one group of criteria"""
    studies = list(zone_maps.keys())
    n = len(studies)
    ncols = min(3, n)
    nrows = int(np.ceil((n + 1) / ncols))  # one slot for the legend
    fig, axes = plt.subplots(nrows, ncols, figsize=(6.4 * ncols, 5.2 * nrows), squeeze=False, layout='constrained')
    axes = axes.ravel()
    letters = 'abcdefghij'
    for k, s in enumerate(studies):
        plot_zone_panel(axes[k], zone_maps[s], masses[s], hatch_crit='delta_e' if group == 'static' else None)
        axes[k].set_title(f"({letters[k]}) {STUDY_INFO[s]['name']} (study {s})", loc='left', pad=12)
    for ax in axes[n:]:
        ax.axis('off')
    handles = [
        mpatches.Patch(facecolor=COLOR_OK, alpha=ZONE_ALPHA, label='Criteria met throughout mission'),
        mpatches.Patch(facecolor=COLOR_FAIL, alpha=ZONE_ALPHA, label='Criteria violated in at least one point'),
    ]
    if group == 'static':
        handles.append(mpatches.Patch(facecolor='none', edgecolor='k', hatch='//', label=f'Untrimmable ($\\delta_e < {DELTA_E_MIN:.0f}^\\circ$)'))
    handles += [mlines.Line2D([], [], color=MP_GREYS[i], linestyle=MP_LINESTYLES[i], linewidth=1.2, label=MP_LABELS[i]) for i in range(6)]
    axes[n].legend(handles=handles, loc='center', frameon=False, title=f"{aircraft.replace('_', ' ')}: {GROUP_TITLES[group].lower()}")
    if SAVE_FIGS:
        fig.savefig(os.path.join(OUT_DIR, f"stability_zones_{group}_{aircraft}.png"), dpi=600)
    return fig

# ---------------------------------------------------------------------
# Export of feasible locations for stability_zone_plotter.py
# ---------------------------------------------------------------------

def _intervals(x, ok):
    """Feasible intervals [[lo, hi], ...] of x along a boolean mask"""
    out, start = [], None
    for k, flag in enumerate(ok):
        if flag and start is None:
            start = k
        if (not flag or k == len(ok) - 1) and start is not None:
            end = k if flag else k - 1
            out.append([float(x[start]), float(x[end])])
            start = None
    return out

def export_zones(aircraft, all_zone_maps, masses):
    """Rows: group (static, longitudinal, lateral), study, location variable, specific power and the feasible
    intervals of the location (fraction of pressure shell length, semi-span or nacelle mass fraction) within the
    swept range [loc_min, loc_max]. Studies
    with a second location variable (5) are not exported. Geometry columns (x_nose, x_shell_fwd, x_shell_aft,
    half_span) are taken from the TASOPT results file if mass_distr_exec.jl wrote them (otherwise nan)."""
    rows = []
    for group, zone_maps in all_zone_maps.items():
        for s, zm in zone_maps.items():
            df_mass = masses[s]
            geom = {c: float(df_mass[c].iloc[0]) if c in df_mass.columns else np.nan for c in ['x_nose', 'x_shell_fwd', 'x_shell_aft', 'half_span']}
            geom.update(loc_min=float(zm['x'][0]), loc_max=float(zm['x'][-1]))  # swept range of the location variable
            if zm['yv'] == 'sigma_fcs':
                for sig in zm['y_samples']:
                    i = int(np.argmin(np.abs(zm['y'] - sig)))
                    rows.append(dict(aircraft=aircraft, group=group, study=s, loc_var=zm['xv'], sigma_fcs=float(sig),
                                     intervals=json.dumps(_intervals(zm['x'], zm['ok_all'][i, :])), **geom))
            elif zm['yv'] is None:
                rows.append(dict(aircraft=aircraft, group=group, study=s, loc_var=zm['xv'], sigma_fcs=float(df_mass['sigma_fcs'].iloc[0]),
                                 intervals=json.dumps(_intervals(zm['x'], zm['ok_all'])), **geom))
    df = pd.DataFrame(rows)
    path = os.path.join(OUT_DIR, f"stability_zones_{aircraft}.csv")
    df.to_csv(path, index=False)
    return path

# =============================================================================
# NILS: used as manual check of below uncommented version of same method
# that logic for choosing most critical flight condition works
# def plot_static_stability_lines_nils(df_mass, margins, study_idx):
    
#     fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
#     print('margins =', margins)
    
#     info = STUDY_INFO[study_idx]
#     xv, yv = info['x'], info['y']
#     x = np.sort(df_mass[xv].dropna().unique())
#     y = np.sort(df_mass[yv].dropna().unique()) if yv is not None else None
#     print('x =', x)
#     print('y =', y)
#     x_f = _fine_axis(x)
#     y_f = _fine_axis(y) if yv is not None else None
#     crit_cols = [c for c in margins.columns if c not in ('case', 'mp')]
#     print('crit_cols =', crit_cols)
#     ok_mp, ok_crit = {}, {c: None for c in crit_cols}
#     for rid in range(1, 7):
#         mrow = margins[margins['mp'] == rid].set_index('case')
#         # print('mrow =', mrow)
#         ok = None
#         for _i, c in enumerate(crit_cols):
#             vals = mrow[c].astype(float)
#             finite = vals.to_numpy()[np.isfinite(vals.to_numpy())]
#             if len(finite) == 0:
#                 fill = -1.0
#             elif finite.min() < 0.0:
#                 fill = finite.min()
#             else:
#                 fill = -0.1 * max(finite.max(), 1e-6)
#             vals = vals.fillna(fill)
            
#             # Arrange values on the (y, x) grid.
#             field = np.full((len(y), len(x)), np.nan)
#             for case, row in df_mass.iterrows():
#                 i = np.searchsorted(y, row[yv])
#                 j = np.searchsorted(x, row[xv])
#                 field[i, j] = vals.get(case, np.nan)
            
#             # Sweep over x at fixed y.
#             # for i, y_value in enumerate(y):
#             ax = ax1 if _i == 0 else ax2
#             for i in sorted({0, len(y) - 1}):
#                 y_value = y[i]
#                 ax.plot(
#                     _scale(xv, x), field[i, :], marker='o',
#                     label=f'mp {rid}, {c}, {yv}={_scale(yv, y_value):g}'
#                 )
            
#             # # Sweep over y at fixed x.
#             # for j, x_value in enumerate(x):
#             #     ax2.plot(
#             #         _scale(yv, y), field[:, j], marker='o',
#             #         label=f'mp {rid}, {c}, {xv}={_scale(xv, x_value):g}'
#             #     )
            
#             # print('vals =', vals.values)
            
#     #         # if c == 'SM':
                
#     #         # elif c == 'delta_e'
#     #         axes[_i].scatter(vals)
            
#     # plt.show()
            
#     ax1.set_xlabel(AXIS_LABELS[xv])
#     ax1.set_ylabel('Static-stability margin')
#     ax1.set_title('Sweep over x at fixed y')
#     # ax1.legend(fontsize=6, ncol=2)
    
#     ax2.set_xlabel(AXIS_LABELS[yv])
#     ax2.set_ylabel('Static-stability margin')
#     ax2.set_title('Sweep over y at fixed x')
#     # ax2.legend(fontsize=6, ncol=2)
    
#     fig.tight_layout()
#     return fig


# NILS: used to generate figures in Chapter 7 of thesis
def plot_static_stability_lines_nils(df_mass, margins, study_idx):

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7, 4), sharey=True)

    info = STUDY_INFO[study_idx]
    xv, yv = info['x'], info['y']
    x = np.sort(df_mass[xv].dropna().unique())
    y = np.sort(df_mass[yv].dropna().unique())

    crit_cols = [c for c in margins.columns if c not in ('case', 'mp')]

    # Track the minimum margin and its corresponding mission point.
    min_fields = {
        c: np.full((len(y), len(x)), np.nan)
        for c in crit_cols
    }
    active_mp = {
        c: np.full((len(y), len(x)), -1, dtype=int)
        for c in crit_cols
    }

    for rid in range(1, 7):
        mrow = margins[margins['mp'] == rid].set_index('case')

        for c in crit_cols:
            vals = mrow[c].astype(float)
            finite = vals.to_numpy()[np.isfinite(vals.to_numpy())]

            if len(finite) == 0:
                fill = -1.0
            elif finite.min() < 0.0:
                fill = finite.min()
            else:
                fill = -0.1 * max(finite.max(), 1e-6)

            vals = vals.fillna(fill)

            # Arrange values on the (y, x) grid.
            field = np.full((len(y), len(x)), np.nan)

            for case, row in df_mass.iterrows():
                if pd.isna(row[xv]) or pd.isna(row[yv]):
                    continue

                i = np.searchsorted(y, row[yv])
                j = np.searchsorted(x, row[xv])
                field[i, j] = vals.get(case, np.nan)

            # Update the minimum and record its active mission point.
            current_min = min_fields[c]
            update = (
                np.isfinite(field)
                & (
                    ~np.isfinite(current_min)
                    | (field < current_min)
                )
            )

            current_min[update] = field[update]
            active_mp[c][update] = rid

    # Plot SM and delta_e on the same axes for each fixed value of y.
    # The margins already contain the differences from their respective limits,
    # so dividing by the limit gives the requested normalisation.
    criterion_specs = [
        (crit_cols[0], SM_min, '$\mathrm{SM_{min}}$'),
        (crit_cols[1], DELTA_E_MIN, '$\delta_\mathrm{e,min}$'),
    ]

    for i in sorted({0, len(y) - 1}):

        ax = ax1 if i == 0 else ax2

        for c, threshold, label in criterion_specs:

            # Normalised criterion margin; the threshold is zero.
            # plot_field = min_fields[c] / threshold
            ###
            plot_field = min_fields[c] / abs(threshold)
            ###

            x_plot = []
            y_plot = []

            for j in range(len(x)):

                x_plot.append(_scale(xv, x[j]))
                y_plot.append(plot_field[i, j])

                # Break the line when the active mission point changes.
                if (
                    j < len(x) - 1
                    and active_mp[c][i, j] != active_mp[c][i, j + 1]
                ):
                    x_plot.append(np.nan)
                    y_plot.append(np.nan)

            ax.plot(
                x_plot,
                y_plot,
                marker='.',
                label=label,
                clip_on=False,
            )

        # # Plot the normalised stability threshold.
        # ax.axhline(
        #     0.0,
        #     color='0.3',
        #     linestyle='--',
        #     linewidth=1.2,
        #     label='Limit = 0',
        #     zorder=3
        # )

        # Ensure that zero is within the plotted y-axis limits.
        ymin, ymax = ax.get_ylim()
        pad = 0.05 * max(ymax - ymin, 1.0)

        if ymin >= 0.0:
            ymin = -pad
        if ymax <= 0.0:
            ymax = pad

        ax.set_ylim(ymin, ymax)

        # Shade x-intervals where both normalised criteria are positive.
        both_positive = np.ones(len(x), dtype=bool)

        for c, threshold, label in criterion_specs:
            plot_field = min_fields[c] / abs(threshold)
            both_positive &= np.isfinite(plot_field[i]) & (plot_field[i] > 0.0)

        x_plot = [_scale(xv, xj) for xj in x]

        for j in range(len(x) - 1):
            if both_positive[j] and both_positive[j + 1]:
                ax.axvspan(
                    x_plot[j], x_plot[j + 1],
                    # color=COLOR_OK,
                    # alpha=ZONE_ALPHA,
                    # zorder=0
                    hatch='///', facecolor='None', edgecolor='lightgrey',
                    clip_on=False,
                )

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

        ax.set_xlabel('Position along fuselage [-]', labelpad=15)
        if i == 0:
            ax.set_ylabel('$(x-x_\mathrm{min})/x_\mathrm{min}$', labelpad=15)
            ax.legend(frameon=False, fontsize=14)
        # ax.set_title(f'Sweep over x at fixed y = {_scale(yv, y[i]):g}')
        ax.set_title('Conservative system' if i == 0 else 'Advanced system', pad=15)
        
    fig.tight_layout()
    return fig
# =============================================================================

# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def run(aircraft_list=None, studies=None, show=True):
    aircraft_list = list(AIRCRAFT.keys()) if aircraft_list is None else aircraft_list
    studies = STUDIES + [s for s in (3, 4) if s not in STUDIES] if studies is None else studies
    figs = {}
    for aircraft in aircraft_list:
        masses, static, dynamic = {}, {}, {}
        for s in studies:
            df_mass = read_mass_csv(aircraft, s)
            df_static = load_static(aircraft, s)
            if df_mass is None or df_static is None:
                continue
            if df_static['case'].nunique() != len(df_mass):
                print(f"{aircraft}, study {s}: {df_static['case'].nunique()} of {len(df_mass)} cases found, skipped")
                continue
            masses[s], static[s] = df_mass, df_static
            dynamic[s] = load_dynamic(aircraft, s)
        if len(masses) == 0:
            print(f"{aircraft}: no complete studies found")
            continue
        all_zone_maps = {}
        # for group in ['static', 'longitudinal', 'lateral']:
        #####
        for group in ['static']:
        #####
            zone_maps = {}
            # print('masses =', masses)
            for s in masses:
                if group == 'static':
                    margins = static_margins(static[s])
                elif dynamic[s] is not None:
                    margins = dynamic_margins(dynamic[s], group)
                else:
                    continue
                zone_maps[s] = zone_map(masses[s], margins, s)
                # =============================================================================
                if s == 7:
                    figs[(aircraft, group)] = plot_static_stability_lines_nils(masses[s], margins, s)
                # =============================================================================
            all_zone_maps[group] = zone_maps
            # figs[(aircraft, group)] = plot_zone_figure(aircraft, group, zone_maps, masses)
        print("Feasible locations written to", export_zones(aircraft, all_zone_maps, masses))
    if show:
        plt.show()
    return figs


if __name__ == "__main__":
    try:
        from nils.matplotlib_custom_settings import *  # same fonts as Figs. 4.9 and 4.10
    except ImportError:
        pass
    run()
