"""
This script implements a physics-based eigenvalue-to-eigenmode
association.
"""

__all__ = ["match_evals_to_emodes"]

import os
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import linear_sum_assignment

from nils.matplotlib_custom_settings import *

# ---------------------------
# Parsers and small helpers
# ---------------------------
def _parse_matrix_file(path):
    lines = open(path, "r", errors="ignore").read().splitlines()
    groups = []
    i, N = 0, len(lines)
    while i < N:
        ln = lines[i].strip().lower()
        if ln.startswith("u w q") or ln.startswith("u  w  q") or ln.startswith("u w q the"):
            i += 1
            mats = []
            for _ in range(6):
                while i < N and lines[i].strip() == "":
                    i += 1
                if i >= N: break
                rows = []
                while i < N and len(rows) < 12:
                    s = lines[i].strip(); i += 1
                    if s == "": continue
                    toks = s.split(); nums = []
                    for t in toks:
                        try: nums.append(float(t))
                        except: pass
                    if len(nums) >= 12:
                        rows.append(nums)
                if len(rows) < 12: break
                arr = np.array(rows, dtype=float)
                A = arr[:, :12]; B = arr[:, 12:] if arr.shape[1] > 12 else None
                mats.append({'A': A, 'B': B})
            if mats:
                groups.append(mats)
        else:
            i += 1
    return groups

def _parse_eigfile_row_per_group(path, row_index=0):
    assert 0 <= row_index <= 5
    lines = open(path, "r", errors="ignore").read().splitlines()
    def _is_header(s): return s.startswith("sigma_fcs") and "Re1" in s and "Im8" in s
    parsed = []
    i, N = 0, len(lines)
    while i < N:
        s = lines[i].strip()
        if _is_header(s):
            numeric_rows = []
            blank_count = 0
            j = i + 1
            while j < N:
                s2 = lines[j].strip()
                if _is_header(s2):
                    break
                if s2 == "":
                    blank_count += 1
                    if blank_count >= 2:
                        j += 1
                        break
                    j += 1
                    continue
                else:
                    blank_count = 0
                toks = s2.split(); nums = []
                for t in toks:
                    try: nums.append(float(t))
                    except: pass
                if len(nums) >= 16:
                    numeric_rows.append(nums)
                j += 1
            if len(numeric_rows) > row_index:
                row = numeric_rows[row_index]
                re_start = len(row) - 16
                parsed.append({'params': row[:re_start], 'raw': row, 'row_index': row_index})
            i = j
        else:
            i += 1
    return parsed

# ---------------------------
# Modal analysis primitives
# ---------------------------
def _eig_from_A(A):
    vals, vecs = np.linalg.eig(A)
    return np.array(vals, dtype=complex), np.array(vecs, dtype=complex)

def _normalize_vec(v):
    v = v.copy().astype(complex)
    m = np.max(np.abs(v))
    if m == 0: return v
    v /= m
    k = int(np.argmax(np.abs(v)))
    phase = np.angle(v[k])
    return v * np.exp(-1j * phase)

def _participation(v):
    vabs = np.abs(v.flatten())
    lon = np.sum(vabs[[0,1,2,3]]); lat = np.sum(vabs[[4,5,6,7]]); oth = np.sum(vabs[8:12])
    s = lon + lat + oth + 1e-12
    return {'lon': lon/s, 'lat': lat/s, 'oth': oth/s, 'comp': vabs}

def _mode_metrics(lam):
    sigma = lam.real; wd = abs(lam.imag)
    wn = math.hypot(sigma, wd)
    zeta = -sigma / wn if wn > 0 else float('nan')
    return {'lambda': lam, 'sigma': sigma, 'wd': wd, 'wn': wn, 'zeta': zeta, 'freq': wd / (2*math.pi)}

def _MAC(a, b):
    num = abs(np.vdot(a, b))**2
    den = (np.vdot(a, a).real) * (np.vdot(b, b).real)
    return float(num/den) if den > 0 else 0.0

# ---------------------------
# Classification and tracking
# ---------------------------
def _classify_group_from_A(A):
    """
    Minimal edit here: skip kinematic integrator eigenpairs (those dominated by x,y,z,psi
    components and numerically zero eigenvalue). This removes the 4 zero eigenpairs that
    are not physical dynamic modes.
    """
    vals, vecs = _eig_from_A(A)
    used = set()
    reps = []
    # threshold for treating eigenvalue as 'zero' and for kinematic-dominance
    ZERO_EIG_TOL = 1e-8
    KIN_FRAC_TOL = 0.7  # if >70% of eigenvector energy in indices 8..11 -> kinematic

    for i, lam in enumerate(vals):
        # check raw eigenvector kinematic fraction before normalizing
        v_raw = vecs[:, i]
        kin_frac = np.sum(np.abs(v_raw[8:12])) / (np.sum(np.abs(v_raw)) + 1e-16)
        if abs(lam) < ZERO_EIG_TOL and kin_frac > KIN_FRAC_TOL:
            # skip kinematic integrator eigenpair entirely
            continue

        if i in used:
            continue

        if abs(lam.imag) < 1e-8:
            v = _normalize_vec(vecs[:, i])
            reps.append((lam, v)); used.add(i)
        else:
            conj = np.conj(lam); j = None
            for k in range(len(vals)):
                if k in used or k == i: continue
                if abs(vals[k] - conj) < 1e-6:
                    # also check and skip if the partner is kinematic-zero dominated
                    v_k_raw = vecs[:, k]
                    kin_frac_k = np.sum(np.abs(v_k_raw[8:12])) / (np.sum(np.abs(v_k_raw)) + 1e-16)
                    if abs(vals[k]) < ZERO_EIG_TOL and kin_frac_k > KIN_FRAC_TOL:
                        continue
                    j = k; break
            if j is None:
                v = _normalize_vec(vecs[:, i])
                reps.append((lam, v)); used.add(i)
            else:
                # choose positive-imag representative
                if lam.imag > 0:
                    reps.append((lam, _normalize_vec(vecs[:, i])))
                else:
                    reps.append((vals[j], _normalize_vec(vecs[:, j])))
                used.add(i); used.add(j)

    classified = []
    for lam, v in reps:
        part = _participation(v)
        met = _mode_metrics(lam)
        if abs(lam.imag) < 1e-8:
            label = 'real'
        else:
            if part['lon'] >= 0.6:
                label = 'longitudinal'
            elif part['lat'] >= 0.5:
                label = 'dutch-roll'
            else:
                label = 'longitudinal' if abs(v[2]) > abs(v[6]) else 'dutch-roll'
        classified.append({'lambda': lam, 'vec': v, 'part': part, 'metrics': met, 'label': label})

    # longitudinal split (unchanged)
    longitudinal = [c for c in classified if c['label'] == 'longitudinal']
    longitudinal = sorted(longitudinal, key=lambda x: x['metrics']['wd'], reverse=True)
    if len(longitudinal) >= 2:
        longitudinal[0]['label'] = 'short-period'
        longitudinal[-1]['label'] = 'phugoid'
        for mid in longitudinal[1:-1]:
            mid['label'] = 'short-period'
    elif len(longitudinal) == 1:
        freq = longitudinal[0]['metrics']['freq']
        longitudinal[0]['label'] = 'short-period' if freq > 0.15 else 'phugoid'

    # reals -> roll/spiral with p-vs-phi heuristic but robust fallback
    reals = [c for c in classified if abs(c['lambda'].imag) < 1e-8]
    if len(reals) >= 2:
        p_vals = np.array([abs(r['vec'][5]) for r in reals])
        phi_vals = np.array([abs(r['vec'][7]) for r in reals])
        p_idx = int(np.argmax(p_vals)); phi_idx = int(np.argmax(phi_vals))
        p_decisive = p_vals[p_idx] > 1.1 * phi_vals[phi_idx]
        phi_decisive = phi_vals[phi_idx] > 1.1 * p_vals[p_idx]
        if (p_idx != phi_idx) and (p_decisive or phi_decisive):
            reals[p_idx]['label'] = 'roll'; reals[phi_idx]['label'] = 'spiral'
            unlabeled = [r for r in reals if r.get('label','real') == 'real']
            if unlabeled:
                unlabeled_sorted = sorted(unlabeled, key=lambda x: x['lambda'].real)
                if len(unlabeled_sorted) == 1:
                    unlabeled_sorted[0]['label'] = 'spiral' if unlabeled_sorted[0]['lambda'].real > reals[p_idx]['lambda'].real else 'roll'
                else:
                    unlabeled_sorted[0]['label'] = 'roll'; unlabeled_sorted[-1]['label'] = 'spiral'
        else:
            reals_sorted = sorted(reals, key=lambda x: x['lambda'].real)
            reals_sorted[0]['label'] = 'roll'; reals_sorted[-1]['label'] = 'spiral'
            if len(reals_sorted) > 2:
                for r in reals_sorted[1:-1]:
                    r['label'] = 'roll' if abs(r['lambda'].real - reals_sorted[0]['lambda'].real) < abs(r['lambda'].real - reals_sorted[-1]['lambda'].real) else 'spiral'
    elif len(reals) == 1:
        r = reals[0]; r['label'] = 'roll' if r['part']['lat'] > r['part']['oth'] and r['part']['lat'] > 0.2 else 'spiral'

    for c in classified:
        if 'label' not in c: c['label'] = 'unknown'
    return classified

def _track_groups(classified_groups):
    canonical = ['roll', 'short-period', 'dutch-roll', 'phugoid', 'spiral']
    n = len(classified_groups)
    tracks = {lab: [None]*n for lab in canonical}
    g0 = classified_groups[0]
    # print('canonical, g0 =', canonical, [m['label'] for m in g0])
    # for lab in canonical:
    #     found = next((m for m in g0 if m['label'] == lab), None)
    #     if found is None:
    #         raise RuntimeError(f"Initial group missing canonical label '{lab}'; cannot proceed.")
    #     tracks[lab][0] = found
    # =============================================================================
    for lab in canonical:
        found = next((m for m in g0 if m['label'] == lab), None)
    
        if found is None:
            # fallback: pick best available candidate of similar type
            if lab in ['roll', 'spiral', 'dutch-roll']:
                candidates = [m for m in g0 if m['label'] in ['roll', 'spiral', 'dutch-roll']]
            elif lab in ['short-period', 'phugoid']:
                candidates = [m for m in g0 if m['label'] in ['short-period', 'phugoid', 'longitudinal']]
            else:
                candidates = g0
    
            if candidates:
                found = candidates[0]  # just take one, tracking will sort it out
            else:
                found = g0[0]  # absolute fallback (never crash)
    
        tracks[lab][0] = found
    # =============================================================================
    for g in range(1, n):
        prevs = [tracks[lab][g-1] for lab in canonical]
        cur = classified_groups[g]
        m = len(canonical); q = len(cur)
        cost = np.full((m, q), 1e6)
        for i, prev in enumerate(prevs):
            for j, cand in enumerate(cur):
                mac = _MAC(prev['vec'], cand['vec'])
                lamdist = abs(prev['lambda'] - cand['lambda']) / max(1.0, abs(prev['lambda']), abs(cand['lambda']))
                penalty = 0.0 if cand['label'] == canonical[i] else 0.2
                cost[i, j] = (1.0 - mac) + 0.5 * lamdist + penalty
        rows, cols = linear_sum_assignment(cost)
        for r, c in zip(rows, cols):
            if cost[r, c] > 1e5:
                raise RuntimeError(f"Assignment cost too large at group {g}: row {r} col {c}")
            tracks[canonical[r]][g] = cur[c]
    return tracks

# ---------------------------
# Main exported function
# ---------------------------
def match_evals_to_emodes(matfile, eigfile, plot_matching=True):
    assert os.path.exists(matfile) and os.path.exists(eigfile), "Input files not found."
    mat_groups = _parse_matrix_file(matfile)
    if len(mat_groups) == 0:
        raise RuntimeError("No matrix groups parsed.")
    total_groups = len(mat_groups)
    # if total_groups != 210:
    #     raise RuntimeError(f"Expected 210 groups; found {total_groups}")
    def _category_for_index(g):
        # =============================================================================
        # 26.03.2026
        if g < 100: return 'fuselage'
        # if g < 100: return 'wing'
        # =============================================================================
        if g < 200: return 'wing'
        return 'nacelle'
    categories = ['fuselage', 'wing', 'nacelle']
    # modes = ['roll', 'short-period', 'dutch-roll', 'phugoid', 'spiral']
    modes = ['short-period', 'phugoid', 'roll', 'spiral', 'dutch-roll']
    result = {cat: {m: {} for m in modes} for cat in categories}
    plot_data = {ridx: {m: {'lon': [], 'lat': [], 'freq': [], 'zeta': [], 'gindex': [], 'params': []} for m in modes} for ridx in range(6)}

    for ridx in range(6):
        eig_rows = _parse_eigfile_row_per_group(eigfile, row_index=ridx)
        if len(eig_rows) != len(mat_groups):
            raise RuntimeError(f"Mismatch groups: mat_groups={len(mat_groups)} vs eig_rows for row {ridx}={len(eig_rows)}")
        n_groups = len(mat_groups)
        classified_groups = []
        for g in range(n_groups):
            A = mat_groups[g][ridx]['A']
            classified_groups.append(_classify_group_from_A(A))
        tracks = _track_groups(classified_groups)
        for g in range(n_groups):
            params = eig_rows[g].get('params', [])
            params_first5 = list(params)[:5] + [np.nan] * max(0, 5 - len(params))
            for m in modes:
                t = tracks[m][g]
                if t is None: continue
                plot_data[ridx][m]['lon'].append(t['part']['lon'])
                plot_data[ridx][m]['lat'].append(t['part']['lat'])
                plot_data[ridx][m]['freq'].append(t['metrics']['freq'])
                plot_data[ridx][m]['zeta'].append(t['metrics']['zeta'])
                plot_data[ridx][m]['gindex'].append(g)
                plot_data[ridx][m]['params'].append(params_first5)
        cols = ['sigma_fcs','span_loc','fcs_loc','wing_frac','nacelle_frac','AoA','Mach','Beta','h','n','W','Re','Im']
        result_data = {cat: {m: [] for m in modes} for cat in categories}
        for g in range(n_groups):
            params = eig_rows[g]['params']
            params_padded = list(params)[:11] + [np.nan] * max(0, 11 - len(params))
            cat = _category_for_index(g)
            for m in modes:
                entry = params_padded.copy()
                assigned = tracks[m][g]
                if assigned is None:
                    entry += [np.nan, np.nan]
                else:
                    lam = assigned['lambda']
                    entry += [float(lam.real), float(lam.imag)]
                result_data[cat][m].append(entry)
        for cat in categories:
            for m in modes:
                df = pd.DataFrame(result_data[cat][m], columns=cols)
                result[cat][m][ridx+1] = df

    # Diagnostic plots (unchanged)
    colors_dict = {
        'roll':colors[0],
        'short-period':colors[1],
        'dutch-roll':colors[2],
        'phugoid':colors[3],
        'spiral':colors[4],
    }
    
    if plot_matching == True:
    
        titles = ['1. Take-off', '2. Climb', '3. Beginning of cruise', '4. End of cruise', '5. Descent', '6. Landing']
        markers = ['o', 's', 'o', 's', 'd']
        marker_colors = [colors[0], colors[0], colors[1], colors[1], colors[1]]
        
        # NILS: plot longitudinal and lateral contributions
        # of various eigenmodes to rows in system matrix A
        
        fig = plt.figure(figsize=(15,9))
        gs = gridspec.GridSpec(
            2, 3,
            figure=fig,
            left=0.075,
            right=0.75,
            bottom=0.075,
            top=0.9,
            hspace=0.1,
            wspace=0.2,
        )
        axes = [fig.add_subplot(gs[_i,_j]) for _i in range(2) for _j in range(3)]
        
        for ridx, ax in enumerate(axes[:6]):
            for mode_counter, m in enumerate(modes):
                ax.scatter(
                    plot_data[ridx][m]['lon'],
                    plot_data[ridx][m]['lat'],
                    marker=markers[mode_counter],
                    alpha=0.8,
                    label=m,
                    s = 40,
                    facecolor=marker_colors[mode_counter],
                    edgecolor='k',
                    linewidth=0.5,
                    clip_on=False,
                )
        
            ax.plot([ax.get_xlim()[0], 1], [ax.get_ylim()[0], 1], color='k', alpha=0.2)
        
            # Axis limits and scales
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_xlim(1e-9, 1)
            ax.set_ylim(1e-7, 1)
            ax.set_aspect('equal')
        
            # Move spines to x=1 and y=1
            ax.spines['bottom'].set_position(('data', 1.0))
            ax.spines['left'].set_position(('data', 1.0))
            ax.spines[['right','top']].set_visible(False)
        
            if (ridx == 0 or ridx == 1 or ridx == 2):
                ax.xaxis.set_label_position('top')
                ax.tick_params(
                    axis='x',
                    which='both',
                    top=True,
                    labeltop=True,
                    bottom=False,
                    labelbottom=False,
                    length=0
                )
            else:
                ax.tick_params(
                    axis='x',
                    which='both',
                    top=False,
                    labeltop=False,
                    bottom=False,
                    labelbottom=False,
                    length=0
                )
        
            if (ridx == 2 or ridx == 5):
                ax.yaxis.set_label_position('right')
                ax.tick_params(
                    axis='y',
                    which='both',
                    right=True,
                    labelright=True,
                    left=False,
                    labelleft=False,
                    length=0
                )
            else:
                ax.tick_params(
                    axis='y',
                    which='both',
                    right=False,
                    labelright=False,
                    left=False,
                    labelleft=False,
                    length=0
                )
            
            if (ridx == 0 or ridx == 1 or ridx == 2):
                ax.set_title(titles[ridx], pad=20, y=-0.3)
            elif (ridx == 3 or ridx == 4 or ridx == 5):
                ax.set_title(titles[ridx], pad=20)
            
        legend_ax = fig.add_axes([0.8, 0.0, 0.2, 1.0])
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
        legend_elements = [
            Line2D([], [], marker='o', markerfacecolor=colors[0], markeredgecolor='k', linestyle='', markeredgewidth=0.5, label='Short period'),
            Line2D([], [], marker='s', markerfacecolor=colors[0], markeredgecolor='k', linestyle='', markeredgewidth=0.5, label='Phugoid'),
            Line2D([], [], marker='o', markerfacecolor=colors[1], markeredgecolor='k', linestyle='', markeredgewidth=0.5, label='Roll subsidence'),
            Line2D([], [], marker='s', markerfacecolor=colors[1], markeredgecolor='k', linestyle='', markeredgewidth=0.5, label='Spiral'),
            Line2D([], [], marker='d', markerfacecolor=colors[1], markeredgecolor='k', linestyle='', markeredgewidth=0.5, label='Dutch roll'),
        ]
        legend_ax.legend(
            handles=legend_elements, loc='center right', bbox_to_anchor=(1.0, 0.5), ncol=1, labelspacing=3.0, frameon=False
        )
        
        fig.text(0.45, 0.05, 'Longitudinal weighting (-)', ha='center', va='bottom')
        fig.text(0.05, 0.5, 'Lateral weighting (-)', ha='right', va='center', rotation='vertical')
        
        # plt.savefig('long_lat_weightings_log.png', format='png', dpi=600)
        plt.show()
        '''
        # NILS: damped frequency or damping ratio for
        # different FCS integration options (nacelle,
        # fuselage, wing) on grid of design parameters
        # (sigma_fcs vs fcs_loc for fuselage, sima_fcs
        # vs span_loc for wing, 1D line for nacelle
        
        # User input
        y_var_plot = 'freq'
        # y_var_plot = 'zeta'

        fig = plt.figure(figsize=(15,9))
        gs = gridspec.GridSpec(
            2, 3,
            figure=fig,
            left=0.125,
            right=0.85,
            bottom=0.075,
            top=0.9,
            hspace=0.3,
            wspace=0.3,
        )
        axes = [fig.add_subplot(gs[_i,_j]) for _i in range(2) for _j in range(3)]
        
        for ridx, ax in enumerate(axes[:6]):
            row, col = divmod(ridx, 3)
            
            gs_inner = gridspec.GridSpecFromSubplotSpec(
                1, 3,
                subplot_spec=gs[row, col],
                hspace=0.1,
            )
            ax_fuse = fig.add_subplot(gs_inner[0])
            ax_wing = fig.add_subplot(gs_inner[1], sharey=ax_fuse)
            ax_nacelle = fig.add_subplot(gs_inner[2], sharey=ax_fuse)
            
            for _ax in [ax_wing, ax_nacelle]:
                _ax.spines[['left', 'right', 'bottom', 'top']].set_visible(False)
                _ax.tick_params(
                    axis='x',
                    which='both',
                    right=False,
                    labelright=False,
                    left=False,
                    labelleft=False,
                    length=0
                )
                _ax.tick_params(
                    axis='y',
                    which='both',
                    right=False,
                    labelright=False,
                    left=False,
                    labelleft=False,
                    length=0
                )
                _ax.set_xticks([])
                _ax.set_yticks([])
            
            ax_fuse.spines[['right', 'bottom', 'top']].set_visible(False)
            ax_fuse.tick_params(
                axis='x',
                which='both',
                right=False,
                labelright=False,
                left=False,
                labelleft=False,
                length=0
            )
            ax_fuse.tick_params(
                axis='y',
                which='both',
                right=False,
                labelright=False,
                left=False,
                labelleft=False,
                length=0
            )
            ax_fuse.set_xticks([])
            
            for m in modes:
                
                x_fuse = np.arange(len(plot_data[ridx][m][y_var_plot]))[:100].reshape(10, 10)
                y_fuse = np.array(plot_data[ridx][m][y_var_plot][:100]).reshape(10, 10)
                
                for i in range(10):
                    ax_fuse.plot(x_fuse[i, :], y_fuse[i, :], color=colors[0], lw=0.6, alpha=0.6)
                    ax_fuse.plot(x_fuse[:, i], y_fuse[:, i], color=colors[0], lw=0.6, alpha=0.6)
                
                
                x_wing = np.arange(len(plot_data[ridx][m][y_var_plot]))[100:200].reshape(10, 10)
                y_wing = np.array(plot_data[ridx][m][y_var_plot][100:200]).reshape(10, 10)
                
                for i in range(10):
                    ax_wing.plot(x_wing[i, :], y_wing[i, :], color=colors[1], lw=0.6, alpha=0.6)
                    ax_wing.plot(x_wing[:, i], y_wing[:, i], color=colors[1], lw=0.6, alpha=0.6)
                
                ax_nacelle.plot(
                    np.arange(len(plot_data[ridx][m][y_var_plot]))[200:],
                    plot_data[ridx][m][y_var_plot][200:],
                    color=colors[2],
                    lw=0.6,
                )
            
            ax.spines[['right', 'bottom', 'top']].set_visible(False)
            ax.tick_params(
                axis='x',
                which='both',
                right=False,
                labelright=False,
                left=False,
                labelleft=False,
                length=0
            )
            ax.tick_params(
                axis='y',
                which='both',
                right=False,
                labelright=False,
                left=True,
                labelleft=True,
                length=0
            )
            ax.set_xticks([])
            ax.set_title(titles[ridx], pad=20)
            
        legend_ax = fig.add_axes([0.9, 0.0, 0.1, 1.0])
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
        legend_elements = [
            Line2D([], [], color=colors[0], linestyle='-', label='Fuselage'),
            Line2D([], [], color=colors[1], linestyle='-', label='Wing'),
            Line2D([], [], color=colors[2], linestyle='-', label='Nacelle'),
        ]
        legend_ax.legend(
            handles=legend_elements, loc='center right', bbox_to_anchor=(1.0, 0.5), ncol=1, labelspacing=3.0, frameon=False
        )
        
        if y_var_plot == 'freq':
            fig.text(0.05, 0.5, 'Damped frequency, $\omega_d$ (1/s)', ha='left', va='center', rotation='vertical')
            # plt.savefig('damped_frequency.png', format='png', dpi=600)
        elif y_var_plot == 'zeta':
            fig.text(0.05, 0.5, 'Damping ratio, $\zeta$  (-)', ha='left', va='center', rotation='vertical')
            # plt.savefig('damping_ratio.png', format='png', dpi=600)
        
        plt.show()
        '''
    return result

#%%

if __name__ == "__main__":

    # If run as script, small demonstration

    matfile = r"C:\Users\nmb48\Documents\GitHub\SUAVE\Tutorials-2.5.0.2\B737_AVL_Tutorial\CUsersnmb48\suave_dynamic_stability_matrix_combined.txt"
    eigfile = r"C:\Users\nmb48\Documents\GitHub\SUAVE\Tutorials-2.5.0.2\B737_AVL_Tutorial\CUsersnmb48\suave_dynamic_stability_outputs_combined.txt"
    res = match_evals_to_emodes(matfile, eigfile, plot_matching=True)
