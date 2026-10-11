r"""
blown_wing_layout_claude.py

CLAUDE (04.10.2026): parametric propeller layout on a FIXED airframe, for the JVL
verification of the blown-wing surrogate (Step 3 of
jvl_blown_wing_verification_plan_041026.md, issue B6).

Takes one geometry .csv written by TASOPT.jl-priv (SuaveGeometryGetter.get_suave_geometry,
e.g. TASOPT.jl-priv\mass_distr_geometry_regional_LH2_041026.csv, which mass_distr_exec.jl
exports from atr72600_fc_retrofit_full_nsub2_3527_0.0012_3.2e-5.jld2) and returns a modified
copy in which only the propeller-related entries are replaced:

    sect_* columns (prop_k_in / prop_k_out / mid_pitch_k_(k+1) sections re-generated),
    nacelle_origin, Dprop, Jgain, dY, dy

so that the result can be passed unchanged to nils.jvl.code.jvl_wrapper.wrap_jvl.

Layout rule (chosen so that eq. (21) of de Vries et al. (2019) holds EXACTLY, i.e. the
x-axis of the target figure is unambiguous, issue B18):

    dy    = dY * b / (N * D) - 1            (spacing as a fraction of D; infeasible if < 0)
    pitch = D * (1 + dy) = dY * b / N
    N/2 propellers per side occupy the band [y0, y0 + dY * b/2], centres at
    y_k = y0 + pitch/2 + k * pitch,  k = 0 .. N/2 - 1,
    y0  = min(body half-width, b/2 * (1 - dY))   (band starts at the fuselage side and is
          moved inboard only if it would otherwise extend beyond the tip; the overlap with
          the fuselage is reported as `y0_inside_fuselage_m`)

Propeller disk centre (nacelle_origin):
    x = x_LE(y_k) - ax_pos * c(y_k)   -> JVL dxdisk = -ax_pos * c   (surrogate: x/c = ax_pos)
    z = z_LE(y_k)                     -> JVL dndisk = 0             (surrogate: no offset, B11)

Section properties are interpolated on the planform sections 'center', 'inboard',
'outboard', 'tip' (wing_<name>_* columns): chord ratio and t/c linearly in eta, dihedral
and quarter-chord sweep piecewise-constant from the inboard planform section, so the
planform itself is unchanged.

Pure numpy/pandas - no SUAVE import - so it can be tested stand-alone.
"""

__all__ = ["make_layout", "base_mass_row", "planform_sections", "check_layout"]

import ast
import math
import re

import numpy as np
import pandas as pd

_PLANFORM = ("center", "inboard", "outboard", "tip")


def planform_sections(df_geom):
    """Planform sections of the main wing as a list of dicts sorted by eta."""
    secs = []
    for name in _PLANFORM:
        dih_key = f"wing_{name}_dihedral_outboard" if name != "tip" else "wing_tip_dihedral_tip"
        secs.append({
            "label": name,
            "eta": float(df_geom[f"wing_{name}_percent_span_location"][0]),
            "rcp": float(df_geom[f"wing_{name}_root_chord_percent"][0]),
            "tc": float(df_geom[f"wing_{name}_thickness_to_chord"][0]),
            "dih": float(df_geom[dih_key][0]),
            "sweep": float(df_geom[f"wing_{name}_sweeps_quarter_chord"][0]),
        })
    return sorted(secs, key=lambda s: s["eta"])


def _inboard_section(secs, eta):
    """Planform section at or immediately inboard of eta."""
    cand = [s for s in secs if s["eta"] <= eta + 1e-12]
    return cand[-1] if cand else secs[0]


def _interp(secs, eta, key):
    etas = [s["eta"] for s in secs]
    vals = [s[key] for s in secs]
    return float(np.interp(eta, etas, vals))


def _le_xz(df_geom, secs, eta):
    """Leading-edge x and z at eta, consistent with SUAVE's segment construction
    (quarter-chord line swept by each segment's quarter-chord sweep, dihedral per segment).
    Sweep and dihedral are taken in radians (SUAVE convention)."""
    b = float(df_geom["wing_spans_projected"][0])
    c_root = float(df_geom["wing_chords_root"][0])
    x0 = float(df_geom["wing_origin_x"][0])
    z0 = float(df_geom["wing_origin_z"][0])
    semi = b / 2.0
    x_qc = x0 + 0.25 * c_root
    z = z0
    for s_in, s_out in zip(secs[:-1], secs[1:]):
        if eta <= s_in["eta"]:
            break
        eta_end = min(eta, s_out["eta"])
        dyy = (eta_end - s_in["eta"]) * semi
        x_qc += dyy * math.tan(s_in["sweep"])
        z += dyy * math.tan(s_in["dih"])
    c = c_root * _interp(secs, eta, "rcp")
    return x_qc - 0.25 * c, z, c


def make_layout(df_geom_base, N, D, dY, ax_pos=0.2, Jgain=1.0, tip_margin=0.0, eta_merge_tol=1e-6,
                scale_nacelles=True):
    """Return (df_geom, info). df_geom is None if the layout is infeasible.

    Parameters
    ----------
    df_geom_base : DataFrame read from a TASOPT geometry .csv (3 rows: nacelle_origin x/y/z)
    N : total number of propellers (even)
    D : propeller diameter (m)
    dY : blown span fraction (band width per side / semi-span), as in de Vries (2019)
    ax_pos : propeller axial position ahead of the LE as a fraction of the local chord
    Jgain : jet velocity RATIO Vjet/Vinf written to the 'Jgain' column (as TASOPT does);
            jvl_patches_claude converts it to the JVL gain (Vjet/Vinf - 1)
    scale_nacelles : scale nacelle length/diameter with D / Dprop_base (only matters for
            full-configuration runs; wing-only runs drop the nacelle bodies)
    """
    df = df_geom_base.copy(deep=True)
    b = float(df["wing_spans_projected"][0])
    semi = b / 2.0
    N = int(N)
    if N % 2 != 0 or N < 2:
        raise ValueError("N must be even and >= 2")
    n_half = N // 2
    dy = dY * b / (N * D) - 1.0
    pitch = dY * b / N
    y_fus = float(df["body_widths_maximum"][0]) / 2.0
    y0 = min(y_fus, semi * (1.0 - dY) - tip_margin)
    info = {
        "N": N, "D": D, "dY": dY, "dy": dy, "dy_m": dy * D, "pitch_m": pitch,
        "y0_m": y0, "y0_inside_fuselage_m": max(0.0, y_fus - y0), "b": b, "feasible": True, "reason": "",
    }
    if dy < 0.0:
        info.update(feasible=False, reason="negative propeller spacing (dy < 0)")
        return None, info
    if y0 < 0.0:
        info.update(feasible=False, reason="band does not fit on the semi-span")
        return None, info

    secs = planform_sections(df)
    centres = [y0 + pitch / 2.0 + k * pitch for k in range(n_half)]
    new_secs = []
    for k, yc in enumerate(centres, start=1):
        new_secs.append((f"prop_{k}_in", (yc - D / 2.0) / semi))
        new_secs.append((f"prop_{k}_out", (yc + D / 2.0) / semi))
    for k in range(1, n_half):
        new_secs.append((f"mid_pitch_{k}_{k + 1}", (y0 + k * pitch) / semi))

    # Drop planform sections that coincide with a generated section (no duplicate y in JVL)
    gen_etas = [e for _, e in new_secs]
    keep_planform = [s for s in secs if all(abs(s["eta"] - e) > eta_merge_tol for e in gen_etas)]
    rows = [{"label": s["label"], "eta": s["eta"]} for s in keep_planform]
    rows += [{"label": lab, "eta": e} for lab, e in new_secs]
    rows.sort(key=lambda r: r["eta"])
    for r in rows:
        if r["eta"] > 1.0 + 1e-9:
            info.update(feasible=False, reason=f"section {r['label']} beyond tip (eta = {r['eta']:.4f})")
            return None, info
        sin = _inboard_section(secs, r["eta"])
        r["rcp"] = _interp(secs, r["eta"], "rcp")
        r["tc"] = _interp(secs, r["eta"], "tc")
        r["dih"] = sin["dih"] if r["label"] != "tip" else secs[-1]["dih"]
        r["sweep"] = sin["sweep"] if r["label"] != "tip" else secs[-1]["sweep"]

    # Replace sect_* columns
    df = df[[c for c in df.columns if not re.match(r"sect_\d+_", c)]].copy()
    n_rows = len(df)
    new_cols = {}
    for i, r in enumerate(rows, start=1):
        new_cols[f"sect_{i}_label"] = [r["label"]] * n_rows
        new_cols[f"sect_{i}_percent_span_location"] = [r["eta"]] * n_rows
        new_cols[f"sect_{i}_root_chord_percent"] = [r["rcp"]] * n_rows
        new_cols[f"sect_{i}_thickness_to_chord"] = [r["tc"]] * n_rows
        new_cols[f"sect_{i}_dihedral_outboard"] = [r["dih"]] * n_rows
        new_cols[f"sect_{i}_sweeps_quarter_chord"] = [r["sweep"]] * n_rows
    df = pd.concat([df, pd.DataFrame(new_cols, index=df.index)], axis=1)

    # Propeller disk centres (one wing half, as in the TASOPT .csv)
    xs, zs = [], []
    for yc in centres:
        x_le, z_le, c = _le_xz(df_geom_base, secs, yc / semi)
        xs.append(x_le - ax_pos * c)
        zs.append(z_le)
    if n_rows != 3:
        raise ValueError("Expected a 3-row TASOPT geometry .csv (nacelle_origin rows x, y, z)")
    df["nacelle_origin"] = [str([float(v) for v in xs]), str([float(v) for v in centres]), str([float(v) for v in zs])]

    D_base = float(df_geom_base["Dprop"][0])
    if scale_nacelles and D_base > 0:
        f = D / D_base
        for col in ("nacelle_length", "nacelle_inlet_diameter", "nacelle_diameter"):
            df[col] = float(df_geom_base[col][0]) * f
        df["nacelle_areas_wetted"] = float(df_geom_base["nacelle_areas_wetted"][0]) * f ** 2

    df["Dprop"] = D
    df["Jgain"] = Jgain
    df["dY"] = dY
    df["dy"] = dy
    info["sections"] = [(r["label"], round(r["eta"] * semi, 4)) for r in rows]
    info["prop_x"] = xs
    info["prop_y"] = centres
    return df, info


def base_mass_row(df_geom):
    """Mass/inertia row in the format wrap_jvl expects (study columns are dummies)."""
    keys = ["x_cg", "y_cg", "z_cg", "mass", "max_takeoff", "takeoff", "max_zero_fuel",
            "Ixx", "Iyy", "Izz", "Ixy", "Iyz", "Izx"]
    row = {k: float(df_geom[k][0]) for k in keys}
    row.update(sigma_fcs=-1.0, span_loc=-1.0, fcs_loc=-1.0, wing_frac=1.0, nacelle_frac=1.0)
    return pd.Series(row)


def check_layout(df_geom, D):
    """Return (n_props_per_side, blown_span_per_side, min_gap) from a laid-out df (for tests)."""
    y = ast.literal_eval(df_geom["nacelle_origin"].iloc[1])
    gaps = np.diff(y) - D if len(y) > 1 else np.array([np.nan])
    return len(y), len(y) * D, float(np.min(gaps))
