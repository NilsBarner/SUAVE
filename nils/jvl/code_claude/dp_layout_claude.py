r"""
dp_layout_claude.py

CLAUDE (10.10.2026): road-map point 2 - "designer's view" propulsor layout. Given the only free
design variable, the number of propulsors N, a placement strategy and a disk-area scaling
strategy, return all dependent layout quantities (D, y_i, pitch, s, s/D, dY, dy, x_p, y_p,
A_tot, Tc, Vjet/V) plus feasibility flags. Pure numpy/pandas (no SUAVE, no flydrogen).

Placement (TASOPT, src/IO/output_plots.jl): N = 2: y = eta_s b/2 (wing break, here the
'outboard' planform section); N > 2: N/2 propellers per side at
    y_i = linspace(b_o/2 + k_root D, f_outer b/2, N/2),   k_root = 2, f_outer = 0.75
(b_o = fuselage width). k_root and f_outer are options.

Disk-area scaling at a fixed mission point (V, rho) and fixed total thrust T:
  'const_area'   N pi D^2/4 = A_ref (= 2 pi D_ref^2/4, D_ref = TASOPT 'Dprop' of the 2-engine
                 aircraft). Equivalent to constant disk loading, constant Pi, constant
                 Tc = T/(q A) and constant (C_T, J) of the propeller (T/N = C_T rho (V/J)^2 D^2
                 => D^2 ~ 1/N), i.e. the user's strategies 1-3 coincide unless a tip-Mach/rpm
                 limit is active.
  'max_diameter' the largest D satisfying (a) no disk overlap, s/D >= s_min (default 0.05,
                 Hawkswell's minimum), (b) fuselage clearance y_1 - D/2 >= b_o/2 + gap_fus,
                 (c) ground clearance of every disk at a bank angle phi:
                     h_wing + y_p c_i + y_i tan(Gamma) - y_i tan(phi) - D/2 >= h_min.
                 Solved by bisection on D.

Chordwise / vertical position (Hawkswell 2024, Ch. 5/8 optimum, p. 126-127 and 185-187 of the PDF):
  x_p = 0.25 c_i upstream of the LE, y_p = -0.05 c_i (folding) or 0 (non-folding).
Optional integration constraint ('--integration motor'): the disk must clear the motor/gearbox
length ahead of the LE, x_p >= l_fwd = l_ref (P_shaft/N / P_ref)^(1/3), so x_p/c = max(0.25, l_fwd/c).
ASSUMPTION (to be replaced by your motor sizing): l_ref = 1.0 m at P_ref = 1 MW.

Ground geometry ASSUMPTIONS (ATR 72-600, to be checked): wing chord line h_wing = 3.4 m above
ground (fuselage underside ~0.6 m + fuselage height 2.77 m), bank angle phi = 5 deg,
h_min = 0.18 m (7 in, CS 25.925(a)), Gamma from the geometry .csv (0 for the TASOPT ATR).

Usage (any env with numpy/pandas):
  python dp_layout_claude.py --geometry <TASOPT geometry .csv> --N 2 4 6 8 10 12 16
prints the layout table and writes dp_layout_table_claude.csv to SUAVE\nils\jvl\data_claude\<ddmmyy>_dp_aspect1_claude\
"""

import argparse
import math
import os
import sys
from datetime import date

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from blown_wing_layout_claude import planform_sections, _interp  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(HERE), "data_claude")
DEFAULT_GEOMETRY = r"C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_geometry_regional_LH2_041026.csv"
GROUND = dict(h_wing=3.4, phi_deg=5.0, h_min=0.18, gap_fus=0.1)


def wing_data(df_geom):
    """Span, area, root chord, fuselage width, break station, dihedral from a TASOPT geometry csv."""
    g = lambda k: float(df_geom[k].iloc[0])  # noqa: E731
    secs = planform_sections(df_geom)
    return {
        "b": g("wing_spans_projected"), "S": g("wing_areas_reference"), "c_root": g("wing_chords_root"),
        "b_o": g("body_widths_maximum"), "D_ref": g("Dprop"), "mass": g("takeoff"),
        "eta_s": [s["eta"] for s in secs if s["label"] == "outboard"][0],
        "gamma": math.radians(g("wing_center_dihedral_outboard")), "secs": secs,
    }


def chord_at(w, y):
    """Local chord (root-chord fraction interpolated on the planform sections)."""
    return w["c_root"] * _interp(w["secs"], y / (w["b"] / 2.0), "rcp")


def stations(w, N, D, k_root=2.0, f_outer=0.75):
    """TASOPT propulsor y positions (one side)."""
    if N == 2:
        return np.array([w["eta_s"] * w["b"] / 2.0])
    return np.linspace(w["b_o"] / 2.0 + k_root * D, f_outer * w["b"] / 2.0, N // 2)


def _constraints(w, N, D, y, s_min, y_p_over_c, ground):
    c = np.array([chord_at(w, yi) for yi in y])
    gaps = np.diff(y) - D if len(y) > 1 else np.array([np.inf])
    clear = (ground["h_wing"] + y_p_over_c * c + y * math.tan(w["gamma"])
             - y * math.tan(math.radians(ground["phi_deg"])) - D / 2.0)
    return {
        "overlap_ok": bool(np.all(gaps >= s_min * D - 1e-9)),
        "fuselage_ok": bool(y[0] - D / 2.0 >= w["b_o"] / 2.0 + ground["gap_fus"] - 1e-9),
        "ground_ok": bool(np.all(clear >= ground["h_min"] - 1e-9)),
        "min_ground_clearance_m": float(clear.min()), "min_gap_over_D": float(gaps.min() / D),
    }


def layout(w, N, strategy, T_total=None, rho=None, V=None, P_shaft=None, s_min=0.05, x_p_over_c=0.25,
           y_p_over_c=-0.05, integration="none", l_ref=1.0, P_ref=1.0e6, k_root=2.0, f_outer=0.75, ground=GROUND):
    """All dependent layout quantities for N propulsors (dict)."""
    if strategy == "const_area":
        D = w["D_ref"] * math.sqrt(2.0 / N)
    elif strategy == "max_diameter":
        lo, hi = 0.05, 10.0
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            k = _constraints(w, N, mid, stations(w, N, mid, k_root, f_outer), s_min, y_p_over_c, ground)
            lo, hi = (mid, hi) if (k["overlap_ok"] and k["fuselage_ok"] and k["ground_ok"]) else (lo, mid)
        D = lo
    else:
        raise ValueError(strategy)
    y = stations(w, N, D, k_root, f_outer)
    k = _constraints(w, N, D, y, s_min, y_p_over_c, ground)
    c = np.array([chord_at(w, yi) for yi in y])
    pitch = float(np.mean(np.diff(y))) if len(y) > 1 else float("nan")
    s = pitch - D if len(y) > 1 else float("nan")
    xpc = x_p_over_c
    l_fwd = 0.0
    if integration == "motor" and P_shaft is not None:
        l_fwd = l_ref * (P_shaft / N / P_ref) ** (1.0 / 3.0)
        xpc = max(x_p_over_c, l_fwd / float(np.mean(c)))
    out = {
        "N": N, "strategy": strategy, "integration": integration, "D": D, "A_tot": N * math.pi * D ** 2 / 4.0,
        "y_inner": float(y[0]), "y_outer": float(y[-1]), "pitch": pitch, "s": s,
        "s_over_D": s / D if len(y) > 1 else float("nan"),
        # de Vries span fraction / spacing: dY b = N D (1 + dy) (eq. 21), blown band = N/2 pitches per side
        "dY": (N * (D + s) if len(y) > 1 else N * D) / w["b"], "dy": s / D if len(y) > 1 else 0.0,
        "blown_frac": N * D / w["b"], "c_mean_props": float(np.mean(c)), "x_p_over_c": xpc, "l_fwd": l_fwd,
        "x_p": xpc * float(np.mean(c)), "y_p_over_c": y_p_over_c, **k,
    }
    out["feasible"] = out["overlap_ok"] and out["fuselage_ok"] and out["ground_ok"]
    if T_total is not None and rho is not None and V is not None:
        q = 0.5 * rho * V ** 2
        out["Tc_disk"] = T_total / (q * out["A_tot"])  # Hawkswell's Tc = T / (q A)
        out["Vjet_over_V"] = math.sqrt(1.0 + out["Tc_disk"])
        out["a_p"] = (out["Vjet_over_V"] - 1.0) / 2.0
    return out


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--geometry", default=DEFAULT_GEOMETRY)
    p.add_argument("--N", type=int, nargs="*", default=[2, 4, 6, 8, 10, 12, 16])
    p.add_argument("--out", default=os.path.join(DATA, f"{date.today().strftime('%d%m%y')}_dp_aspect1_claude"))
    args = p.parse_args(argv)
    w = wing_data(pd.read_csv(args.geometry))
    rows = [layout(w, N, st) for st in ("const_area", "max_diameter") for N in args.N]
    df = pd.DataFrame(rows)
    os.makedirs(args.out, exist_ok=True)
    df.to_csv(os.path.join(args.out, "dp_layout_table_claude.csv"), index=False)
    with pd.option_context("display.width", 220, "display.max_columns", 30):
        print(df[["strategy", "N", "D", "A_tot", "s_over_D", "dY", "blown_frac", "c_mean_props",
                  "min_ground_clearance_m", "feasible"]].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
