r"""
jvl_patterson_2d_claude.py

CLAUDE (05.10.2026): JVL against the 2D Euler (OVERFLOW) actuator-disk results that underlie
the beta surrogate of Patterson (2016), PhD thesis, Georgia Tech, Sec. 3.2.1:
  - Table 1 (p. 95): six validation cases (R/c, u/c, Vj/V, alpha, cl, beta), symmetric
    (modified NACA 0012) airfoil, M = 0.2;
  - the fitted surrogate (eqs. 56-59, R^2 = 0.99 to the OVERFLOW data) on its fit range
    R/c = 0.125 ... 3, u/c = 0.25 ... 1.5, Vj/V = 1.25 ... 2.0.

Patterson's definitions: 2D disk of height 2R placed u upstream of the LE, slipstream far-wake
velocity Vj; lift multiplier K_L = cl_blown / cl_isolated and (eq. 54)
    K_L = (beta (Vj - V) + V)^2 / V^2   ->   beta = (sqrt(K_L) - 1) / (Vj/V - 1).
Note: beta multiplies the FAR-WAKE excess velocity Vj - V = 2 a_p V, whereas eq. (28) of de Vries
et al. (and analyse_blown_wing) multiply it with a_w, the excess velocity at the wing LE.

JVL model (uniformly blown AR-100 wing as in jvl_calibration_2d_claude.py, whose
geometry_text / run_jvl are imported unchanged): hdisk = 2R/c (2D disk height), fh = 1,
JVL gain 1, DVjet = Vj/V - 1, dxdisk = -u/c, flat camber line, run-case Mach = 0.2.
beta_JVL = (sqrt(CLff_on / CLff_off) - 1) / (Vj/V - 1) (Trefftz-plane lift; AR effects cancel
to first order in the ratio).

Run (suave or py311b):
  python C:\Users\nmb48\Documents\GitHub\SUAVE\nils\jvl\code_claude\jvl_patterson_2d_claude.py
Output: SUAVE\nils\jvl\data_claude\<ddmmyy>_patterson_2d\patterson_2d_results.csv (+ case folders)
"""

import argparse
import math
import os
import sys
from datetime import date

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import jvl_calibration_2d_claude as cal  # noqa: E402

# Patterson (2016), eq. (58)
K = np.array([
    [0.378269, 0.748135, -0.179986, -0.056464, -0.146746, -0.015255],
    [3.071020, -1.769885, 0.436595, 0.148643, -0.989332, 0.197940],
    [-2.827730, 2.054064, -0.467410, -0.277325, 0.698981, -0.008226],
    [0.997936, -0.916118, 0.199829, 0.157810, -0.143368, -0.057385],
    [-0.127645, 0.135543, -0.028919, -0.026546, 0.010470, 0.012221],
])
# Patterson (2016), Table 1: R/c, u/c, Vj/V, alpha (deg), cl (OVERFLOW), beta (OVERFLOW)
TABLE1 = [
    (0.35, 0.35, 2.0, 3.0, 0.9542, 0.6114),
    (0.40, 0.80, 2.0, 3.0, 1.0822, 0.7161),
    (0.60, 0.60, 2.0, 3.0, 1.2032, 0.8095),
    (0.80, 0.80, 2.0, 3.0, 1.3474, 0.9149),
    (0.70, 0.75, 2.0, 5.0, 2.1279, 0.8646),
    (0.50, 0.75, 2.0, 1.0, 0.3850, 0.7721),
]


def beta_surrogate(Rc, uc, r):
    """Eq. (59) with X = [1, u/c, (u/c)^2, (u/c) Vj/V, Vj/V, (Vj/V)^2]^T."""
    X = np.array([1.0, uc, uc ** 2, uc * r, r, r ** 2])
    return float(sum((K[i] @ X) * Rc ** i for i in range(5)))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--jvl", default=cal.DEFAULT_JVL)
    p.add_argument("--out-root", default=cal.DEFAULT_OUT_ROOT)
    p.add_argument("--mach", type=float, default=0.2)
    p.add_argument("--alpha", type=float, default=3.0, help="alpha of the surrogate-grid cases (deg)")
    p.add_argument("--Rc", type=float, nargs="*", default=[0.125, 0.25, 0.35, 0.5, 0.75, 1.0, 1.5])
    p.add_argument("--uc", type=float, nargs="*", default=[0.25, 0.5, 1.0])
    p.add_argument("--r", type=float, nargs="*", default=[1.25, 1.5, 2.0], help="Vj / V")
    p.add_argument("--nspan", type=int, default=40)
    p.add_argument("--nchord", type=int, default=8)
    args = p.parse_args(argv)

    out_dir = os.path.join(args.out_root, f"{date.today().strftime('%d%m%y')}_patterson_2d")
    os.makedirs(out_dir, exist_ok=True)
    cal.RUN_CASE = cal.RUN_CASE.replace(" Mach      =   0.0\n", f" Mach      =   {args.mach}\n")  # run-case Mach (in memory only)

    # cases: (R/c, u/c, alpha) -> list of Vj/V; Table 1 first
    cases = {}
    for Rc, uc, r, a, cl, b in TABLE1:
        cases.setdefault((Rc, uc, a), set()).add(r)
    for Rc in args.Rc:
        for uc in args.uc:
            cases.setdefault((Rc, uc, args.alpha), set()).update(args.r)

    rows = []
    for (Rc, uc, a), rs in sorted(cases.items()):
        rs = sorted(rs)
        folder = os.path.join(out_dir, f"Rc{Rc:g}_uc{uc:g}_a{a:g}")
        geom = cal.geometry_text("full", 2.0 * Rc, 0.0, args.nspan, args.nchord, -uc, 0.0)
        res = cal.run_jvl(args.jvl, folder, geom, [0.0] + [r - 1.0 for r in rs], a)
        off = res[0]
        for r, on in zip(rs, res[1:]):
            KL = on["CLff"] / off["CLff"]
            row = {"Rc": Rc, "uc": uc, "alpha_deg": a, "Vj_over_V": r, "hdisk_over_c": 2.0 * Rc,
                   "CLff_off": off["CLff"], "CLff_on": on["CLff"], "KL_jvl": KL,
                   "beta_jvl": (math.sqrt(KL) - 1.0) / (r - 1.0),
                   "beta_surrogate": beta_surrogate(Rc, uc, r),
                   "CT_jvl": on.get("CTtot", np.nan)}
            t1 = [t for t in TABLE1 if t[0] == Rc and t[1] == uc and t[2] == r and t[3] == a]
            row["beta_overflow_table1"] = t1[0][5] if t1 else np.nan
            row["cl_overflow_table1"] = t1[0][4] if t1 else np.nan
            rows.append(row)
        print(f"R/c = {Rc:g}, u/c = {uc:g}, alpha = {a:g}: done")
    df = pd.DataFrame(rows)
    df["ratio_jvl_surrogate"] = df["beta_jvl"] / df["beta_surrogate"]
    out_csv = os.path.join(out_dir, "patterson_2d_results.csv")
    df.to_csv(out_csv, index=False)
    with pd.option_context("display.width", 200, "display.max_rows", 200):
        print(df[["Rc", "uc", "alpha_deg", "Vj_over_V", "beta_jvl", "beta_surrogate", "beta_overflow_table1"]].round(3).to_string(index=False))
    print("Written to", out_csv)


if __name__ == "__main__":
    main()
