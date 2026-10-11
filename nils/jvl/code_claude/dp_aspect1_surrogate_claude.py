r"""
dp_aspect1_surrogate_claude.py

CLAUDE (10.10.2026): road-map point 3 - Aspect 1 results with the blown-wing surrogate
(analyse_blown_wing, imported UNCHANGED) as a function of the only free design variable N,
for the layouts of dp_layout_claude.py (TASOPT placement; 'const_area' and 'max_diameter'
disk scaling; x_p/c = 0.25 and y_p/c = -0.05 after Hawkswell, optionally with the motor-length
integration constraint).

Boundary conditions per mission point (fixed, from blown_wing_reference_claude.operating_point,
imported unchanged): altitude, V, T/W (total thrust, carried entirely by the N DP propellers),
CL_airframe (-> alpha_w by eq. 29). Airframe: TASOPT ATR 72-600 FC retrofit geometry .csv.

  A1  dCL vs N for takeoff (CL 2.2, alpha_w ~ 23 deg: outside the surrogate's validity, shown for
      reference), takeoff_linear (CL 0.8) and cruise (CL 0.6). Shaded band: surrogate accuracy
      from V1 (dCL_true = dCL / ratio, ratio in [min, max] of dp_verification_stats.csv,
      default 0.86-1.14). Dotted: beta applied to the far-wake excess (Patterson eq. 54).
  A2  cruise L/D vs N: (i) DP propellers producing the cruise thrust (blown), (ii) DP propellers
      folded (no blowing; cruise thrust from elsewhere), both with the nacelle-count penalty
          dCD_nac(N) = CD_nac,ref [(N/2)^(1/3) - 1]
      (total nacelle wetted area ~ N (P/N)^(2/3); CD_nac,ref = 0.00098, the TASOPT 2-nacelle value of
      the sigma = 3528 W/kg ATR retrofit). Baseline CD0 = 0.0226 (TASOPT ATR start-of-cruise value),
      span efficiency e from the flydrogen CONSTRAINT used inside analyse_blown_wing. ASSUMPTIONS,
      override on the command line.

Run in py311b:
  (py311b) C:\Users\nmb48\Documents\GitHub\Flydrogen> python C:\Users\nmb48\Documents\GitHub\SUAVE\nils\jvl\code_claude\dp_aspect1_surrogate_claude.py
Outputs: SUAVE\nils\jvl\data_claude\<ddmmyy>_dp_aspect1_claude\
  fig_A1_dCL_vs_N, fig_A2_LD_vs_N (.png/.svg), dp_aspect1_results_claude.csv, dp_layout_table_claude.csv
"""

import argparse
import glob
import math
import os
import sys
from datetime import date

import numpy as np
import pandas as pd
from ambiance import Atmosphere

from flydrogen.utilities.matplotlib_custom_settings import *  # noqa: F401,F403
from flydrogen.config import CONFIG
from flydrogen.methods import PerfectState
from flydrogen.inputs import Inputs
from flydrogen.analyses.aero_analysis import analyse_blown_wing, estimate_slipstream_velocity
from flydrogen.analyses.tests_claude.blown_wing_reference_claude import operating_point, _beta
from flydrogen.analyses.tests_claude.devries_fig5_surrogate_check_claude import alpha_w_eq29

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dp_layout_claude import wing_data, layout, DEFAULT_GEOMETRY, DATA  # noqa: E402

G = 9.81
OPS = ["takeoff", "takeoff_linear", "cruise"]
STRATS = ["const_area", "max_diameter"]


def evaluate(w, lay, op, eta_p, cd0, cd_nac_ref):
    """Surrogate quantities for one layout at one operating point (dict)."""
    h, T0, p0, V, ToverW, CL = operating_point(op)
    air = PerfectState.from_fluid_name("Air").set_p_T(p0, T0).set_V(V)
    W = w["mass"] * G
    INPUTS = Inputs()
    INPUTS.AC.WoverS, INPUTS.AC.ToverW, INPUTS.AC.CL_airframe = W / w["S"], ToverW, CL
    AR = w["b"] ** 2 / w["S"]
    N, D, xpc = lay["N"], lay["D"], lay["x_p_over_c"]
    dCL, dCD0, dCDi, V_w = analyse_blown_wing(air, INPUTS, AR=AR, lambda_cover2=0.0, ax_pos=xpc, dY=lay["dY"],
                                              dy=lay["dy"], N_prop=N, alpha_p=0.0, D_p=D, S=w["S"], b=w["b"])
    _, a_w, a_p, Rpoverc = estimate_slipstream_velocity(air, xpc, INPUTS, lay["dY"], lay["dy"], N, AR, D_p=D, S=w["S"])
    beta = _beta(a_p, xpc, Rpoverc)
    al = alpha_w_eq29(CL, AR, air.Ma)
    frac = N * D / w["b"]
    dCL_far = 2 * np.pi * np.sin(al) * ((1 + 2 * a_p * beta) ** 2 - 1) * frac  # Patterson eq. (54) basis
    e = getattr(INPUTS.OPS.CONSTRAINTS, CONFIG.constr_type).e
    cdi0 = CL ** 2 / (math.pi * AR * e)
    dcd_nac = cd_nac_ref * ((N / 2.0) ** (1.0 / 3.0) - 1.0)
    CD_blown = cd0 + dcd_nac + dCD0 + cdi0 + dCDi
    CD_folded = cd0 + dcd_nac + cdi0
    return {
        "op": op, "V": V, "Mach": air.Ma, "T_total": ToverW * W, "P_shaft": ToverW * W * V / eta_p,
        "CL_airframe": CL, "alpha_w_deg": math.degrees(al), "a_p": a_p, "a_w": a_w, "beta": beta, "Rp_over_c": Rpoverc,
        "Vjet_over_V": 1 + 2 * a_p, "dCL": dCL, "dCL_far_wake": dCL_far, "dCD0": dCD0, "dCDi": dCDi,
        "dCD_nacelles": dcd_nac, "CDi0": cdi0, "e": e,
        "LD_blown": (CL + dCL) / CD_blown, "LD_folded": CL / CD_folded, "LD_ref_N2_noblow": CL / (cd0 + cdi0),
    }


def _save(fig, out, name):
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(out, f"{name}.{ext}"), dpi=200)
    plt.close(fig)  # noqa: F405


def accuracy_band(path):
    """(min, max) of dCL_surrogate / dCL_RANS from the newest dp_verification_stats.csv."""
    files = sorted(glob.glob(path))
    if not files:
        return 0.86, 1.14
    st = pd.read_csv(files[-1])
    r = st[st["model"] == "surrogate"].iloc[0]
    return float(r["min_ratio"]), float(r["max_ratio"])


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--geometry", default=DEFAULT_GEOMETRY)
    p.add_argument("--N", type=int, nargs="*", default=[2, 4, 6, 8, 10, 12, 14, 16])
    p.add_argument("--eta-p", type=float, default=0.8, help="propeller efficiency for the shaft power (integration constraint)")
    p.add_argument("--cd0", type=float, default=0.0226)
    p.add_argument("--cd-nac-ref", type=float, default=0.00098)
    p.add_argument("--l-ref", type=float, default=1.0, help="motor + gearbox length ahead of the LE at P_ref (m), ASSUMPTION")
    p.add_argument("--P-ref", type=float, default=1.0e6)
    p.add_argument("--out", default=os.path.join(DATA, f"{date.today().strftime('%d%m%y')}_dp_aspect1_claude"))
    args = p.parse_args(argv)
    os.makedirs(args.out, exist_ok=True)
    w = wing_data(pd.read_csv(args.geometry))
    lo, hi = accuracy_band(os.path.join(DATA, "*_dp_verification_claude", "dp_verification_stats.csv"))

    rows, lays = [], []
    for op in OPS:
        h, T0, p0, V, ToverW, CL = operating_point(op)
        rho = p0 / (287.05 * T0)
        T = ToverW * w["mass"] * G
        P = T * V / args.eta_p
        for st in STRATS:
            for integ in ("none", "motor"):
                for N in args.N:
                    lay = layout(w, N, st, T_total=T, rho=rho, V=V, P_shaft=P, integration=integ,
                                 l_ref=args.l_ref, P_ref=args.P_ref)
                    res = evaluate(w, lay, op, args.eta_p, args.cd0, args.cd_nac_ref)
                    rows.append({**lay, **res})
                    if op == OPS[0] and integ == "none":
                        lays.append(lay)
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(args.out, "dp_aspect1_results_claude.csv"), index=False)
    pd.DataFrame(lays).to_csv(os.path.join(args.out, "dp_layout_table_claude.csv"), index=False)

    cs = {st: colors[k] for k, st in enumerate(STRATS)}  # noqa: F405
    titles = {"takeoff": r"take-off, $C_L$ = 2.2 ($\alpha_w \approx 23^\circ$, outside validity)",
              "takeoff_linear": r"take-off thrust, $C_L$ = 0.8", "cruise": r"cruise, $C_L$ = 0.6"}
    # A1
    fig, axs = plt.subplots(1, 3, figsize=(16, 5.2))  # noqa: F405
    for ax, op in zip(axs, OPS):
        for st in STRATS:
            for integ, ls in (("none", "-"), ("motor", "--")):
                d = df[(df["op"] == op) & (df["strategy"] == st) & (df["integration"] == integ)].sort_values("N")
                ax.plot(d["N"], d["dCL"], ls, color=cs[st], lw=1.3,
                        label=f"{st}, " + (r"$x_p/c$ = 0.25" if integ == "none" else "motor-length constraint"))
                f = d["feasible"].astype(bool)
                ax.scatter(d["N"][f], d["dCL"][f], s=18, color=cs[st])
                ax.scatter(d["N"][~f], d["dCL"][~f], s=24, facecolors="none", edgecolors=cs[st])
                if integ == "none":
                    ax.fill_between(d["N"], d["dCL"] / hi, d["dCL"] / lo, color=cs[st], alpha=0.12, lw=0)
                    ax.plot(d["N"], d["dCL_far_wake"], ":", color=cs[st], lw=1.0)
        ax.set_title(titles[op], fontsize=10)
        ax.set_xlabel("Number of DP propulsors $N$")
        ax.set_ylabel(r"$\Delta C_L$ (surrogate)")
        ax.set_ylim(bottom=0)
    axs[0].plot([], [], "k:", label=r"$\beta$ on far-wake excess")
    axs[0].scatter([], [], facecolors="none", edgecolors="k", label="infeasible layout (disks overlap)")
    axs[0].fill_between([], [], [], color="grey", alpha=0.2, label=f"surrogate accuracy (V1: {lo:.2f}-{hi:.2f})")
    axs[0].legend(frameon=False, fontsize=7)
    _save(fig, args.out, "fig_A1_dCL_vs_N")

    # A2
    fig, ax = plt.subplots(figsize=(7.5, 5.5))  # noqa: F405
    d0 = df[(df["op"] == "cruise")]
    for st in STRATS:
        d = d0[(d0["strategy"] == st) & (d0["integration"] == "none")].sort_values("N")
        ax.plot(d["N"], d["LD_blown"], "-o", ms=3, color=cs[st], label=f"{st}, DP propellers blowing (cruise thrust)")
        ax.plot(d["N"], d["LD_folded"], "--", color=cs[st], label=f"{st}, DP propellers folded")
    ax.axhline(float(d0["LD_ref_N2_noblow"].iloc[0]), color="grey", lw=0.8, ls=":", label="N = 2 nacelles, no blowing")
    ax.set_xlabel("Number of DP propulsors $N$")
    ax.set_ylabel(r"Cruise $L/D$ at fixed $\alpha_w$ (wing not resized)")
    ax.legend(frameon=False, fontsize=8)
    _save(fig, args.out, "fig_A2_LD_vs_N")

    cols = ["op", "strategy", "integration", "N", "D", "s_over_D", "x_p_over_c", "feasible", "Vjet_over_V",
            "beta", "dCL", "dCL_far_wake", "LD_blown", "LD_folded"]
    with pd.option_context("display.width", 220, "display.max_rows", 300):
        print(df[df["integration"] == "none"][cols].round(4).to_string(index=False))
    print("Written to", args.out)


if __name__ == "__main__":
    main()
