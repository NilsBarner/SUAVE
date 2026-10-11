r"""
dp_verification_figures_claude.py

CLAUDE (10.10.2026): road-map point 1 - "designer's view" verification figures V1-V3 comparing
Hawkswell's (2024) RANS, JVL and the blown-wing surrogate for the spanwise-periodic ATR-42-section
cases (NACA 43018, D/c = 0.5, x_p/c = -0.33, y_p/c = -0.05, uniform disks, no flap).

  V1  parity: dCL_model vs dCL_CFD for all CFD points (alpha 1-12 deg, Tc = 1/2, s/D = 0-2);
      mean ratio and RMS relative error per model (alpha >= 3 deg) in the legend
  V2  constant TOTAL thrust: dCL vs wing thrust coefficient C_T,w = T_tot/(q S)
      = Tc (pi D/4c) / (1 + s/D) at alpha = 6 deg. The Tc = 1 and Tc = 2 families overlap for
      C_T,w = 0.26 ... 0.37: at equal C_T,w the gap between the families is the effect of
      distributing the same thrust over fewer, more highly loaded disks (an N-sweep at fixed
      thrust). The x-pair of Hawkswell Fig. 6.15 (Tc = 1, s/D = 0.05 vs Tc = 2, s/D = 1.1) is marked.
  V3  ratio dCL_model / dCL_CFD vs alpha (operating condition) and vs s/D (design choice)

Increments dCL = CL_blown - CL_unblown incl. T sin(alpha) (thesis eq. 1): CFD w.r.t. the CFD clean
wing (Fig. 6.4a), JVL w.r.t. its own unblown wing (AR -> infinity), surrogate as in
jvl_hawkswell_postprocess_claude.surrogate_dCL (imported unchanged; beta applied to a_w, as in
analyse_blown_wing). Variant "surrogate (far wake)": Patterson's own definition (eq. 54), i.e.
beta applied to the far-wake excess 2 a_p instead of a_w.

Run in py311b (imports flydrogen):
  (py311b) C:\Users\nmb48\Documents\GitHub\Flydrogen> python C:\Users\nmb48\Documents\GitHub\SUAVE\nils\jvl\code_claude\dp_verification_figures_claude.py ^
        [--jvl <...>\hawkswell_jvl_results.csv]
Default JVL results: SUAVE\nils\jvl\data_claude\051026_linux_preview_claude\hawkswell_jvl_results.csv
Outputs: SUAVE\nils\jvl\data_claude\<ddmmyy>_dp_verification_claude\
  fig_V1_parity, fig_V2_constant_thrust, fig_V3_ratios (.png/.svg), dp_verification_points.csv,
  dp_verification_stats.csv
"""

import argparse
import os
from datetime import date

import numpy as np
import pandas as pd

from flydrogen.utilities.matplotlib_custom_settings import *  # noqa: F401,F403
from flydrogen.analyses.tests_claude.jvl_hawkswell_postprocess_claude import cfd_clean, surrogate_dCL, D_OVER_C

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(HERE), "data_claude")
MODELS = ["JVL", "surrogate", "surrogate (far wake)"]


def ct_wing(Tc, sD):
    """Wing thrust coefficient T_tot / (q S) of a periodic row (per pitch D + s, chord c)."""
    return Tc * np.pi * D_OVER_C / 4.0 / (1.0 + sD)


def surrogate_far_wake(Tc, alpha_deg, sD):
    """Surrogate increment with beta applied to the far-wake excess 2 a_p (Patterson eq. 54)."""
    d_aw, beta, a_w = surrogate_dCL(Tc, alpha_deg, sD)
    a_p = (np.sqrt(1.0 + Tc) - 1.0) / 2.0
    s = np.sin(np.radians(alpha_deg))
    frac = 1.0 / (1.0 + sD)
    return d_aw + 2 * np.pi * s * (((1 + 2 * a_p * beta) ** 2 - 1) - ((1 + a_w * beta) ** 2 - 1)) * frac


def cfd_points(ref):
    """All CFD lift points with Tc > 0 as (Tc, s/D, alpha, dCL_cfd); duplicates averaged."""
    c = ref[(ref["quantity"] == "CL") & (ref["Tc"] > 0)].copy()
    c["s_over_D"] = c["s_over_D"].fillna(0.0)
    c = c.groupby(["Tc", "s_over_D", "alpha_deg"], as_index=False)["value"].mean()
    c["dCL_cfd"] = c["value"] - cfd_clean(ref, c["alpha_deg"].values)
    return c[c["alpha_deg"] > 0].rename(columns={"value": "CL_cfd"})


def jvl_dCL(jvl, Tc, sD, alpha):
    """JVL increment at (Tc, s/D), linearly interpolated in alpha on JVL's alpha grid."""
    on = jvl[(jvl["Tc"] == Tc) & np.isclose(jvl["s_over_D"], sD)].sort_values("alpha_deg")
    off = jvl[(jvl["Tc"] == 0.0) & np.isclose(jvl["s_over_D"], sD)].sort_values("alpha_deg")
    if on.empty or off.empty:
        return np.nan
    d = on["CL_hawkswell_inf"].values - off["CL_hawkswell_inf"].values
    return float(np.interp(alpha, on["alpha_deg"].values, d))


def _save(fig, out, name):
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(out, f"{name}.{ext}"), dpi=200)
    plt.close(fig)  # noqa: F405


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--jvl", default=os.path.join(DATA, "051026_linux_preview_claude", "hawkswell_jvl_results.csv"))
    p.add_argument("--reference", default=os.path.join(DATA, "hawkswell_2024_digitised_claude.csv"))
    p.add_argument("--out", default=os.path.join(DATA, f"{date.today().strftime('%d%m%y')}_dp_verification_claude"))
    args = p.parse_args(argv)
    os.makedirs(args.out, exist_ok=True)
    jvl = pd.read_csv(args.jvl)
    ref = pd.read_csv(args.reference, comment="#")

    pts = cfd_points(ref)
    pts["C_T_w"] = ct_wing(pts["Tc"], pts["s_over_D"])
    pts["JVL"] = [jvl_dCL(jvl, t, s, a) for t, s, a in zip(pts["Tc"], pts["s_over_D"], pts["alpha_deg"])]
    pts["surrogate"] = [surrogate_dCL(t, a, s)[0] for t, s, a in zip(pts["Tc"], pts["s_over_D"], pts["alpha_deg"])]
    pts["surrogate (far wake)"] = [surrogate_far_wake(t, a, s) for t, s, a in zip(pts["Tc"], pts["s_over_D"], pts["alpha_deg"])]
    for m in MODELS:
        pts["ratio " + m] = pts[m] / pts["dCL_cfd"]
    pts.to_csv(os.path.join(args.out, "dp_verification_points.csv"), index=False)

    rows = []
    for m in MODELS:
        r = pts.loc[pts["alpha_deg"] >= 3.0, "ratio " + m].dropna()  # alpha < 3 deg: dCL_CFD < 0.2 (Fig. 6.4a), ratios ill-conditioned
        rows.append({"model": m, "n": len(r), "mean_ratio": r.mean(), "min_ratio": r.min(), "max_ratio": r.max(),
                     "rms_rel_err": np.sqrt(((r - 1) ** 2).mean())})
    stats = pd.DataFrame(rows)
    stats.to_csv(os.path.join(args.out, "dp_verification_stats.csv"), index=False)
    cm = {m: colors[k] for k, m in enumerate(MODELS)}  # noqa: F405
    mk = {3.0: "o", 6.0: "s", 9.0: "^", 12.0: "D"}

    # V1 parity
    fig, ax = plt.subplots(figsize=(6.5, 6.0))  # noqa: F405
    lim = [0, 1.05 * float(pts[["dCL_cfd"] + MODELS].max().max())]
    ax.plot(lim, lim, color="k", lw=0.8)
    ax.fill_between(lim, [0.9 * v for v in lim], [1.1 * v for v in lim], color="grey", alpha=0.15, lw=0)
    for m in MODELS:
        st = stats[stats["model"] == m].iloc[0]
        for a, g in pts.groupby("alpha_deg"):
            ax.scatter(g["dCL_cfd"], g[m], s=16, marker=mk.get(a, "."), color=cm[m],
                       facecolors="none" if "far" in m else cm[m], alpha=0.85)
        ax.scatter([], [], color=cm[m], label=f"{m}: mean {st.mean_ratio:.2f}, RMS err {100 * st.rms_rel_err:.0f} %")
    for a in mk:
        ax.scatter([], [], marker=mk[a], color="grey", label=rf"$\alpha$ = {a:g}$^\circ$")
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel(r"$\Delta C_L$ RANS (Hawkswell 2024)")
    ax.set_ylabel(r"$\Delta C_L$ model")
    ax.legend(frameon=False, fontsize=7, loc="upper left")
    _save(fig, args.out, "fig_V1_parity")

    # V2 constant total thrust, alpha = 6 deg
    fig, ax = plt.subplots(figsize=(7.5, 5.5))  # noqa: F405
    a6 = pts[pts["alpha_deg"] == 6.0]
    sD_fine = np.linspace(0.0, 2.0, 81)
    for tc, ls, mkr in ((1.0, "-", "o"), (2.0, "--", "s")):
        g = a6[a6["Tc"] == tc].sort_values("C_T_w")
        ax.scatter(g["C_T_w"], g["dCL_cfd"], s=30, marker=mkr, facecolors="none", edgecolors="k",
                   label=rf"RANS, $T_c$ = {tc:g}")
        x = ct_wing(tc, sD_fine)
        ax.plot(x, [surrogate_dCL(tc, 6.0, s)[0] for s in sD_fine], ls, color=cm["surrogate"], lw=1.2)
        ax.plot(x, [surrogate_far_wake(tc, 6.0, s) for s in sD_fine], ls, color=cm["surrogate (far wake)"], lw=1.0)
        j = jvl[(jvl["Tc"] == tc) & (jvl["alpha_deg"] == 6.0)].sort_values("s_over_D")
        jx = ct_wing(tc, j["s_over_D"].values)
        ax.plot(jx, [jvl_dCL(jvl, tc, s, 6.0) for s in j["s_over_D"]], ls, color=cm["JVL"], lw=1.2)
    for m in MODELS:
        ax.plot([], [], color=cm[m], label=m)
    ax.plot([], [], "k-", label=r"$T_c$ = 1 family")
    ax.plot([], [], "k--", label=r"$T_c$ = 2 family")
    ax.axvspan(ct_wing(2.0, 2.0), ct_wing(1.0, 0.0), color="grey", alpha=0.08, lw=0)  # C_T,w range covered by both families
    x_pair = ct_wing(1.0, 0.05)
    ax.axvline(x_pair, color="grey", lw=0.6, ls=":")
    ax.annotate(r"Fig. 6.15 pair: $T_c$=1, $s/D$=0.05 vs $T_c$=2, $s/D$=1.1", (x_pair, 0.02), fontsize=7, rotation=90,
                va="bottom", ha="right")
    ax.set_xlabel(r"Wing thrust coefficient $C_{T,w} = T_\mathrm{tot}/(qS)$")
    ax.set_ylabel(r"$\Delta C_L$ at $\alpha = 6^\circ$ (incl. $T\sin\alpha$)")
    ax.set_title("Same total thrust (vertical lines), different distribution (families)", fontsize=10)
    ax.legend(frameon=False, fontsize=7, loc="upper left")
    _save(fig, args.out, "fig_V2_constant_thrust")

    # V3 ratios
    fig, axs = plt.subplots(1, 2, figsize=(13, 5.0))  # noqa: F405
    for k, m in enumerate(MODELS):
        g = pts[pts["Tc"] == 1.0]
        jit = (k - 1) * 0.25
        axs[0].scatter(g["alpha_deg"] + jit, g["ratio " + m], s=12, color=cm[m], alpha=0.6)
        mm = g.groupby("alpha_deg")["ratio " + m].mean()
        axs[0].plot(mm.index + jit, mm.values, "-", color=cm[m], label=m)
        for tc, ls in ((1.0, "-"), (2.0, "--")):
            h = a6[a6["Tc"] == tc].sort_values("s_over_D")
            axs[1].plot(h["s_over_D"], h["ratio " + m], ls, marker="o", ms=3, color=cm[m],
                        label=f"{m}, $T_c$ = {tc:g}")
    for ax in axs:
        ax.axhline(1.0, color="k", lw=0.8)
        ax.axhspan(0.9, 1.1, color="grey", alpha=0.12, lw=0)
        ax.set_ylabel(r"$\Delta C_{L,\mathrm{model}} / \Delta C_{L,\mathrm{RANS}}$")
    axs[0].set_xlabel(r"$\alpha$ (deg), $T_c$ = 1, all $s/D$ (line: mean)")
    axs[1].set_xlabel(r"Spacing $s/D$ ($\alpha = 6^\circ$)")
    axs[0].legend(frameon=False, fontsize=8)
    axs[1].legend(frameon=False, fontsize=7, ncol=2)
    _save(fig, args.out, "fig_V3_ratios")

    with pd.option_context("display.width", 200):
        print(stats.round(3).to_string(index=False))
    print("Written to", args.out)


if __name__ == "__main__":
    main()
