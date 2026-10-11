r"""
jvl_hawkswell_periodic_claude.py

CLAUDE (05.10.2026): JVL verification against the 3D/2D actuator-disk RANS results of
Hawkswell (2024), PhD thesis, Cambridge, Ch. 6 (Figs. 6.2, 6.4a, 6.13, 6.15): a straight,
untapered NACA 43018 wing (ATR-42 root section) blown by a spanwise-periodic row of uniform
actuator disks (no swirl, no nacelles), no flap, M = 0.18.

  D/c = 0.5, x_p/c = -0.33, y_p/c = -0.05, Tc = T / (q pi D^2/4) = 1 and 2,
  s/D = 0.05 ... 2 (tip-to-tip gap s), alpha = 0 ... 12 deg, plus the continuous ("2D") jet.

Hawkswell's 3D CFD is spanwise periodic (no tips). JVL has no periodic boundary condition, so
each case is run on two finite wings made of N2 and N1 = N2/2 pitches per semispan
(identical strip layout, aspect ratios AR2 = 2 AR1) and extrapolated to AR -> infinity
assuming CL(AR) = CL_inf - k/AR (lifting-line leading term):

    CL_inf = CL2 + (CL2 - CL1) * AR1 / (AR2 - AR1)

JVL model of one propeller: a jet band of width D between two SECTIONs carrying JETCONTROL,
hdisk = pi D / 4 (Hawkswell's 2D mean jet height, same mass and momentum flux as the disk),
fh = 1, Vjet/V = sqrt(1 + Tc) (actuator disk), JVL gain 1 and DVjet = Vjet/V - 1 (fix of
issue B1). Each gap carries one SECTION without JETCONTROL at mid-gap, so that no jet is
active in the gap (jets act only between two JETCONTROL sections; issue B10, checked in the
2D calibration). The root (y = 0) and the tip are at mid-gap. dxdisk = x_p/c and
dndisk = y_p/c are written for completeness (JVL's thin-jet model has no jet position, so
Hawkswell's x_p and y_p effects cannot be reproduced).

Re-used code (imported unchanged):
  jvl_calibration_2d_claude.run_jvl / RUN_CASE   writes .run/.deck, calls jvl.exe, parses 'st'
  jvl_patches_claude.postprocess_jvl_text         per-SECTION Nspanwise + total strip cap
The run-case Mach of RUN_CASE (0.0) is replaced in memory by --mach (0.18).

Recorded per (s/D, N, alpha, Tc): CLff, CLtot, CLffc, CLffj, CDff, CTtot, DCQtot (Trefftz-plane
CLff is the reference lift metric, see Appendix C of the plan). Hawkswell's lift includes the
thrust component normal to the freestream (thesis eq. 1: CL = (L + T sin(alpha)) / (q S)) with
the disk axis along the chord; the post-processing therefore compares CLff + CTtot sin(alpha).

Run in the "suave" (or py311b) environment:

  python C:\Users\nmb48\Documents\GitHub\SUAVE\nils\jvl\code_claude\jvl_hawkswell_periodic_claude.py --smoke
  python C:\Users\nmb48\Documents\GitHub\SUAVE\nils\jvl\code_claude\jvl_hawkswell_periodic_claude.py

Outputs: SUAVE\nils\jvl\data_claude\<ddmmyy>_hawkswell_periodic[_<tag>]\
  hawkswell_jvl_raw.csv        one row per JVL run case (finite wings)
  hawkswell_jvl_results.csv    AR -> infinity extrapolation per (s/D, alpha, Tc)
  naca43018_claude.dat, case folders with cal.jvl / cal.run / st_*.txt / jvl_log.txt
Then plot with (py311b):
  python -m flydrogen.analyses.tests_claude.jvl_hawkswell_postprocess_claude --jvl <...>\hawkswell_jvl_results.csv
"""

import argparse
import math
import os
import sys
import time
from datetime import date

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import jvl_calibration_2d_claude as cal  # noqa: E402
from jvl_patches_claude import postprocess_jvl_text  # noqa: E402

DEFAULT_S_OVER_D = [0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.1, 1.5, 2.0]  # 0.0 = continuous jet ("2D")
DEFAULT_ALPHA = [0.0, 3.0, 6.0, 9.0, 12.0]
DEFAULT_TC = [1.0, 2.0]
AFILE_NAME = "naca43018_claude.dat"


def naca5_coordinates(L=4, P=3, TT=18, n=81):
    """NACA LPQTT (Q = 0, standard camber) coordinates in Selig order (TE-upper-LE-lower-TE).

    Camber line of the 5-digit family with p = P/20 (r, k1 tabulated for 230: r = 0.2025,
    k1 = 15.957 at cl_i = 0.3; k1 scales with the design lift coefficient cl_i = 0.15 L).
    NACA 43018: cl_i = 0.6, max camber at 15 % chord, t/c = 0.18.
    """
    if P != 3:
        raise ValueError("only P = 3 (p = 0.15) is tabulated here")
    r, k1 = 0.2025, 15.957 * (0.15 * L) / 0.3
    t = TT / 100.0
    beta = np.linspace(0.0, np.pi, n)
    x = 0.5 * (1.0 - np.cos(beta))
    yc = np.where(x < r, k1 / 6.0 * (x ** 3 - 3 * r * x ** 2 + r ** 2 * (3 - r) * x), k1 * r ** 3 / 6.0 * (1 - x))
    dyc = np.where(x < r, k1 / 6.0 * (3 * x ** 2 - 6 * r * x + r ** 2 * (3 - r)), -k1 * r ** 3 / 6.0)
    yt = 5 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x ** 2 + 0.2843 * x ** 3 - 0.1036 * x ** 4)  # closed TE
    th = np.arctan(dyc)
    xu, yu = x - yt * np.sin(th), yc + yt * np.cos(th)
    xl, yl = x + yt * np.sin(th), yc - yt * np.cos(th)
    X = np.concatenate([xu[::-1], xl[1:]])
    Y = np.concatenate([yu[::-1], yl[1:]])
    return X, Y


def write_afile(path):
    X, Y = naca5_coordinates()
    with open(path, "w") as f:
        f.write("NACA 43018 (CLAUDE, 05.10.2026)\n")
        for x, y in zip(X, Y):
            f.write(f"{x:.6f}  {y:.6f}\n")


def section_text(y, jet, hdisk, dxdisk, dndisk, afile):
    """One SECTION (chord 1, LE at x = 0), optionally with AFILE and JETCONTROL/JETPARAM.
    JETCONTROL / JETPARAM lines as in jvl_calibration_2d_claude.section_text."""
    txt = (
        "SECTION\n"
        "#Xle    Yle      Zle      Chord     Ainc\n"
        f"0.0  {y:.6f}  0.0  1.0  0.0\n"
    )
    if afile:
        txt += f"AFILE\n{afile}\n"
    if jet:
        txt += (
            "JETCONTROL\n"
            "#Jname   Jgain    SgnDup\n"
            "DVjet    1.0      1.0\n"
            "JETPARAM\n"
            "#hdisk   fh      djet0   djet1   djet3       dxdisk  dndisk\n"
            f"{hdisk:.6f}    1.0     0.0    0.0   0.0    {dxdisk:.6f}   {dndisk:.6f}\n"
        )
    return txt + "\n"


def section_layout(s_over_D, D, n_pitch):
    """[(y, has_jet)] for one semispan of n_pitch periodic pitches; root and tip at mid-gap."""
    if s_over_D <= 0.0:  # continuous jet over the whole semispan
        return [(0.0, True), (n_pitch * D, True)]
    s = s_over_D * D
    p = D + s
    secs = [(0.0, False)]
    for k in range(n_pitch):
        y0 = k * p + 0.5 * s
        secs += [(y0, True), (y0 + D, True), ((k + 1) * p, False)]
    return secs


def geometry_text(s_over_D, D, n_pitch, mach, xp, yp, afile, nchord):
    secs = section_layout(s_over_D, D, n_pitch)
    semispan = secs[-1][0]
    Sref = 2.0 * semispan * 1.0
    txt = (
        "hawkswell_periodic_claude\n"
        f"#Mach\n{mach}\n"
        "#IYsym IZsym Zsym\n0 0 0.0\n"
        "#Sref   Cref   Bref\n"
        f"{Sref:.6f} 1.0 {2.0 * semispan:.6f}\n"
        "#Xref   Yref   Zref\n0.25 0.0 0.0\n"
        "#CDp\n0.0\n"
        "SURFACE\nmain_wing\n"
        "#Nchordwise  Cspace   Nspanwise  Sspace\n"
        f"{nchord}  1.0\n"
        "YDUPLICATE\n0.0\n\n"
    )
    hdisk = math.pi * D / 4.0
    for y, jet in secs:
        txt += section_text(y, jet, hdisk, xp, yp, afile)
    return txt, semispan, Sref


def strips_per_pitch(s_over_D, D, strips_per_D):
    rate = strips_per_D / D
    if s_over_D <= 0.0:
        return int(math.ceil(rate * D))
    return int(math.ceil(rate * D)) + 2 * int(math.ceil(rate * 0.5 * s_over_D * D))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--jvl", default=cal.DEFAULT_JVL)
    p.add_argument("--out-root", default=cal.DEFAULT_OUT_ROOT)
    p.add_argument("--tag", default="")
    p.add_argument("--s-over-D", type=float, nargs="*", default=DEFAULT_S_OVER_D)
    p.add_argument("--alpha", type=float, nargs="*", default=DEFAULT_ALPHA)
    p.add_argument("--Tc", type=float, nargs="*", default=DEFAULT_TC, help="T / (q pi D^2/4) per disk")
    p.add_argument("--D", type=float, default=0.5, help="propeller diameter / chord")
    p.add_argument("--xp", type=float, default=-0.33, help="disk x / chord (written as dxdisk)")
    p.add_argument("--yp", type=float, default=-0.05, help="disk vertical offset / chord (written as dndisk)")
    p.add_argument("--mach", type=float, default=0.18)
    p.add_argument("--camber", choices=["naca43018", "flat"], default="naca43018")
    p.add_argument("--strips-per-D", type=int, default=6)
    p.add_argument("--nchord", type=int, default=8)
    p.add_argument("--max-wing-strips", type=int, default=240, help="both halves; 240 x (Nchord + Nwjet = 2 x 8) = 3840 vortices < nvmax = 5000 of the JVL 2.16 source")
    p.add_argument("--max-pitches", type=int, default=40)
    p.add_argument("--smoke", action="store_true", help="s/D = 0, 0.05, 1.1 at alpha = 6 only")
    args = p.parse_args(argv)
    if args.smoke:
        args.s_over_D, args.alpha = [0.0, 0.05, 1.1], [6.0]

    name = f"{date.today().strftime('%d%m%y')}_hawkswell_periodic" + (f"_{args.tag}" if args.tag else "") + ("_smoke" if args.smoke else "")
    out_dir = os.path.join(args.out_root, name)
    os.makedirs(out_dir, exist_ok=True)
    afile = None
    if args.camber == "naca43018":
        write_afile(os.path.join(out_dir, AFILE_NAME))
        afile = "../" + AFILE_NAME  # relative to the case folder (JVL runs with cwd = case folder)

    cal.RUN_CASE = cal.RUN_CASE.replace(" Mach      =   0.0\n", f" Mach      =   {args.mach}\n")  # run-case Mach (in memory only)
    r_values = [math.sqrt(1.0 + tc) for tc in args.Tc]
    dvjets = [0.0] + [r - 1.0 for r in r_values]

    rows = []
    t0 = time.time()
    for sD in args.s_over_D:
        spp = strips_per_pitch(sD, args.D, args.strips_per_D)
        n2 = max(2, min(args.max_pitches, (args.max_wing_strips // 2) // spp))
        n2 -= n2 % 2
        for n_pitch in (n2 // 2, n2):
            geom, semispan, Sref = geometry_text(sD, args.D, n_pitch, args.mach, args.xp, args.yp, afile, args.nchord)
            geom, info = postprocess_jvl_text(geom, fix_gain=False, D=args.D, strips_per_D=args.strips_per_D,
                                              n_span_min=1, nchord=args.nchord, max_wing_strips=args.max_wing_strips)
            AR = (2.0 * semispan) ** 2 / Sref
            for a in args.alpha:
                folder = os.path.join(out_dir, f"sD{sD:g}_N{n_pitch}_a{a:g}")
                res = cal.run_jvl(args.jvl, folder, geom, dvjets, a)
                for k, st in enumerate(res):
                    rows.append({
                        "s_over_D": sD, "n_pitch": n_pitch, "semispan": semispan, "AR": AR, "alpha_deg": a,
                        "Tc": 0.0 if k == 0 else args.Tc[k - 1], "Vjet_over_V": 1.0 if k == 0 else r_values[k - 1],
                        "strips_wing": info["nspan"][0][2] * 2 if info["nspan"] else np.nan,
                        "strip_cap_factor": info["strip_cap_factor"],
                        **{key: st.get(key, np.nan) for key in ("CLff", "CLtot", "CLffc", "CLffj", "CLcir", "CLjet",
                                                                "CDff", "CDtot", "CTtot", "DCQtot", "e")},
                    })
            print(f"s/D = {sD:g}: N = {n_pitch} pitches/semispan, AR = {AR:.1f}, strips (both halves) = "
                  f"{info['nspan'][0][2] * 2 if info['nspan'] else '?'}, cap factor = {info['strip_cap_factor']}  "
                  f"[{time.time() - t0:.0f} s]")

    raw = pd.DataFrame(rows)
    blown_frac = 1.0 / (1.0 + raw["s_over_D"])
    hdisk = math.pi * args.D / 4.0
    raw["CT_expected"] = hdisk * blown_frac * (raw["Vjet_over_V"] ** 2 - 1.0)  # = Tc * (pi D^2/4) / (c (D + s))
    raw["CT_rel_err"] = raw["CTtot"] / raw["CT_expected"] - 1.0
    raw.to_csv(os.path.join(out_dir, "hawkswell_jvl_raw.csv"), index=False)

    # AR -> infinity (Richardson, leading 1/AR term)
    out = []
    keys = ["CLff", "CLtot", "CLffc", "CLffj", "CDff", "CTtot"]
    for (sD, a, tc), g in raw.groupby(["s_over_D", "alpha_deg", "Tc"]):
        g = g.sort_values("AR")
        if len(g) != 2:
            continue
        lo, hi = g.iloc[0], g.iloc[1]
        row = {"s_over_D": sD, "alpha_deg": a, "Tc": tc, "Vjet_over_V": hi["Vjet_over_V"], "AR1": lo["AR"], "AR2": hi["AR"]}
        for k in keys:
            row[k + "_AR2"] = hi[k]
            row[k + "_inf"] = hi[k] + (hi[k] - lo[k]) * lo["AR"] / (hi["AR"] - lo["AR"])
        out.append(row)
    res = pd.DataFrame(out)
    sa = np.sin(np.radians(res["alpha_deg"]))
    ca = np.cos(np.radians(res["alpha_deg"]))
    res["CL_hawkswell_inf"] = res["CLff_inf"] + res["CTtot_inf"] * sa  # thesis eq. (1): L + T sin(alpha)
    res["Cx_hawkswell_inf"] = res["CDff_inf"] - res["CTtot_inf"] * ca  # thesis eq. (3), inviscid
    res.to_csv(os.path.join(out_dir, "hawkswell_jvl_results.csv"), index=False)
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(res[["s_over_D", "alpha_deg", "Tc", "AR1", "AR2", "CLff_AR2", "CLff_inf", "CTtot_inf", "CL_hawkswell_inf"]].round(4))
        print("CT check (JVL / expected - 1): min %.4f, max %.4f" % (raw.loc[raw.Tc > 0, "CT_rel_err"].min(), raw.loc[raw.Tc > 0, "CT_rel_err"].max()))
    print("Written to", out_dir)


if __name__ == "__main__":
    main()
