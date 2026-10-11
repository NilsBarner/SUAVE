r"""
jvl_calibration_2d_claude.py

CLAUDE (04.10.2026): JVL calibration on a 2D-like (AR = 100), flat-plate, uniformly or
band-blown wing - Step 1 of jvl_blown_wing_verification_plan_041026.md. No SUAVE, no TASOPT:
it writes .jvl/.run/.deck files itself (format copied from
SUAVE\nils\jvl\jvl2.16\JVL2.16\bin\jvl_calibrate.jvl and the SUAVE-generated batch files),
calls jvl.exe and parses the outputs with jvl_patches_claude.

What it checks
--------------
(1) Jet strength convention: with gain = 1 and DVjet = r - 1, JVL must give
        DCQ = (r - 1) * hjet * L / Sref,   CT = 2 * hjet * L / Sref * r * (r - 1),
    hjet = hdisk * 0.5 * (1 + 1/r) for fh = 1, L = blown span (both wing halves).
    This validates the fix of issue B1 (JVL gain = Vjet/V - 1).
(2) Jet spanwise confinement (issue B10): for a band-blown wing (JETCONTROL only on the two
    sections bounding the band) the same DCQ formula must hold with L = band width,
    i.e. no jet ramps into the neighbouring unblown interval.
(3) "Effective beta" of JVL's thin-jet model: for the uniformly blown wing,
        beta_JVL = a_eff / a,  with  dcl_JVL = 2 pi sin(alpha) [(1 + a_eff)^2 - 1],  a = r - 1,
    i.e. the value of the surrogate's finite-slipstream factor beta (Patterson) that would make
    the surrogate's 2D lift increment (alpha_p = 0) equal JVL's, as a function of hjet/c.

Run (any Python with numpy + pandas, e.g. the "suave" or "py311b" environment):

  python C:\Users\nmb48\Documents\GitHub\SUAVE\nils\jvl\code_claude\jvl_calibration_2d_claude.py

Outputs: SUAVE\nils\jvl\data_claude\<ddmmyy>_calibration_2d\calibration_2d_results.csv
"""

import argparse
import math
import os
import shutil
import subprocess
import sys
from datetime import date

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from jvl_patches_claude import parse_stability_file, parse_surface_forces_file  # noqa: E402

DEFAULT_JVL = r"C:\Users\nmb48\Documents\GitHub\SUAVE\nils\jvl\jvl2.16\jvl.exe"
DEFAULT_OUT_ROOT = r"C:\Users\nmb48\Documents\GitHub\SUAVE\nils\jvl\data_claude"

SEMI_SPAN = 50.0  # chord = 1 -> AR = 100  # CLAUDE (04.10.2026, fix 1): was 250 (AR = 500), which stopped JVL with "SDUPL: Vortex array overflow"
CHORD = 1.0


def section_text(y, jet, hdisk, dxdisk):
    txt = (
        "SECTION\n"
        "#Xle    Yle      Zle      Chord     Ainc\n"
        f"0.0  {y:.6f}  0.0  {CHORD}  0.0\n"
    )
    if jet:
        txt += (
            "JETCONTROL\n"
            "#Jname   Jgain    SgnDup\n"
            "DVjet    1.0      1.0\n"
            "JETPARAM\n"
            "#hdisk   fh      djet0   djet1   djet3       dxdisk  dndisk\n"
            f"{hdisk:.6f}    1.0     0.0    0.0   0.0    {dxdisk:.6f}   0.0\n"
        )
    return txt + "\n"


def geometry_text(layout, hdisk, band_half_width, nspan, nchord, dxdisk, eps):
    Sref = 2.0 * SEMI_SPAN * CHORD
    txt = (
        "calibration_claude\n"
        "#Mach\n0.0\n"
        "#IYsym IZsym Zsym\n0 0 0.0\n"
        "#Sref   Cref   Bref\n"
        f"{Sref} {CHORD} {2.0 * SEMI_SPAN}\n"
        "#Xref   Yref   Zref\n0.25 0.0 0.0\n"
        "#CDp\n0.0\n"
        "SURFACE\nWing\n"
        "#Nchordwise  Cspace   Nspanwise  Sspace\n"
        f"{nchord}  1.0  {nspan}  0.0  0  1.0  {nchord}  1.0\n"  # CLAUDE (04.10.2026, fix 1): explicit Nujet = 0, Nwjet = Nchord
        "YDUPLICATE\n0.0\n\n"
    )
    if layout == "full":
        txt += section_text(0.0, True, hdisk, dxdisk)
        txt += section_text(SEMI_SPAN, True, hdisk, dxdisk)
    elif layout == "band":
        txt += section_text(0.0, True, hdisk, dxdisk)
        txt += section_text(band_half_width, True, hdisk, dxdisk)
        if eps > 0:
            txt += section_text(band_half_width + eps, False, hdisk, dxdisk)
        txt += section_text(SEMI_SPAN, False, hdisk, dxdisk)
    else:
        raise ValueError(layout)
    return txt


RUN_CASE = """
 ---------------------------------------------
 Run case  {idx}:   case_{idx:02d}

 alpha        ->  alpha       =   {alpha:.5f}
 beta         ->  beta        =   0.0
 pb/2V        ->  pb/2V       =   0.0
 qc/2V        ->  qc/2V       =   0.0
 rb/2V        ->  rb/2V       =   0.00000
 DVjet        ->  DVjet       =   {dvjet:.5f}

 alpha     =   0.00000     deg
 beta      =   0.00000     deg
 pb/2V     =   0.0
 qc/2V     =   0.0
 rb/2V     =   0.00000
 CL        =   0.00000
 CDo       =   0.0
 bank      =   0.00000     deg
 elevation =   0.00000     deg
 heading   =   0.00000     deg
 Mach      =   0.0
 velocity  =   0.0     m/s
 density   =   1.225     kg/m^3
 grav.acc. =   9.81     m/s^2
 turn_rad. =   0.00000     m
 load_fac. =   1.0
 X_cg      =   0.25     m
 Y_cg      =   0.0     m
 Z_cg      =   0.0     m
 mass      =   1.0     kg
 Ixx       =   1.0     kg-m^2
 Iyy       =   1.0     kg-m^2
 Izz       =   1.0     kg-m^2
 Ixy       =   0.0     kg-m^2
 Iyz       =   0.0     kg-m^2
 Izx       =   0.0     kg-m^2
"""


def run_jvl(jvl_exe, folder, geom_txt, dvjets, alpha_deg, timeout=120):
    if os.path.exists(folder):
        shutil.rmtree(folder)
    os.makedirs(folder)
    with open(os.path.join(folder, "cal.jvl"), "w") as f:
        f.write(geom_txt)
    with open(os.path.join(folder, "cal.run"), "w") as f:
        for i, dv in enumerate(dvjets, start=1):
            f.write(RUN_CASE.format(idx=i, alpha=alpha_deg, dvjet=dv))
    deck = "CASE cal.run\nOPER\n"
    for i in range(1, len(dvjets) + 1):
        deck += f"{i}\nx\nst\nst_{i:02d}.txt\nfn\nfn_{i:02d}.txt\n"
    deck += "\nQUIT\n"
    with open(os.path.join(folder, "commands.deck"), "w") as f:
        f.write(deck)
    with open(os.path.join(folder, "jvl_log.txt"), "w") as log:
        proc = subprocess.Popen([jvl_exe, "cal.jvl"], cwd=folder, stdin=subprocess.PIPE, stdout=log, stderr=log)
        try:
            proc.communicate(deck.encode("utf-8"), timeout=timeout)
        except subprocess.TimeoutExpired:
            proc.kill()
            raise RuntimeError(f"JVL timed out in {folder}")
    out = []
    for i in range(1, len(dvjets) + 1):
        p = os.path.join(folder, f"st_{i:02d}.txt")
        if not os.path.exists(p):
            # CLAUDE (04.10.2026, fix 1): report why JVL stopped (e.g. "Vortex array overflow") instead of only the missing file
            with open(os.path.join(folder, "jvl_log.txt"), "r") as log:
                tail = "".join(log.readlines()[-8:])
            raise FileNotFoundError(f"{p}\n--- last lines of jvl_log.txt ---\n{tail}")
        out.append(parse_stability_file(p))
    return out


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--jvl", default=DEFAULT_JVL)
    p.add_argument("--out-root", default=DEFAULT_OUT_ROOT)
    p.add_argument("--alpha", type=float, default=5.0, help="angle of attack (deg)")
    p.add_argument("--hdisk", type=float, nargs="*", default=[0.1, 0.25, 0.5, 1.0], help="hdisk / chord")
    p.add_argument("--r", type=float, nargs="*", default=[1.05, 1.1, 1.2, 1.4], help="Vjet / Vinf")
    p.add_argument("--band-half-width", type=float, default=10.0, help="blown half-band width (chords)")  # CLAUDE (04.10.2026, fix 1): was 25
    p.add_argument("--nspan", type=int, default=40)  # CLAUDE (04.10.2026, fix 1): was 100
    p.add_argument("--nchord", type=int, default=8)  # CLAUDE (04.10.2026, fix 1): was 12
    p.add_argument("--dxdisk", type=float, default=-0.2, help="propulsor x offset from LE (chords)")
    args = p.parse_args(argv)

    out_dir = os.path.join(args.out_root, f"{date.today().strftime('%d%m%y')}_calibration_2d")
    os.makedirs(out_dir, exist_ok=True)
    Sref = 2.0 * SEMI_SPAN * CHORD
    a_rad = math.radians(args.alpha)
    rows = []
    dvjets = [0.0] + [r - 1.0 for r in args.r]
    for layout, eps in (("full", 0.0), ("band", 0.0), ("band", 2.5)):  # eps in chords  # CLAUDE (04.10.2026, fix 3): was 0.5 (< one strip of 50/40 = 1.25 c, so JVL could not "adjust spanwise spacing"); 2.5 c = 2 strips, on the strip grid
        for h in args.hdisk:
            folder = os.path.join(out_dir, f"{layout}_eps{eps:g}_h{h:g}")
            res = run_jvl(args.jvl, folder,
                          geometry_text(layout, h, args.band_half_width, args.nspan, args.nchord, args.dxdisk, eps),
                          dvjets, args.alpha)
            L = 2.0 * SEMI_SPAN if layout == "full" else 2.0 * args.band_half_width
            off = res[0]
            for r, on in zip(args.r, res[1:]):
                hjet = h * 0.5 * (1.0 + 1.0 / r)
                dcl = on["CLtot"] - off["CLtot"]
                ref_beta1 = 2 * math.pi * math.sin(a_rad) * (r ** 2 - 1.0)  # surrogate 2D formula, beta = 1, a = r - 1
                # effective a (hence beta) reproducing JVL's increment: 2 pi sin(a) [(1+a_eff)^2 - 1] = dcl_2D
                dcl_2d = dcl * (2.0 * SEMI_SPAN) / L  # per unit blown span (exact for 'full')
                a_eff = math.sqrt(1.0 + dcl_2d / (2 * math.pi * math.sin(a_rad))) - 1.0 if dcl_2d > -2 * math.pi * math.sin(a_rad) else float("nan")
                rows.append({
                    "layout": layout, "eps": eps, "hdisk_over_c": h, "Vjet_over_V": r, "alpha_deg": args.alpha,
                    "hjet_over_c": hjet, "blown_span_over_c": L,
                    "DCQ_jvl": on.get("DCQtot", np.nan), "DCQ_expected": (r - 1.0) * hjet * L / Sref,
                    "CT_jvl": on.get("CTtot", np.nan), "CT_expected": 2.0 * hjet * L / Sref * r * (r - 1.0),
                    "CL_off": off["CLtot"], "CL_on": on["CLtot"], "delta_CL": dcl,
                    "CLcir_on": on.get("CLcir", np.nan), "CLjet_on": on.get("CLjet", np.nan),
                    "dcl_2D_jvl": dcl_2d, "dcl_2D_surrogate_beta1": ref_beta1,
                    "ratio_jvl_to_beta1": dcl_2d / ref_beta1, "beta_eff_jvl": a_eff / (r - 1.0),
                })
            print(f"{layout} eps={eps:g} hdisk/c={h:g}: done")
    df = pd.DataFrame(rows)
    df["DCQ_rel_err"] = df["DCQ_jvl"] / df["DCQ_expected"] - 1.0
    df["CT_rel_err"] = df["CT_jvl"] / df["CT_expected"] - 1.0
    out_csv = os.path.join(out_dir, "calibration_2d_results.csv")
    df.to_csv(out_csv, index=False)
    with pd.option_context("display.width", 200, "display.max_columns", 30):
        print(df[["layout", "eps", "hdisk_over_c", "Vjet_over_V", "DCQ_rel_err", "CT_rel_err",
                  "dcl_2D_jvl", "dcl_2D_surrogate_beta1", "beta_eff_jvl"]].round(4))
    print("Results written to", out_csv)


if __name__ == "__main__":
    main()
