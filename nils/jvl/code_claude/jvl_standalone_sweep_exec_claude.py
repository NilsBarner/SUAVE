r"""
jvl_standalone_sweep_exec_claude.py

CLAUDE (04.10.2026): stand-alone JVL sweep on a FIXED airframe for the verification of the
blown-wing surrogate (Steps 4-6 of jvl_blown_wing_verification_plan_041026.md).

For every (N, dY, Pi) point of a surrogate reference .csv written by
Flydrogen\flydrogen\analyses\tests_claude\blown_wing_reference_claude.py (Step 0):

  1. lay out N propellers of the surrogate's diameter D on the fixed TASOPT wing
     (blown_wing_layout_claude.make_layout; spacing satisfies de Vries eq. (21) exactly),
  2. set the jet velocity ratio Vjet/V from the surrogate (mode 'V': far wake 1 + 2 a_p,
     which for hdisk = pi D/4 and fh = 1 also reproduces the surrogate thrust exactly;
     mode 'L': leading-edge value V_w/V, sensitivity only),
  3. run JVL once with two run cases at identical alpha (= surrogate alpha_w, eq. (29)):
     DVjet = 0 (unblown) and DVjet = 1 (blown), by calling the EXISTING
     nils.jvl.code.jvl_wrapper.wrap_jvl through the in-memory patches of
     jvl_patches_claude (no existing file is modified),
  4. append one row per point to a results .csv (resumable).

Run from the SUAVE root in the "suave" environment, e.g.

  (suave) C:\Users\nmb48\Documents\GitHub\SUAVE> python -m nils.jvl.code_claude.jvl_standalone_sweep_exec_claude ^
        --reference C:\Users\nmb48\Documents\GitHub\Flydrogen\flydrogen\analyses\tests_claude\data_claude\blown_wing_surrogate_reference_takeoff_linear.csv ^
        --smoke

See README_claude.md in this folder for all options and the recommended run order.
"""

__all__ = []

import argparse
import os
import sys
import time
from datetime import date

import numpy as np
import pandas as pd

import SUAVE
assert SUAVE.__version__ == '2.5.0', 'These tutorials only work with the SUAVE 2.5.0 release'
from SUAVE.Core import Units

from nils.jvl.code.jvl_wrapper import wrap_jvl  # existing wrapper, reused unchanged
from nils.jvl.code_claude.jvl_patches_claude import apply_patches, set_overrides
from nils.jvl.code_claude.blown_wing_layout_claude import make_layout, base_mass_row

DEFAULT_GEOMETRY = r"C:\Users\nmb48\Documents\GitHub\TASOPT.jl-priv\mass_distr_geometry_regional_LH2_041026.csv"
DEFAULT_OUT_ROOT = r"C:\Users\nmb48\Documents\GitHub\SUAVE\nils\jvl\data_claude"


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--reference", required=True, help="surrogate reference .csv from Step 0")
    p.add_argument("--geometry", default=DEFAULT_GEOMETRY, help="baseline TASOPT geometry .csv")
    p.add_argument("--mode", choices=["V", "L"], default="V",
                   help="jet velocity: V = far wake 1+2a_p (primary, thrust-consistent); L = LE value V_w/V")
    p.add_argument("--N", type=int, nargs="*", default=None, help="subset of N (default: all in reference)")
    p.add_argument("--dY", type=float, nargs="*", default=None, help="subset of dY (default: all in reference)")
    p.add_argument("--pi-stride", type=int, default=7,
                   help="use every k-th Pi of the reference grid (50 Pi values -> 8 with k = 7)")
    p.add_argument("--full-config", action="store_true",
                   help="keep tails, fuselage and nacelle bodies (default: wing-only model, issue B14)")
    p.add_argument("--cambered", action="store_true",
                   help="keep the B737a AFILE camber (default: flat camber line, issue B12)")
    p.add_argument("--strips-per-D", type=float, default=12.0, help="spanwise strips per propeller diameter")
    p.add_argument("--n-span-min", type=int, default=20, help="minimum spanwise strips per main-wing surface")
    p.add_argument("--sspace", type=float, default=0.0, help="spanwise spacing parameter (0 = uniform)")
    p.add_argument("--nchord", type=int, default=10, help="chordwise vortices on the main wing")
    p.add_argument("--max-wing-strips", type=int, default=360,
                   help="cap on total main-wing strips (avoids JVL 'Vortex array overflow')")  # CLAUDE (04.10.2026, fix 2)
    p.add_argument("--no-fix-gain", action="store_true",
                   help="write the JETCONTROL gain as is (reproduces issue B1; for comparison only)")
    p.add_argument("--out-root", default=DEFAULT_OUT_ROOT)
    p.add_argument("--tag", default="", help="extra tag appended to the output folder name")
    p.add_argument("--smoke", action="store_true", help="3-point smoke test: N = 8, dY = 0.6, Pi = min/mid/max")
    p.add_argument("--limit", type=int, default=None, help="stop after this many JVL runs")
    p.add_argument("--print-output", action="store_true", help="show JVL console output")
    return p.parse_args(argv)


def select_points(ref, args):
    ref = ref.copy()
    if args.smoke:
        ref = ref[(ref["N"] == 8) & (np.isclose(ref["dY"], 0.6)) & (ref["dy"] >= 0.0)]  # CLAUDE (04.10.2026, fix 4): feasible points only
        idx = sorted(ref["pi_index"].unique())
        keep = [idx[0], idx[len(idx) // 2], idx[-1]]
        return ref[ref["pi_index"].isin(keep)]
    if args.N:
        ref = ref[ref["N"].isin(args.N)]
    if args.dY:
        ref = ref[np.any([np.isclose(ref["dY"], v) for v in args.dY], axis=0)]
    n_pi = ref["pi_index"].max() + 1
    keep = [i for i in range(n_pi) if i % args.pi_stride == 0]
    if (n_pi - 1) not in keep:
        keep.append(n_pi - 1)  # always include the largest Pi
    return ref[ref["pi_index"].isin(keep)]


def main(argv=None):
    args = parse_args(argv)
    ref = pd.read_csv(args.reference)
    op = str(ref["op_point"].iloc[0])
    df_base = pd.read_csv(args.geometry)
    mass_row = base_mass_row(df_base)

    model = "full" if args.full_config else "wing"
    camber = "camb" if args.cambered else "flat"
    run_name = f"{date.today().strftime('%d%m%y')}_{op}_{args.mode}_{model}_{camber}_spd{args.strips_per_D:g}"
    if args.no_fix_gain:
        run_name += "_nofixgain"
    if args.smoke:
        run_name += "_smoke"
    if args.tag:
        run_name += "_" + args.tag
    out_dir = os.path.join(args.out_root, run_name)
    runs_dir = os.path.join(out_dir, "runs")
    os.makedirs(runs_dir, exist_ok=True)
    out_csv = os.path.join(out_dir, "jvl_standalone_sweep_results.csv")
    print("Output folder:", out_dir)

    done = set()
    if os.path.exists(out_csv):
        done = set(pd.read_csv(out_csv)["run_tag"].astype(str))

    apply_patches()
    pts = select_points(ref, args)
    print(f"{len(pts)} points selected ({op}, mode {args.mode}, {model}, {camber})")

    n_run = 0
    cwd0 = os.getcwd()
    for _, r in pts.iterrows():
        N, dY, Pi, pi_index = int(r["N"]), float(r["dY"]), float(r["Pi"]), int(r["pi_index"])
        run_tag = f"N{N}_dY{dY:.3f}_pi{pi_index:02d}"
        if run_tag in done:
            continue
        rec = {
            "run_tag": run_tag, "op_point": op, "mode": args.mode, "model": model, "camber": camber,
            "N": N, "dY": dY, "Pi": Pi, "pi_index": pi_index, "D_prop": r["D_prop"], "dy": r["dy"],
            "dy_m": r["dy_m"], "alpha_deg": r["alpha_w_deg"], "Mach": r["Mach"], "altitude_m": r["altitude_m"],
            "delta_CL_surrogate": r["delta_CL"], "CT_target": r["CT"],
        }
        Vjet_over_V = float(r["Vjet_far_over_V"]) if args.mode == "V" else float(r["V_w_over_V"])
        rec["Vjet_over_V"] = Vjet_over_V

        df_geom, info = make_layout(df_base, N, float(r["D_prop"]), dY, ax_pos=float(r["ax_pos"]), Jgain=Vjet_over_V)
        rec["y0_inside_fuselage_m"] = info["y0_inside_fuselage_m"]
        if df_geom is None:
            rec["status"] = "infeasible: " + info["reason"]
            _append(out_csv, rec)
            continue

        alpha = float(r["alpha_w_deg"])
        set_overrides(
            Prop_PR_des=Pi,
            training_angle_of_attack=np.array([alpha, alpha]) * Units.degrees,  # two run cases: unblown, blown
            training_Mach=np.array([float(r["Mach"])]),
            training_side_slip_angle=np.array([0.0]) * Units.degrees,
            training_altitude=np.array([float(r["altitude_m"])]),
            training_load_factor=np.array([1.0]),
            training_mass=np.array([float(mass_row["takeoff"])]),
            print_output=True,  # CLAUDE (04.10.2026, fix 2): JVL console output -> avl_log.txt in each run folder (False sends it to os.devnull)
            jvl_options=dict(
                DVjet_cases=[0.0, 1.0],
                fix_gain=not args.no_fix_gain,
                wing_only=not args.full_config,
                symmetric_airfoil=not args.cambered,
                D=float(r["D_prop"]),
                strips_per_D=args.strips_per_D,
                n_span_min=args.n_span_min,
                sspace=args.sspace,
                nchord=args.nchord,
                max_wing_strips=args.max_wing_strips,  # CLAUDE (04.10.2026, fix 2)
            ),
        )

        t0 = time.time()
        os.chdir(runs_dir)  # sample_training creates its run folder relative to the cwd
        try:
            results_list = wrap_jvl(df_geom, 0, mass_row, True, run_tag)
            res = results_list[0]
            off = res["case_01_01"].aerodynamics
            on = res["case_01_02"].aerodynamics
            rec.update(
                status="ok",
                alpha_jvl=on.get("alpha", np.nan),
                CL_off=off.CLtot, CL_on=on.CLtot, delta_CL_total=on.CLtot - off.CLtot,
                CLwing_off=off.CL_main_wing, CLwing_on=on.CL_main_wing,
                delta_CL_wing=on.CL_main_wing - off.CL_main_wing,
                CLcir_off=off.get("CLcir", np.nan), CLcir_on=on.get("CLcir", np.nan),
                delta_CLcir=on.get("CLcir", np.nan) - off.get("CLcir", np.nan),
                CLjet_on=on.get("CLjet", np.nan),
                # CLAUDE (04.10.2026, fix 5): Trefftz-plane lift - grid-converged (0.2039/0.2020/0.2008 at 6/12/24
                # strips per D), whereas the surface-force CLtot is inconsistent with it once the jets are on
                CL_ff_off=off.get("CLff", np.nan), CL_ff_on=on.get("CLff", np.nan),
                delta_CL_ff=on.get("CLff", np.nan) - off.get("CLff", np.nan),
                CL_ffc_on=on.get("CLffc", np.nan),
                delta_CL_ffc=on.get("CLffc", np.nan) - off.get("CLffc", np.nan),
                CLffj_on=on.get("CLffj", np.nan),
                CDind_off=off.get("CDind", np.nan), CDind_on=on.get("CDind", np.nan),
                CDtot_on=on.get("CDtot", np.nan), CDff_on=on.get("CDff", np.nan),
                DCQ_on=on.get("DCQtot", np.nan), DCJ_on=on.get("DCJtot", np.nan),
                CT_jvl=on.get("CTtot", np.nan), CP_jvl=on.get("CPtot", np.nan),
            )
        except (FileNotFoundError, RuntimeError, KeyError, ValueError, IndexError) as e:
            rec["status"] = f"failed: {type(e).__name__}: {e}"
        finally:
            os.chdir(cwd0)
        rec["runtime_s"] = round(time.time() - t0, 2)
        _append(out_csv, rec)
        print(run_tag, rec["status"], rec.get("delta_CL_wing"), "(surrogate", rec["delta_CL_surrogate"], ")")
        n_run += 1
        if args.limit is not None and n_run >= args.limit:
            break

    print("Done. Results in", out_csv)


COLUMNS = [
    "run_tag", "status", "op_point", "mode", "model", "camber", "N", "dY", "Pi", "pi_index", "D_prop", "dy", "dy_m",
    "y0_inside_fuselage_m", "alpha_deg", "alpha_jvl", "Mach", "altitude_m", "Vjet_over_V", "CT_target", "CT_jvl", "CP_jvl",
    "DCQ_on", "DCJ_on", "delta_CL_surrogate",
    "delta_CL_ff", "delta_CL_ffc", "CLffj_on", "CL_ff_off", "CL_ff_on", "CL_ffc_on",  # CLAUDE (04.10.2026, fix 5)
    "delta_CL_wing", "delta_CL_total", "delta_CLcir", "CLjet_on",
    "CL_off", "CL_on", "CLwing_off", "CLwing_on", "CLcir_off", "CLcir_on", "CDind_off", "CDind_on", "CDtot_on",
    "CDff_on", "runtime_s",
]


def _append(path, rec):
    df = pd.DataFrame([rec]).reindex(columns=COLUMNS)  # fixed column order, missing entries -> NaN
    if os.path.exists(path) and list(pd.read_csv(path, nrows=0).columns) != COLUMNS:
        # CLAUDE (04.10.2026, fix 5): existing file written with other columns -> rewrite with the current ones
        df = pd.concat([pd.read_csv(path).reindex(columns=COLUMNS), df], ignore_index=True)
        df.to_csv(path, index=False)
        return
    df.to_csv(path, mode="a", header=not os.path.exists(path), index=False)


if __name__ == "__main__":
    main()
