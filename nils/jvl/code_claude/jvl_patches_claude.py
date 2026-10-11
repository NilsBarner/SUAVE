r"""
jvl_patches_claude.py

CLAUDE (04.10.2026): run-time patches and result parsers for the existing SUAVE-JVL
path, written so that NO existing file in SUAVE has to be modified. See Section 3 and
Step 2 of jvl_blown_wing_verification_plan_041026.md.

How it works
------------
`apply_patches()` replaces, *in memory only* (for the current Python process), four
names that `SUAVE.Analyses.Stability.AVL` looks up at call time:

  AVL.sample_training                     -> _sample_training_claude
      Calls the original method unchanged, after replacing selected keyword
      arguments with the values set via `set_overrides()` (operating point,
      Prop_PR_des, ...). This lets the existing `nils.jvl.code.jvl_wrapper.wrap_jvl`
      be reused verbatim, although it hard-codes alpha = 3 deg, M = 0.15 and
      Prop_PR_des = 1.005 (issues B3, B5, B7).

  <AVL module>.write_geometry             -> _write_geometry_claude
      Calls the original write_geometry (incl. adjust_wing_text_for_jvl) and then
      post-processes the written .jvl file (`postprocess_jvl_text`):
        - JETCONTROL gain g -> g - 1, because JVL defines Vjet/Vinf = 1 + gain*DVjet
          whereas TASOPT/wrap_jvl pass the velocity ratio Jgain = Vjet/Vinf (issue B1)
        - optional wing-only model: drops all non-main-wing SURFACEs and all BODYs (B14)
        - optional symmetric (flat camber line) sections: drops AFILE lines (B12)
        - optional spanwise resolution: Nspanwise per main-wing surface such that each
          propeller disk is covered by >= `strips_per_D` strips, with spacing `sspace` (B13)

  <AVL module>.write_run_cases            -> _write_run_cases_claude
      Calls the original write_run_cases and then sets DVjet per run case from
      `DVjet_cases` (e.g. [0.0, 1.0] -> unblown / blown at identical alpha) (B4).

  <AVL module>.run_analysis               -> _run_analysis_claude
      Calls the original `call_avl` and parses the JVL result files with regular
      expressions instead of the fixed-width slicing in read_results.py (B17); adds
      CLcir, the main-wing CL from the surface-forces file, etc.

The patches are only active in a process that calls `apply_patches()`. All other
scripts (incl. jvl_point_exec.py / jvl_sweep_exec.py) are unaffected.

The text-processing functions (`postprocess_jvl_text`, `set_dvjet_in_run_text`,
`parse_stability_file`, `parse_surface_forces_file`) do not import SUAVE and can be
tested stand-alone.
"""

__all__ = [
    "apply_patches",
    "set_overrides",
    "postprocess_jvl_text",
    "set_dvjet_in_run_text",
    "parse_stability_file",
    "parse_surface_forces_file",
]

import math
import os
import re
import sys

# -----------------------------------------------------------------------------
# Global state (process-local)
# -----------------------------------------------------------------------------

_OVERRIDES = {}  # keyword arguments of AVL.sample_training to replace, plus 'jvl_options'
_ORIG = {}  # original SUAVE callables, filled by apply_patches()


def set_overrides(**kwargs):
    """Set the sample_training keyword overrides used by the next wrap_jvl call.

    Recognised keys: any keyword argument of AVL.sample_training (e.g. Prop_PR_des,
    training_angle_of_attack, training_Mach, training_altitude, training_mass, ...)
    and `jvl_options` (dict, see postprocess_jvl_text and _write_run_cases_claude).
    """
    _OVERRIDES.clear()
    _OVERRIDES.update(kwargs)


def apply_patches():
    """Install the in-memory patches (idempotent)."""
    if _ORIG:
        return
    import SUAVE  # noqa: F401
    from SUAVE.Analyses.Stability import AVL as AVLClass  # the package __init__ re-binds AVL to the class
    from SUAVE.Methods.Aerodynamics.AVL.run_analysis import call_avl

    avl_module = sys.modules["SUAVE.Analyses.Stability.AVL"]

    _ORIG["AVLClass"] = AVLClass
    _ORIG["sample_training"] = AVLClass.sample_training
    _ORIG["write_geometry"] = avl_module.write_geometry
    _ORIG["write_run_cases"] = avl_module.write_run_cases
    _ORIG["run_analysis"] = avl_module.run_analysis
    _ORIG["call_avl"] = call_avl

    AVLClass.sample_training = _sample_training_claude
    avl_module.write_geometry = _write_geometry_claude
    avl_module.write_run_cases = _write_run_cases_claude
    avl_module.run_analysis = _run_analysis_claude


def remove_patches():
    """Restore the original SUAVE callables."""
    if not _ORIG:
        return
    avl_module = sys.modules["SUAVE.Analyses.Stability.AVL"]
    _ORIG["AVLClass"].sample_training = _ORIG["sample_training"]
    avl_module.write_geometry = _ORIG["write_geometry"]
    avl_module.write_run_cases = _ORIG["write_run_cases"]
    avl_module.run_analysis = _ORIG["run_analysis"]
    _ORIG.clear()


# -----------------------------------------------------------------------------
# Patched SUAVE callables
# -----------------------------------------------------------------------------

def _sample_training_claude(self, **kwargs):
    jvl_options = _OVERRIDES.get("jvl_options", None)
    for key, value in _OVERRIDES.items():
        if key in ("jvl_options", "print_output"):  # CLAUDE (04.10.2026, fix 4): not sample_training arguments
            continue
        kwargs[key] = value
    # Options read by _write_geometry_claude / _write_run_cases_claude
    self.settings.claude_options = jvl_options
    print_output = _OVERRIDES.get("print_output", None)
    if print_output is not None:
        self.settings.print_output = print_output
    return _ORIG["sample_training"](self, **kwargs)


def _write_geometry_claude(avl_object, run_script_path):
    _ORIG["write_geometry"](avl_object, run_script_path)
    opts = getattr(avl_object.settings, "claude_options", None)
    if not opts:
        return
    geometry_file = avl_object.settings.filenames.features
    with open(geometry_file, "r") as f:
        text = f.read()
    keys = ("fix_gain", "wing_only", "symmetric_airfoil", "D", "strips_per_D", "n_span_min", "sspace", "nchord",
            "max_wing_strips")  # CLAUDE (04.10.2026, fix 2): added max_wing_strips
    new_text, info = postprocess_jvl_text(text, **{k: opts[k] for k in keys if k in opts})
    with open(geometry_file, "w") as f:
        f.write(new_text)
    with open("claude_jvl_postprocess_info.txt", "w") as f:  # written into the run folder for traceability
        for k, v in info.items():
            f.write(f"{k} = {v}\n")


def _write_run_cases_claude(avl_object, trim_aircraft, backend="AVL"):
    _ORIG["write_run_cases"](avl_object, trim_aircraft, backend)
    opts = getattr(avl_object.settings, "claude_options", None)
    if not opts or opts.get("DVjet_cases") is None:
        return
    batch_file = avl_object.current_status.batch_file
    with open(batch_file, "r") as f:
        text = f.read()
    with open(batch_file, "w") as f:
        f.write(set_dvjet_in_run_text(text, opts["DVjet_cases"]))


def _run_analysis_claude(avl_object, print_output, backend="AVL"):
    from SUAVE.Core import Data

    _ORIG["call_avl"](avl_object, print_output)
    results = Data()
    for case in avl_object.current_status.cases:
        if not os.path.exists(case.aero_result_filename_1):
            # CLAUDE (04.10.2026, fix 2): include why JVL stopped (e.g. "Vortex array overflow") in the error
            tail = ""
            log_file = avl_object.settings.filenames.log_filename
            if isinstance(log_file, str) and os.path.exists(log_file):
                with open(log_file, "r") as log:
                    tail = "".join(log.readlines()[-8:])
            raise FileNotFoundError(
                f"JVL did not write {case.aero_result_filename_1} in {os.getcwd()} - last lines of {log_file}:\n{tail}"
            )
        tot = parse_stability_file(case.aero_result_filename_1)
        surf = (
            parse_surface_forces_file(case.aero_result_filename_2)
            if os.path.exists(case.aero_result_filename_2) else []
        )
        case_res = Data()
        case_res.tag = case.tag
        case_res.aerodynamics = Data()
        aero = case_res.aerodynamics
        for k, v in tot.items():
            aero[k] = v
        # names used elsewhere in the existing code (sample_training prints these)
        aero.total_lift_coefficient = tot.get("CLtot", float("nan"))
        aero.total_drag_coefficient = tot.get("CDtot", float("nan"))
        aero.induced_drag_coefficient = tot.get("CDind", float("nan"))
        aero.AoA = tot.get("alpha", float("nan"))
        aero.CL_main_wing = sum(s["CL"] for s in surf if s["name"].startswith("main_wing"))
        aero.CDi_main_wing = sum(s["CDi"] for s in surf if s["name"].startswith("main_wing"))
        aero.surfaces = surf
        results[case.tag] = case_res
    return results


# -----------------------------------------------------------------------------
# Text processing (SUAVE-independent)
# -----------------------------------------------------------------------------

_NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][-+]?\d+)?"


def _split_blocks(text):
    """Split a .jvl file into a header and SURFACE/BODY blocks.

    Returns (header_lines, blocks), each block a dict with keys 'kind', 'name', 'lines'.
    Lines are kept verbatim (without trailing newline).
    """
    lines = text.splitlines()
    starts = [i for i, ln in enumerate(lines) if ln.strip() in ("SURFACE", "BODY")]
    if not starts:
        return lines, []
    header = lines[: starts[0]]
    blocks = []
    for j, s in enumerate(starts):
        e = starts[j + 1] if j + 1 < len(starts) else len(lines)
        blk = lines[s:e]
        name = ""
        for ln in blk[1:]:
            if ln.strip() and not ln.strip().startswith(("#", "!")):
                name = ln.strip()
                break
        blocks.append({"kind": blk[0].strip(), "name": name, "lines": blk})
    return header, blocks


def _section_y_values(block_lines):
    """Yle of every SECTION in a SURFACE block (first numeric line after 'SECTION')."""
    ys = []
    i = 0
    while i < len(block_lines):
        if block_lines[i].strip() == "SECTION":
            k = i + 1
            while k < len(block_lines):
                s = block_lines[k].strip()
                if s and not s.startswith(("#", "!")):
                    nums = re.findall(_NUM, s)
                    if len(nums) >= 2:
                        ys.append(float(nums[1].replace("D", "E").replace("d", "e")))
                    break
                k += 1
            i = k
        i += 1
    return ys


def _interval_strips(ys, per_length):
    """CLAUDE (04.10.2026, fix 3): strips per section interval (>= 1 each), so that every
    section - incl. closely spaced prop edges, mid-pitch and planform breaks - bounds at least
    one strip. Avoids JVL's "Cannot adjust spanwise spacing at SECTION ... Insufficient number
    of spanwise vortices", which a single surface-level Nspanwise triggers."""
    return [max(1, int(math.ceil(per_length * (y1 - y0)))) for y0, y1 in zip(ys[:-1], ys[1:])]


def postprocess_jvl_text(
    text,
    fix_gain=True,
    wing_only=False,
    symmetric_airfoil=False,
    D=None,
    strips_per_D=None,
    n_span_min=8,
    sspace=0.0,
    nchord=None,
    max_wing_strips=360,
):
    """Post-process a .jvl geometry written by SUAVE's write_geometry (JVL backend).

    Parameters
    ----------
    fix_gain : subtract 1 from every JETCONTROL gain (Jgain = Vjet/V -> JVL gain = Vjet/V - 1)
    wing_only : keep only SURFACEs whose name starts with 'main_wing' (drops tails,
        fuselage and nacelle surfaces/bodies)
    symmetric_airfoil : remove AFILE entries of the main wing (flat camber line)
    D, strips_per_D : if both given, specify the spanwise vortices per SECTION interval of
        each main-wing SURFACE (>= 1 per interval) at a density of
        max(strips_per_D / D, n_span_min / surface_width) strips per unit span
        (CLAUDE, 04.10.2026, fix 3: was one Nspanwise per SURFACE)
    sspace : spanwise spacing parameter written when Nspanwise is reset (0 = uniform)
    nchord : chordwise vortices written when Nspanwise is reset (None = keep)
    max_wing_strips : cap on the total number of main-wing strips (both halves); all
        Nspanwise are scaled down proportionally if exceeded. 360 strips x (Nchord + Nwjet
        = 2 x 10) = 7200 wing vortices, below the 9364/9680 vortices of the 26.01.2026 runs
        that JVL handled, and far below its "Vortex array overflow" (CLAUDE, 04.10.2026, fix 2)

    Returns (new_text, info_dict).
    """
    header, blocks = _split_blocks(text)
    info = {"n_jetcontrol": 0, "gains_in": [], "gains_out": [], "dropped": [], "nspan": []}
    out_blocks = []

    # CLAUDE (04.10.2026, fix 2): scale factor that keeps the total main-wing strips <= max_wing_strips
    f_cap = 1.0
    if D is not None and strips_per_D is not None and max_wing_strips is not None:
        n_req = 0
        for blk in blocks:
            if blk["kind"] == "SURFACE" and blk["name"].startswith("main_wing"):
                ys = _section_y_values(blk["lines"])
                width = (max(ys) - min(ys)) if ys else 0.0
                rate = max(strips_per_D / D, n_span_min / width) if width > 0 else 0.0
                n_req += 2 * sum(_interval_strips(ys, rate))  # x2: YDUPLICATE  # CLAUDE (04.10.2026, fix 3): per interval
        f_cap = min(1.0, max_wing_strips / n_req) if n_req else 1.0
    info["strip_cap_factor"] = round(f_cap, 4)

    for blk in blocks:
        is_wing = blk["kind"] == "SURFACE" and blk["name"].startswith("main_wing")
        if wing_only and not is_wing:
            info["dropped"].append(f'{blk["kind"]}:{blk["name"]}')
            continue
        lines = list(blk["lines"])
        if is_wing:
            # 1) JETCONTROL gain (data line follows 'JETCONTROL' and an optional comment line)
            if fix_gain:
                for i, ln in enumerate(lines):
                    if ln.strip() == "JETCONTROL":
                        k = i + 1
                        while k < len(lines) and lines[k].strip().startswith(("#", "!")):
                            k += 1
                        m = re.match(r"^(\s*\S+\s+)(" + _NUM + r")(.*)$", lines[k])
                        if m:
                            g_in = float(m.group(2))
                            g_out = g_in - 1.0
                            lines[k] = f"{m.group(1)}{g_out:.6f}{m.group(3)}"
                            info["n_jetcontrol"] += 1
                            info["gains_in"].append(round(g_in, 6))
                            info["gains_out"].append(round(g_out, 6))
            # 2) flat camber line
            if symmetric_airfoil:
                new_lines, skip = [], False
                for ln in lines:
                    if skip:
                        skip = False
                        continue
                    if ln.strip() == "AFILE":
                        skip = True
                        continue
                    new_lines.append(ln)
                lines = new_lines
            # 3) spanwise resolution
            # CLAUDE (04.10.2026, fix 3): Nspanwise/Sspace per SECTION interval instead of one value per
            # SURFACE (surface header keeps only Nchord, Cspace, so JVL reads them from the SECTION lines)
            if D is not None and strips_per_D is not None:
                ys = _section_y_values(lines)
                width = (max(ys) - min(ys)) if ys else 0.0
                rate = max(strips_per_D / D, n_span_min / width) if width > 0 else 0.0
                n_int = _interval_strips(ys, rate * f_cap)
                for i, ln in enumerate(lines):
                    if ln.strip().startswith("#Nchordwise"):
                        old = re.findall(_NUM, lines[i + 1])
                        nc = int(nchord) if nchord is not None else int(float(old[0]))
                        lines[i] = "#Nchordwise  Cspace"
                        lines[i + 1] = f"{nc}         1.0"
                        break
                j = 0
                for i, ln in enumerate(lines):
                    if ln.strip() == "SECTION" and j <= len(n_int):  # CLAUDE (04.10.2026, fix 4): incl. last SECTION
                        k = i + 1
                        while lines[k].strip().startswith(("#", "!")) or not lines[k].strip():
                            k += 1
                        lines[k] = lines[k].rstrip() + f"   {n_int[min(j, len(n_int) - 1)]}   {sspace}"  # Xle Yle Zle Chord Ainc Nspanwise Sspace (last SECTION: value unused, but present)
                        j += 1
                info["nspan"].append((blk["name"], round(width, 4), sum(n_int), min(n_int)))
        out_blocks.append(lines)

    new_lines = list(header)
    for lines in out_blocks:
        new_lines.extend(lines)
    return "\n".join(new_lines) + "\n", info


def set_dvjet_in_run_text(text, dvjet_cases):
    """Set the DVjet constraint value of the k-th run case to dvjet_cases[k]."""
    pattern = re.compile(r"(DVjet\s+->\s+DVjet\s+=\s+)(" + _NUM + r")")
    matches = list(pattern.finditer(text))
    if len(matches) != len(dvjet_cases):
        raise ValueError(
            f"Found {len(matches)} DVjet constraints in the .run file but {len(dvjet_cases)} values were given"
        )
    out, last = [], 0
    for m, val in zip(matches, dvjet_cases):
        out.append(text[last:m.start(2)])
        out.append(f"{float(val):.5f}")
        last = m.end(2)
    out.append(text[last:])
    return "".join(out)


_STAB_KEYS = (
    "alpha", "beta", "Mach",
    "CXtot", "CYtot", "CZtot", "Cltot", "Cmtot", "Cntot",
    "CLtot", "CLff", "CDtot", "CDff", "e",
    "CLcir", "CLffc", "CLjet", "CLffj",
    "CDind", "CDffi", "e_vec", "CDjet", "CDffj", "CDvis", "CDffv",
    "DCQtot", "DCJtot", "DCEtot", "CTtot", "CPtot",
    "Sref", "Cref", "Bref",
)


def parse_stability_file(path):
    """Parse the 'Total Forces' part of a JVL 'st' output file into a dict.

    Only the first occurrence of each key is kept; 'NaN' values are returned as nan.
    """
    out = {}
    with open(path, "r") as f:
        for ln in f:
            if "derivatives" in ln.lower():  # stop at the derivative tables
                break
            for key, val in re.findall(r"([A-Za-z_]\w*)\s*=\s*(" + _NUM + r"|NaN|nan)", ln):
                if key in _STAB_KEYS and key not in out:
                    try:
                        out[key] = float(val.replace("D", "E").replace("d", "e"))
                    except ValueError:
                        out[key] = float("nan")
    return out


def parse_surface_forces_file(path):
    """Parse the first table of a JVL 'fn' output (forces referred to Sref, Cref, Bref).

    Returns a list of dicts with keys n, Area, CL, CD, Cm, CY, Cn, Cl, CDi, CDv, name.
    """
    rows = []
    with open(path, "r") as f:
        lines = f.readlines()
    i0 = None
    for i, ln in enumerate(lines):
        if re.match(r"^\s*n\s+Area\s+CL\s+CD", ln):
            i0 = i + 1
            break
    if i0 is None:
        return rows
    for ln in lines[i0:]:
        if not ln.strip():
            break
        toks = ln.split()
        try:
            vals = [float(t) for t in toks[:10]]
        except ValueError:
            break
        rows.append({
            "n": int(vals[0]), "Area": vals[1], "CL": vals[2], "CD": vals[3], "Cm": vals[4],
            "CY": vals[5], "Cn": vals[6], "Cl": vals[7], "CDi": vals[8], "CDv": vals[9],
            "name": " ".join(toks[10:]),
        })
    return rows
