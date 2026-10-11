# `nils/jvl/code_claude` — JVL verification of the blown-wing surrogate

Created by Claude on 04.10.2026. This implements Steps 1–6 of `jvl_blown_wing_verification_plan_041026.md` (saved in the Aircraft Modelling project).

**No existing file in SUAVE, TASOPT.jl-priv or Flydrogen was modified.** Changes to SUAVE behaviour are applied *in memory* by `jvl_patches_claude.apply_patches()`, and only in the process that calls it.

| File | Purpose | Environment |
|---|---|---|
| `jvl_patches_claude.py` | In-memory patches of `SUAVE.Analyses.Stability.AVL` (sample_training, write_geometry, write_run_cases, run_analysis), the `.jvl`/`.run` post-processing and the JVL result parsers | suave (parsers are pure Python) |
| `blown_wing_layout_claude.py` | Parametric propeller layout on the fixed TASOPT wing; returns a geometry DataFrame for the unchanged `wrap_jvl` | any (numpy/pandas) |
| `jvl_calibration_2d_claude.py` | Step 1: 2D-like (AR 100) blown flat plate, calling `jvl.exe` directly | any (numpy/pandas) |
| `jvl_standalone_sweep_exec_claude.py` | Steps 4–6: the sweep over (N, ΔY, Π) with two run cases per geometry (unblown/blown) | suave |

Outputs go to `SUAVE\nils\jvl\data_claude\<ddmmyy>_<op>_<mode>_<model>_<camber>_spd<k>[...]\`, containing `jvl_standalone_sweep_results.csv` and the JVL run folders in `runs\`.

## What the patches change compared with `jvl_sweep_exec.py`

| Issue (plan §3) | Original behaviour | With the patches |
|---|---|---|
| B1 | JETCONTROL gain = `Jgain` = Vjet/V, giving Vjet/V ≈ 2 | gain = `Jgain − 1` |
| B3, B5 | α = 3°, M = 0.15, SL, from `wrap_jvl` | α, M and altitude from the surrogate reference (α = α_w, eq. 29) |
| B4 | only blown case | two run cases at the same α: DVjet = 0 and 1 |
| B7 | `Prop_PR_des = 1.005` in the run-folder name | true Π, plus a unique `run_tag` (N, ΔY, Π index) |
| B11 | dxdisk/dndisk from the TASOPT nacelle position | dxdisk = −ax_pos·c, dndisk = 0 (set by the layout) |
| B12 | B737a camber | flat camber line (option `--cambered` restores B737a) |
| B13 | Nspanwise fixed per surface, cosine spacing | ≥ `--strips-per-D` strips per disk, uniform spacing, total main-wing strips ≤ `--max-wing-strips` (default 360) |
| B14 | full configuration | wing-only (option `--full-config`) |
| B17 | fixed-width parsing | regex parsing; adds CLcir and the main-wing CL from the `fn` file |

### Checked offline against your 26.01.2026 JVL files

- `parse_stability_file` and `parse_surface_forces_file` reproduce all printed values.
- The N = 4, Π = 1.005 output confirms B1: gain 1.0131 → Vjet/V = 2.013, CTtot = 2.357.
- The same output also settles B10. JVL's DCQtot = 0.58549 equals `(r − 1)·hjet·ΣD/Sref` with ΣD = 4 × 4.0457 m exactly. The jets are therefore confined to the intervals between two JETCONTROL sections, and no extra ε-sections are needed.

### Jet strength (modes)

Mode `V` (default) uses `Vjet/V = 1 + 2·a_p` (far wake). For hdisk = πD/4 and fh = 1, JVL's jet thrust is `CT = A_disk·(r² − 1)/Sref = T/(q·S)` exactly. Mode V is therefore also the thrust-matched mode T of the plan, and no JVL thrust constraint is needed. The `CT_jvl` against `CT_target` columns check this per point.

Mode `L` uses the surrogate's leading-edge value `V_w/V`, for sensitivity only.

## How to run (suave environment, from the SUAVE root)

```bat
conda activate suave
cd C:\Users\nmb48\Documents\GitHub\SUAVE
```

**Step 1 — calibration** (≈ 12 JVL runs, needs only `jvl.exe`):

```bat
python nils\jvl\code_claude\jvl_calibration_2d_claude.py
```

- Check that `DCQ_rel_err` and `CT_rel_err` are ≈ 0 for all rows. This validates the gain convention and the jet confinement.
- Read `beta_eff_jvl` against `hdisk_over_c` and compare with the surrogate's `beta` column from Step 0.

**Step 5/6 — smoke test** (3 geometries: N = 8, ΔY = 0.6, Π min/mid/max):

```bat
python -m nils.jvl.code_claude.jvl_standalone_sweep_exec_claude --reference C:\Users\nmb48\Documents\GitHub\Flydrogen\flydrogen\analyses\tests_claude\data_claude\blown_wing_surrogate_reference_takeoff_linear.csv --smoke
```

- Inspect `runs\avl_files_N8_*\vehicle.jvl`: main wing only, gains ≈ 2·a_p, no AFILE.
- Also inspect `claude_jvl_postprocess_info.txt` in the same folder.
- In the results CSV, check that `CT_jvl ≈ CT_target`.

**Step 6 — convergence** (worst cases: N = 16 at the largest Π, N = 4 at the smallest Π). Run the same points with different `--strips-per-D` values and compare `delta_CL_wing`:

```bat
for %k in (6 12 24 48) do python -m nils.jvl.code_claude.jvl_standalone_sweep_exec_claude --reference <...>_takeoff_linear.csv --N 4 16 --dY 0.6 --pi-stride 49 --strips-per-D %k
```

Optional sensitivities on the same points: `--nchord 16`, `--sspace 1.0`, `--cambered`, `--full-config`, `--mode L`.

**Step 5 — full sweep** (4 N × 5 ΔY × 8 Π ≤ 160 geometries, infeasible ones skipped). Resumable: re-running skips points already in the results CSV.

```bat
python -m nils.jvl.code_claude.jvl_standalone_sweep_exec_claude --reference <...>_takeoff_linear.csv
python -m nils.jvl.code_claude.jvl_standalone_sweep_exec_claude --reference <...>_cruise.csv
python -m nils.jvl.code_claude.jvl_standalone_sweep_exec_claude --reference <...>_takeoff.csv
```

The last of these is the slide's take-off point, α_w ≈ 23°, outside VLM validity; report it with that caveat. To parallelise, start several processes with disjoint `--N` subsets and the same other options; they write to the same results CSV.

**Step 7 — figures:** run `jvl_standalone_sweep_postprocess_claude.py` in py311b (see the README in `Flydrogen\flydrogen\analyses\tests_claude`).

## Notes and limitations

- **JVL array limits (fix of 04.10.2026).** The calibration first stopped with `SDUPL: Vortex array overflow. Increase nvmax.`
  - The source copy `jvl2.16/JVL2.16/src/jdimen.inc` has `nvmax = 5000` and `nsmax = 500`. Your compiled `jvl.exe` evidently has larger limits: the 26.01.2026 runs used 9680 vortices and 844 strips. The exact values are unknown.
  - The calibration wing is therefore AR 100 with 40 × 8 vortices per half.
  - In the sweep, `--max-wing-strips` (default 360, i.e. about 7200 wing vortices) scales Nspanwise down when needed. Each run folder's `claude_jvl_postprocess_info.txt` reports `strip_cap_factor`; a value < 1 means the cap was active.
  - The JVL console output is now written to `avl_log.txt` in each run folder, and its last lines are included in any error.
- **Spanwise spacing per section (fix of 04.10.2026).** Closely spaced sections (propeller edges, mid-pitch, planform breaks) made JVL stop with `Cannot adjust spanwise spacing at SECTION ... Insufficient number of spanwise vortices`. The main-wing surfaces now give Nspanwise/Sspace per SECTION interval (at least 1 each) instead of one value per surface. `claude_jvl_postprocess_info.txt` lists (surface, width, strips, smallest interval count).

- The layout band starts at the fuselage side (body half-width 1.385 m). For ΔY = 0.9 it is moved 5.7 cm inboard so that it ends at the tip; the `y0_inside_fuselage_m` column records this.
- In wing-only mode the nacelle bodies are still created by `wrap_jvl`, because SUAVE's jet placement uses their positions, but they are removed from the `.jvl` file before JVL runs.
- `wrap_jvl` still prints `CTtot`/`CLtot` of case 1, the unblown case. The blown values are in the results CSV.
- Re-running with `--no-fix-gain` reproduces the old (B1) behaviour, for comparison only.

## JVL verification against Hawkswell (2024) (05.10.2026)

`jvl_hawkswell_periodic_claude.py` reproduces the spanwise-periodic actuator-disk RANS cases of Hawkswell's thesis, Ch. 6:
- NACA 43018 (the ATR-42 root section), no flap, M = 0.18.
- D/c = 0.5, x_p/c = −0.33, y_p/c = −0.05.
- Tc = 1 and 2, with Tc = T / (q πD²/4).
- s/D from 0 (a continuous jet, i.e. 2D) to 2, at α = 0–12°.

JVL has no periodic boundary condition, so each case is run on two wings, N/2 and N pitches per semispan, and extrapolated to AR → ∞. It re-uses `run_jvl` and `RUN_CASE` from `jvl_calibration_2d_claude.py` (the run-case Mach is set in memory) and `postprocess_jvl_text` from `jvl_patches_claude.py`. The reference data are in `data_claude\hawkswell_2024_digitised_claude.csv`.

```bat
python C:\Users\nmb48\Documents\GitHub\SUAVE\nils\jvl\code_claude\jvl_hawkswell_periodic_claude.py --smoke
python C:\Users\nmb48\Documents\GitHub\SUAVE\nils\jvl\code_claude\jvl_hawkswell_periodic_claude.py
```

The full run takes about 1 h. `--max-wing-strips` defaults to 240, which keeps the vortex count below nvmax = 5000 of the source. Then plot in py311b with `jvl_hawkswell_postprocess_claude.py` (see the README in `Flydrogen\flydrogen\analyses\tests_claude`).

Notes:
- **Near-field lift with jets.** A JVL 2.16 build made from the same source on Linux gives identical Trefftz-plane results (CLff, CDff, CTtot, DCQ) to your `jvl.exe`. Its near-field `CLjet`, however, is 0.073 where yours is 0.000, with a non-zero rolling moment on a symmetric wing. Use the Trefftz-plane `CLff` only.
- **What JVL cannot represent.** JVL's thin-jet model has no jet position, so the effects of x_p and y_p are not reproduced. It also cannot represent airfoil thickness or viscous de-cambering, so compare the increments ΔCL = CL_blown − CL_unblown.

## Distributed propulsion, road-map points 1–3 (10.10.2026)

Run all three in py311b. The verification and Aspect 1 scripts import flydrogen; the layout script is pure numpy/pandas.

| Script | Purpose | Output (`data_claude\<ddmmyy>_…`) |
|---|---|---|
| `dp_verification_figures_claude.py` | V1 parity, V2 constant total thrust (ΔC_L against C_T,w at α = 6°), V3 ratios against α and s/D, for RANS (Hawkswell), JVL and the surrogate (β on a_w, and on the far-wake excess) | `_dp_verification_claude` |
| `dp_layout_claude.py` | N → D, stations, s/D, ΔY, x_p, y_p, feasibility. Uses TASOPT placement and the `const_area` / `max_diameter` scaling | `_dp_aspect1_claude\dp_layout_table_claude.csv` |
| `dp_aspect1_surrogate_claude.py` | A1: ΔC_L against N (take-off, take-off_linear, cruise) with the V1 accuracy band. A2: cruise L/D against N, blowing or folded, including nacelle drag | `_dp_aspect1_claude` |

```bat
python C:\Users\nmb48\Documents\GitHub\SUAVE\nils\jvl\code_claude\dp_verification_figures_claude.py --jvl <PC run>\hawkswell_jvl_results.csv
python C:\Users\nmb48\Documents\GitHub\SUAVE\nils\jvl\code_claude\dp_aspect1_surrogate_claude.py
```

Assumptions to check are listed in each script's docstring:
- ground geometry: h_wing = 3.4 m, bank angle 5°, minimum clearance 0.18 m;
- motor length 1.0 m at 1 MW, scaled with P^(1/3);
- CD0 = 0.0226 and CD_nac,ref = 0.00098;
- all thrust is carried by the N DP propellers.
