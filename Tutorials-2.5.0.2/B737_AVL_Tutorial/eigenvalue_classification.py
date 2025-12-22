# python code to parse the uploaded eigenvalue file, associate 8 eigenvalues to 5 aircraft modes,
# track modes across parameter combinations, and produce visual + quantitative checks.
# This code is self-contained. It will:
# 1. Read '/mnt/data/suave_dynamic_stability_outputs_combined.txt'
# 2. Parse groups (parameter combinations) separated by blank lines.
# 3. Expect 6 rows per group (six flight conditions). For each row: extract Re1..Re8, Im1..Im8 -> eigenvalues.
# 4. For each flight-condition index (0..5) track modes across groups using a heuristic classifier + Hungarian matching.
# 5. Produce plots: eigenvalue trajectories in complex plane (colored by mode), frequency & damping vs parameter index.
# 6. Print summary and ambiguous-case flags.
#
# If scipy is not available, the code falls back to a simple greedy assignment.
# Author: Assistant (adapted to user's data).

import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import math, sys, os

filepath = r"C:/Users/nmb48/suave_dynamic_stability_outputs_combined.txt"
if not os.path.exists(filepath):
    raise FileNotFoundError(f"Expected file at {filepath} not found.")

# Try to import Hungarian assignment; fallback to greedy if unavailable
use_scipy = True
try:
    from scipy.optimize import linear_sum_assignment
except Exception as e:
    print("scipy.linear_sum_assignment not available, will use greedy matching fallback.", file=sys.stderr)
    use_scipy = False

def parse_groups(path):
    text = open(path, "r", errors="ignore").read()
    # Split into groups separated by blank lines (two or more newlines)
    raw_groups = [g.strip() for g in text.split("\n\n") if g.strip()!='']
    groups = []
    for g in raw_groups:
        lines = [ln.strip() for ln in g.splitlines() if ln.strip()!='']
        # Find header line (contains 'Re1' or 'sigma_fcs'); collect numeric lines after header
        header_idx = None
        for i,ln in enumerate(lines):
            if 'Re1' in ln or 'sigma_fcs' in ln:
                header_idx = i
                break
        # Numeric lines are those that start with digit or '-' or '.'
        numeric_lines = []
        for ln in lines:
            if len(ln)==0: continue
            first = ln[0]
            if first.isdigit() or first in "-.":
                numeric_lines.append(ln)
            else:
                # attempt to detect numeric line by splitting and checking first token
                toks = ln.split()
                if len(toks)>0 and (toks[0][0].isdigit() or toks[0][0] in "-."):
                    numeric_lines.append(ln)
        # If numeric_lines empty, try all lines except header
        if not numeric_lines and header_idx is not None:
            numeric_lines = lines[header_idx+1:]
        # Convert numeric lines to lists of floats (tokens separated by whitespace)
        parsed = []
        for ln in numeric_lines:
            toks = ln.split()
            try:
                vals = [float(t) for t in toks]
                parsed.append(vals)
            except:
                # skip unparsable lines
                continue
        if parsed:
            groups.append((lines, parsed))
    return groups

groups = parse_groups(filepath)
n_groups = len(groups)
if n_groups == 0:
    raise RuntimeError("No numeric groups were parsed - file format unexpected.")

print(f"Parsed {n_groups} parameter-combination groups.")

# Expecting 6 rows per group (user indicated 6 flight conditions across each parameter combination)
rows_per_group = [len(parsed) for (_,parsed) in groups]
unique_counts = sorted(set(rows_per_group))
print(f"Rows per group counts observed (unique): {unique_counts}")
# Choose the most common as expected rows per group
from collections import Counter
cnt = Counter(rows_per_group)
expected_rows = cnt.most_common(1)[0][0]
print(f"Using expected rows per group = {expected_rows} (most common)")

# Build structured dataset: groups_data[g][r] -> dict with params and eigenvalues (list of 8 complex)
groups_data = []
for idx,(raw_lines, parsed) in enumerate(groups):
    if len(parsed) < expected_rows:
        # skip groups with fewer rows (if any)
        continue
    # inspect header (if included) to find column indices
    header_tokens = None
    for ln in raw_lines:
        if 'Re1' in ln and 'Im1' in ln:
            header_tokens = ln.split()
            break
    # If no header, assume Re1..Re8 start after the 11th number as in example: W is 11th index
    # We'll attempt to find Re1/Im1 by matching token counts:
    sample = parsed[0]
    ncols = len(sample)
    # Heuristic: last 16 columns are Re1..Re8 Im1..Im8
    if ncols >= 16:
        re_start = ncols - 16  # index of Re1
        im_start = re_start + 8
    else:
        raise RuntimeError("Unexpected number of columns in numeric line; can't locate Re/Im columns.")
    rows = []
    for vals in parsed[:expected_rows]:
        Re = vals[re_start:re_start+8]
        Im = vals[im_start:im_start+8]
        eigs = np.array([complex(Re[i], Im[i]) for i in range(8)])
        # Also capture a few identifying params (first 11 columns as in header)
        params = vals[:min(11, len(vals))]
        rows.append({"params": params, "eigs": eigs})
    groups_data.append(rows)

n_groups = len(groups_data)
print(f"Using {n_groups} clean groups (each with {expected_rows} rows). Total flight-condition sequences = {expected_rows}.")

# Helper functions for pairing and classification
def pair_eigenvalues(eigs, imag_tol=1e-6):
    """
    Given 8 eigenvalues, pair conjugate complex eigenvalues and leave real ones.
    Returns:
      complex_reps: list of representative eigenvalues for conjugate pairs (positive-imag chosen)
      real_vals: list of real eigenvalues (imag ~ 0)
      pairs_full: list of tuples (pos_im, neg_im) for complex pairs
    """
    eigs = np.array(eigs)
    reals = []
    complex_pos = []
    complex_neg = []
    for lam in eigs:
        if abs(lam.imag) <= imag_tol:
            reals.append(np.real(lam))
        else:
            if lam.imag > 0:
                complex_pos.append(lam)
            else:
                complex_neg.append(lam)
    # Pair each positive imag with the closest negative imag by distance
    pairs = []
    used_neg = set()
    for p in complex_pos:
        if not complex_neg:
            break
        dists = [abs(p - n) for n in complex_neg]
        j = int(np.argmin(dists))
        pairs.append((p, complex_neg[j]))
        used_neg.add(j)
    # if any negative left unpaired, try to pair remaining negatives to positives (unlikely)
    remaining_neg = [n for i,n in enumerate(complex_neg) if i not in used_neg]
    for n in remaining_neg:
        if complex_pos:
            dists = [abs(n - p) for p in complex_pos]
            j = int(np.argmin(dists))
            pairs.append((complex_pos[j], n))
    # Complex reps choose the positive-imag eigenvalue as representative
    complex_reps = [p for p,n in pairs if p.imag>0]
    # Convert real list to floats
    real_vals = [float(r) for r in reals]
    # If due to numeric issues we end up with wrong counts, attempt fallback pairing by magnitude sorting
    if len(complex_reps) + len(real_vals) != 5:
        # fallback: sort eigs by imaginary absolute descending and make first 3 as complex pairs and remaining 2 as reals
        sorted_by_im = sorted(eigs, key=lambda x: abs(x.imag), reverse=True)
        complex_reps = []
        real_vals = []
        taken = set()
        for lam in sorted_by_im:
            if len(complex_reps) < 3 and abs(lam.imag) > imag_tol:
                if lam.imag > 0:
                    complex_reps.append(lam)
                elif lam.imag < 0:
                    # will be ignored as representative
                    continue
            else:
                if abs(lam.imag) <= imag_tol:
                    real_vals.append(float(np.real(lam)))
        # reduce lists to expected sizes
        complex_reps = complex_reps[:3]
        real_vals = real_vals[:2]
    return complex_reps, real_vals, pairs

def compute_mode_descriptors(rep):
    """Compute descriptors for a representative eigenvalue (complex for oscillatory, real for non-oscillatory)."""
    if isinstance(rep, complex) or np.iscomplexobj(rep):
        sigma = rep.real
        wd = abs(rep.imag)
        wn = math.hypot(sigma, wd)
        zeta = -sigma / wn if wn>0 else np.nan
        freq_hz = wd / (2*math.pi)
        return {"lambda": rep, "sigma": sigma, "wd": wd, "wn": wn, "zeta": zeta, "freq_hz": freq_hz, "osc": True}
    else:
        sigma = float(rep)
        t_half = math.log(2)/abs(sigma) if sigma!=0 else np.inf
        return {"lambda": sigma, "sigma": sigma, "wd": 0.0, "wn": abs(sigma), "zeta": np.nan, "freq_hz": 0.0, "t_half": t_half, "osc": False}

# Classification per group-row: assign the 5 modes using heuristic:
def classify_modes_from_eigs(eigs):
    complex_reps, real_vals, pairs = pair_eigenvalues(eigs)
    # descriptors for complex reps: keep positive-imag representatives
    complex_desc = [compute_mode_descriptors(rep) for rep in complex_reps]
    # sort complexes by damped frequency (wd) descending
    complex_desc_sorted = sorted(complex_desc, key=lambda d: d["wd"], reverse=True)
    # Assign ordering: highest wd -> short-period, next -> dutch-roll, lowest -> phugoid
    mode_names = []
    mode_map = {}
    if len(complex_desc_sorted) >= 3:
        labels = ["short-period", "dutch-roll", "phugoid"]
        for lab,d in zip(labels, complex_desc_sorted[:3]):
            mode_map[lab] = d
            mode_names.append(lab)
    else:
        # fallback: whatever we have
        labels = ["short-period", "dutch-roll", "phugoid"]
        for i,d in enumerate(complex_desc_sorted):
            lab = labels[i]
            mode_map[lab] = d
            mode_names.append(lab)
    # For the reals: more negative sigma -> roll (fast), less negative -> spiral (slow)
    real_sorted = sorted(real_vals)
    # real_sorted ascending (most negative first)
    if len(real_sorted) >= 2:
        roll = real_sorted[0]
        spiral = real_sorted[1]
        mode_map["roll"] = compute_mode_descriptors(roll)
        mode_map["spiral"] = compute_mode_descriptors(spiral)
        mode_names += ["roll","spiral"]
    else:
        # fallback: assign any leftover
        if len(real_sorted)==1:
            mode_map["roll"] = compute_mode_descriptors(real_sorted[0])
            mode_map["spiral"] = compute_mode_descriptors(real_sorted[0])
            mode_names += ["roll","spiral"]
    return mode_map

# Build per-flight-condition sequences across groups
n_conditions = expected_rows
sequences = [ [] for _ in range(n_conditions) ]  # sequences[c][g] -> dict of modes for group g at condition c
for gidx,rows in enumerate(groups_data):
    for cidx in range(n_conditions):
        row = rows[cidx]
        eigs = row["eigs"]
        mode_map = classify_modes_from_eigs(eigs)
        sequences[cidx].append(mode_map)

# Now do tracking across groups (parameter index) for each sequence (flight condition)
# We'll produce tracked trajectories for 5 mode labels for each condition.
tracked = []  # list per condition: dict label->list of descriptors across groups
ambiguous_flags = defaultdict(list)  # condition -> list of group indices flagged ambiguous

for cidx,seq in enumerate(sequences):
    # seq is list of mode_maps (length n_groups)
    # initialize tracking with the first group's labels and representatives
    labels = list(seq[0].keys())
    # Ensure consistent label ordering: prefer ['roll','short-period','dutch-roll','phugoid','spiral']
    preferred_order = ["roll","short-period","dutch-roll","phugoid","spiral"]
    # Create current_order: intersection preserving preferred order
    current_order = [lab for lab in preferred_order if lab in labels] + [lab for lab in labels if lab not in preferred_order]
    # tracked_data: dict label -> list of descriptors per group
    tracked_data = {lab: [seq[0].get(lab, None)] for lab in current_order}
    # track across subsequent groups using assignment
    for g in range(1, n_groups):
        prev_labels = list(tracked_data.keys())
        prev_descs = [tracked_data[lab][-1] for lab in prev_labels]
        next_map = seq[g]
        next_labels = list(next_map.keys())
        next_descs = [next_map[lab] for lab in next_labels]
        # Build cost matrix between prev_labels and next_labels based on complex-plane distance of representative eigenvalues
        m = len(prev_labels); n = len(next_labels)
        cost = np.zeros((m,n))
        for i, pd in enumerate(prev_descs):
            for j, nd in enumerate(next_descs):
                if pd is None or nd is None:
                    cost[i,j] = 1e6
                else:
                    # Use complex distance normalized by magnitude
                    lam_p = pd["lambda"]
                    lam_n = nd["lambda"]
                    dist = abs(lam_p - lam_n)
                    norm = max(1.0, abs(lam_p), abs(lam_n))
                    # also include damping difference term (zeta), scaled
                    zeta_p = pd.get("zeta", 0.0) if pd.get("zeta", None) is not None else 0.0
                    zeta_n = nd.get("zeta", 0.0) if nd.get("zeta", None) is not None else 0.0
                    dz = abs((zeta_p or 0.0) - (zeta_n or 0.0))
                    cost[i,j] = dist / norm + 0.5 * dz
        # Solve assignment
        assignment = []
        if use_scipy:
            try:
                row_ind, col_ind = linear_sum_assignment(cost)
                assignment = list(zip(row_ind, col_ind))
            except Exception as e:
                use_scipy = False
        if (not use_scipy) or (len(assignment)==0):
            # greedy fallback
            assignment = []
            cost_cp = cost.copy()
            m,n = cost_cp.shape
            assigned_rows = set(); assigned_cols = set()
            while True:
                if cost_cp.size==0:
                    break
                i,j = np.unravel_index(np.argmin(cost_cp), cost_cp.shape)
                if cost_cp[i,j] > 1e5:
                    break
                assignment.append((i,j))
                # invalidate row i and col j
                cost_cp[i,:] = 1e6
                cost_cp[:,j] = 1e6
                assigned_rows.add(i); assigned_cols.add(j)
                if len(assigned_rows) >= m or len(assigned_cols) >= n:
                    break
        # Build new tracked_data by reordering next labels to match prev_labels according to assignment
        new_tracked = {}
        assigned_next = set()
        for i,j in assignment:
            prev_lab = prev_labels[i]
            next_lab = next_labels[j]
            new_tracked[prev_lab] = tracked_data[prev_lab] + [ next_map[next_lab] ]
            assigned_next.add(next_lab)
        # For any prev_label unassigned, append None
        for pl in prev_labels:
            if pl not in new_tracked:
                new_tracked[pl] = tracked_data[pl] + [ None ]
        # For any next_label unassigned (new emergent), create a series with None for previous groups
        for nl in next_labels:
            if nl not in assigned_next:
                new_tracked[nl] = [ None ]*(g) + [ next_map[nl] ]
        tracked_data = new_tracked
        # detect ambiguous cases: if any two representatives are very close in the next_map
        # compute pairwise distances
        reps = [d["lambda"] for d in next_descs]
        for i in range(len(reps)):
            for j in range(i+1, len(reps)):
                if abs(reps[i] - reps[j]) / max(1.0, abs(reps[i]), abs(reps[j])) < 0.02:
                    ambiguous_flags[cidx].append((g, i, j, next_labels[i], next_labels[j], float(abs(reps[i]-reps[j]))))
    tracked.append(tracked_data)

# Visualization: for each flight condition, plot eigenvalue trajectories (complex plane) with one color per tracked label
colors = {
    "roll":"tab:blue", "short-period":"tab:orange", "dutch-roll":"tab:green", "phugoid":"tab:red", "spiral":"tab:purple"
}
for cidx,tracked_data in enumerate(tracked):
    plt.figure(figsize=(8,6))
    ax = plt.gca()
    for lab, series in tracked_data.items():
        # series is list of descriptors length n_groups
        xs = []
        ys = []
        for d in series:
            if d is None:
                xs.append(np.nan); ys.append(np.nan)
            else:
                lam = d["lambda"]
                xs.append(lam.real); ys.append(lam.imag)
        ax.scatter(xs, ys, marker='.', label=lab, color=colors.get(lab,None))
    ax.set_xlabel("Real (rad/s)")
    ax.set_ylabel("Imag (rad/s)")
    ax.set_title(f"Condition {cidx+1}: Eigenvalue trajectories in complex plane across {n_groups} parameter combos")
    ax.grid(True)
    ax.legend()
    plt.tight_layout()
    plt.show()

    # Plot frequency and damping vs parameter index
    plt.figure(figsize=(10,4))
    ax1 = plt.subplot(1,2,1)
    for lab, series in tracked_data.items():
        freqs = [ (d["wd"]/(2*np.pi) if d is not None else np.nan) for d in series ]
        ax1.scatter(range(n_groups), freqs, marker='.', label=lab, color=colors.get(lab,None))
    ax1.set_xlabel("Parameter combo index")
    ax1.set_ylabel("Damped freq (Hz)")
    ax1.set_title("Frequency vs parameter index")
    ax1.grid(True)
    ax1.legend()

    ax2 = plt.subplot(1,2,2)
    for lab, series in tracked_data.items():
        zetas = [ (d["zeta"] if d is not None else np.nan) for d in series ]
        ax2.scatter(range(n_groups), zetas, marker='.', label=lab, color=colors.get(lab,None))
    ax2.set_xlabel("Parameter combo index")
    ax2.set_ylabel("Damping ratio zeta")
    ax2.set_title("Damping ratio vs parameter index")
    ax2.grid(True)
    plt.tight_layout()
    plt.show()

import sys
sys.exit()

# Print a compact summary and ambiguous flags
print("\nSummary per flight condition (first few groups shown):\n")
for cidx,tracked_data in enumerate(tracked[:min(6,len(tracked))]):
    print(f"Condition {cidx+1}:")
    for lab, series in tracked_data.items():
        first = series[0]
        if first is not None:
            lam = first["lambda"]
            print(f"  {lab:12s} repr lambda = {lam.real:+8.4f}{lam.imag:+8.4f}j, zeta={first.get('zeta',np.nan):6.3f}, f_d={first.get('wd',0)/ (2*np.pi):6.3f} Hz")
        else:
            print(f"  {lab:12s} repr = None")
    print("")

# Ambiguity report
total_amb = sum(len(v) for v in ambiguous_flags.values())
print(f"\nAmbiguity flags detected: {total_amb} events across conditions. Showing up to 20:")
count = 0
for cidx, events in ambiguous_flags.items():
    for ev in events:
        g, i, j, lab_i, lab_j, dist = ev
        print(f" Condition {cidx+1}, group idx {g}: close pair {lab_i} vs {lab_j} (abs dist = {dist:.4e})")
        count += 1
        if count>=20:
            break
    if count>=20:
        break

print("\nDone. The figures above show eigenvalue trajectories and mode metrics. \nIf you want, I can (a) save these plots to disk, (b) output a CSV with the tracked modal data, or (c) adjust thresholds/weights used in the matching and ambiguity detection. Which would you like?")
