import os
import re
import glob
import numpy as np
import pandas as pd

# -------------------------
# utils: build the dataset
# -------------------------
def build_dyn_dataset(root_dir):
    """
    Scan root_dir for files named
      suave_dynamic_stability_outputs_<x>_<y>.txt
    Read each file, annotate with x (int), file_y (int) and row (int),
    and return a single DataFrame with a MultiIndex:
    (x, sigma_fcs, span_loc, fcs_loc, wing_frac, nacelle_frac, row)
    """
    pattern = os.path.join(root_dir, "suave_dynamic_stability_outputs_*.txt")
    files = sorted(glob.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No matching files in {root_dir}")

    rows = []
    filename_re = re.compile(r"suave_dynamic_stability_outputs_(\d+)_(\d+)\.txt$")
    for fn in files:
        m = filename_re.search(os.path.basename(fn))
        if not m:
            continue
        x = int(m.group(1))
        y = int(m.group(2))
        # read as whitespace-delimited table
        df = pd.read_csv(fn, sep=r'\s+', engine='python', header=0)
        # annotate meta columns
        df = df.copy()
        df['x'] = x
        df['file_y'] = y
        df['row'] = np.arange(len(df))  # row number within this file
        rows.append(df)

    df_all = pd.concat(rows, ignore_index=True, sort=False)

    # ensure numeric columns have numeric dtype
    numeric_cols = [c for c in df_all.columns if c not in ('x','file_y','row')]
    df_all[numeric_cols] = df_all[numeric_cols].apply(pd.to_numeric, errors='coerce')

    # create the MultiIndex requested:
    key_cols = ['x', 'sigma_fcs', 'span_loc', 'fcs_loc', 'wing_frac', 'nacelle_frac', 'row']
    df_all.set_index(key_cols, inplace=True)

    return df_all


def _expand_key(x, first5, row):
    """
    Replace None with slice(None) for MultiIndex lookup.
    """
    xk = slice(None) if x is None else x
    f5 = tuple(slice(None) if v is None else v for v in first5)
    rk = slice(None) if row is None else row
    return (xk, *f5, rk)


# -------------------------
# small selection helpers
# -------------------------
RE_COLS = [f"Re{i}" for i in range(1,9)]
IM_COLS = [f"Im{i}" for i in range(1,9)]

def get_eigen_parts(df_all, x, first5, row, part='Re', idx=None, exact=True, tol=1e-9):
    assert part in ('Re','Im'), "part must be 'Re' or 'Im'"
    cols = RE_COLS if part == 'Re' else IM_COLS

    key = _expand_key(x, first5, row)

    # ---- wildcard path (fast, native MultiIndex slicing) ----
    if x is None or None in first5 or row is None:
        out = df_all.loc[key, cols]
        vals = out.to_numpy(dtype=float)
        if idx is None:
            return vals
        else:
            assert 1 <= int(idx) <= 8
            return vals[:, int(idx)-1]

    # ---- original exact / fuzzy logic (unchanged) ----
    try:
        ser = df_all.loc[key]
    except KeyError:
        if exact:
            raise

        x_mask = (df_all.index.get_level_values('x') == x)
        sigma = df_all.index.get_level_values('sigma_fcs').astype(float)
        span  = df_all.index.get_level_values('span_loc').astype(float)
        fcs   = df_all.index.get_level_values('fcs_loc').astype(float)
        wing  = df_all.index.get_level_values('wing_frac').astype(float)
        nac   = df_all.index.get_level_values('nacelle_frac').astype(float)
        row_vals = df_all.index.get_level_values('row').astype(int)

        cond = (
            x_mask &
            (np.abs(sigma - first5[0]) <= tol) &
            (np.abs(span  - first5[1]) <= tol) &
            (np.abs(fcs   - first5[2]) <= tol) &
            (np.abs(wing  - first5[3]) <= tol) &
            (np.abs(nac   - first5[4]) <= tol) &
            (row_vals == row)
        )
        matches = df_all[cond]
        if len(matches) == 0:
            raise KeyError("No matching row found (even fuzzy).")
        ser = matches.iloc[0]

    vals = ser[cols].to_numpy(dtype=float)
    return vals if idx is None else float(vals[int(idx)-1])

    
    
def get_flight_condition(df_all, x, first5, row, exact=True, tol=1e-9):
    key = _expand_key(x, first5, row)

    # ---- wildcard path ----
    if x is None or None in first5 or row is None:
        cols = ['AoA','Mach','Beta','h','n','W']
        out = df_all.loc[key, cols]
        return out.to_numpy(dtype=float)

    # ---- original exact / fuzzy logic ----
    try:
        ser = df_all.loc[key]
    except KeyError:
        if exact:
            raise

        x_mask = (df_all.index.get_level_values('x') == x)
        sigma = df_all.index.get_level_values('sigma_fcs').astype(float)
        span  = df_all.index.get_level_values('span_loc').astype(float)
        fcs   = df_all.index.get_level_values('fcs_loc').astype(float)
        wing  = df_all.index.get_level_values('wing_frac').astype(float)
        nac   = df_all.index.get_level_values('nacelle_frac').astype(float)
        row_vals = df_all.index.get_level_values('row').astype(int)

        cond = (
            x_mask &
            (np.abs(sigma - first5[0]) <= tol) &
            (np.abs(span  - first5[1]) <= tol) &
            (np.abs(fcs   - first5[2]) <= tol) &
            (np.abs(wing  - first5[3]) <= tol) &
            (np.abs(nac   - first5[4]) <= tol) &
            (row_vals == row)
        )
        matches = df_all[cond]
        if len(matches) == 0:
            raise KeyError("No matching row found (even fuzzy).")
        ser = matches.iloc[0]

    return {
        'AoA':  float(ser['AoA']),
        'Mach': float(ser['Mach']),
        'Beta': float(ser['Beta']),
        'h':    float(ser['h']),
        'n':    float(ser['n']),
        'W':    float(ser['W']),
    }


def get_mass_parameters(df_all, x, first5, row, exact=True, tol=1e-9):
    """
    Retrieve mass magnitude/location parameters:
    (sigma_fcs, span_loc, fcs_loc, wing_frac, nacelle_frac)

    Supports wildcards (None) for x, first5 entries, and row.
    Returns an (N,5) numpy array.
    """

    # build boolean mask explicitly (robust to partial slicing)
    mask = np.ones(len(df_all), dtype=bool)

    if x is not None:
        mask &= (df_all.index.get_level_values('x') == x)

    names = ['sigma_fcs','span_loc','fcs_loc','wing_frac','nacelle_frac']
    for name, val in zip(names, first5):
        if val is not None:
            vals = df_all.index.get_level_values(name).astype(float)
            mask &= (np.abs(vals - val) <= tol)

    if row is not None:
        mask &= (df_all.index.get_level_values('row') == row)

    sel = df_all[mask]
    if len(sel) == 0:
        raise KeyError("No matching mass-parameter rows found.")

    idx = sel.index

    return np.column_stack([
        idx.get_level_values('sigma_fcs').astype(float),
        idx.get_level_values('span_loc').astype(float),
        idx.get_level_values('fcs_loc').astype(float),
        idx.get_level_values('wing_frac').astype(float),
        idx.get_level_values('nacelle_frac').astype(float),
    ])


# -------------------------
# Example usage
# -------------------------
if __name__ == "__main__":
    root_dir = r"C:/Users/nmb48/"   # change to folder containing your .txt files
    df_all = build_dyn_dataset(root_dir)

    first5 = (4000.0, None, None, None, None)
    x = None
    row = None
    
    Re = get_eigen_parts(df_all, x, first5, row, part='Re')
    Im = get_eigen_parts(df_all, x, first5, row, part='Im')
    fc = get_flight_condition(df_all, x, first5, row)
    mp = get_mass_parameters(
        df_all,
        x=4,
        first5=(4000.0, None, -1.0, 1.0, 0.2),
        row=2
    )


    evs = get_eigen_parts(df_all, x=3, first5=(4000.0, None, None, None, None), row=None)
    
    evs = get_eigen_parts(df_all, 4, (4000,0.99,-1,1,0.2), row=2)
    
    print('Mass parameters:', mp)
    print("Flight condition:", fc)
    print("Re eigenvalues:", Re)
    print("Im eigenvalues:", Im)

