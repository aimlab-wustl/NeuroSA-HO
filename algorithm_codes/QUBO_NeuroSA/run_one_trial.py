#!/usr/bin/env python3
import argparse, math, csv
from pathlib import Path
import numpy as np
from decimal import Decimal

FOLDER_MAP = {'uf20':'Q_matrix_uf20','uf50':'Q_matrix_uf50','uf75':'Q_matrix_uf75','uf100':'Q_matrix_uf100','uf175':'Q_matrix_uf175','uf250':'Q_matrix_uf250'}
CNF_FOLDER_MAP = {
    'uf20': 'uf20-91',
    'uf50': 'uf50-218',
    'uf75': 'uf75-325',
    'uf100': 'uf100-430',
    'uf175': 'uf175-753',
    'uf250': 'uf250-1065'
}
C_MAP      = {'uf20':91,'uf50':218,'uf75':325,'uf100':430,'uf175':753,'uf250':1065}
ANNEAL_COEFF_MAP = {'uf20':88,'uf50':96,'uf75':96,'uf100':96,'uf175':104,'uf250':104}


def delta_dir_name_from_token(tok: str, max_frac: int = 5) -> str:
    """
    Build 'delta_<mantissa>e<±NN>' from the original CLI token 'tok' without rounding.
    - One digit before '.', up to 'max_frac' after (no padding zeros).
    - Exponent always has sign and at least two digits (e.g., +00, -02, +03).
    """
    d = Decimal(tok)
    if d.is_zero():
        return "delta_0e+00"

    neg = d.is_signed()
    a = d.copy_abs()
    e = a.adjusted()
    m = a.scaleb(-e)

    s = format(m.normalize(), 'f')
    if '.' in s:
        intp, frac = s.split('.', 1)
        frac = frac[:max_frac].rstrip('0')
        s = intp + ('.' + frac if frac else '')

    sign = '-' if neg else ''
    exp = f"{e:+03d}" if -99 <= e <= 99 else f"{e:+d}"
    return f"delta_{sign}{s}e{exp}"


def load_npy_any(path: str) -> np.ndarray:
    """
    Load a .npy file whether it was saved normally or with pickled/object payloads.
    - Tries allow_pickle=False first.
    - Falls back to allow_pickle=True only if necessary.
    - If it's a 0-D object array containing an ndarray, unwraps it.
    - If it's an object array of a single element, unwraps that too.
    """
    try:
        arr = np.load(path, allow_pickle=False)
        return arr
    except ValueError as e:
        msg = str(e).lower()
        if ("pickle" not in msg) and ("object array" not in msg) and ("object arrays" not in msg):
            raise

    arr = np.load(path, allow_pickle=True)

    if isinstance(arr, np.ndarray) and arr.dtype == object:
        if arr.shape == () or arr.size == 1:
            obj = arr.reshape(()).item()
            if isinstance(obj, np.ndarray):
                return obj
            try:
                return np.array(obj)
            except Exception:
                pass

        flat = list(arr.flat)
        if flat and all(isinstance(x, np.ndarray) for x in flat):
            shapes = {x.shape for x in flat}
            if len(shapes) == 1:
                return np.stack(flat).squeeze()

        raise TypeError(f"{path}: object-array payload; expected a single ndarray, got shape {arr.shape}")

    return arr


def parse_dimacs(filename):
    clauses = []
    num_vars = None
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('c'):
                continue
            if line in ['%', '0']:
                continue
            if line.startswith('p'):
                parts = line.split()
                if len(parts) >= 4:
                    num_vars = int(parts[2])
                else:
                    raise ValueError("Invalid problem line: " + line)
            else:
                tokens = line.split()
                clause = [int(tok) for tok in tokens if int(tok) != 0]
                if clause:
                    clauses.append(clause)
    if num_vars is None:
        raise ValueError("Number of variables not specified in file.")
    return num_vars, clauses


def run_single_trial(Q, Q_diag, C, num_vars, clauses_array, anneal_coeff, thld_delta, MAX_ITER):
    pos = clauses_array > 0
    indices = np.abs(clauses_array) - 1

    (N, M) = Q.shape

    beta_end = 2.0
    thld = 1.0
    thld_max = anneal_coeff * beta_end

    X_init = np.random.randint(0, 2, N)
    s_p = X_init
    s_n = (1 - X_init)
    vmem = np.matmul(Q, (s_p - s_n))

    sat = np.where(pos, s_p[indices] == 1, s_p[indices] == 0)
    clause_sat = np.any(sat, axis=1)
    satisfied_count = np.sum(clause_sat)
    best_cost = satisfied_count
    best_iter = 0
    reached_opt = 0
    iter_count = 1

    while iter_count < MAX_ITER:
        p = np.random.randint(0, N)
        s_p_p = s_p[p].item()
        s_n_p = s_n[p].item()
        vmem_p = vmem[p].item()

        noisethld = anneal_coeff * np.log(np.random.random() + 1e-6) / (thld_max * np.log(1 + thld / thld_max))
        temp = -((s_p_p - s_n_p) * (vmem_p + Q_diag[p]))
        spike = 1 if 0.5 * noisethld + temp < 0 else 0

        if spike:
            new_ds_p_p = s_n_p > 0
            new_ds_n_p = s_p_p > 0
            s_p[p] = 0 if s_p_p else 1
            s_n[p] = 0 if s_n_p else 1

            sat = np.where(pos, s_p[indices] == 1, s_p[indices] == 0)
            clause_sat = np.any(sat, axis=1)
            satisfied_count = np.sum(clause_sat)

            vmem += 2 * (new_ds_p_p - new_ds_n_p) * Q[p, :]

            if satisfied_count > best_cost:
                best_cost = satisfied_count
                best_iter = iter_count
                if best_cost == C and not reached_opt:
                    reached_opt = 1
                    return reached_opt, int(best_cost), int(best_iter), MAX_ITER

        thld += thld_delta
        iter_count += 1

    return reached_opt, int(best_cost), np.nan, MAX_ITER


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Run a single trial using repo-relative inputs/outputs.")
    ap.add_argument("--prefix", required=True, choices=FOLDER_MAP.keys())
    ap.add_argument("--instance", type=int, required=True)
    ap.add_argument("--thld-delta", type=str, required=True)
    ap.add_argument("--trial-id", type=int, default=0)
    ap.add_argument("--max-iter", type=int, default=None,
                    help="If not provided, computed from (1 + (anneal_coeff*beta_end - 1)/thld_delta).")
    args = ap.parse_args()

    base_dir = Path(__file__).resolve().parent

    prefix = args.prefix
    file_num = "0" + str(args.instance)

    Q_path = base_dir / "Q_matrices" / FOLDER_MAP[prefix] / f"Q_matrix_{prefix}_{file_num}.npy"
    cnf_filename = base_dir / "SATLIB_DATA" / CNF_FOLDER_MAP[prefix] / f"{prefix}-{file_num}.cnf"

    Q = load_npy_any(str(Q_path))
    Q = -Q
    Q_diag = np.diag(Q).copy()
    np.fill_diagonal(Q, 0)

    num_vars, clauses_array = parse_dimacs(str(cnf_filename))
    clauses_array = np.array(clauses_array)

    C = C_MAP[prefix]
    anneal_coeff = ANNEAL_COEFF_MAP[prefix]

    thld_delta_token = args.thld_delta
    thld_delta_num = float(Decimal(thld_delta_token))

    if args.max_iter is None:
        beta_end = 2.0
        MAX_ITER = int(math.ceil(1 + (anneal_coeff * beta_end - 1) / thld_delta_num))
    else:
        MAX_ITER = int(args.max_iter)

    reached_opt, best_cost, best_iter, MAX_ITER = \
        run_single_trial(Q, Q_diag, C, num_vars, clauses_array, anneal_coeff, thld_delta_num, MAX_ITER)

    out_dir = (
        base_dir / "results" / prefix / f"inst_{args.instance:02d}" /
        delta_dir_name_from_token(args.thld_delta, max_frac=5)
    )
    out_dir.mkdir(parents=True, exist_ok=True)

    out_csv = out_dir / f"trial_{args.trial_id:05d}.csv"
    with out_csv.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["prefix","instance","thld_delta","trial_id","C","MAX_ITER",
                    "best_cost","reached_opt","best_iter"])
        w.writerow([prefix,args.instance,thld_delta_num,args.trial_id,C,MAX_ITER,
                    best_cost,reached_opt,
                    ("" if isinstance(best_iter,float) and np.isnan(best_iter) else int(best_iter))])

    print(f"Wrote: {out_csv}")
