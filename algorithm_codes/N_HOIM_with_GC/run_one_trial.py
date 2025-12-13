#!/usr/bin/env python3
import argparse, math, csv
from pathlib import Path
import numpy as np
from scipy.sparse import load_npz, csr_matrix, issparse
from decimal import Decimal

FOLDER_MAP = {'uf20':'H_matrix_uf20','uf50':'H_matrix_uf50','uf75':'H_matrix_uf75','uf100':'H_matrix_uf100','uf175':'H_matrix_uf175','uf250':'H_matrix_uf250'}
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
    e = a.adjusted()             # exponent for scientific notation
    m = a.scaleb(-e)             # mantissa in [1,10)

    s = format(m.normalize(), 'f')
    if '.' in s:
        intp, frac = s.split('.', 1)
        frac = frac[:max_frac].rstrip('0')
        s = intp + ('.' + frac if frac else '')

    sign = '-' if neg else ''
    exp = f"{e:+03d}" if -99 <= e <= 99 else f"{e:+d}"
    return f"delta_{sign}{s}e{exp}"


def load_sparse_any(path: str):
    """Load a SciPy sparse matrix from a .npz written by save_npz or from a pickled/legacy .npz."""
    try:
        return load_npz(path).tocsr()
    except Exception as e:
        msg = str(e)
        if "pickled data" not in msg and "Invalid NPY" not in msg:
            raise

    with np.load(path, allow_pickle=True) as z:
        files = set(z.files)

        if {"data", "indices", "indptr", "shape"} <= files:
            shape_arr = z["shape"]
            if hasattr(shape_arr, "tolist"):
                shape = tuple(int(x) for x in np.array(shape_arr).ravel().tolist())
            else:
                shape = tuple(int(x) for x in shape_arr)
            H = csr_matrix((z["data"], z["indices"], z["indptr"]), shape=shape)

        else:
            H = None
            for k in z.files:
                if z[k].dtype == object:
                    cand = z[k].item()
                    if issparse(cand):
                        H = cand
                        break
            if H is None:
                raise TypeError(f"{path}: object inside npz is not a SciPy sparse matrix")

    return H.tocsr()


def read_colormap(filename):
    with open(filename, "r") as f:
        return [int(line.strip()) for line in f if line.strip()]


def form_groups(colormap):
    uniq = sorted(set(colormap))
    groups=[[] for _ in range(len(uniq))]
    for i,c in enumerate(colormap): groups[c-1].append(i)
    return [np.array(g,dtype=int) for g in groups]


def calculate_clause_products(H, variables_pm1):
    out = np.empty(H.shape[0], dtype=int)
    for i in range(H.shape[0]):
        idx = H[i].indices
        out[i] = np.prod(variables_pm1[idx])
    return out


def run_single_trial(H_csr, H_csc, groups, C, anneal_coeff, thld_delta, MAX_ITER):
    num_clauses, num_variables = H_csr.shape
    column_sums = np.asarray(H_csr.sum(axis=0)).ravel()
    W = np.array([row.data[0] if row.nnz>0 else 0 for row in H_csr])

    required_colour = len(groups)
    group_masks = [np.zeros(num_variables, bool) for _ in groups]
    for m,g in zip(group_masks,groups): m[g]=True

    row_to_var = np.full((required_colour, num_clauses), -1, int)
    var_colour = np.empty(num_variables, int)
    for g,m in enumerate(group_masks): var_colour[m]=g
    for j in range(num_variables):
        g = var_colour[j]
        for i in H_csc.indices[H_csc.indptr[j]:H_csc.indptr[j+1]]:
            row_to_var[g, i] = j
    group_clause_idxs=[]; group_var_idxs=[]
    for g in range(required_colour):
        ci = np.where(row_to_var[g]>=0)[0]
        group_clause_idxs.append(ci)
        group_var_idxs.append(row_to_var[g, ci])

    beta_end = 2.0
    thld = 1.0
    thld_max = anneal_coeff * beta_end

    X_init = np.random.randint(0,2,size=num_variables)
    spins  = 2*X_init - 1
    T      = calculate_clause_products(H_csr, spins)
    T_bin  = (T==1)

    best_cost  = (np.sum(T*W) + 7*C)//8
    sat_clause = best_cost
    best_iter  = 0

    q = np.zeros(num_variables, bool)
    q_cal = np.empty(num_variables, int)
    sigma_bool = np.zeros(num_clauses, bool)

    iter_count = 1
    reached_opt = 0

    while iter_count < MAX_ITER:
        for g in range(required_colour):
            cur = groups[g]
            noisethld = (-anneal_coeff*np.log(np.random.rand(len(cur))+1e-6) /
                     (thld_max*np.log(1+thld/thld_max)))

            q_cal = H_csc.T.tocsr().dot(T_bin); q_cal *= 2; q_cal -= column_sums
            q_current = q_cal[cur] < 0.5*noisethld

            flip_vars = cur[q_current]
            q[cur]=False
            if flip_vars.size: q[flip_vars]=True
            sum_q = int(q_cal[flip_vars].sum()) if flip_vars.size else 0

            sigma_bool.fill(False)
            ci = group_clause_idxs[g]; vj = group_var_idxs[g]
            if vj.size: sigma_bool[ci] = q[vj]
            T_bin ^= sigma_bool

            sat_clause -= sum_q//4
            if sat_clause > best_cost:
                best_cost = sat_clause; best_iter = iter_count
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

    H_path  = base_dir / "H_matrices" / FOLDER_MAP[prefix] / f"H_matrix_{prefix}_{file_num}.npz"
    cm_path = base_dir / "Colormap" / f"Colormap_{prefix}" / f"Colormap_{prefix}_{file_num}.txt"

    H = load_sparse_any(str(H_path))
    H_csc = H.tocsc()
    groups = form_groups(read_colormap(str(cm_path)))

    C = C_MAP[prefix]
    anneal_coeff = ANNEAL_COEFF_MAP[prefix]

    thld_delta_token = args.thld_delta
    thld_delta_num   = float(Decimal(thld_delta_token))

    if args.max_iter is None:
        beta_end = 2.0
        MAX_ITER = int(math.ceil(1 + (anneal_coeff * beta_end - 1) / thld_delta_num))
    else:
        MAX_ITER = int(args.max_iter)

    reached_opt, best_cost, best_iter, MAX_ITER = \
        run_single_trial(H, H_csc, groups, C, anneal_coeff, thld_delta_num, MAX_ITER)

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
