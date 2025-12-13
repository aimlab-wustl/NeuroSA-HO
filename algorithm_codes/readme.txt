Inside each of these folder, there is the script: run_one_trial.py 
# Single-trial runner (run_one_trial.py)

This script runs ONE randomized trial for a selected SATLIB/UF benchmark prefix + instance, and writes ONE CSV result file under the local results/ folder.

---

## How to run (Command line)

From the folder containing run_one_trial.py:

1. Minimal run (uses default MAX_ITER)
   python run_one_trial.py --prefix uf50 --instance 1 --thld-delta 0.00333333

2. Override MAX_ITER
   python run_one_trial.py --prefix uf50 --instance 1 --thld-delta 0.00333333 --max-iter 5000

3. Set a custom trial id (only affects output filename)
   python run_one_trial.py --prefix uf50 --instance 1 --thld-delta 0.00333333 --trial-id 7

To run multiple trials, repeat the command with different --trial-id values.

---

## How to run (example_trial_run.py)

Open example_trial_run.py, edit the variables at the top (PREFIX, INSTANCE, THLD_DELTA, MAX_ITER, TRIAL_IDS), then run the file.
This calls run_one_trial.py internally once per trial id in TRIAL_IDS, producing multiple CSV outputs without typing terminal arguments.

---

## Arguments

Required arguments:

--prefix
Benchmark family identifier. Must be one of:
uf20 uf50 uf75 uf100 uf175 uf250

--instance
Integer instance id (e.g., 1, 2, 3, ...). This selects which input instance file to load.

--thld-delta
Threshold increment token provided as a STRING. Examples:
0.00333333
1e-3
2.5E-4
The token is parsed using Python Decimal for stable interpretation and is also used to form the output subfolder name.

Optional arguments:

--trial-id
Default: 0
Integer identifier used ONLY to name the output CSV file. It does not change the algorithm.

--max-iter
Default: computed automatically
Maximum number of iterations used in the run.
If not provided, the script computes:
MAX_ITER = ceil( 1 + (anneal_coeff * beta_end - 1) / thld_delta )
where:
beta_end = 2.0
anneal_coeff is determined by the chosen prefix (fixed mapping inside the script)
thld_delta is the numeric value of --thld-delta

---

## Output (CSV)

Each run writes one CSV file containing:

prefix, instance, thld_delta, trial_id, C, MAX_ITER, best_cost, reached_opt, best_iter

Meaning:

* reached_opt = 1 means the run reached the target optimum (best_cost == C) and returned immediately on the first hit.
* best_iter is the iteration index when the best_cost was first achieved; it is blank/NaN if the optimum was never hit.
