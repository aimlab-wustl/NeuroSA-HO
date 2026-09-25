# NeuroSA-HO: Higher-Order Neuromorphic Ising Machines

This repository contains data, software, and hardware code for the paper:

**Higher-order neuromorphic Ising machines—autoencoders and Fowler-Nordheim annealers are all you need for scalability**  
Faiek Ahsan *et al.*, *Nature Communications*, 2026.

The code implements higher-order neuromorphic Ising-machine solvers based on the autoencoder formulation described in the paper. The main objective has the form

```math
\min_{\mathbf{s}\in\{\pm 1\}^{N}} E(\mathbf{s})
=
-\sum_{k=1}^{M} J_k \prod_{i=1}^{N} s_i^{H_{k,i}}
```

or equivalently, the solver maximizes

```math
F(\mathbf{s})
=
\sum_{k=1}^{M} J_k \prod_{i=1}^{N} s_i^{H_{k,i}}
```

Here:

- `N` is the number of spin variables.
- `M` is the number of higher-order interaction terms / monomials / clauses.
- `H` is an `M x N` sparse interconnection matrix.
- `H[k, i] = 1` means spin `i` participates in interaction term `k`.
- `J[k]` is the coefficient / weight of term `k`.
- The product term `T[k] = prod_i s_i^H[k,i]` is the clause or monomial output.

The key idea is that the solver works in the sparse clause-output space rather than explicitly recomputing higher-order products from scratch at every step. Fowler-Nordheim annealing is used to generate the time-varying noisy thresholds used by the latent spin neurons.

---

## Repository structure

```text
NeuroSA-HO/
├── algorithm_codes/
│   ├── N_HOIM_single_sipin_flip/
│   ├── N_HOIM_with_GC/
│   └── QUBO_NeuroSA/
├── hardware_codes/
├── TTS_Data_and_Figure_generation_codes/
└── README.md
```

### `algorithm_codes/N_HOIM_single_sipin_flip`

CPU simulation of the higher-order Ising machine without graph coloring. At each iteration, all candidate latent-neuron events are tested, but a global arbiter selects only one accepted spin flip. This preserves an asynchronous single-spin-flip update.

Use this folder when you want the simplest higher-order solver implementation.

### `algorithm_codes/N_HOIM_with_GC`

CPU simulation of the graph-colored higher-order Ising machine. Variables with the same color are conditionally independent, so they can be updated in parallel. This is the preferred implementation when a valid graph coloring is available.

Use this folder when you want the higher-throughput version of the higher-order solver.

### `algorithm_codes/QUBO_NeuroSA`

CPU simulation of the second-order NeuroSA/QUBO solver used for the quadratized comparison experiments.

Use this folder when comparing the direct higher-order formulation against a second-order quadratized formulation.

### `hardware_codes`

Scripts and RTL for generating FPGA solver modules from `H` matrices and color maps. This folder was used for the FPGA experiments reported in the paper.

### `TTS_Data_and_Figure_generation_codes`

Recorded CPU/FPGA data and MATLAB scripts used to generate the time-to-solution and iteration-to-solution figures in the paper.

---

## Software requirements

The CPU simulation scripts require Python with NumPy and SciPy.

A minimal setup is:

```bash
python -m pip install numpy scipy
```

The figure-generation scripts require MATLAB.

The FPGA code requires an FPGA/RTL tool flow compatible with the provided Verilog modules. The hardware implementation in the paper was developed around the RFSoC 4x2 / Zynq UltraScale+ workflow.

---

## Running the existing benchmark instances

The main CPU runner in each algorithm folder is:

```text
run_one_trial.py
```

Each run performs one randomized trial and writes one CSV result file under the local `results/` directory.

### Higher-order solver without graph coloring

```bash
cd algorithm_codes/N_HOIM_single_sipin_flip

python run_one_trial.py --prefix uf50 --instance 1 --thld-delta 0.00333333
```

To set the maximum number of iterations manually:

```bash
python run_one_trial.py --prefix uf50 --instance 1 --thld-delta 0.00333333 --max-iter 5000
```

To label the output file with a custom trial ID:

```bash
python run_one_trial.py --prefix uf50 --instance 1 --thld-delta 0.00333333 --trial-id 7
```

### Higher-order solver with graph coloring

```bash
cd algorithm_codes/N_HOIM_with_GC

python run_one_trial.py --prefix uf50 --instance 1 --thld-delta 0.00333333
```

This version also requires a color-map file for the selected instance.

### Second-order QUBO / NeuroSA solver

```bash
cd algorithm_codes/QUBO_NeuroSA

python run_one_trial.py --prefix uf50 --instance 1 --thld-delta 0.00333333
```

This version requires the corresponding quadratized `Q` matrix and SATLIB CNF file.

---

## Command-line arguments

Required arguments:

```text
--prefix
```

Benchmark family identifier. In the current scripts this must be one of:

```text
uf20, uf50, uf75, uf100, uf175, uf250
```

These names are hard-coded in `FOLDER_MAP`, `C_MAP`, and `ANNEAL_COEFF_MAP` inside each `run_one_trial.py`.

```text
--instance
```

Integer instance ID. For example:

```bash
--instance 1
```

The scripts internally form the file number as:

```python
file_num = "0" + str(instance)
```

So instance `1` is read as `01`.

```text
--thld-delta
```

Sampling period of the discretized Fowler-Nordheim annealing schedule. This controls the annealing speed.

Examples:

```bash
--thld-delta 0.00333333
--thld-delta 1e-3
--thld-delta 2e-4
```

Smaller `thld_delta` gives slower annealing and usually higher reliability, but it increases runtime. Larger `thld_delta` gives faster annealing, but can reduce success probability.

Optional arguments:

```text
--trial-id
```

Integer used only for naming the output CSV file. It does not change the algorithm.

```text
--max-iter
```

Maximum number of iterations. If this is not provided, the script computes

```text
MAX_ITER = ceil(1 + (anneal_coeff * beta_end - 1) / thld_delta)
```

with

```text
beta_end = 2.0
```

---

## Output format

Each run writes one CSV file:

```text
results/<prefix>/inst_<instance>/delta_<delta-token>/trial_<trial-id>.csv
```

The CSV columns are:

```text
prefix, instance, thld_delta, trial_id, C, MAX_ITER, best_cost, reached_opt, best_iter
```

where:

- `best_cost` is the best score reached during the run.
- `reached_opt = 1` means the run reached the target score `C`.
- `best_iter` is the iteration where the best score was first reached.
- `C` is the target score used by the benchmark script.

For MAX-3SAT benchmark instances, `C` is the number of clauses, so `best_cost == C` means all clauses were satisfied.

---

## Running multiple trials

Each algorithm folder also contains:

```text
example_trial_run.py
```

Open this file and edit the variables near the top, such as:

```python
PREFIX = "uf50"
INSTANCE = 1
THLD_DELTA = "0.00333333"
MAX_ITER = None
TRIAL_IDS = range(10)
```

Then run:

```bash
python example_trial_run.py
```

This repeatedly calls `run_one_trial.py` and creates one CSV output per trial.

---

# Using the code for a new higher-order Ising problem

The current scripts are benchmark runners. To solve a new problem, the main task is to convert the problem into the sparse higher-order Ising form

```math
F(\mathbf{s})
=
\sum_{k=1}^{M} J_k \prod_{i=1}^{N} s_i^{H_{k,i}}
```

Then you provide the corresponding `H` matrix, target score, and annealing parameters.

---

## 1. Build the higher-order term list

Represent the objective as a list of monomials. Each monomial should contain:

```text
variables_in_term, coefficient
```

For example, a third-order term involving variables `0`, `2`, and `5` with coefficient `+1` is represented conceptually as:

```text
([0, 2, 5], +1)
```

Important conventions:

- Variable indices are zero-based in the code.
- Drop additive constants; they do not affect the optimizer.
- Combine duplicate monomials when possible.
- If a coefficient becomes zero after combining terms, remove that term.
- For a term involving variables `[i, j, k]`, the monomial is `s[i] * s[j] * s[k]`.
- A first-order term is allowed.
- A second-order term is allowed.
- A higher-order term of any order is allowed, as long as it can be represented as a sparse row in `H`.

---

## 2. Save the `H` matrix

The higher-order scripts expect a SciPy sparse `.npz` matrix.

Create a sparse matrix `H` with shape:

```text
number_of_terms x number_of_variables
```

Each row corresponds to one higher-order Ising term. For a term `k`, place nonzero entries in the columns corresponding to the variables that participate in that term.

In the current benchmark scripts, the nonzero entries of each row store the coefficient `J[k]`. Therefore, for row `k`, every participating variable column should contain the same coefficient `J[k]`.

For example, if term `k` is:

```text
J[k] * s[0] * s[2] * s[5]
```

then row `k` of `H` should have the value `J[k]` in columns `0`, `2`, and `5`, and zeros elsewhere.

Save the matrix using SciPy's sparse `.npz` format. For a custom problem named `custom`, use a filename such as:

```text
H_matrix_custom_01.npz
```

For the non-graph-colored solver, place the file in a folder such as:

```text
algorithm_codes/N_HOIM_single_sipin_flip/H_matrices/H_matrix_custom/H_matrix_custom_01.npz
```

For the graph-colored solver, place the file in:

```text
algorithm_codes/N_HOIM_with_GC/H_matrices/H_matrix_custom/H_matrix_custom_01.npz
```

Then update the corresponding `FOLDER_MAP` in `run_one_trial.py`:

```python
FOLDER_MAP = {
    "custom": "H_matrix_custom",
}
```

You may also add `"custom"` to the existing dictionary instead of replacing the benchmark entries.

---

## 3. Set the target score `C`

The current scripts use `C_MAP[prefix]` as the target score for early stopping.

For MAX-SAT benchmarks, `C` is the number of clauses.

For a new problem, set:

```python
C_MAP = {
    "custom": target_score,
}
```

where `target_score` is the score at which you want the script to stop.

If the exact optimum is unknown, you can either:

1. set `C` to the best-known target score, or
2. disable early stopping and simply run until `MAX_ITER`.

For a general weighted Ising objective, the current benchmark scoring lines may need to be changed. The present MAX-SAT code computes a SAT-style score using

```python
best_cost = (np.sum(T * W) + 7 * C) // 8
```

This formula is specific to the MAX-3SAT polynomial used in the paper. For a generic higher-order Ising objective, replace this score with

```python
best_cost = np.sum(T * W)
```

and update the later `sat_clause` / `best_cost` logic consistently.

---

## 4. Add graph coloring if using `N_HOIM_with_GC`

The graph-colored solver requires a color map file.

Two variables must not have the same color if they appear together in any higher-order term. Equivalently, build a conflict graph where variables are connected if they share at least one monomial, then color that graph.

The color map can be generated during preprocessing using DSATUR, a fast greedy graph-coloring heuristic introduced by Brélaz (1979). Optimal coloring is not required for this algorithm. The purpose of coloring is to identify conditionally independent latent neurons that can be updated in parallel. As discussed in the paper, computing an optimal coloring is NP-hard, and a fast greedy heuristic is sufficient during preprocessing.

The color map file contains one integer color per variable, one per line. Colors are one-indexed.

Example for six variables:

```text
1
2
1
3
2
3
```

For the current path convention, save it as:

```text
algorithm_codes/N_HOIM_with_GC/Colormap/Colormap_custom/Colormap_custom_01.txt
```

Then update `FOLDER_MAP`, `C_MAP`, and `ANNEAL_COEFF_MAP` in `algorithm_codes/N_HOIM_with_GC/run_one_trial.py`.

It is not necessary for the coloring to use the minimum possible number of colors, although fewer colors usually means more variables can be updated in parallel per color cycle.

---

# Mapping common problems to `H`

## General higher-order Ising / PUBO

Given

```math
F(\mathbf{s})
=
\sum_k J_k \prod_{i\in S_k} s_i
```

create one row of `H` for each term `k`.

- Put nonzero entries in the columns corresponding to variables in `S_k`.
- Store `J_k` as the nonzero value in that row.
- Drop constant terms.

---

## MAX-CUT

For an edge `(i, j)` with weight `w_ij`, the MAX-CUT objective can be written, up to an additive constant, using a second-order Ising term.

One convenient convention for this code is to add one row with variables `[i, j]` and coefficient

```text
J = -w_ij / 2
```

because maximizing the cut corresponds to favoring `s_i s_j = -1`.

For unweighted MAX-CUT, use `J = -0.5` for every edge.

For a generic MAX-CUT run, you should also replace the MAX-3SAT-specific score computation with the appropriate cut-value computation or with the generic Ising objective score.

---

## 3R-3X / planted XOR-SAT

For each XOR equation involving variables `(i, j, k)` with right-hand side bit `b`, use the spin form

```math
s_i s_j s_k = (-1)^b
```

Add one row with variables `[i, j, k]` and coefficient

```text
J = (-1)**b
```

Maximizing the sum of these weighted products maximizes the number of satisfied XOR equations.

---

## MAX-SAT

For a clause with `p` literals, the higher-order Ising expansion used in the paper is

```math
\Phi_k =
\sum_{r=1}^{p} (-1)^{r-1}
\sum_{1 \le i_1 < \cdots < i_r \le p}
s_{k,i_1}s_{k,i_2}\cdots s_{k,i_r}
```

The literal-to-spin convention is:

```text
positive literal x_i      ->  +s_i
negative literal not x_i  ->  -s_i
```

For MAX-3SAT, each clause produces seven monomial terms before collecting duplicates:

```text
+s1 +s2 +s3
-s1*s2 -s1*s3 -s2*s3
+s1*s2*s3
```

After accounting for negated literals, add each resulting monomial as a row in `H`. If multiple clauses produce the same monomial, their coefficients can be summed into a single row.

For the MAX-3SAT benchmark scripts in this repository, the SAT score is recovered from the Ising polynomial using

```text
number_satisfied = (sum_k J_k T_k + 7*C) / 8
```

where `C` is the number of clauses.

---

# Choosing the Fowler-Nordheim annealing parameters

The main annealing parameters in the CPU scripts are:

```python
anneal_coeff
thld_delta
MAX_ITER
beta_end = 2.0
```

The code sets:

```python
thld = 1.0
thld_max = anneal_coeff * beta_end
MAX_ITER = ceil(1 + (anneal_coeff * beta_end - 1) / thld_delta)
```

Important naming note: in the benchmark code, `C_MAP` is the target score, while `thld_max` is the FN annealing scale. These are different quantities.

---

## `anneal_coeff`

`anneal_coeff` sets the energy scale of the FN threshold.

A practical rule is:

```text
Choose anneal_coeff slightly larger than the largest local flip-energy scale observed in preliminary runs.
```

For a general higher-order Ising objective, inspect the local field

```math
\sum_k J_k H_{k,i} T_k
```

or the corresponding spin-flip energy scale over random initial states and early trial runs, then choose `anneal_coeff` slightly above the largest typical value.

If `anneal_coeff` is too small, the threshold scale is too low and the system may become too greedy too early.

If `anneal_coeff` is too large, the system can remain too noisy for too long.

Values used in the paper/code for MAX-3SAT are:

```python
ANNEAL_COEFF_MAP = {
    "uf20": 88,
    "uf50": 96,
    "uf75": 96,
    "uf100": 96,
    "uf175": 104,
    "uf250": 104,
}
```

Other values used in the paper include:

```text
3R3X: A = 10
5R5X: A = 14
7R7X: A = 14
```

For MAX-CUT G-set experiments, the paper used graph-dependent values such as:

```text
G4:  A = 60
G11: A = 40
G15: A = 60
G22: A = 60
```

For a new problem, do not assume these values are automatically optimal. Use them only as scale references.

---

## `thld_delta`

`thld_delta` is the discrete sampling period of the FN annealer.

Smaller values mean slower annealing:

```text
smaller thld_delta -> slower schedule -> higher reliability -> longer runtime
```

Larger values mean faster annealing:

```text
larger thld_delta -> faster schedule -> lower runtime -> possibly lower success probability
```

For larger and harder problem instances, the annealing schedule should often be made deliberately slower by decreasing `thld_delta`. This gives the dynamics more time to approach high-quality or state-of-the-art solutions instead of freezing too quickly into suboptimal states. In practice, this slower schedule is especially important when the goal is not merely to get a reasonable solution, but to reliably converge to the best-known or SOTA-quality solution.

For reproducibility or high-success runs, start with a slow schedule such as:

```text
2e-3, 1e-3, 5e-4, or 2e-4
```

For quick debugging, a faster schedule can be used, but the solution quality should not be judged from debugging runs.

A useful starting grid is problem-dependent, but the following values are reasonable for initial sweeps:

```text
1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 3.333e-3, 5e-3, 7e-3, 1e-2
```

Then refine around the best-performing value.

---

## `MAX_ITER`

If `MAX_ITER` is not provided, the script computes it from `anneal_coeff`, `beta_end`, and `thld_delta`.

For quick debugging, manually set a small value:

```bash
--max-iter 1000
```

For real optimization runs, either omit `--max-iter` or set a sufficiently large value.

If the solver is not reaching the target score:

1. decrease `thld_delta`,
2. increase `MAX_ITER`,
3. retune `anneal_coeff`,
4. run more independent trials.

---

## TTS / ITS parameter selection

For time-to-solution (TTS) or iteration-to-solution (ITS) calculations, do not choose `thld_delta` from a single run. Instead, sweep a grid of `thld_delta` values and compute the success probability and TTS/ITS for each schedule.

The grid search is intended to select the annealing speed automatically:

```text
fast schedules    -> lower runtime per trial, but lower success probability
slow schedules    -> higher runtime per trial, but higher success probability
optimal TTS/ITS   -> best trade-off between runtime and success probability
```

Thus, when reporting TTS/ITS, run multiple independent trials for each `thld_delta`, estimate the success probability for each schedule, and report the best TTS/ITS obtained over the grid.

This is particularly important for larger problem sizes. As the problem becomes harder, the best TTS/ITS point often shifts toward slower annealing schedules, because the solver needs more gradual annealing to reliably reach SOTA-quality solutions.

---

# Getting the final spin assignment

The current benchmark scripts mainly save the best score and hit time. They do not save the best spin assignment by default.

For a new optimization problem, you may want to save the best spin vector. The simplest approach is to maintain and store the spin state during the run.

For `N_HOIM_single_sipin_flip`, after a selected spin flip, add:

```python
spins[chosen_index] *= -1
```

and whenever the best score improves, store:

```python
best_spins = spins.copy()
```

For `N_HOIM_with_GC`, after computing `flip_vars`, add:

```python
spins[flip_vars] *= -1
```

and similarly store:

```python
best_spins = spins.copy()
```

Then save `best_spins` at the end of the run using, for example:

```python
np.save(out_dir / f"best_spins_{args.trial_id:05d}.npy", best_spins)
```

For general problems, maintaining the spin vector directly is the safest readout method.

---

# Hardware code

The hardware generator is in:

```text
hardware_codes/generate_verilog_scripts.py
```

It reads:

```text
hardware_codes/H_matrices/
hardware_codes/Colormaps/
```

and generates instance-specific Verilog modules under:

```text
hardware_codes/Verilog_files/solver_modules/
```

The generated files include:

```text
solver_<instance>.v
solver_readout_<instance>.v
```

Shared RTL modules include:

```text
adder_tree.v
top_solver_axi.v
```

To run the generator from the `hardware_codes` folder:

```bash
cd hardware_codes
python generate_verilog_scripts.py
```

No command-line arguments are required. Input/output paths and the problem-size-to-parameter mapping are defined inside the script.

The hardware code is mainly intended to reproduce the FPGA-style solver generation used in the paper. For new problems, first verify the CPU `N_HOIM_with_GC` implementation, then generate RTL using the corresponding `H` matrix and color map.

---

# TTS/ITS data and figure generation

The folder

```text
TTS_Data_and_Figure_generation_codes/
```

contains recorded data and MATLAB scripts used to generate time-to-solution and iteration-to-solution figures.

The main recorded-data CSV files are:

```text
coarse_summary.csv
coarse_instance_tts.csv
```

`coarse_summary.csv` contains summary-level TTS/ITS information, including the problem prefix, threshold delta, selected runtime, and median TTS.

`coarse_instance_tts.csv` contains instance-level success probabilities and TTS estimates.

---

# Practical workflow for a new problem

A recommended workflow is:

1. Convert the problem into the higher-order Ising form  
   `F(s) = sum_k J_k prod_i s_i^H[k,i]`.

2. Build and save the sparse `H` matrix as `.npz`.

3. If using graph coloring, build a conflict graph from `H` and generate a valid greedy coloring, for example using DSATUR.

4. Add a new prefix to:
   - `FOLDER_MAP`
   - `C_MAP`
   - `ANNEAL_COEFF_MAP`

5. Replace the MAX-3SAT-specific scoring rule if the new problem is not MAX-3SAT.

6. Run a few short debugging trials with a small `MAX_ITER`.

7. Tune `anneal_coeff` based on the local flip-energy scale.

8. Run a grid over `thld_delta`.

9. For larger problems or SOTA-quality targets, deliberately include slower schedules in the `thld_delta` grid.

10. Run many independent trials and estimate success probability / TTS.

11. Save the best spin assignment if the application requires the actual solution, not only the objective value.

---

# References

## Paper

If you use this repository, please cite:

```bibtex
@Article{Ahsan2026,
author={Ahsan, Faiek
and Maiti, Saptarshi
and Chen, Zihao
and Kaiser, Jakob
and Nandi, Ankita
and Srivatsav, Madhuvanthi
and Schemmel, Johannes
and Andreou, Andreas G.
and Eshraghian, Jason
and Thakur, Chetan Singh
and Chakrabartty, Shantanu},
title={Higher-order neuromorphic Ising machines---autoencoders and Fowler-Nordheim annealers are all you need for scalability},
journal={Nature Communications},
year={2026},
month={Apr},
day={16},
volume={17},
number={1},
pages={5293},
issn={2041-1723},
doi={10.1038/s41467-026-71937-4},
url={https://doi.org/10.1038/s41467-026-71937-4}
}
```

## License

Except where otherwise noted, original material in this repository is licensed under the [Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0) License](https://creativecommons.org/licenses/by-nc/4.0/). See the [LICENSE](LICENSE) file for the full license text.

You are free to share and adapt the licensed material for non-commercial purposes, provided that appropriate credit is given and any changes are indicated. Third-party materials remain subject to their respective licenses.
