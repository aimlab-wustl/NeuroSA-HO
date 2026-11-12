This folder contains the scripts and RTL used to generate synthesizable Verilog solver modules for the benchmark instances used in the paper.

---

## Quick summary

* **Generator**: `generate_verilog_scripts.py` — reads H-matrix CSV files and colormap files and generates Verilog solver modules.
* **Inputs**: `H_matrices/` (H_matrix_<instance>.csv) and `Colormaps/` (colormap_<instance>.txt).
* **Outputs**: `Verilog_files/solver_modules/<instance>/solver_<instance>.v` and `solver_readout_<instance>.v`.
* **Shared RTL**: `Verilog_files/other_modules/adder_tree.v` and `Verilog_files/other_modules/top_solver_axi.v`.

---

## How to run

From the repository root:

```
python generate_verilog_scripts.py
```

* No command-line arguments are required.
* Input/output paths and the problem-size-to-`C_VAL` mapping are defined at the top of the script and can be modified if needed.

---

## Required inputs (per instance)

1. **H-matrix CSV file**

   `./H_matrices/<subfolder>/H_matrix_<instance>.csv`

2. **Colormap file**

   `./Colormaps/colormap_<instance>.txt`

---

## Output (generated Verilog)

For each problem instance, the script generates:

```
Verilog_files/solver_modules/<instance>/
    solver_<instance>.v
    solver_readout_<instance>.v
```

* `solver_<instance>.v` contains the full solver datapath and control logic derived from the H-matrix and colormap.
* `solver_readout_<instance>.v` implements solution readout and interfaces to higher-level modules.

---

## Additional RTL modules

The following RTL modules are required for integration and synthesis:

* `adder_tree.v`: shared arithmetic module used by all solver instances.
* `top_solver_axi.v`: top-level hardware wrapper providing AXI4-Stream interfaces and connectivity to the processing system.

---

## Noise scaling: software rule vs. RTL comparator

Software (reference) rule: a flip is accepted when ΔE < noise. With
`q_cal = 2a − colsum` (the local field; `a = Σ H·T_bin`), a flip changes the energy by
ΔE = 2·q_cal, so the rule is

    q_cal < 0.5 · w_int16,     w_int16 = trunc(v)

where `v = −A·ln(U + ε) / (C·ln(1 + t/C))`.

RTL: `generate_verilog_scripts.py` emits `(a <<< 1) < noisethld_adj`, with
`noisethld_adj = w + colsum`. This is `2a − colsum < w`, i.e. `q_cal < w`.
(The Verilog's `q_cal` is the raw adder-tree sum `a`, not `2a − colsum`.)

Consequence: streaming the same integer `w_int16` to the RTL gives an effective
temperature twice as high as the software rule.

Required host-side conversion: to make the RTL implement ΔE < noise, the host streams

    w_host = 0.5 · w_int16 = ceil(trunc(v) / 2) = floor(v / 2 + 1/2)

(round-half-up). A plain right shift, `trunc(v/2)`, is not equivalent: for odd
`trunc(v)` it rejects equal-energy moves (`q_cal = 0`) that the software rule accepts.
