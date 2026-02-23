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