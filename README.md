# NeuroSA-HO: Data and code for Higher-Order Neuromorphic Ising Machines

This repository contains the **data and code** used in several experiments from the paper:

**“Higher-Order Neuromorphic Ising Machines — Autoencoders and Fowler-Nordheim Annealers are all you need for Scalability.”**

## What’s included

- **CPU simulation code** for two variants of the *neuromorphic higher-order Ising machine* discussed in the paper.
- **NeuroSA solver code** for solving **quadratized SAT** problems.
- A **data** folder containing recorded results from:
  - **FPGA experiments**
  - **CPU experiments**
  including **TTS/ITS** measurements.
- The corresponding **MATLAB scripts** used to generate the **TTS/ITS figures** reported in the paper.

For detailed instructions and folder-specific descriptions, please refer to the `readme.txt` inside each folder.

## Citation

If you use this repository, please cite the preprint:

```bibtex
@misc{ahsan2025higherorderneuromorphicisingmachines,
      title={Higher-Order Neuromorphic Ising Machines -- Autoencoders and Fowler-Nordheim Annealers are all you need for Scalability},
      author={Faiek Ahsan and Saptarshi Maiti and Zihao Chen and Jakob Kaiser and Ankita Nandi and Madhuvanthi Srivatsav and Johannes Schemmel and Andreas G. Andreou and Jason Eshraghian and Chetan Singh Thakur and Shantanu Chakrabartty},
      year={2025},
      eprint={2506.19964},
      archivePrefix={arXiv},
      primaryClass={cs.NE},
      url={https://arxiv.org/abs/2506.19964},
}
