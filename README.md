# Best-Worst Method (BWM) in Julia: Euclidean and Linear Models

This repository provides an efficient and robust Julia implementation of the **Best-Worst Method (BWM)** for multi-criteria decision-making (MCDM), solved using `JuMP` and `Ipopt`.

The algorithm computes optimal criteria weights, consistency indices, and total deviations using two distinct formulations:
1. **Euclidean BWM:** A non-linear formulation based on minimizing the sum of squared Euclidean deviations.
2. **Linear Chebyshev BWM:** A linear programming formulation minimizing the maximum absolute deviation.

---

## Installation

Ensure Julia is installed on your system. Install the required optimization dependencies via the Julia REPL:

```julia
using Pkg
Pkg.add(["JuMP", "Ipopt"])
```

---

## Quick Start

```julia

# Best-to-Others (BO) and Others-to-Worst (OW) comparison vectors
BO = [3, 1, 5, 3, 7]
OW = [5, 7, 3, 5, 1]

# Solve both Euclidean and Linear BWM models
BWM(BO, OW)
```

---

## Citation

If you use this repository or any part of these models in your research, please cite the corresponding original papers:

### 1. For the Euclidean BWM Model
If your study utilizes the **Euclidean distance-based model**, please cite:

> Kocak, H., Caglar, A., & Oztas, G. Z. (2018). Euclidean best–worst method and its application. *International Journal of Information Technology & Decision Making*, 17(05), 1587–1605. https://doi.org/10.1142/S021962201850029X

```bibtex
@article{kocak2018euclidean,
  title     = {Euclidean best--worst method and its application},
  author    = {Kocak, Huseyin and Caglar, Ali and Oztas, G. Z.},
  journal   = {International Journal of Information Technology \& Decision Making},
  volume    = {17},
  number    = {05},
  pages     = {1587--1605},
  year      = {2018},
  publisher = {World Scientific},
  doi       = {10.1142/S021962201850029X}
}
```

### 2. For the Linear Chebyshev BWM Model
If your study utilizes the **Linear BWM model**, please cite:

> Rezaei, J. (2016). Best-worst multi-criteria decision-making method: Some properties and a linear model. *Omega*, 64, 126–130. https://doi.org/10.1016/j.omega.2015.12.001

```bibtex
@article{rezaei2016best,
  title     = {Best-worst multi-criteria decision-making method: Some properties and a linear model},
  author    = {Rezaei, Jafar},
  journal   = {Omega},
  volume    = {64},
  pages     = {126--130},
  year      = {2016},
  publisher = {Elsevier},
  doi       = {10.1016/j.omega.2015.12.001}
}
```

---

## License
This project is open-source and available under the [MIT License](LICENSE).
