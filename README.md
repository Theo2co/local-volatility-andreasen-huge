# Local Volatility Calibration via the Andreasen–Huge Scheme

Calibration of a local-volatility surface on **Euro Stoxx 50 (SX5E)** option quotes using the fully-implicit one-step PDE scheme of **Andreasen & Huge (2010)**, plus extrapolation of the surface to new maturities.

> Submitted as Assignment 3 of *Advanced Derivatives* (Prof. Elena Perazzi, EPFL — Fall 2025).

---

## Overview

Standard Dupire-style local-volatility calibration differentiates the implied-volatility surface twice in strike and once in maturity, which amplifies any market noise. **Andreasen & Huge (2010)** sidestep this by propagating call prices forward one maturity at a time with a **fully-implicit finite-difference step**, and fitting a piecewise-constant local-volatility parameter per maturity layer by least squares. The result is a smooth, arbitrage-friendly surface obtained from a sequence of tiny convex problems instead of a single ill-posed differentiation.

This project implements the full pipeline:

1. Read the quoted IV grid and convert to Black–Scholes call prices on the $(K, T)$ nodes.
2. For each maturity step, calibrate the piecewise-constant parameter $\theta(K)$ by least squares against quoted prices and propagate with the implicit AH step.
3. Invert Black–Scholes at every node to recover the calibrated implied-volatility surface $\hat\sigma(K, T)$.
4. **Extrapolate** the surface to $T = 1$ and $T = 1.5$ by propagating from the nearest quoted expiry with no refit.
5. Visualise the resulting price and IV surfaces and the strike cross-sections.

---

## Method

### Setup

Spot $S_0 = 100$, rates $r = q = 0$. Strikes are absolute, $K_i = S_0 \cdot k_i$ with $k_i$ the relative-strike grid (40% to 200% of spot). Only cells with a quote are usable, so a binary mask `mask_valid = (implied_vols > 0)` is kept throughout to avoid polluting the least-squares objective with missing data.

### From IVs to call prices

At each quoted node, IVs are converted to Black–Scholes prices

$$C^\text{mkt}(K_i, T_j) = C_\text{BS}(S_0, K_i, T_j ; \sigma_{i,j}, r, q),$$

and the $t = 0$ column is set to the intrinsic value $\max(S_0 - K, 0)$.

### The Andreasen–Huge step

Between two consecutive expiries $T_j \to T_{j+1}$, the AH fully-implicit scheme propagates one time step:

$$A_{j+1}(v) \, C_{j+1} = C_j,$$

where $A_{j+1}(v)$ is the full one-step matrix (identity already included) parameterised by a piecewise-constant local-volatility vector $v$, with the partition defined by the midpoints of the quoted strike indices. Letting

$$z_i = \frac{\Delta t}{2\, \Delta k^2}\, \theta_{i+1}^2 \qquad (i = 1, \ldots, N-2),$$

the interior rows of $A_{j+1}(v)$ form the tridiagonal

$$A_{j+1}(v) = \text{tridiag}\big( -z_i,\; 1 + 2z_i,\; -z_i \big),$$

with identity boundary rows. The linear system is solved with a single tridiagonal pass.

### Per-layer calibration

For each maturity step, $v$ is chosen by least squares on the available quoted strikes:

$$\min_v \sum_{i \in \text{quoted}(j)} \left( C^\text{model}_{i,j+1}(v) - C^\text{mkt}_{i,j+1} \right)^2,$$

with initial guess $v^{(0)} = 0.3\, S_0$ (per parameter) and bounds $v \in [0, 2 S_0]$.

### Implied-volatility inversion

Once the price grid is built, Black–Scholes is inverted at each node by Nelder–Mead from a seed of $0.5$,

$$\hat\sigma(K_i, T_j) = \arg\min_{\sigma \geq 0} \big( C_\text{BS}(S_0, K_i, T_j ; \sigma) - C_{i,j} \big)^2,$$

with $|\hat\sigma|$ recorded to guard against tiny negative steps from the unconstrained search.

### Extrapolation

For target maturities $T_\text{new} \in \{1, 1.5\}$, the price surface is propagated from the nearest quoted expiry $T_j < T_\text{new}$ using the parameters already calibrated on that interval. **No refit is done.** BS is then re-inverted to obtain $\hat\sigma(K, T_\text{new})$.

---

## Results

### On the original grid

The interpolated **call-price surface** is monotone in maturity and decreasing in strike, as expected. The **implied-volatility surface** reproduces the pronounced left skew at short maturities and flattens for longer tenors — the SX5E equity-index pattern.

### Extrapolation to T = 1 and T = 1.5

The price surface remains monotone in $T$ (no calendar arbitrage). The extrapolated IVs extend smoothly from the nearest quoted terms and no spurious ripples appear at the junctions.

### Strike cross-sections

The shortest maturity has the largest curvature and the steepest left-skew, with smiles flattening as $T$ grows. The extrapolated smiles at $T = 1$ and $T = 1.5$ sit naturally between their neighbouring quoted terms and preserve the overall term structure.

All figures are reproduced in [`andreasen_huge_report.pdf`](andreasen_huge_report.pdf) and rendered inline in the notebook.

---

## Implementation notes

- **Tridiagonal solve.** Each AH step amounts to one banded linear solve. The implementation embeds the band in a sparse/dense system; for a production-grade implementation, a dedicated Thomas algorithm would be the next step.
- **Per-layer least squares.** Calibration is local in maturity, which keeps each subproblem small and convex-friendly. Calendar arbitrage is structurally avoided by the forward-only propagation: every $C_{j+1}$ is derived from $C_j$ through a positive monotone operator.
- **Unconstrained IV inversion with absolute value.** A simple guard rather than a full constrained optimiser, justified by the small grid size. A constrained Brent or bracketed bisection would be more robust at scale.

### Limitations and natural extensions

- No hard no-arbitrage constraints beyond the implicit scheme itself and visual checks.
- The IV inversion uses unconstrained Nelder–Mead; bracketed bisection would be safer.
- **Possible improvements:** (i) warm-starting $v$ from the previous tenor, (ii) adding light Tikhonov regularisation across $K$ to smooth the local-vol profile, (iii) testing alternative seeds/optimisers to quantify sensitivity.

---

## Repository structure

```
.
├── andreasen_huge_calibration.ipynb   # Jupyter notebook: full pipeline
├── andreasen_huge_report.pdf          # Method, figures and discussion
├── requirements.txt                   # Python dependencies
├── .gitignore
└── README.md
```

---

## How to run

```bash
git clone https://github.com/Theo2co/local-volatility-andreasen-huge.git
cd local-volatility-andreasen-huge
python -m venv .venv
source .venv/bin/activate            # Windows: .venv\Scripts\activate
pip install -r requirements.txt
jupyter lab                          # then open andreasen_huge_calibration.ipynb
```

Requires Python 3.10+.

---

## Reference

Andreasen, J. and Huge, B. (2010). *Volatility Interpolation*. Danske Markets working paper. SSRN: https://ssrn.com/abstract=1694972

---

## Authors

Group work submitted for *Advanced Derivatives* at EPFL:

- **Théodore Decaux** — [@Theo2co](https://github.com/Theo2co) (repository maintainer)
- **Elias Bourgon**
- **Jason Santangelo**

---

## License

Released under the [MIT License](LICENSE).
