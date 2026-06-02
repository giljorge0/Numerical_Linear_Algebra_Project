

---

## Overview

This project implements four core numerical linear algebra algorithms in MATLAB, exploring their behaviour on classical test matrices (Vandermonde, tridiagonal, Hilbert, Lehmer) and a practical image-compression application via SVD.

All code is contained in a single file (`projeto_ALN_2026.m`) structured as a MATLAB Live Script, with sections corresponding to each question of the project statement.



---

## Topics Covered

### Question 1 — Newton-Schulz Iteration

Matrix inversion via the iterative scheme:

$$X_{k+1} = 2X_k - X_k A X_k$$

- **Theoretical proof** that the residuals satisfy $R_{k+1} = R_k^2$, guaranteeing quadratic convergence when $\rho(I - AX_0) < 1$.
- **Proof of convergence** for the initial guess $X_0 = A^\top / \|A\|_F^2$.
- **Implementation** with convergence criteria on $\|X_{k+1} - X_k\|_F$ and $\|R_k\|_F$.
- **Tests** on Vandermonde $V([1,3,5,7,9])$, tridiagonal $A_n$, Hilbert $H_n$, and Lehmer $L_n$ matrices.
- **Computational convergence order** computed as:

$$p \approx \frac{\ln(\|R_{k+1}\|_F / \|R_k\|_F)}{\ln(\|R_k\|_F / \|R_{k-1}\|_F)}$$

**Key observations:** Quadratic convergence (~order 2) is confirmed for well-conditioned matrices. The method diverges for Hilbert matrices with $n \geq 8$ due to extreme ill-conditioning.

---

### Question 2 — Householder QR Decomposition

QR factorisation $A = QR$ via Householder reflections:

$$v = x + \text{sign}(x_1)\|x\| e_1, \quad H = I - 2vv^\top / \|v\|^2$$

- Applied to Lehmer and Hilbert matrices $L_{100 \times n}$, $H_{100 \times n}$ for $n = 10, 20, \ldots, 100$.
- Residuals $\|QR - A\|_F$ and orthogonality error $\|Q^\top Q - I\|_F$ reported for all sizes.

**Key observations:** Both errors remain at machine-precision level ($\sim 10^{-14}$–$10^{-15}$) across all tested sizes, confirming the numerical stability of the Householder approach even for the severely ill-conditioned Hilbert matrix.

---

### Question 3 — Image Compression via SVD

A rank-$p$ approximation of an image matrix:

$$A_p = \sum_{i=1}^{p} \sigma_i u_i v_i^\top$$

with quality metric:

$$q = \frac{\sum_{i=1}^{p} \sigma_i^2}{\sum_{i=1}^{r} \sigma_i^2} \times 100\%$$

**3a — Grayscale image** (`cameraman.tif`, 256×256):

| % SVs kept | p | Quality % | Storage ratio |
|:---:|:---:|:---:|:---:|
| 1% | 3 | 94.08% | 0.023 |
| 5% | 13 | 97.88% | 0.102 |
| 10% | 26 | 99.04% | 0.204 |
| 25% | 64 | 99.78% | 0.501 |
| 50% | 128 | 99.98% | 1.002 |

**3b — Colour image** (`peppers.png`): SVD applied independently to each RGB channel, with per-channel and mean quality reported.

---

### Question 4 — Condition Number via the Power Method

The spectral condition number $\text{cond}_2(A) = |\lambda|_{\max} / |\lambda|_{\min}$ is estimated without forming $A^{-1}$:

- **Power method** → dominant eigenvalue $|\lambda|_{\max}$.
- **Inverse power method** (via LU factorisation) → smallest eigenvalue $|\lambda|_{\min}$.

**Hilbert matrices** $H_n$, $n = 5, \ldots, 14$:

| n | cond₂ (power method) | cond₂ (MATLAB) | Relative error |
|:---:|:---:|:---:|:---:|
| 5 | 4.766 × 10⁵ | 4.766 × 10⁵ | ~10⁻¹¹ |
| 10 | 1.602 × 10¹³ | 1.602 × 10¹³ | ~10⁻⁵ |
| 13 | 1.965 × 10¹⁷ | 1.326 × 10¹⁹ | ~0.99 |

Accuracy degrades for $n \geq 12$ as the matrix approaches numerical singularity.

**Lehmer matrices** $L_n$, $n = 10, 100, \ldots, 500$: the power method matches MATLAB's `cond` accurately up to $n = 200$; the inverse power method requires significantly more iterations as $n$ grows (up to 800 iterations for $n = 500$).

---
