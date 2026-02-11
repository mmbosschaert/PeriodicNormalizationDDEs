# 🌀 PeriodicNormalizationDDEs.jl

> **Unveiling the Critical Dynamics of Limit Cycles in Delay Differential Equations.**

[![arXiv](https://img.shields.io/badge/arXiv-2505.19786-b31b1b.svg)](https://arxiv.org/abs/2505.19786)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Julia](https://img.shields.io/badge/Julia-1.9+-9558B2.svg?logo=julia)](https://julialang.org)
[![BifurcationKit](https://img.shields.io/badge/Plugin-BifurcationKit.jl-blue.svg)](https://github.com/bifurcationkit/BifurcationKit.jl)

**PeriodicNormalizationDDEs.jl** is a specialized Julia package for computing the **critical normal form coefficients** of codimension-1 bifurcations of limit cycles in Delay Differential Equations (DDEs). 

Designed to work seamlessly as a plugin for the powerful [BifurcationKit.jl](https://github.com/bifurcationkit/BifurcationKit.jl) ecosystem, this tool bridges the gap between abstract functional analysis and high-performance numerics. It implements the explicit computational formulas derived in the 2025 breakthrough work using the **Sun-Star Calculus** framework.

## 🚀 Key Features

* **BifurcationKit Integration:** Directly consumes bifurcation points detected by `BifurcationKit.jl`, streamlining your workflow from continuation to criticality analysis.
* **Codimension-1 Support:** Handles **Fold** (Limit Point), **Flip** (Period-Doubling), and **Neimark-Sacker** (Torus) bifurcations of limit cycles.
* **Criticality Detection:** Explicitly computes the coefficients needed to distinguish between **nondegenerate, subcritical, and supercritical** bifurcations.
* **Numerical Robustness:** Introduces the **Characteristic Operator** to transform abstract operator equations into boundary-value problems (BVPs) solvable via **Orthogonal Collocation**.
* **Sun-Star Powered:** Built on the rigorous "Sun-Star" dual semigroup perturbation framework, ensuring theoretical soundness for delay equations.

---

## 📚 The Science Behind the Code

This software is the numerical realization of a trilogy of theoretical advancements in the bifurcation theory of DDEs by **Bram Lentjes, Maikel M. Bosschaert, Len Spek, and Yuri A. Kuznetsov**.

### 1. The Numerical Breakthrough: Periodic Normalization
**[Numerical Periodic Normalization at Codim 1 Bifurcations of Limit Cycles in DDEs](https://arxiv.org/abs/2505.19786)** This preprint provides the explicit computational formulas implemented here. It introduces the **Characteristic Operator**, a crucial innovation that allows the abstract homological equations of the normal form theory to be solved using robust, standard BVP solvers rather than unstable shooting methods.

### 2. The Foundation: Periodic Center Manifolds
**[Periodic Center Manifolds for DDEs in the Light of Suns and Stars](https://link.springer.com/article/10.1007/s10884-023-10289-9)** *Journal of Dynamics and Differential Equations (2023)* Establishes the existence of smooth periodic center manifolds for non-hyperbolic cycles in DDEs using the sun-star calculus, providing the rigorous domain where our computations take place.

### 3. The Framework: Periodic Normal Forms
**[Periodic Normal Forms for Bifurcations of Limit Cycles in DDEs](https://www.sciencedirect.com/science/article/abs/pii/S0022039625000725)** *Journal of Differential Equations (2025)* Provides the rigorous theoretical construction of the periodic normalization framework. It proves the existence of a special coordinate system on the center manifold, constructed via **time-periodic smooth Jordan chains** for the original and adjoint operators, which allows the local dynamics to be described by canonical periodic normal forms.

---

## 🛠️ How It Works

The package automates the complex "Periodic Normalization" pipeline:

1.  **Input:** Takes a periodic orbit bifurcation point detected by `BifurcationKit.jl`.
2.  **Operator Construction:** Constructs the *Characteristic Operator* associated with the linearization around the cycle.
3.  **BVP Solving:** Instead of tackling the infinite-dimensional DDE operator directly, we solve a sequence of Linear Boundary Value Problems (LBVPs) on the interval $[0, T]$ using **Orthogonal Collocation**.
4.  **Coefficient Extraction:** Integral formulas involving the solutions of these LBVPs (and the adjoint eigenfunctions) are evaluated to yield the critical normal form coefficient (e.g., the first Lyapunov coefficient).
5.  **Output:** A numerical value determining the bifurcation's criticality (e.g., `c < 0` $\to$ Supercritical/Stable, `c > 0` $\to$ Subcritical/Unstable).

---

## 💻 Installation

To install the package, use the Julia package manager:

```julia
using Pkg
Pkg.add(url="[https://github.com/mmbosschaert/PeriodicNormalizationDDEs.jl](https://github.com/mmbosschaert/PeriodicNormalizationDDEs.jl)")
