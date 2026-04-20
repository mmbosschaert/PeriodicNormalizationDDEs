# 🌀 PeriodicNormalizationDDEs.jl

> **Bifurcations of Limit Cycles in Delay Differential Equations: Theory & Software.**

[![Publication](https://img.shields.io/badge/DOI-10.1137%2F25M1763573-blue)](https://doi.org/10.1137/25M1763573)
[![Julia](https://img.shields.io/badge/Julia-1.9+-9558B2.svg?logo=julia)](https://julialang.org)
[![BifurcationKit](https://img.shields.io/badge/Plugin-BifurcationKit.jl-blue.svg)](https://bifurcationkit.github.io/BifurcationKitDocs.jl/stable/)

**PeriodicNormalizationDDEs.jl** is a specialized Julia package for computing the **critical normal form coefficients** of codimension one bifurcations of limit cycles in Delay Differential Equations (DDEs). 

Designed to work seamlessly as a plugin for the powerful [BifurcationKit.jl](https://github.com/bifurcationkit/BifurcationKit.jl) ecosystem, this tool bridges the gap between abstract functional analysis and high-performance numerics. It implements the explicit computational formulas derived in the 2026 breakthrough work using the **Numerical Periodic Normalization** and **Sun-Star Calculus** framework.

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
**[Numerical Periodic Normalization at Codim 1 Bifurcations of Limit Cycles in DDEs](https://bifurcationkit.github.io/BifurcationKitDocs.jl/stable/)** *SIAM Journal on Applied Dyanamical Systems (2026)* This paper bridges the gap between theory and code, providing the explicit formulas we implement. A major challenge in DDEs is solving the "homological equations" for normal form coefficients. This work introduces the **Characteristic Operator**, a crucial innovation that bypasses shooting methods entirely. It reformulates the abstract operator equations into Periodic Linear Boundary Value Problems. This allows us to use Orthogonal Collocation, solving for the solution over the entire period simultaneously, to compute the coefficients with high precision and stability.

### 2. The Framework: Periodic Normal Forms
**[Periodic Normal Forms for Bifurcations of Limit Cycles in DDEs](https://www.sciencedirect.com/science/article/abs/pii/S0022039625000725)** *Journal of Differential Equations (2025)* This paper provides the rigorous theoretical construction of the periodic normalization framework. It establishes a coordinate system on the periodic center manifold using periodic smooth Jordan chains for the original and adjoint operators. This framework allows us to study the local dynamics on the periodic center manifold in terms of periodic normal forms.

### 3. The Foundation: Periodic Center Manifolds
**[Periodic Center Manifolds for DDEs in the Light of Suns and Stars](https://link.springer.com/article/10.1007/s10884-023-10289-9)** *Journal of Dynamics and Differential Equations (2023)* This paper lays the necessary mathematical groundwork for the entire project. In Delay Differential Equations, the state space is infinite-dimensional because the system's future depends on its history, not just its current state. This "memory" makes standard geometric tools difficult to apply. This work overcomes those hurdles using the **Sun-Star Calculus** framework. It proves that despite the infinite-dimensional nature of the delays, there exists a smooth, finite-dimensional "center manifold" effectively capturing the critical dynamics. This method allows us to reduce the complex, infinite-dimensional DDE problem into a manageable, finite-dimensional form that our software can solve numerically.

---

## 🛠️ How It Works

The package automates the complex "Periodic Normalization" pipeline:

1.  **Input:** Takes a periodic orbit bifurcation point detected by `BifurcationKit.jl` or `DDE-BifTool.m`.
2.  **Operator Construction:** Constructs the *Characteristic Operator* associated with the linearization around the cycle.
3.  **BVP Solving:** Instead of tackling the infinite-dimensional DDE operator directly, we solve a sequence of Periodic Linear Boundary Value Problems on the interval $[0, T]$ using **Orthogonal Collocation**.
4.  **Coefficient Extraction:** Integral formulas involving the solutions of these periodic BVPs (and the adjoint eigenfunctions) are evaluated to yield the critical normal form coefficient.
5.  **Output:** A numerical value determining the bifurcation's criticality (e.g., `c < 0` $\to$ Supercritical/Stable, `c > 0` $\to$ Subcritical/Unstable).

---

## 💻 Installation

To install the package, use the Julia package manager:

```julia
using Pkg
Pkg.add(url="[https://github.com/mmbosschaert/PeriodicNormalizationDDEs.jl](https://github.com/mmbosschaert/PeriodicNormalizationDDEs.jl)")
