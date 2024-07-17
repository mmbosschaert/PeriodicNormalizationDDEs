using Revise
using DDEBifTool
using GLMakie
using LinearAlgebra
using SparseArrays
using Symbolics
using Setfield
using Infiltrator
using DifferentialEquations
using JLD
using CairoMakie

GLMakie.activate!()

function bifurcationvalues(a, b, c, d, chi, r, S)
  xstar = 1 / (3 * a) * (b - d + sqrt((b - d)^2 - 3 * a * c * S))
  Iapp = xstar^2 * (a * xstar - b + d) + c * (S * (xstar - chi) - 1)

  # Auxiliary variables to calculate tau
  A = xstar^2 * ((3 * a * xstar + 2 * d)^2 - 4 * b^2) - 2 * r * xstar * (2 * d * xstar - 1) * (3 * a * xstar - 2 * b + 2 * d) + r^2 * (4 * d * xstar * (-2 * b * xstar + d * xstar - 1) + 1)
  B = 9 * a^2 * xstar^4 + 2 * r * xstar * (3 * a * xstar - 2 * b + 2 * d) - 4 * b^2 * xstar^2 - 4 * d * xstar + r^2 + 1
  D = -1 / 2 * B - 1 / 2 * sqrt(B^2 - 4 * A)
  omega = if D < 0
    sqrt(-1 / 2 * B + 1 / 2 * sqrt(B^2 - 4 * A))
  else
    sqrt(-1 / 2 * B - 1 / 2 * sqrt(B^2 - 4 * A))
  end

  Y = omega / (2 * b) * (-1 / xstar + 2 * d / (1 + omega^2) + r * (2 * b - 2 * d - 3 * a * xstar) / (r^2 + omega^2))
  tau = 1 / omega * (asin(Y) + 2 * π)

  return xstar, Iapp, tau
end

# define constants and inline functions
const dims = 3
const a = 1.0
const b = 3.0
const c = 1.0
const d = 5.0
const χ = -1.6

# define model
function RoseHindmarsh(u, p; r=1.4)
  Iapp, S = p
  x, y, z = u[:, 1]
  xτ, yτ, zτ = u[:, 2]
  du = similar(u[:, 1])
  du[1] = y - a * x^3 + b * xτ^2 - c * z + Iapp
  du[2] = c - d * x^2 - y
  du[3] = r * (S * (x - χ) - z)
  du
end

S = -8.0
xstar, Iapp, τ = bifurcationvalues(a, b, c, d, χ, 1.4, S)
τs = [0.0 τ]

# calculate multi-linear forms
jet = DDEBifTool.getJet(RoseHindmarsh, dims, τs)

# calculate stability of equilibrium
ystar = c - d * xstar^2
zstar = S * (xstar - χ)
p = [Iapp, S]

# check for quilibrium
u = repeat([xstar, ystar, zstar], 1, 2)
RoseHindmarsh(u, p) ≈ zeros(3)

stst1 = stst(u[:,1], p)
# The target σ is the center around which eiganvalues are computed. The default value of 0.0 causes a problem.
stst1 = stability(jet, stst1, τs; σ = -1.0)

# calculate normal form coefficients
zeho1 = point_to_zeho(jet, stst1, τs)
zeho1 = normalform(jet, zeho1, τs)
