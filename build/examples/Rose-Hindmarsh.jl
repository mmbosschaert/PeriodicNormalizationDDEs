
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

parameterbounds = (min=[-19.0; -8.05], max=[-18.70; -7.96])

# check for quilibrium
u = repeat([xstar, ystar, zstar], 1, 2)
RoseHindmarsh(u, p) ≈ zeros(3)

stst1 = stst(u[:, 1], p)
# The target σ is the center around which eiganvalues are computed. The default value of 0.0 causes a problem.
stst1 = stability(jet, stst1, τs; σ=-1.0)

# calculate normal form coefficients
zeho1 = point_to_zeho(jet, stst1, τs)
zeho1 = normalform(jet, zeho1, τs)

zeho1.nmfm.h200.(-τs)

# ϵ = 1.0e-4
# ntst = 20
# ncol = 4
# psol_guess = DDEBifTool.foldHopfToNS(jet, zeho1, ϵ, ntst, ncol, τs)
# fig = Figure()
# ax1 = Axis(fig[1, 1])
# lines!([p[1] for p in psol_guess.profile])
# lines!(repeat([zeho1.coords[1]], 60))
# ax2 = Axis(fig[1, 2])
# lines!([p[2] for p in psol_guess.profile])
# lines!(repeat([zeho1.coords[2]], 60))
# ax3 = Axis(fig[1, 3])
# lines!([p[3] for p in psol_guess.profile])
# lines!(repeat([zeho1.coords[3]], 60))

# scatter(vcat([[real(z) imag(z)] for z in psol_guess.stability]...))
# lines!(cos.(range(0, 2π, 100)), sin.(range(0, 2π, 100)))


# parameter predictions
psol_guess = [DDEBifTool.foldHopfToNS(jet, zeho1, ϵ, ntst, ncol, τs) for ϵ in 0:0.0001:0.01]
fig = Figure()
ax1 = Axis(fig[1, 1])
scatter!(ax1, [zeho1.parameters[1]], [zeho1.parameters[2]])
scatter!(hcat([psol_guess[i].parameters for i in 1:100]...))

ϵ = 0.001
ns_branch = SetupNSBranch(jet, zeho1, ϵ, ntst, ncol, τs; parameterbounds=parameterbounds, δmax=0.1);
continue!(ns_branch)

# calculate normal form coefficients
ap = [1, 2]
normal_form_coefficients!(jet, ns_branch, τs, ap)

# to continue Hopf curve emanating from zeho1
hopf = DDEBifTool.zeho_to_hopf(zeho1)
hopf_br = SetupHopfBranch(jet, hopf, τs; parameterbounds=parameterbounds, δ=0.04, δmax=0.04, MaxNumberofSteps=250, NumberOfFails=4)

continue!(hopf_br)
reverse_branch!(hopf_br)
continue!(hopf_br)

fig = Figure()
ax1 = Axis(fig[1, 1])
scatter!(ax1, [zeho1.parameters[1]], [zeho1.parameters[2]], color=:black, markersize=12, label="Zero-Hopf point")
lines!(ax1, get_params(hopf_br.points), label="Hopf branch", color=Cycled(2))
lines!(ax1, get_params(ns_branch.points), label="Neimark-Sacker branch", color=Cycled(3))
scatter!(ax1, get_params(ns_branch.points[real(get_nmfm_coefficients(ns_branch)).<0.0]), color=:red, label=L"Neimark-Sacker normal form coefficient $\text{Re}(d)<0.0$")
# plot zero Hopf point again to get on top
scatter!(ax1, [zeho1.parameters[1]], [zeho1.parameters[2]], color=:black, markersize=12)
ax1.xlabel = "Iapp"
ax1.ylabel = "S"
axislegend(ax1, position=:lt)






# interactive stability plot of Hopf branch
stability(jet.Ms, hopf_br.points, τs, σ=-1.0)

fig = Figure()
set_theme!(theme_light())
ax1 = Axis(fig[1, 1])
lines!(ax1, get_params(hopf_br.points), color=:royalblue)
ax2 = Axis(fig[1, 2])
vlines!(ax2, 0.0, color=:red, linewidth=2)
slider = Slider(fig[2, 1:2], range=1:length(hopf_br.points))
ind = Observable(1)
connect!(ind, slider.value)
pars = @lift hopf_br.points[$ind].parameters
scatter!(ax1, @lift(Point2f($pars)), markersize=20)
μ = @lift hopf_br.points[$ind].stability
μ_to_points = @lift(Point2f.(vec(real($μ)), vec(imag($μ))))
scatter!(ax2, μ_to_points)
