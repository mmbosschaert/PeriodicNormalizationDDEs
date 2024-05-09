using Revise

using DDEBifTool
using GLMakie
using LinearAlgebra
using DataFrames
using GLM
using Infiltrator

# define constants, and active parameters (ap)
const b = 0.9
const ϵ₀ = 0.08
const c = 2.0528
const d = -3.2135
const ap = [1, 2]
const dims = 2
const τs = [0.0; 1.7722]

# define model
function fhn(u, p)
  β, α = p
  x, y = u[:, 1]
  xτ, _ = u[:, 2]
  du = similar(u[:, 1])

  du[1] = -x^3 / 3 + (c + α) * x^2 + d * x - y + 2 * β * tanh(xτ)
  du[2] = ϵ₀ * (x - b * y)

  du
end

# calculate multi-linear forms
jet = getJet(fhn, dims, τs);

# define equilibria near Hopf point
# parameter names to index
par_indx = (β=1, α=2)

# set parameters for generalized Hopf point 
β = 1.899972376314736
α = -1.042940901468666
params = [β; α]

# set parameter bounds
# parameterbounds = (min=[-0.5; 5.0], max=[0.5; 15.0])

# calculate stability of equilibrium
stst1 = stst(zeros(dims), params)
stst1 = stability(jet, stst1, τs)

# calculate normal form coefficients
genh1 = point_to_genhopf(jet, stst1, τs)
genh1 = normalform(jet, genh1, τs)

ϵ = 0.01
ntst = 20
ncol = 3
lpc_guess = generalizedHopfToPsol(jet, genh1, ϵ, ntst, ncol, τs)
lpc_brI = SetupLPCBranch(jet, lpc_guess, τs, δ=0.001, δmin=1e-06, δmax=0.01, MaxNumberofSteps=2000);
continue!(lpc_brI)
get_params(br) = hcat([point.parameters for point in br]...)
lines(get_params(lpc_brI.points))

genh1_higherorder = DDEBifTool.normalform_beta(jet, genh1, τs)

@unpack K10, K01, K02, K11, c₂, c₃, a3201 = genh1_higherorder.nmfm

β₁ = ϵ -> real(c₂) * ϵ^4 + 2(real(c₃) - a3201 * real(c₂)) * ϵ^6
β₂ = ϵ -> -2real(c₂) * ϵ^2 + (4a3201 * real(c₂) - 3 * real(c₃)) * ϵ^4

# continue Hopf branch
hopf1 = point_to_hopf(jet, stst1, τs)

hopf_branchI = SetupHopfBranch(jet, hopf1, τs, MaxNumberofSteps=250, δ=0.01, δmin=1e-08, δmax=0.1)
continue!(hopf_branchI)
reverse_branch!(hopf_branchI)
continue!(hopf_branchI)

fig = Figure()
ax = Axis(fig[1, 1])
lines!(get_params(hopf_branchI.points), label="Hopf branch")
scatter!(get_params(lpc_brI.points[1:10:300]), label="Computed Limit point cycle branch")
lines!(hcat([params + K10 * β₁(ϵ) + K01 * β₂(ϵ) + 0.5K02 * β₂(ϵ)^2 + K11 * β₁(ϵ) * β₂(ϵ) for ϵ in 0.001:0.001:0.07]...),
      label="Predicted Limit point cycle branch")
scatter!([genh1.parameters[1]], [genh1.parameters[2]], color=:black, label="Generalized Hopf point", markersize=10)
xlims!(ax, 1.89985, 1.90005)
ylims!(ax, -1.06, -1.02)
axislegend(ax, position=:lb)
