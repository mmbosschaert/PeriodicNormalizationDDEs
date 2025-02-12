using Revise
using PeriodicNormalizationDDEs
const PN = PeriodicNormalizationDDEs
using GLMakie
using LinearAlgebra
using DataFrames
using GLM

# define constants, and active parameters (ap)
const gᵤ = 0.1;
const gᵥ = 0.52;
const β = 0.1;
const ap = [1, 2];
const dims = 2;
const τs = [0.0; 1.0];

# define model
function acs(u, p)
  ζ, τ = p
  x, y = u[:, 1]
  xτ, yτ = u[:, 2]
  du = similar(u[:, 1])
  du[1] = τ * y
  du[2] = τ * (-x - gᵤ * xτ - 2 * ζ * y - gᵥ * yτ + β * xτ^3)
  du
end

# calculate multi-linear forms
jet = getJet(acs, dims, τs);

# define equilibria near Hopf point
# parameter names to index
par_indx = (ζ=1, τ=2)

# set parameters for generalized Hopf point 
ζ = -0.2632524975782207
τ = 6.157979479822465
params = [ζ; τ]

# set parameter bounds
parameterbounds = (min=[-0.5; 5.0], max=[0.5; 15.0])

# calculate stability of equilibrium
stst1 = stst(zeros(dims), params)
stst1 = stability(jet, stst1, τs)

# calculate normal form coefficients
genh1 = point_to_genhopf(jet, stst1, τs)
genh1 = normalform(jet, genh1, τs)

# Create convergence plot for second order
ϵ = 0.01
ntst = 20
ncol = 3

δps1 = exp10.(LinRange(-1.2, -0.9, 20))
relative_errors_second_order = zeros(length(δps1))
for (i, δp) in enumerate(δps1)
  lpc_guess = generalizedHopfToPsol(jet, genh1, δp, ntst, ncol, τs)
  lpc_brI = SetupLPCBranch(jet, lpc_guess, τs)
  lpc_corrected = lpc_brI.points[1]
  relative_errors_second_order[i] = norm([vcat(lpc_corrected.profile...); lpc_corrected.parameters; lpc_corrected.period] -
                                         [vcat(lpc_guess.profile...); lpc_guess.parameters; lpc_guess.period], Inf) /
                                    norm([vcat(lpc_corrected.profile...); lpc_corrected.parameters; lpc_corrected.period], Inf)
end

df_second_order = DataFrame([log10.(collect(δps1)) log10.(relative_errors_second_order)], [:ρ, :relative_error])
linearRegressor_second_order = lm(@formula(relative_error ~ ρ), df_second_order)
coeffs_second_order = coef(linearRegressor_second_order)

fig = Figure()
ax = Axis(fig[1, 1], xlabel="δp", ylabel="Relative error", title="Convergence plot for Active control system")
scatter!(log10.(collect(δps1)), log10.(relative_errors_second_order), label="second order")
lines!(log10.(collect(δps1)), coeffs_second_order[1] .+ log10.(collect(δps1 .^ coeffs_second_order[2])))

# Create convergence plot for mixed order
genh1 = point_to_genhopf(jet, stst1, τs)
genh1 = PN.normalform_beta(jet, genh1, τs)

δps2 = exp10.(LinRange(-1.2, -0.9, 20))
relative_errors_higher_order = zeros(length(δps2))
for (i, δp) in enumerate(δps2)
  lpc_guess = PN.generalizedHopfToPsolHigherOrder(jet, genh1, δp, ntst, ncol, τs)
  lpc_brI = SetupLPCBranch(jet, lpc_guess, τs)
  lpc_corrected = lpc_brI.points[1]
  relative_errors_higher_order[i] = norm([vcat(lpc_corrected.profile...); lpc_corrected.parameters; lpc_corrected.period] -
                                         [vcat(lpc_guess.profile...); lpc_guess.parameters; lpc_guess.period], Inf) /
                                    norm([vcat(lpc_corrected.profile...); lpc_corrected.parameters; lpc_corrected.period], Inf)
end

df_higher_order = DataFrame([log10.(collect(δps2)) log10.(relative_errors_higher_order)], [:ρ, :relative_error])
linearRegressor_higher_order = lm(@formula(relative_error ~ ρ), df_higher_order)
coeffs_higher_order = coef(linearRegressor_higher_order)

fig = Figure()
ax = Axis(fig[1, 1], xlabel="δp", ylabel="Relative error", title="Convergence plot for Active control system")
scatter!(log10.(collect(δps1)), log10.(relative_errors_second_order), label="second order")
lines!(log10.(collect(δps1)), coeffs_second_order[1] .+ log10.(collect(δps1 .^ coeffs_second_order[2])))
scatter!(log10.(collect(δps2)), log10.(relative_errors_higher_order), label="higher order")
lines!(log10.(collect(δps2)), coeffs_higher_order[1] .+ log10.(collect(δps2 .^ coeffs_higher_order[2])))
axislegend(position = :rb)
