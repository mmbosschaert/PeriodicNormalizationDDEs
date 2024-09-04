using Revise

using DDEBifTool
using GLMakie
using LinearAlgebra
using DataFrames
using GLM
using Infiltrator
using UnPack

# define constants, and active parameters (ap)
const b = 0.9;
const ϵ₀ = 0.08;
const c = 2.0528;
const d = -3.2135;
const ap = [1, 2];
const dims = 2;
const τs = [0.0; 1.7722];

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
lpc_guess = generalizedHopfToPsol(jet, genh1, ϵ, ntst, ncol, τs);
lpc_brI = SetupLPCBranch(jet, lpc_guess, τs, δ=0.001, δmin=1e-06, δmax=0.01, MaxNumberofSteps=2000);
continue!(lpc_brI);
get_params(br) = hcat([point.parameters for point in br]...);
lines(get_params(lpc_brI.points));

genh1_higherorder = DDEBifTool.normalform_beta(jet, genh1, τs)

@unpack K10, K01, K02, K11, K03, c₂, c₃, a3201 = genh1_higherorder.nmfm

β₁ = ϵ -> real(c₂) * ϵ^4 + 2(real(c₃) - a3201 * real(c₂)) * ϵ^6
β₂ = ϵ -> -2real(c₂) * ϵ^2 + (4a3201 * real(c₂) - 3 * real(c₃)) * ϵ^4

# continue Hopf branch
hopf1 = point_to_hopf(jet, stst1, τs)

hopf_branchI = SetupHopfBranch(jet, hopf1, τs, MaxNumberofSteps=250, δ=0.01, δmin=1e-08, δmax=0.1)
continue!(hopf_branchI)
reverse_branch!(hopf_branchI)
continue!(hopf_branchI)

fig = Figure(size = (750, 800),dpi=300)
#fig = Figure(dpi=300)
ax = Axis(fig[2, 1], xlabel = "β", ylabel = "α")
lines!(get_params(hopf_branchI.points), label="subcritical Hopf branch", color=Cycled(2))
lines!(get_params(hopf_branchI.points)[:,1:77], label="supercritical Hopf branch", color=Cycled(6))
scatter!(get_params(lpc_brI.points[1:10:300]), label="Computed Limit point cycle branch")
lines!(hcat([params  + K10 * β₁(ϵ) +  K01 * β₂(ϵ)+ 0.5K02 * β₂(ϵ)^2 + K11 * β₁(ϵ) * β₂(ϵ) + (1/6) * K03 * β₂(ϵ)^3  for ϵ in 0.001:0.001:0.07]...),
       label="Higher order LPC predictor", color=Cycled(1))
scatter!([genh1.parameters[1]], [genh1.parameters[2]], color=:black, label="Generalized Hopf point", markersize=10)
xlims!(ax, 1.89985, 1.90005)
ylims!(ax, -1.06, -1.02)
hidedecorations!(ax, ticklabels=false, label=false, ticks=false)
#axislegend(ax, position=:lb)

ax2 = Axis(fig[1, 1], title = "FitzHugh-Nagumo model", xlabel = "β", ylabel = "α")
lines!(get_params(hopf_branchI.points), label="subcritical Hopf branch", color=Cycled(2))
lines!(get_params(hopf_branchI.points)[:,1:77], label="supercritical Hopf branch", color=Cycled(6))
scatter!(get_params(lpc_brI.points[1:10:300]), label="Computed Limit point cycle branch")
lines!(hcat([params  - 2 * real(c₂) * K01 * ϵ^2 for ϵ in 0.001:0.001:0.05]...),
       label="Linear LPC predictor", color=Cycled(1),linestyle=:dash)
lines!(hcat([params  + K10 * β₁(ϵ) +  K01 * β₂(ϵ) + 0.5K02 * β₂(ϵ)^2 + K11 * β₁(ϵ) * β₂(ϵ) + (1/6) * K03 * β₂(ϵ)^3 for ϵ in 0.001:0.001:0.01]...),
       label="Higher order LPC predictor", color=Cycled(1))
scatter!([genh1.parameters[1]], [genh1.parameters[2]], color=:black, label="Generalized Hopf point", markersize=10)
xlims!(ax2, 1.89985, 1.90005)
ylims!(ax2, -1.06, -1.02)
hidedecorations!(ax2, ticklabels=false, label=false, ticks=false)
#axislegend(ax2, position=:lb)
fig[1,2] = Legend(fig, ax2, "", framevisible=false)


# Create convergence plot for second order
genh1 = point_to_genhopf(jet, stst1, τs)
genh1 = normalform(jet, genh1, τs)

ϵ = 0.01
ntst = 20
ncol = 3

δps1 = exp10.(LinRange(-2.5, -0.8, 40))
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
ax = Axis(fig[1, 1], xlabel="δp", ylabel="Relative error", title="Convergence plot for FHN")
scatter!(log10.(collect(δps1)), log10.(relative_errors_second_order), label="second order")
lines!(log10.(collect(δps1)), coeffs_second_order[1] .+ log10.(collect(δps1 .^ coeffs_second_order[2])))


# Create convergence plot for mixed order
genh1 = point_to_genhopf(jet, stst1, τs)
genh1 = DDEBifTool.normalform_beta(jet, genh1, τs)

δps2 = exp10.(LinRange(-2.5, -0.8, 40))
relative_errors_higher_order = zeros(length(δps2))
for (i, δp) in enumerate(δps2)
  lpc_guess = DDEBifTool.generalizedHopfToPsolHigherOrder(jet, genh1, δp, ntst, ncol, τs)
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
ax = Axis(fig[1, 1], xlabel="δp", ylabel="Relative error", title="Convergence plot for FHN")
scatter!(log10.(collect(δps1)), log10.(relative_errors_second_order), label="second order")
#lines!(log10.(collect(δps1)), coeffs_second_order[1] .+ log10.(collect(δps1 .^ coeffs_second_order[2])))
scatter!(log10.(collect(δps2)), log10.(relative_errors_higher_order), label="higher order")
#lines!(log10.(collect(δps2)), coeffs_higher_order[1] .+ log10.(collect(δps2 .^ coeffs_higher_order[2])))
axislegend(position = :rb)
hidedecorations!(ax, ticklabels=false, label=false, ticks=false)
