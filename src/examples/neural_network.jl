
using Revise
using PeriodicNormalizationDDEs
using PeriodicNormalizationDDEs:getJet, LocateDoubleHopf, delete_one_from_approx_equal
using PeriodicNormalizationDDEs:determine_unfolding_case
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

# define constants and inline functions
const dims = 2
const b = 2.0
const d = 1.2
const τ₁ = 11.6
const τ₂ = 20.3

@inline g(x) = (tanh(x - 1.0) + tanh(1.0)) * cosh(1.0)^2;

# define model
function neuralMassModel(u, p)
  a, c = p
  x₁, x₂ = u[:, 1]
  x₁τ₁, x₂τ₁ = u[:, 2]
  x₁τ₂, x₂τ₂ = u[:, 3]
  du = similar(u[:, 1])
  du[1] = -x₁ - a * g(b * x₁τ₁) + c * g(d * x₂τ₂)
  du[2] = -x₂ - a * g(b * x₂τ₁) + c * g(d * x₁τ₂)
  du
end

# define model for simulation
function neuralMassModel!(du, u, h, p, t)
  a, c = p
  x₁, x₂ = u
  x₁τ₁, x₂τ₁ = h(p, t - τ₁)
  x₁τ₂, x₂τ₂ = h(p, t - τ₂)
  du[1] = -x₁ - a * g(b * x₁τ₁) + c * g(d * x₂τ₂)
  du[2] = -x₂ - a * g(b * x₂τ₁) + c * g(d * x₁τ₂)
end

# set delays
τs = [0.0; τ₁; τ₂]

# calculate multi-linear forms
jet = getJet(neuralMassModel, dims, τs)

## continue equilibria
# parameter names to index
par_indx = (a=1, c=2)

# continuation parameter for equilibria
con_par = par_indx.c

# parameters
a = 0.069
c = 0.4
params = [a; c]

# set parameter bounds
parameterbounds = (min=[-0.1; -0.1], max=[1.2; 1.0])

stst_branch = SetupStstBranch(neuralMassModel, con_par, zeros(dims), params, τs, dims; parameterbounds=parameterbounds, δmax=0.01)
continue!(stst_branch)

# calculate stability of the equilibria points
stability(jet.Ms, stst_branch.points, τs)

# Function to reverse and concatenate parameters
reverse_params(p) = reduce(hcat, [p.parameters[end:-1:1] for p in p])

fig = Figure(; fontsize=18)
ax1 = Axis(fig[1, 1])
ax2 = Axis(fig[1, 2])
lines!(ax1, reverse_params(stst_branch.points))
vlines!(ax2, [0.0], color=:red, linewidth=6)
ax1.xlabel = "c"
ax1.ylabel = "a"
active_branch = stst_branch;
slider = Slider(fig[2, 1:2], range=1:length(stst_branch.points))
ind = Observable(1)
connect!(ind, slider.value)
pars = @lift active_branch.points[$ind].parameters[end:-1:1]
scatter!(ax1, @lift(Point2f($pars)), markersize=20)
λ = @lift active_branch.points[$ind].stability
eigenvalues_stst = @lift(Point2f.(vec(real($λ)), vec(imag($λ))))
scatter!(ax2, eigenvalues_stst)
xlims!(ax1, (0.0, 1.2))
ylims!(ax1, (-1.0, 1.0))

# detect Hopf points
hopf_points = LocateHopfPoints(jet, stst_branch, τs)

# add Hopf points to the bifurcation diagram
for p in hopf_points
  scatter!(ax1, [p.parameters[2]], [p.parameters[1]], color=:red)
end

hopf_branches = []
for p in hopf_points
  hopf_branch = SetupHopfBranch(jet, p, τs; parameterbounds=parameterbounds)
  continue!(hopf_branch)
  reverse_branch!(hopf_branch)
  continue!(hopf_branch)
  # calculate stability along the Hopf curve
  stability(jet.Ms, hopf_branch.points, τs)
  push!(hopf_branches, hopf_branch)
end

fig = Figure()
set_theme!(theme_light())
ax1 = Axis(fig[1, 1])
lines!(ax1, reverse_params(stst_branch.points))
for hopf_branch in hopf_branches
  lines!(ax1, reverse_params(hopf_branch.points), color=:royalblue)
end
ax2 = Axis(fig[1, 2])
lines!(ax2, hcat([zeros(100), range(-2.0, 2.0, 100)]...), color=:red, linewidth=6)
active_branch = hopf_branches[1]
slider = Slider(fig[2, 1:2], range=1:length(active_branch.points))
#
ind = Observable(1)
connect!(ind, slider.value)
#
pars = @lift active_branch.points[$ind].parameters[end:-1:1]
scatter!(ax1, @lift(Point2f($pars)), markersize=20)
μ = @lift active_branch.points[$ind].stability
μ_to_points = @lift(Point2f.(vec(real($μ)), vec(imag($μ))))
scatter!(ax2, μ_to_points)

# detect and locate double hopf points
double_hopf_points = []
for hopf_branch in hopf_branches
  ind_double_hopf = locate_double_hopf(hopf_branch.points)
  for p in hopf_branch.points[ind_double_hopf]
    p_corrected = LocateDoubleHopf(jet, p, τs; MaxIter=100, tol=1e-12)
    push!(double_hopf_points, p_corrected)
    # scatter!(ax1, [p_corrected.parameters[2]], [p_corrected.parameters[1]], color=:black)
    # put text near the double hopf point
    # text!(p_corrected.parameters[2], p_corrected.parameters[1], text="DH")
  end
end

# remove duplicates
delete_one_from_approx_equal!(double_hopf_points)

# plot double hopf points
for (i, p) in enumerate(double_hopf_points)
  scatter!(ax1, [p.parameters[2]], [p.parameters[1]], color=:black)
  text!(ax1, p.parameters[2], p.parameters[1], text="DH$i")
end

# add normal form coefficients
# ugly workaround to avoid error
function compute_nmfm_coefficients(double_hopf_points)
  new_double_hopf_points = []
  for hoho in double_hopf_points
    try
      push!(new_double_hopf_points, normalform(jet, hoho, τs))
    catch
      println("Failed to calculate normal form coefficients")
    end
  end
  new_double_hopf_points
end

# add normal form coefficients
double_hopf_points = compute_nmfm_coefficients(double_hopf_points)

determine_unfolding_case(double_hopf_points)

# detect zero hopf points 
ind_zero_hopf1 = locate_zero_hopf(hopf_branches[1].points)
ind_zero_hopf2 = locate_zero_hopf(hopf_branches[2].points)
zero_hopf_points1 = hopf_branches[1].points[ind_zero_hopf1]
zero_hopf_points2 = hopf_branches[2].points[ind_zero_hopf2]
scatter!(ax1, reverse_params(hopf_branches[1].points[ind_zero_hopf1]), color=:red)
scatter!(ax1, reverse_params(hopf_branches[2].points[ind_zero_hopf2]), color=:red)

zeho1 = point_to_zeho(jet, zero_hopf_points1[1], τs)

normalform(jet, zeho1, τs)


# detect generalized Hopf points
genh_points = []
for (i, hopf_branch) in enumerate(hopf_branches)
  ind_genh = detect_genh(jet, hopf_branch.points, τs)
  p_correct = [locate_genh(jet, hopf_branch, indx, par_indx.c, τs; MaxIter=100, tol=1e-10) for indx in ind_genh]
  p_correct = [p for p in p_correct if p !== nothing]
  push!(genh_points, p_correct)
  # if genh_points !== []
  #   for p in genh_points
  #     scatter!(ax1, , color=:red)
  #   end
  # end
end
genh_points = vcat(genh_points...)

delete_one_from_approx_equal!(genh_points)

for (i, p) in enumerate(genh_points)
  scatter!(ax1, [p.parameters[2]], [p.parameters[1]], color=:red)
  text!(ax1, p.parameters[2], p.parameters[1], text="GH$i")
end

hoho1 = double_hopf_points1[1] # OK

determine_unfolding_case(double_hopf_points)

hoho1 = double_hopf_points[3]

# plot predicteors in parameter plane
# fig = Figure()
# set_theme!(theme_light())
# ax1 = Axis(fig[1, 1])
# lines!(ax1, reverse_params(stst_branch.points))
# for hopf_branch in hopf_branches
#   lines!(ax1, reverse_params(hopf_branch.points), color=:royalblue)
# end
# scatter!(ax1, [hoho1.parameters[2]], [hoho1.parameters[1]], color=:black, label="Double Hopf point")
# scatter!(ax1, hcat([hoho1.parameters + hoho1.nmfm.K * [-real(hoho1.nmfm.g1011); -real(hoho1.nmfm.g0021)] .* ϵ^2 for ϵ ∈ range(0, 0.1, 100)]...)[end:-1:1, :], color=:royalblue, label="Predictor Neimark-Sacker curve I")
# scatter!(ax1, hcat([hoho1.parameters + hoho1.nmfm.K * [-real(hoho1.nmfm.g2100); -real(hoho1.nmfm.g1110)] .* ϵ^2 for ϵ ∈ range(0, 0.1, 100)]...)[end:-1:1, :], color=:orangered, label="Predictor Neimark-Sacker curve II")

# continuation of Neimark-Sacker curves emanating from a double Hopf point
ϵ₁ = 0.01
ϵ₂ = 0.01
ntst = 20
ncol = 3
nsbranch = SetupNSBranch(jet, hoho1, ϵ₁, ϵ₂, ntst, ncol, τs;
  parameterbounds=parameterbounds, δ=0.001, δmin=1e-06, δmax=[0.1; 0.10],
  MaxNumberofSteps=[200; 200]);

continue!(nsbranch[1])
reverse_branch!(nsbranch[1])
continue!(nsbranch[1])
continue!(nsbranch[2])
# reverse_branch!(nsbranch[2])
# continue!(nsbranch[2])

# [p.period for p in nsbranch[2].points]

# indx = findall(p -> maximum(vcat(p.profile...)) - minimum(vcat(p.profile...)) > 1e-02, ns_brII.points)
# ns_brII = @set ns_brII.points = ns_brII.points[indx]

# compute normalform coefficients
ap = [1, 2]
normal_form_coefficients!(jet, nsbranch[1], τs, ap)
normal_form_coefficients!(jet, nsbranch[2], τs, ap)

# plot Neimark-Sacker branches
fig = Figure()
set_theme!(theme_light())
ax1 = Axis(fig[1, 1])
for hopf_branch in hopf_branches
  lines!(ax1, reverse_params(hopf_branch.points), color=:royalblue)
end
lines!(ax1, reverse_params(nsbranch[1].points))
lines!(ax1, reverse_params(nsbranch[2].points))
scatter!(ax1, reverse_params(nsbranch[1].points[real(get_nmfm_coefficients(nsbranch[1])).<0.0]), color=:orangered)
scatter!(ax1, reverse_params(nsbranch[1].points[real(get_nmfm_coefficients(nsbranch[1])).>0.0]), color=:royalblue)
scatter!(ax1, reverse_params(nsbranch[2].points[real(get_nmfm_coefficients(nsbranch[2])).<0.0]), color=:orangered)
scatter!(ax1, reverse_params(nsbranch[2].points[real(get_nmfm_coefficients(nsbranch[2])).>0.0]), color=:royalblue)

multipliers!(jet, nsbranch[1], τs)
multipliers!(jet, nsbranch[2], τs)

# interactive plot
set_theme!(theme_dark())
fig = Figure()
ax1 = Axis(fig[1, 1])
ax2 = Axis(fig[1, 2])
ax3 = Axis(fig[2, 2])
ax4 = Axis(fig[2, 1])
ax5 = Axis3(fig[1, 3])
lines!(ax1, reverse_params(hopf_brI.points), color=Hopf_color)
lines!(ax1, reverse_params(hopf_brII.points), color=Hopf_color)
# lines!(ax1,reverse_params(psol_brI.points))
lines!(ax2, cos.(2 * pi * range(0, 1, 100)), sin.(2 * pi * range(0, 1, 100)))
active_branch = nsbranch[1].points
scatter!(ax1, reverse_params(nsbranch[2].points[real(get_nmfm_coefficients(nsbranch[2])).<0.0]), color=:red, label="NS curve")
scatter!(ax1, reverse_params(nsbranch[2].points[real(get_nmfm_coefficients(nsbranch[2])).>0.0]), color=:blue, label="NS curve")
scatter!(ax1, reverse_params(nsbranch[1].points[real(get_nmfm_coefficients(nsbranch[1])).<0.0]), color=:red, label="NS curve")
scatter!(ax1, reverse_params(nsbranch[1].points[real(get_nmfm_coefficients(nsbranch[1])).>0.0]), color=:blue, label="NS curve")
lines!(ax1, reverse_params(nsbranch[2].points), color=:green)
lines!(ax1, reverse_params(lpc_brI.points), color=:purple)
slider = Slider(fig[3, :], range=1:length(active_branch))
ind = Observable(1)
connect!(ind, slider.value)
pars = @lift active_branch[$ind].parameters[end:-1:1]
μ = @lift active_branch[$ind].stability
scatter!(ax1, @lift(Point2f($pars)), markersize=14, color=:purple)
μ_points = @lift(Point2f.(vec(real($μ)), vec(imag($μ))))
scatter!(ax2, μ_points)
# plot profile
ts = @lift active_branch[$ind].mesh
profile = @lift active_branch[$ind].profile
data1 = @lift(Point2f.($ts[:], hcat($profile...)[1, :]))
data2 = @lift(Point2f.($ts[:], hcat($profile...)[2, :]))
lines!(ax3, data1)
lines!(ax3, data2)

points = Observable(Point2f[])
scatter!(ax1, points, color=:green, marker=:xcross, markersize=10, label="selected points")

on(events(ax1.scene).mousebutton) do event
  if event.button == Mouse.left && is_mouseinside(ax1.scene)
    if event.action == Mouse.press
      mp = mouseposition(ax1.scene)
      points[] = [mp]
      notify(points)
    end
  end
end

condition(u, t, integrator) = u[1] > 40.0
affect!(integrator) = terminate!(integrator)
cb = DiscreteCallback(condition, affect!)
prob(h, p, tspan) = DDEProblem(neuralMassModel!, h(p, 0), h, tspan, p; constant_lags=[τ₁, τ₂], callback=cb)
alg = MethodOfSteps(Tsit5())

tspan = (0.0, 1e5)
tspan = (0.0, 5e4)

# h(out, p,t) = (out .= [0.3; 0.0; 0.0])

period = active_branch[ind[]].period

active_branch[ind[]].profile[1][1]

dt = 1.0e-4

# h(p, t) = 0.07ones(dims)
# h(p, t) = [active_branch[ind[]].profile[10][1]*cos(2pi * t / period);active_branch[ind[]].profile[10][1]*sin(2pi * t / period)]
# h(p, t) = [active_branch[ind[]].profile[10][1]*cos(2pi * t / period);active_branch[ind[]].profile[10][1]*sin(2pi * t / period)]

function h(p, t)
  γ = active_branch[ind[]].profile
  mesh = active_branch[ind[]].mesh
  ncol = active_branch[ind[]].ncol
  T = active_branch[ind[]].period
  # interpolate(t/period,0.0,γ,mesh,ncol)
  0.99 * interpolate(0.0, t, γ, T, mesh, ncol)
end

# scatter!(ax3, hcat(h.(Ref(nothing), range(0,period,100))...)[1,:],hcat(h.(Ref(nothing), range(0,period,100))...)[2,:])
# scatter!(ax3, hcat(h.(Ref(nothing), range(0,period,100))...)[1,:])

button = Button(fig, label="Simulate")

on(button.clicks) do n

  tspan = (0.0, 1e5)
  c, a = points[][1].data
  # sol = solve(prob(h,[a; c], tspan), alg, reltol=1e-08, abstol=1e-06, dtmax=0.1) 
  sol = solve(prob(h, [a; c], tspan), alg)
  empty!(ax3)
  lines!(ax3, sol.t, sol[1, :])
  scatter!(ax3, range(0, period, 100), hcat(h.(Ref(nothing), range(0, period, 100))...)[2, :])

  # empty!(ax3)
  # lines!(ax3, sol[1, :], sol[2, :])
  scatter!(ax3, range(0, period, 100), hcat(h.(Ref(nothing), range(0, period, 100))...)[2, :])

  tend = sol.t[end]
  # method I
  sol_u = hcat(sol(range(tend - tend / 10, tend, step=dt)).u...)
  xdelayed = hcat(sol(range(tend - tend / 10 - 1.0, tend - 1.0, step=dt)).u...)[1, :]
  # Calculate the sign difference between adjacent elements
  sign_diff = diff(sign.(xdelayed))
  # Find the indices where the sign changes (zero crossings)
  zero_crossings = findall(sign_diff .!= 0)
  empty!(ax4)
  scatter!(ax4, sol_u[:, zero_crossings], color=:orangered)
  # scatter!(ax4, sol_u)
  delete!(ax5)
  ax5 = Axis3(fig[1, 3])
  lines!(ax5, [xdelayed'; sol_u])
  scatter!(ax5, [zeros(length(zero_crossings))'; sol_u[:, zero_crossings]], color=:orangered)
end

set_theme!(theme_light())

# Create 3d plot
fig = Figure()
ax1 = Axis3(fig[1, 1], xlabel="c", ylabel="a", zlabel="normal form coefficient real(d)")
# lines!(get_params(lpc_branchI.points), color=:green)
lines!([reverse_params(nsbranch[2].points); 0.0 * real(get_nmfm_coefficients(nsbranch[2])')])
scatter!([reverse_params(nsbranch[2].points); real(get_nmfm_coefficients(nsbranch[2])')])
zlims!(ax1, (-0.001, 0.001))












# select Hopf point near c ≈ 0.7 hopf_branchI
indx = findfirst(p -> 0.68 < p.parameters[2] < 0.7, hopf_branches[1].points)
hopf_point = hopf_branches[1].points[indx]
scatter!(ax1, reverse_params([hopf_branches[1].points[indx]]), color=:magenta)

# continu psol branch emanating from hopf point
psol_branch = SetupPsolBranch(jet, con_par, hopf_point, τs; parameterbounds=parameterbounds, MaxNumberofSteps=120, ntst=20, ncol=4, δ=0.01, δmax=0.5)
reverse_branch!(psol_branch)
continue!(psol_branch)

# Plot parameters vs max x1 of profile
psol_param_max_profile_data = hcat(
  [p.parameters[2] for p in psol_branch.points],
  [maximum(reduce(hcat, p.profile)[1, :]) for p in psol_branch.points])
lines(psol_param_max_profile_data)
scatter!(psol_param_max_profile_data, markersize=6)

# # Plot parameters vs amplitude
# psol_param_amplitude_profile_data = hcat(
#     [p.parameters[2] for p in psol_branch.points], 
#     [sum(abs.(extrema(reduce(hcat,p.profile)[1,:]))) for p in psol_branch.points])
# 
# lines(psol_param_amplitude_profile_data)
# scatter!(psol_param_amplitude_profile_data)

multipliers!(jet, psol_branch::NamedTuple, τs)

psol_branch = detectSpecialPoints(psol_branch)

# add bifucation points to the plot
fig = Figure()
ax1 = Axis(fig[1, 1])
lines!(ax1, psol_param_max_profile_data)
scatter!(ax1, [p.parameters[2] for p in psol_branch.points[psol_branch.specialpoints.indx_fold]],
  [maximum(reduce(hcat, p.profile)[1, :]) for p in psol_branch.points[psol_branch.specialpoints.indx_fold]],
  color=:orange, marker=:square, markersize=14, label="fold")
scatter!(ax1, [p.parameters[2] for p in psol_branch.points[psol_branch.specialpoints.indx_pd]],
  [maximum(reduce(hcat, p.profile)[1, :]) for p in psol_branch.points[psol_branch.specialpoints.indx_pd]],
  color=:magenta, marker=:xcross, markersize=14, label="period doubling")
scatter!(ax1, [p.parameters[2] for p in psol_branch.points[psol_branch.specialpoints.indx_ns]],
  [maximum(reduce(hcat, p.profile)[1, :]) for p in psol_branch.points[psol_branch.specialpoints.indx_ns]],
  color=:royalblue, marker=:diamond, markersize=14, label="Neimark-Sacker")
axislegend(ax1)

# continuation of lpc
lpc_guess = psol_branch.points[psol_branch.specialpoints.indx_fold[1]]
lpc_branchI = SetupLPCBranch(jet, lpc_guess, τs; parameterbounds=parameterbounds, δ=0.1, δmax=0.3)
continue!(lpc_branchI)
reverse_branch!(lpc_branchI)
continue!(lpc_branchI)

# # plot parameters
lines(reverse_params(lpc_branchI.points))

using BifurcationKit
lpc_poI = load("/home/maikel/.julia/dev/DDEBifTool/src/examples/lpc_poI.jld2", "lpc_poI")
lpc_poI_40 = load("/home/maikel/.julia/dev/DDEBifTool/src/examples/lpc_poI_40.jld2", "lpc_poI_40")
lpc_poI_N13 = load("/home/maikel/.julia/dev/DDEBifTool/src/examples/lpc_poI_N13.jld2", "lpc_poI")
hp_codim2_1 = load("/home/maikel/.julia/dev/DDEBifTool/src/examples/hp_codim2_1.jld2", "hp_codim2_1")
hp_codim2_2 = load("/home/maikel/.julia/dev/DDEBifTool/src/examples/hp_codim2_2.jld2", "hp_codim2_2")

# plot parameters from pseudo-spectral metohod
fig = Figure()
ax1 = Axis(fig[1, 1])
lines!(reverse_params(lpc_branchI.points))
lines!(lpc_poI.c, lpc_poI.a, linewidth=2)
scatter!(lpc_poI.c, lpc_poI.a)

# 
# # compute periodic normal form coefficients
# ap = [1, 2]
# # lpc_br1_coefficients = [PeriodicNormalizationDDEs.fold_coefficient_julia(jet,jet.M,psol,τs,dims,3,20,ap) for psol ∈ lpc_branchI.points]
# normal_form_coefficients!(jet, lpc_branchI, τs, ap)
# get_nmfm_coefficients(branch) = [p.nmfm for p in branch.points]
# 
# 
# # using GaussQuadrature
# # nodes = first(legendre(3))
# # lpc_br1_coefficients = [fold_coefficient(jet,jet.M,psol.profile,
# #     psol.period*vec(psol.mesh),psol.parameters,τs,psol.period,
# #     hcat([psol.period*vec(psol.mesh)[i*3+1] .+ (psol.period*vec(psol.mesh)[(i+1)*3+1] - psol.period*vec(psol.mesh)[i*3+1])/2 .* (nodes .+ 1) for i in 0:20-1]...),
# #     2,3,psol.period*vec(psol.mesh)[1:3:end],20) for psol ∈ lpc_branchI.points]
# 
# scatter!(reverse_params(lpc_branchI.points[get_nmfm_coefficients(lpc_branchI).<0.0]))
# scatter!(reverse_params(lpc_branchI.points[get_nmfm_coefficients(lpc_branchI).>0.0]))
# 
# # plot both coordinate of the periodic orbit against parameter c
# fig = Figure()
# ax1 = Axis3(fig[1, 1])
# # lines!(reverse_params(hopf_branchI.points), color=:royalblue)
# # lines!(reverse_params(hopf_branchII.points), color=:royalblue)
# lines!(reverse_params(lpc_branchI.points), color=:green)
# lines!([reverse_params(lpc_branchI.points); get_nmfm_coefficients(lpc_branchI)'])
# zlims!(ax1, (-4, 4))
# 
# # compute multipliers
# multipliers!(jet, lpc_branchI, τs)
# 
# # continue ns point
# ns_guess = psol_branch.points[psol_branch.specialpoints.indx_ns][1]
# ns_guess = psol_to_ns(ns_guess)
# 
# ns_branchI = SetupNSBranch(jet, ns_guess, τs; parameterbounds=parameterbounds, δ=0.002, δmax=0.1, MaxNumberofSteps=100)
# continue!(ns_branchI)
# reverse_branch!(ns_branchI)
# continue!(ns_branchI)
# 
# # compute normalform coefficients
# normal_form_coefficients!(jet, ns_branchI, τs, ap)
# 
# # compute multipliers
# multipliers!(jet, lpc_branchI, τs)
# 
# fig = Figure()
# ax1 = Axis(fig[1:2, 1])
# ax2 = Axis(fig[1, 2])
# ax3 = Axis(fig[2, 2])
# lines!(ax1, reverse_params(lpc_branchI.points))
# scatter!(ax1, reverse_params(lpc_branchI.points[get_nmfm_coefficients(lpc_branchI).<0.0]))
# scatter!(ax1, reverse_params(lpc_branchI.points[get_nmfm_coefficients(lpc_branchI).>0.0]))
# scatter!(ax1, reverse_params(ns_branchI.points[real(get_nmfm_coefficients(ns_branchI)).<0.0]))
# scatter!(ax1, reverse_params(ns_branchI.points[real(get_nmfm_coefficients(ns_branchI)).>0.0]))
# scatter!(ax1, hcat([p.parameters[[2, 1]] for p in ns_branchI.points]...))
# 
# # 3d figure showing normal form coefficient along the parameter bifurcation curve
# fig = Figure()
# ax1 = Axis3(fig[1, 1])
# # lines!(reverse_params(hopf_branchI.points), color=:royalblue)
# # lines!(reverse_params(hopf_branchII.points), color=:royalblue)
# lines!(reverse_params(lpc_branchI.points), color=:green)
# lines!(reverse_params(ns_branchI.points), color=:magenta)
# lines!([reverse_params(ns_branchI.points); real(get_nmfm_coefficients(ns_branchI))'])
# zlims!(ax1, (-0.01, 0.01))



# continue pd as ns with ω = π
pd_guess = psol_branch.points[psol_branch.specialpoints.indx_pd][2]
pd_branchI = SetupPDBranch(jet, pd_guess, τs; parameterbounds=parameterbounds, δ=0.001, δmax=0.1, MaxNumberofSteps=500)
continue!(pd_branchI)
reverse_branch!(pd_branchI)
continue!(pd_branchI)

# compute normalform coefficients
normal_form_coefficients!(jet, pd_branchI, τs, ap)

# compute multipliers
multipliers!(jet, pd_branchI, τs)

# 3d figure showing normal form coefficient along the parameter bifurcation curve
fig = Figure(fontsize=18)
ax1 = Axis3(fig[1, 1], xlabel="c", ylabel="a", zlabel="normal form coefficient c", title="Period doubling branch with normal coefficients")
# lines!(reverse_params(lpc_branchI.points), color=:green)
# lines!(reverse_params(ns_branchI.points), color=:magenta)
# lines!(reverse_params(pd_branchI.points), color=:red)
# lines!([reverse_params(ns_branchI.points); real(get_nmfm_coefficients(ns_branchI))'])
lines!([reverse_params(pd_branchI.points); real(get_nmfm_coefficients(pd_branchI))'])
scatter!(ax1, reverse_params(pd_branchI.points[real(get_nmfm_coefficients(pd_branchI)).<0.0]))
scatter!(ax1, reverse_params(pd_branchI.points[real(get_nmfm_coefficients(pd_branchI)).>0.0]))
zlims!(ax1, (-0.01, 0.01))

## continue lpc curse with double period near GPD point
# detect GPD points
pd_branchI = LocateGPDPoints(pd_branchI)

# filter GPD with parameter c > 0.5
gpd_indx = pd_branchI.specialpoints.gpd_indx
gpd_indx = [i for (i, p) in enumerate(pd_branchI.points[gpd_indx]) if p.parameters[end] > 0.55]
gpd_indx = pd_branchI.specialpoints.gpd_indx[gpd_indx]
gpd_points = pd_branchI.points[gpd_indx.+1]

# branch off near GPD2 to detect LPC curve with double period
psol1 = [pd_branchI.points[gpd_indx[2]-10]]

scatter!(reverse_params(gpd_points), color=:red)
scatter!(reverse_params(psol1), color=:magenta)

psol_branchII = SetupPsolBranch(jet, con_par, psol1[1], τs;
  parameterbounds=parameterbounds, MaxNumberofSteps=100, δ=0.01, δmax=0.01)
# reverse_branch!(psol_branchII)
continue!(psol_branchII)

# plot
lines([p.parameters[2] for p in psol_branchII.points], [maximum(reduce(hcat, p.profile)[1, :]) for p in psol_branchII.points])
scatter!([p.parameters[2] for p in psol_branchII.points], [maximum(reduce(hcat, p.profile)[1, :]) for p in psol_branchII.points])

multipliers!(jet, psol_branchII, τs)

psol_branchII = detectSpecialPoints(psol_branchII)

# plot in two parameters
lines(reverse_params(pd_branchI.points), color=:royalblue)
lines!(reverse_params(psol_branchII.points), color=:royalblue)
scatter!(reverse_params(gpd_points), color=:red)
scatter!(reverse_params(psol1), color=:lightgreen)
scatter!(reverse_params([psol_branchII.points[psol_branchII.specialpoints.indx_fold[1]]]), color=:purple)

# continuation of lpc
lpc_guess = psol_branchII.points[psol_branchII.specialpoints.indx_fold[2]]
lpc_branchII = SetupLPCBranch(jet, lpc_guess, τs; parameterbounds=parameterbounds, δ=0.01, δmax=0.1, MaxNumberofSteps=900)
continue!(lpc_branchII)
reverse_branch!(lpc_branchII)
continue!(lpc_branchII)

# length(lpc_branchII.points)
# # remove points with profiles that are too small
# indx = findall(p -> maximum(vcat(p.profile...)) - minimum(vcat(p.profile...)) > 1e-02, lpc_branchII.points)
# lpc_branchII = @set lpc_branchII.points = lpc_branchII.points[indx]
# length(lpc_branchII.points)

normal_form_coefficients!(jet, lpc_branchII, τs, ap)

# plot parameters
fig = Figure(fontsize=24)
ax1 = Axis3(fig[1, 1], xlabel="c", ylabel="a", zlabel="normal form coefficient c", title="Period doubling branch with normal coefficients")
# lines!(reverse_params(lpc_branchI.points), color=:green)
# lines!(reverse_params(ns_branchI.points), color=:magenta)
# lines!(reverse_params(pd_branchI.points), color=:red)
# lines!([reverse_params(ns_branchI.points); real(get_nmfm_coefficients(ns_branchI))'])
lines!(ax1, [reverse_params(pd_branchI.points); real(get_nmfm_coefficients(pd_branchI))'], linewidth=2)
scatter!(ax1, reverse_params(pd_branchI.points[real(get_nmfm_coefficients(pd_branchI)).<0.0]))
scatter!(ax1, reverse_params(pd_branchI.points[real(get_nmfm_coefficients(pd_branchI)).>0.0]))
zlims!(ax1, (-0.01, 0.01))
ax2 = Axis(fig[1, 2], title="LPC branch emanting from, and connecting, the two GPD point")
lines!(ax2, reverse_params(pd_branchI.points), linewidth=2)
# scatter!(reverse_params(pd_branchI.points[real(get_nmfm_coefficients(pd_branchI)).<0.0]))
# scatter!(reverse_params(pd_branchI.points[real(get_nmfm_coefficients(pd_branchI)).>0.0]))
scatter!(ax2, reverse_params(lpc_branchII.points[get_nmfm_coefficients(lpc_branchII).<0.0]))
scatter!(ax2, reverse_params(lpc_branchII.points[get_nmfm_coefficients(lpc_branchII).>0.0]))
scatter!(ax2, reverse_params(gpd_points), color=:red, markersize=18)
ax2.xlabel = "c"
ax2.ylabel = "a"
aCairoMakie.activate!()
CairoMakie.save("neural_mass_model_period_doubling.pdf", fig)
GLMakie.activate!()

# find sign change index of get_nmfm_coefficients(lpc_branchII)
sign_diff = diff(sign.(get_nmfm_coefficients(lpc_branchII)))
# extract non zero indices
zero_crossings = findall(sign_diff .!= 0)
# plot the sign change
lines(reverse_params(lpc_branchII.points), linewidth=2, label="Period doubling branch")
scatter!(reverse_params(lpc_branchII.points[zero_crossings]))

# plot parameters
fig = Figure(fontsize=25)
ax1 = Axis(fig[1, 1], title="LPC branch emanting from, and connecting, the two GPD points")
# lines!(ax1, reverse_params(pd_branchI.points), linewidth=2, label="Period doubling branch")
scatter!(reverse_params(pd_branchI.points[real(get_nmfm_coefficients(pd_branchI)).<0.0]), label="Period Doubling curve: c < 0")
scatter!(reverse_params(pd_branchI.points[real(get_nmfm_coefficients(pd_branchI)).>0.0]), label="Period Doubling curve: c > 0")
scatter!(ax1, reverse_params(lpc_branchII.points[get_nmfm_coefficients(lpc_branchII).<0.0]), color=Cycled(5), label="LPC curve: b < 0")
scatter!(ax1, reverse_params(lpc_branchII.points[get_nmfm_coefficients(lpc_branchII).>0.0]), color=Cycled(6), label="LPC curve: b > 0")
# plot GH
for (i, p) in enumerate(gpd_points)
  scatter!(ax1, [p.parameters[2]], [p.parameters[1]], color=:black, markersize=12)
  text!(ax1, p.parameters[2], p.parameters[1], text="GPD$i", align=(:right, :top))
end
ax1.xlabel = "c"
ax1.ylabel = "a"
axislegend(ax1, position=(0.72, 0.85))
# axislegend(ax1, orientation = :horizontal)
# scatter!(ax1, reverse_params(gpd_points), color=:black, markersize=12)
ax2 = Axis(fig[1, 2], title="Zoom near Cusp")
scatter!(reverse_params(pd_branchI.points[real(get_nmfm_coefficients(pd_branchI)).<0.0]))
scatter!(reverse_params(pd_branchI.points[real(get_nmfm_coefficients(pd_branchI)).>0.0]))
scatter!(ax2, reverse_params(lpc_branchII.points[get_nmfm_coefficients(lpc_branchII).<0.0]), color=Cycled(5))
scatter!(ax2, reverse_params(lpc_branchII.points[get_nmfm_coefficients(lpc_branchII).>0.0]), color=Cycled(6))
ax2.xlabel = "c"
ax2.ylabel = "a"
xlims!(ax2, (0.5014001278354234, 0.5181064070294856))
ylims!(ax2, (0.14164443894351003, 0.15094391363995052))
# fig[2, :] = Legend(fig, ax1, "", orientation=:horizontal)
CairoMakie.activate!()
CairoMakie.save("neural_mass_model_period_doubling_zoom_near_cusp.pdf", fig)
GLMakie.activate!()

# plot both coordinate of the periodic orbit against parameter c
fig = Figure()
ax1 = Axis3(fig[1, 1])
lines!(reverse_params(lpc_branchII.points))
lines!([reverse_params(lpc_branchII.points); get_nmfm_coefficients(lpc_branchII)'])
zlims!(ax1, (-1, 1))

# compute multipliers
multipliers!(jet, psol_branchII, τs)

# select Hopf point near c ≈ 0.7 hopf_branchII
indx = findfirst(p -> 0.68 < p.parameters[2] < 0.7, hopf_branchII.points)
hopf_pointII = hopf_branchII.points[indx]

# continu psol branch emanating from hopf point
psol_branchIII = SetupPsolBranch(jet, con_par, hopf_pointII, τs; parameterbounds=parameterbounds, MaxNumberofSteps=120, ntst=40, ncol=4, δ=0.01, δmax=0.5)
reverse_branch!(psol_branchIII)
continue!(psol_branchIII)

# Plot parameters vs max x1 of profile
psol_param_max_profile_data = hcat(
  [p.parameters[2] for p in psol_branchIII.points],
  [maximum(reduce(hcat, p.profile)[1, :]) for p in psol_branchIII.points])
lines(psol_param_max_profile_data)
scatter!(psol_param_max_profile_data, markersize=6)

# # Plot parameters vs amplitude
# psol_param_amplitude_profile_data = hcat(
#     [p.parameters[2] for p in psol_branchIII.points], 
#     [sum(abs.(extrema(reduce(hcat,p.profile)[1,:]))) for p in psol_branchIII.points])
# 
# lines(psol_param_amplitude_profile_data)
# scatter!(psol_param_amplitude_profile_data)

multipliers!(jet, psol_branchIII, τs)

psol_branchIII = detectSpecialPoints(psol_branchIII)

# add bifucation points to the plot
scatter!([p.parameters[2] for p in psol_branchIII.points[psol_branchIII.specialpoints.indx_fold]],
  [maximum(reduce(hcat, p.profile)[1, :]) for p in psol_branchIII.points[psol_branchIII.specialpoints.indx_fold]],
  color=:orange, marker=:square, markersize=14)

# continuation of lpc
lpc_guess = psol_branchIII.points[psol_branchIII.specialpoints.indx_fold[1]]
lpc_branchIII = SetupLPCBranch(jet, lpc_guess, τs; parameterbounds=parameterbounds, δ=0.1, δmax=2, MaxNumberofSteps=1000)
continue!(lpc_branchIII)
reverse_branch!(lpc_branchIII)
lpc_branchIII = @set lpc_branchIII.MaxNumberofSteps = 200
continue!(lpc_branchIII)

# plot parameters
lines(reverse_params(lpc_branchIII.points))
scatter!(reverse_params(lpc_branchIII.points[end-1:end]))

# compute periodic normal form coefficients
normal_form_coefficients!(jet, lpc_branchIII, τs, ap)

# compute stability
multipliers!(jet, lpc_branchIII, τs)

lines(reverse_params(lpc_branchI.points))
scatter!(reverse_params(lpc_branchI.points[get_nmfm_coefficients(lpc_branchI).<0.0]))
scatter!(reverse_params(lpc_branchI.points[get_nmfm_coefficients(lpc_branchI).>0.0]))
lines!(reverse_params(lpc_branchII.points))
scatter!(reverse_params(lpc_branchII.points[get_nmfm_coefficients(lpc_branchII).<0.0]))
scatter!(reverse_params(lpc_branchII.points[get_nmfm_coefficients(lpc_branchII).>0.0]))
lines!(reverse_params(lpc_branchIII.points))
scatter!(reverse_params(lpc_branchIII.points[get_nmfm_coefficients(lpc_branchIII).<0.0]))
scatter!(reverse_params(lpc_branchIII.points[get_nmfm_coefficients(lpc_branchIII).>0.0]))

# 3d figure showing normal form coefficient along the parameter bifurcation curve
fig = Figure()
ax1 = Axis3(fig[1, 1])
coeffs = real(get_nmfm_coefficients(lpc_branchIII))
coeffs /= first(findmax(abs.(coeffs)))
lines!([reverse_params(lpc_branchIII.points); coeffs'], markersize=2000)
scatter!(ax1, reverse_params(lpc_branchIII.points[real(get_nmfm_coefficients(lpc_branchIII)).<0.0]))
scatter!(ax1, reverse_params(lpc_branchIII.points[real(get_nmfm_coefficients(lpc_branchIII)).>0.0]))
zlims!(ax1, (-0.1, 0.1))

# interactive plot
fig = Figure()
ax1 = Axis(fig[1:2, 1])
ax2 = Axis(fig[1, 2])
ax3 = Axis(fig[2, 2])
lines!(ax1, reverse_params(lpc_branchIII.points))
scatter!(ax1, reverse_params(lpc_branchIII.points[get_nmfm_coefficients(lpc_branchIII).<0.0]))
scatter!(ax1, reverse_params(lpc_branchIII.points[get_nmfm_coefficients(lpc_branchIII).>0.0]))
lines!(ax2, cos.(2 * pi * range(0, 1, 100)), sin.(2 * pi * range(0, 1, 100)))
active_branch = lpc_branchIII.points
slider = Slider(fig[3, :], range=1:length(active_branch))
ind = Observable(1)
connect!(ind, slider.value)
pars = @lift active_branch[$ind].parameters[end:-1:1]
μ = @lift active_branch[$ind].stability
scatter!(ax1, @lift(Point2f($pars)), markersize=12, color=:red)
μ_points = @lift(Point2f.(vec(real($μ)), vec(imag($μ))))
scatter!(ax2, μ_points)
# plot profile
ts = @lift active_branch[$ind].mesh
profile = @lift active_branch[$ind].profile
data = @lift(Point2f.($ts[:], hcat($profile...)[2, :]))
lines!(ax3, data)
























# Continue limit point curve emanating from the third generalized Hopf point
ϵ = 0.1
ntst = [20 40 20]
ncol = 4
lpc_branches = []
for (i, genh_point) in enumerate(genh_points)
  lpc_guess = generalizedHopfToPsol(jet, genh_point, ϵ, ntst[i], ncol, τs)
  lpc_branch = SetupLPCBranch(jet, lpc_guess, τs;
    parameterbounds=parameterbounds, δ=0.01, δmax=0.1, MaxNumberofSteps=5000, NumberOfFails=2)
  @time continue!(lpc_branch)
  push!(lpc_branches, lpc_branch)
end

lpc_guess = generalizedHopfToPsol(jet, genh_points[2], ϵ, 40, ncol, τs)
lpc_branch_II_ntst20 = SetupLPCBranch(jet, lpc_guess, τs;
  parameterbounds=parameterbounds, δ=0.01, δmax=0.1, MaxNumberofSteps=5000, NumberOfFails=2)
@time continue!(lpc_branch_II_ntst20)

fig = Figure()
set_theme!(theme_light())
ax1 = Axis(fig[1, 1])
for hopf_branch in hopf_branches
  lines!(ax1, reverse_params(hopf_branch.points), color=:royalblue)
end
for lpc_branch in lpc_branches
  lines!(ax1, reverse_params(lpc_branch.points), color=:green)
end

# compute normal form coefficients
for lpc_branch in lpc_branches
  normal_form_coefficients!(jet, lpc_branch, τs, ap)
end

#$ # add normal form coefficients to the plot
#$ for lpc_branch in lpc_branches
#$   scatter!(ax1, reverse_params(lpc_branch.points[get_nmfm_coefficients(lpc_branch).<0.0]), color=:red)
#$   scatter!(ax1, reverse_params(lpc_branch.points[get_nmfm_coefficients(lpc_branch).>0.0]), color=:blue)
#$ end

lpc_branch = lpc_branches[3]
scatter!(ax1, reverse_params(lpc_branch.points[get_nmfm_coefficients(lpc_branch).<0.0]), color=:red)
scatter!(ax1, reverse_params(lpc_branch.points[get_nmfm_coefficients(lpc_branch).>0.0]), color=:blue)

# plot pseudo-spectral method
lines!(ax1, lpc_poI.c, lpc_poI.a, linewidth=2)
scatter!(ax1, lpc_poI_40.c, lpc_poI_40.a)
lines!(ax1, hp_codim2_1.c, hp_codim2_1.a, linewidth=2)
lines!(ax1, hp_codim2_2.c, hp_codim2_2.a, linewidth=2)

using BifurcationKit
lpc_poI = load("/home/maikel/.julia/dev/DDEBifTool/src/examples/lpc_poI.jld2", "lpc_poI")
lpc_poI_40 = load("/home/maikel/.julia/dev/DDEBifTool/src/examples/lpc_poI_40.jld2", "lpc_poI_40")
lpc_poI_N13 = load("/home/maikel/.julia/dev/DDEBifTool/src/examples/lpc_poI_N13.jld2", "lpc_poI")
hp_codim2_1 = load("/home/maikel/.julia/dev/DDEBifTool/src/examples/hp_codim2_1.jld2", "hp_codim2_1")
hp_codim2_2 = load("/home/maikel/.julia/dev/DDEBifTool/src/examples/hp_codim2_2.jld2", "hp_codim2_2")

# create figure for Odo, just comparing one of the LPC curves with the pseudo-spectral method
fig = Figure(fontsize=22, linewidth=3)
fig[0, 1:3] = Label(fig, "Comparing LPC curve approximation via direct discretization of the DDE or via the pseudo-spectral method (PSM)")
ax1 = Axis(fig[1, 1], xlabel="c", ylabel="a", title="Direct discretization of the DDE")
lines!(ax1, reverse_params(hopf_branches[2].points), label="Hopf branch")
scatter!(ax1, reverse_params([genh_points[2]]), color=:red, label="Generalized Hopf point")
lines!(ax1, reverse_params(lpc_branch_II_ntst20.points), label="LPC curve with ntst=40 (5 min)")
axislegend(ax1, position=:lt)
# lines!(ax1, reverse_params(lpc_branches[2].points))
xlims!(ax1, (0.40877622991118046, 0.8))
ax2 = Axis(fig[1, 2], xlabel="c", ylabel="a", title="Comparison with pseudo-spectral method")
lines!(ax2, reverse_params(hopf_branches[2].points))
scatter!(ax2, reverse_params([genh_points[2]]), color=:red)
lines!(ax2, reverse_params(lpc_branches[2].points), label="LPC curve with ntst=20 (5 min)")
lines!(ax2, lpc_poI.c, lpc_poI.a, label="LPC PSM with N=12, ntst=20 (113 minutes)")
lines!(ax2, lpc_poI_40.c, lpc_poI_40.a, label="LPC PSM with N=12, ntst=40 (6.6 hours)", linestyle=:dot)
lines!(ax2, lpc_poI_N13.c, lpc_poI_N13.a, label="LPC PSM with N=13, ntst=20 (130 minutes)")
xlims!(ax2, (0.40877622991118046, 0.8))
axislegend(ax2, position=:lt)
ax3 = Axis(fig[1, 3], xlabel="c", ylabel="a", title="Zoom near left most CUSP point ")
lines!(ax3, reverse_params(hopf_branches[2].points), label="Hopf branch")
lines!(ax3, reverse_params(lpc_branches[2].points))
lines!(ax3, lpc_poI.c, lpc_poI.a)
lines!(ax3, lpc_poI_40.c, lpc_poI_40.a, linestyle=:dot)
lines!(ax3, lpc_poI_N13.c, lpc_poI_N13.a)
xlims!(ax3, (0.4315667234279448, 0.5586676164360006))
ylims!(ax3, (-0.025404869154176224, 0.18664204132014747))
# fig[2, :] = Legend(fig, ax2, "", orientation=:horizontal)
CairoMakie.activate!()
CairoMakie.save("LPC_curve_neural_mass_model.pdf", fig)
GLMakie.activate!()

