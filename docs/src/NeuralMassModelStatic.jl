# # Neural Mass Model
#
# This tutorial demonstrates the correctness of our derivations for critical periodic normal form coefficients associated with the fold, period-doubling, and Neimark-Sacker bifurcations. For the periodic fold bifurcation, we follow a branch of limit points of cycles, showing that the critical normal form coefficient's sign correctly vanishes when passing a cusp point. Similarly, for the period-doubling bifurcation, we continue a curve of period-doubling bifurcations and identify points where the critical normal form coefficient vanishes, indicating a generalized period-doubling point. Normal form analysis predicts a curve of limit points of cycles with double the period emanating from these points, and we will demonstrate this prediction. Lastly, for the Neimark-Sacker bifurcation, we confirm that the normal form coefficients of the Neimark-Sacker curves emanating from a double Hopf point align with the predictions derived from the normal form coefficients at the double Hopf point.

# The model is defined by a system of coupled delay differential equations with two time delays, representing the interaction of neural populations in a simplified framework:


##!CODEBOUNDARY math
\begin{align}
  \dot{x}_1 &= - x_1 - a \cdot g(b \cdot x_{1, \tau_1}) + c \cdot g(d \cdot x_{2, \tau_2}) \\
  \dot{x}_2 &= - x_2 - a \cdot g(b \cdot x_{2, \tau_1}) + c \cdot g(d \cdot x_{1, \tau_2})
\end{align}
##!CODEBOUNDARY

#
# This tutorial guides readers through the process of verifying these derivations by integrating symbolic calculations with numerical continuation implemented in Julia. Along the way, we demonstrate how to set up the model, compute periodic solutions, and analyze bifurcation diagrams. By providing a detailed, step-by-step framework, this tutorial not only validates the correctness of our theoretical derivations but also equips researchers with the tools to investigate other models that incorporate delays.
#
# ## Load package
#
# We start by loading the necessary packages.
#

##!CODEBOUNDARY @example NeuralMassModel
using Revise
using PeriodicNormalizationDDEs
const PN = PeriodicNormalizationDDEs
using LinearAlgebra
using SparseArrays
using Symbolics
using Setfield
using Infiltrator
using DifferentialEquations
using JLD
using CairoMakie
##!CODEBOUNDARY

#
# ## Define constants
#
# Here, we define the constants used throughout the model. The `dims` variable represents the number of components. The delays, `τ₁` and `τ₂`, are provided in a vector `τs`, which includes a zero delay for convenience in calculations.
#

##!CODEBOUNDARY @example NeuralMassModel
const dims = 2
const b = 2.0
const d = 1.2
const τ₁ = 11.6
const τ₂ = 20.3
τs = [0.0; τ₁; τ₂]
##!CODEBOUNDARY

#
# ## Model definition
#
# The neural mass model is defined as a system of delay differential equations. The function `g(x)` acts as a nonlinear transfer function, using a hyperbolic tangent for sigmoid-like behavior. The main function, `neuralMassModel`, defines the system's dynamics in terms of delayed variables, parameters `a` and `c`, and the nonlinear interactions between two coupled units.
#

##!CODEBOUNDARY @example NeuralMassModel
@inline g(x) = (tanh(x - 1.0) + tanh(1.0)) * cosh(1.0)^2;
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
##!CODEBOUNDARY

#
# ## Model definition for simulation
#
# For simulation purposes, the model is adapted to use a standard form compatible with solvers for delay differential equations. The function `neuralMassModel!` computes the derivatives in place, given the current state `u`, the delayed state `h`, parameters `p`, and the current time `t`.
#

##!CODEBOUNDARY @example NeuralMassModel
function neuralMassModel!(du, u, h, p, t)
  a, c = p
  x₁, x₂ = u
  x₁τ₁, x₂τ₁ = h(p, t - τ₁)
  x₁τ₂, x₂τ₂ = h(p, t - τ₂)
  du[1] = -x₁ - a * g(b * x₁τ₁) + c * g(d * x₂τ₂)
  du[2] = -x₂ - a * g(b * x₂τ₁) + c * g(d * x₁τ₂)
end
##!CODEBOUNDARY

#
# ## Generate multi-linear forms
#
# We use the `getJet` function to create a structure containing the multi-linear forms of the neural mass model. These forms are used in calculations of normal form coefficients, continuation of equilibria and periodic orbits, and their stability analysis. By calculating the multi-linear forms symbolically, we ensure precise and efficient computation for subsequent analysis steps.
#

##!CODEBOUNDARY @example NeuralMassModel
# Calculate multi-linear forms
jet = getJet(neuralMassModel, dims, τs)
##!CODEBOUNDARY

#
# ## Equilibria continuation
#
# Continuing the tutorial, we now perform a continuation of equilibria for the
# neural mass model. This involves tracking steady states as a chosen parameter
# (`c`) varies. First, we define a named tuple `par_indx` for easy access to
# parameters and set `c` as the continuation parameter. Initial values for `a`
# and `c` are stored in `params`, and parameter bounds are defined with
# `parameterbounds` to constrain the values during continuation. The steady-state
# branch is initialized using `SetupStstBranch`, which prepares the model,
# parameter settings, and state variables for continuation. The process is then
# executed with `continue!`, calculating equilibria over the specified range of
# the continuation parameter. This step allows us to explore how the system's
# steady states respond to changes in the parameter `c`.
#

##!CODEBOUNDARY @example NeuralMassModel
# Create a named tuple for accessing parameters
par_indx = (a=1, c=2)

# continuation parameter for equilibria
con_par = par_indx.c

# parameters
a = 0.069
c = 0.4
params = [a; c]

# set parameter bounds
parameterbounds = (min=[0.0; 0.0], max=[1.2; 0.93])

stst_branch = SetupStstBranch(neuralMassModel, con_par, zeros(dims), params, τs, dims; parameterbounds=parameterbounds, δmax=0.01)
continue!(stst_branch)
##!CODEBOUNDARY

#
# ## Stability Calculation of Equilibria
#
# Next, we compute the stability of the steady-state branch. The function `stability` calculates the eigenvalues of the characteristics matrix at each point on the branch.
#

##!CODEBOUNDARY @example NeuralMassModel
stability(jet, stst_branch, τs)
##!CODEBOUNDARY

#
# ## Parameter Manipulation
#
# In this section, we define a convenience function to extract parameters from branch points, enabling easy access and manipulation of the parameter values. To align with the continuation setup, the function reverses the order of the parameters, placing `c` as the first coordinate.
#

##!CODEBOUNDARY @example NeuralMassModel
extract_pars(p) = reduce(hcat, [p.parameters[end:-1:1] for p in p])
##!CODEBOUNDARY

#
# ## Detect Hopf points and add stability information
#

##!CODEBOUNDARY @example NeuralMassModel
hopf_points = LocateHopfPoints(jet, stst_branch, τs)
stability(jet, hopf_points, τs)
##!CODEBOUNDARY

#
# ## Inspect Eigenvalues on steady-state branch
#
# Here we plot the steady-state branch and the eigenvalues at the first Hopf point. In the plot to the right, We see that are indeed two complex eigenvalues on the imaginary axis.
#

##!CODEBOUNDARY @example NeuralMassModel
fig = Figure();
ax1 = Axis(fig[1, 1], title = "Bifurcation diagram")
lines!(ax1, extract_pars(stst_branch.points), label="Steady-state branch")
xlims!(ax1, (0.0, 1.2))
ylims!(ax1, (0.0, 1.0))
ax1.xlabel = "c"
ax1.ylabel = "a"
scatter!(ax1, extract_pars(hopf_points), color=:red, label="Hopf points")
ax2 = Axis(fig[1, 2], title = "Eigenvalues at first Hopf point")
vlines!(ax2, [0.0], color=:red, linewidth=4)
λ = hopf_points[1].stability
eigenvalues_stst = Point2f.(vec(real(λ)), vec(imag(λ)))
scatter!(ax2, eigenvalues_stst)
axislegend(ax1, position=:lt)
fig
##!CODEBOUNDARY

#
# ## Continue Hopf points
#
# We continue the Hopf points in both parameters and in both directions.
#

##!CODEBOUNDARY @example NeuralMassModel
hopf_branches = PN.Branch[]
for hopf_point in hopf_points
  hopf_branch = SetupHopfBranch(jet, hopf_point, τs; parameterbounds=parameterbounds)
  continue!(hopf_branch)
  reverse_branch!(hopf_branch)
  continue!(hopf_branch)
  push!(hopf_branches, hopf_branch)
end
##!CODEBOUNDARY

#
# ## Calculate stability along the Hopf curves
#
# Here we calculate the stability of the Hopf points on the Hopf branches. Note that we
# need to increase the number of iterations to ensure convergence.
#

##!CODEBOUNDARY @example NeuralMassModel
## calculate stability along the Hopf curves
for hopf_branch in hopf_branches
  stability(jet, hopf_branch, τs, maxit=200)
end
##!CODEBOUNDARY

#
# ## Detect special points
#
# Using the spectral information derived in the previous step we use the function 
# `find_corrected_double_hopf_points` to locates double Hopf points present on the
# Hopf branches. We set the maximum number of iterations to 100 and the tolerance to 1e-15.
# After this we order the double hopf points such the order remains
# the same for each execution of the code.

##!CODEBOUNDARY @example NeuralMassModel
double_hopf_points = PN.find_corrected_double_hopf_points(
  hopf_branches, jet, τs; 
  max_iter=100, tol=1e-15)
double_hopf_points = sort(double_hopf_points, by = x -> x.parameters[2], rev=true)
##!CODEBOUNDARY

#
# ## Detect generalized Hopf points
#
# Similar we use the function `find_genh_points` to locate generalized Hopf points present on the Hopf branches. We set the maximum number of iterations to 100 and the tolerance to 1e-10.

##!CODEBOUNDARY @example NeuralMassModel
genh_points = PN.find_genh_points(
  hopf_branches, jet, τs, par_indx.c;
  max_iter=100, tol=1e-10)
##!CODEBOUNDARY

#
# ## Visualizing Hopf branches with bifurcation points
#
# We creat a figure with the continued Hopf branches together with the several double Hopf and generalized Hopf pionts

##!CODEBOUNDARY @example NeuralMassModel
figHopf = Figure();
ax1 = Axis(figHopf[1, 1], title = "Bifurcation diagram")
lines!(ax1, extract_pars(stst_branch.points), label="Steady-state branch", color=:black)
for (i,hopf_branch) in enumerate(hopf_branches)
  lines!(ax1, extract_pars(hopf_branch.points), label="Hopf branch", color=:royalblue)
end
for (i, p) in enumerate(double_hopf_points)
  scatter!(ax1, [p.parameters[2]], [p.parameters[1]], color=:orangered)
  text!(ax1, p.parameters[2], p.parameters[1], text="HH$i")
end
for (i,p) in enumerate(genh_points)
  scatter!(ax1, [p.parameters[2]], [p.parameters[1]], color=:purple)
  text!(ax1, p.parameters[2], p.parameters[1], text="GH$i")
end
ax1.xlabel = "c"
ax1.ylabel = "a"
axislegend(ax1, merge=true)
figHopf
##!CODEBOUNDARY

#
# ## Continue limit point of cycles emanating from the generalized Hopf points
#
# We continue the limit point of cycles emanating from the generalized Hopf points. We create three branches, one for each generalized Hopf point. We use the function `generalizedHopfToPsol` to create an initial guess for the limit point of cycles branch. We then set up the limit point of cycles branch using the function `SetupLPCBranch` and continue the branch using the function `continue!`.
#

##!CODEBOUNDARY @example NeuralMassModel
ϵ = 0.01
ntst = [20, 40, 20]
ncol = 4
MaxNumberofSteps = 10^6
δmax = 0.1

lpc_branches = map(1:3) do i
    # Create the initial guess from the i-th generalized Hopf point
    lpc_guess = generalizedHopfToPsol(jet, genh_points[i], ϵ, ntst[i], ncol, τs)
    
    # Set up the LPC branch
    branch = SetupLPCBranch(
        jet, lpc_guess, τs;
        parameterbounds = parameterbounds,
        δ = 0.01,
        δmax = δmax,
        MaxNumberofSteps = MaxNumberofSteps,
        NumberOfFails = 2,
        con_par = [par_indx.a; par_indx.c]
    )
    
    # Only reverse the first branch if needed
    if i == 1
        reverse_branch!(branch)
    end
    
    @time continue!(branch)
    branch
end
##!CODEBOUNDARY

#
# ## Plot limit point of cycles branches
#
# Here we plot the three continued limit point of cycles branches. We see that there are multiple cusp point where we should expect the normal form coefficients to change sign.
#

##!CODEBOUNDARY @example NeuralMassModel
fig = Figure();
set_theme!(theme_light())
ax1 = Axis(fig[1, 1])
for hopf_branch in hopf_branches
  lines!(ax1, extract_pars(hopf_branch.points), color=:royalblue, label="Hopf branch")
end
for lpc_branch in lpc_branches
  lines!(ax1, extract_pars(lpc_branch.points), color=:green, label="LPC branch")
end
for (i,p) in enumerate(genh_points)
  scatter!(ax1, [p.parameters[2]], [p.parameters[1]], color=:purple, label="GH")
  text!(ax1, p.parameters[2], p.parameters[1], text="GH$i")
end
axislegend(ax1, merge=:true)
fig
##!CODEBOUNDARY

#
# ## Calculate normal form coefficients
#
# We calculate the normal form coefficients for each limit point of cycles branch.
#

##!CODEBOUNDARY @example NeuralMassModel
## calculate normal form coefficients
for branch in lpc_branches
  normal_form_coefficients!(jet, branch, τs)
end
##!CODEBOUNDARY

#
# ## Plot sign of normal form coefficients along the branches
#
# In this visualization, we compare the signs of the normal form coefficients along each limit-point cycle (LPC) branch, coloring positive coefficients in blue and negative coefficients in orange-red to draw attention to sign changes. We see that the three cusp points at the top of the figure indeed exhibit sign flips, consistent with the presence of cusps. However, there are additional sign changes that do not align with known cusp points. In the subsequent figure, we will show that at those locations, the normal form coefficients pass through infinity, indicating a different type of discontinuity or singularity in the normal form expression.
#

##!CODEBOUNDARY @example NeuralMassModel
fig = Figure();
set_theme!(theme_light())
ax1 = Axis(fig[1, 1], title="Normal form coefficients along the LPC branches")
for hopf_branch in hopf_branches
  lines!(ax1, extract_pars(hopf_branch.points), color=:lightgray, label="Hopf branch")
end
for lpc_branch in lpc_branches
  scatter!(
    ax1,
    extract_pars(lpc_branch.points[get_nmfm_coefficients(lpc_branch).<0.0]),
    markersize=5,
    color=:orangered,
    label="LPC branch negative sign"
  )
  scatter!(
    ax1,
    extract_pars(lpc_branch.points[get_nmfm_coefficients(lpc_branch).>0.0]),
    markersize=5,
    color=:royalblue,
    label="LPC branch positive sign"
  )
end
for (i,p) in enumerate(genh_points)
  scatter!(ax1, [p.parameters[2]], [p.parameters[1]], 
           color=:purple,
           markersize=5,
           label="Generalized Hopf point",)
  text!(ax1, p.parameters[2], p.parameters[1], text="GH$i", color=:black)
end
axislegend(ax1, merge=:true, position=:lb)
fig
##!CODEBOUNDARY

#
# ## Closer inspection of the sign changes
#
# We first collect the normal form coefficients from each limit-point cycle (LPC) branch and apply a signed logarithmic transformation (sign(c)*log(1+|c|)) to make both the sign and magnitude of large coefficients more visually discernible. We then plot these transformed coefficients as colored markers along the parameter axis. By using a diverging color map, values near zero appear in a neutral shade, while large positive or negative values are highlighted in contrasting colors. A color bar indicates the range of transformed coefficient values.
#
# In addition to plotting the Hopf branches as black lines for context, we look for abrupt sign changes in the coefficients—points where they switch from significantly positive to significantly negative (or vice versa). These sign changes can flag potential dynamical transitions in the system. Whenever two consecutive points in a branch satisfy a defined “blow-up” criterion (e.g., opposite signs with large magnitudes), those points are marked with a star symbol to emphasize them. This targeted highlighting helps pinpoint where the model’s dynamics may require closer examination for bifurcation behavior.
#

##!CODEBOUNDARY @example NeuralMassModel
using Statistics  # for reduce/quantile if needed

fig = Figure();
ax = Axis(fig[1, 1], title="Inspection of sign changes on the lpc branches")

# 1) Maybe plot Hopf branches in a simple color:
for hopf_branch in hopf_branches
  lines!(
    ax,
    extract_pars(hopf_branch.points),
    color=:lightgray,
    label="Hopf branch")
end

# 2) Gather all coefficients from your LPC branches
all_coefs = reduce(vcat, [get_nmfm_coefficients(b) for b in lpc_branches])

# Transform them for color scaling
trans_all = sign.(all_coefs) .* log.(1 .+ abs.(all_coefs))

# We'll pick symmetrical color limits around min and max of 'trans_all'
cmin, cmax = minimum(trans_all), maximum(trans_all)
colorlims = (cmin, cmax)

# A diverging colormap from red->white->blue or whichever you like:
my_cmap = cgrad([:orangered, :lightgray, :royalblue], 256)

for lpc_branch in lpc_branches
    coeffs = get_nmfm_coefficients(lpc_branch)
    # Transform each coefficient
    trans = sign.(coeffs) .* log.(1 .+ abs.(coeffs))
    
    scatter!(
        ax,
        extract_pars(lpc_branch.points),
        color = trans,
        colormap = my_cmap,
        colorrange = colorlims,
        markersize = 5
    )
end

Colorbar(fig[1, 2],
         colormap  = my_cmap,
         limits    = colorlims,
         label     = "sign(c)*log(1+|c|)")


for lpc_branch in lpc_branches
    pnts = lpc_branch.points
    coefs = get_nmfm_coefficients(lpc_branch)
    
    for i in 1:(length(pnts)-1)
        c1, c2 = coefs[i], coefs[i+1]
        # Suppose your criterion for "blowing up sign change" is that
        # both are large in absolute value but opposite signs
        if (sign(c1) != sign(c2)) && (abs(c1) > 2) && (abs(c2) > 100)
            # Mark with a star marker
            scatter!(ax, extract_pars([pnts[i+1]]),
                     marker=:star5,
                     color=:purple,
                     markersize=7)
        end
    end
end

for (i,p) in enumerate(genh_points)
  scatter!(ax, [p.parameters[2]], [p.parameters[1]], 
           color=:purple,
           markersize=6,
           label="Generalized Hopf point",)
  text!(ax, p.parameters[2], p.parameters[1], text="GH$i", color=:black)
end
axislegend(ax, merge=:true, position=:lc)

fig
##!CODEBOUNDARY

#
# ## Zoom in on the second limit point of cycles branch
#
# By zooming in on the bottom part we obtain a clearer view of the sign changes in the normal form coefficients.
#

##!CODEBOUNDARY @example NeuralMassModel
xlims!(ax, (0.45, 0.54))
ylims!(ax, (-0.01, 0.14))
fig
##!CODEBOUNDARY

#
# ## Compute normal form coefficients at the double Hopf point
#

##!CODEBOUNDARY @example NeuralMassModel
double_hopf_points = PN.compute_nmfm_coefficients(double_hopf_points, jet, τs)
##!CODEBOUNDARY

#
# ## Continue Neimark-Sacker branches emanating from the first double Hopf points
#

##!CODEBOUNDARY @example NeuralMassModel
ϵ = 0.01
ntst = 30
ncol = 3
nsbranchI = SetupNSBranch(
  jet, double_hopf_points[1], ϵ, ϵ, ntst, ncol, τs;
  parameterbounds=parameterbounds, δ=0.01, δmin=1e-06, δmax=[0.1; 0.1],
  MaxNumberofSteps=[300; 300], con_par = [par_indx.a; par_indx.c]
);
# continue!(nsbranchI[1])
continue!(nsbranchI[2])

ns_guess, psol_branch = PN.hopf_to_psol_to_ns(hopf_branches[2],
                                              double_hopf_points[1], jet, τs,
                                              distance_from_hh = 0.001)
nsbranchI[1] = SetupNSBranch(jet, ns_guess, τs; 
                parameterbounds=parameterbounds, δ=0.1, δmax=0.1,
                            con_par = [par_indx.a; par_indx.c],
                            MaxNumberofSteps=300)
continue!(nsbranchI[1])
##!CODEBOUNDARY

#

##!CODEBOUNDARY @example NeuralMassModel
ϵ1 = 0.01
ϵ2 = 0.01
nsbranchII = SetupNSBranch(
  jet, double_hopf_points[2], ϵ1, ϵ2, ntst, ncol, τs;
  parameterbounds=parameterbounds, δ=0.01, δmin=1e-06, δmax=[0.1; 0.1],
  MaxNumberofSteps=[300; 300], con_par = [par_indx.a; par_indx.c]
);
# continue!(nsbranchII[1])
continue!(nsbranchII[2])

ns_guess, psol_branch = PN.hopf_to_psol_to_ns(hopf_branches[2], double_hopf_points[2],
                                              jet, τs, distance_from_hh = 0.02)
nsbranchII[1] = SetupNSBranch(jet, ns_guess, τs; 
                parameterbounds=parameterbounds, δ=0.1, δmax=0.1,
                            con_par = [par_indx.a; par_indx.c],
                            MaxNumberofSteps=1000)
# reverse_branch!(nsbranchII[1])
continue!(nsbranchII[1])
##!CODEBOUNDARY

#

##!CODEBOUNDARY @example NeuralMassModel
ϵ1 = 0.01
nsbranchIII = SetupNSBranch(
  jet, double_hopf_points[3], ϵ1, ϵ2, ntst, ncol, τs;
  parameterbounds=parameterbounds, δ=0.01, δmin=1e-06, δmax=[0.1; 0.1],
  MaxNumberofSteps=[800; 800], con_par = [par_indx.a; par_indx.c]
);
# continue!(nsbranchIII[1])
continue!(nsbranchIII[2])

ns_guess, psol_branch = PN.hopf_to_psol_to_ns(hopf_branches[2], double_hopf_points[3], jet, τs,
  distance_from_hh = 0.02, sign=+,
  MaxNumberofSteps=800, δ=0.001, δmax=0.001, par_index=2, 
  parameter_before=false)
nsbranchIII[1] = SetupNSBranch(jet, ns_guess, τs;
                parameterbounds=parameterbounds, δ=0.1, δmax=0.1,
                con_par = [par_indx.a; par_indx.c], MaxNumberofSteps=800)
reverse_branch!(nsbranchIII[1])
continue!(nsbranchIII[1])
##!CODEBOUNDARY

#

##!CODEBOUNDARY @example NeuralMassModel
# ## Plot the periodic orbit branch
fig = Figure();
# ax1 = Axis(fig[1,1], title="Periodic orbit branch with bifurcation points")
# for hopf_branch in hopf_branches
#   lines!(ax1, extract_pars(hopf_branch.points), color=:lightgray, label="Hopf branch")
# end
# lines!(ax1, extract_pars(psol_branch.points), label="Periodic orbit branch", color=:green)
# lines!(ax1, extract_pars(nsbranchI[2].points), color=:gray, label="NS branch")
# # add special points
# scatter!(ax1, extract_pars(psol_branch.points[psol_branch.specialpoints.indx_pd]),
#          color=:magenta, marker=:xcross, markersize=5, label="Period doubling")
# scatter!(ax1, extract_pars(psol_branch.points[psol_branch.specialpoints.indx_ns]),
#          color=:purple, marker=:diamond, markersize=5, label="Neimark-Sacker")
# scatter!(ax1, extract_pars(psol_branch.points[psol_branch.specialpoints.indx_fold]),
#          color=:orange, marker=:circle, markersize=5, label="Fold")
# axislegend(ax1, position=:rt, labelsize=8)
# ax1.xlabel = "c"
# ax1.ylabel = "a"
fig
##!CODEBOUNDARY

#
# ## Plot Neimarck-Sacker branches
#

##!CODEBOUNDARY @example NeuralMassModel
fig = Figure();
set_theme!(theme_light())
ax1 = Axis(fig[1, 1], title="Neimark-Sacker branches")
for hopf_branch in hopf_branches
  lines!(ax1, extract_pars(hopf_branch.points), color=:lightgray, label="Hopf branch")
end
for nsbranch in [nsbranchI; nsbranchII; nsbranchIII]
  lines!(ax1, extract_pars(nsbranch.points), color=:green, label="NS branch")
end
for (i,p) in enumerate(double_hopf_points)
  scatter!(ax1, [p.parameters[2]], [p.parameters[1]], color=:purple, label="HH")
  text!(ax1, p.parameters[2], p.parameters[1], text="HH$i")
end
ax1.xlabel = "c"
ax1.ylabel = "a"
axislegend(ax1, merge=:true)
fig
##!CODEBOUNDARY

#
# ## Compute normal form coefficients for the Neimark-Sacker branches
#

##!CODEBOUNDARY @example NeuralMassModel
## calculate normal form coefficients
for branch in [nsbranchI; nsbranchII; nsbranchIII]
  normal_form_coefficients!(jet, branch, τs)
end
##!CODEBOUNDARY

#
# ## Show normal form coefficients at the double Hopf point
#
# Based on the critical normal forms of the double Hopf points, we determine the relevant unfolding case. We find that all of our identified points fall into the simple case I (following the classification in [Yuri]). However, in the first and third points, a time-reversal element must be considered, indicating that the Neimark–Sacker branches emanating from these points should exhibit stable tori.
#
# !!! info
#
# It is important to note that the stability of these tori pertains only to the three-dimensional center manifold of the Neimark-Sacker points. If there are unstable eigenvalues at the double Hopf point, leading to unstable multipliers at the emanating Neimark-Sacker points, the tori emanating from the Neimark-Sacker points will be unstable within the four-dimensional center manifold of the double Hopf point. As a result, forward integration of the full system will not allow simulation of these tori.
#

##!CODEBOUNDARY @example NeuralMassModel
## calculate normal form coefficients
PN.determine_unfolding_case(double_hopf_points)
##!CODEBOUNDARY

#
# ## Plot sign of normal form coefficients along the branches
#
# In this visualization, we compare the signs of the real part of critical normal form coefficients along each Neimark-Sacker (NS) branch, coloring positive coefficients in blue and negative coefficients in orange-red to draw attention to sign changes. We see that the first and third have negative normal form coefficients along the emanating Neimark-Sacker curves, while the fifth double Hopf point has positive normal form coefficients along the emanating Neimark-Sacker curves. This indeed corresponds to the analysis as presented in [Yuri], taking into account the time-reversal for the first and third double Hopf point. At the fifth double Hopf point there is no time-reversal to take into account, therefore, we expect positive critical normal form coefficients along the emanating Neimark-Sacker curvers. We see below that this is indeed the case.
#
# Lastly, we also see various sign changes along the Neimark-Sacker branches. We will return to these points in further research when analysis higher codimension periodic bifurcations.
#

##!CODEBOUNDARY @example NeuralMassModel
fig = Figure();
set_theme!(theme_light())
ax1 = Axis(fig[1, 1], title="Normal form coefficients along the Neimark-Sacker branches")
for hopf_branch in hopf_branches
  lines!(ax1, extract_pars(hopf_branch.points), color=:lightgray, label="Hopf branch")
end
for branch in [nsbranchI; nsbranchII; nsbranchIII]
  if ~isempty(branch.points[real.(get_nmfm_coefficients(branch)).<0.0])
    scatter!(
      ax1,
      extract_pars(branch.points[real.(get_nmfm_coefficients(branch)).<0.0]),
      markersize=4,
      color=:orangered,
      label="NS branch negative sign"
    )
  end
  if ~isempty(branch.points[real.(get_nmfm_coefficients(branch)).>0.0])
    scatter!(
      ax1,
      extract_pars(branch.points[real.(get_nmfm_coefficients(branch)).>0.0]),
      markersize=4,
      color=:royalblue,
      label="NS branch positive sign"
    )
  end
end
for (i,p) in enumerate(double_hopf_points)
  scatter!(ax1, [p.parameters[2]], [p.parameters[1]], 
           color=:purple,
           markersize=5,
           label="double Hopf point",)
  text!(ax1, p.parameters[2], p.parameters[1], text="HH$i", color=:black)
end
ax1.xlabel = "c"
ax1.ylabel = "a"
axislegend(ax1, merge=:true, position=:lb)
fig
##!CODEBOUNDARY

#
# ## Compute multipliers on the Neimark-Sacker branches
#
# First we calculate the multipliers along the branch with the function `multiples!`. Then we detect the special points with the function `detectSpecialPoints`.
#

##!CODEBOUNDARY @example NeuralMassModel
for nsbranch in [nsbranchI; nsbranchII]
  multipliers!(jet, nsbranch, τs)
end
##!CODEBOUNDARY

#

##!CODEBOUNDARY @setup NeuralMassModel
if ENV["INTERACTIVE"] == "true"
  using WGLMakie
  using Bonito, Markdown
  using Bonito:Slider
  App() do session::Session
      index_slider = Slider(1:length(nsbranchII[1].points))
      fig = Figure(; size=(800, 600))
      ax1 = Axis(fig[1, 1], title="LPC branch")
      ax2 = Axis(fig[1, 2], title="Multipliers")
      ax3 = Axis(fig[2, 1], title="Profile")
      ax4 = Axis(fig[2, 2], title="Profile")
      psol_param_max_profile_data = hcat(
        [p.parameters[2] for p in nsbranchII[1].points],
        [maximum(reduce(hcat, p.profile)[1, :]) for p in nsbranchII[1].points])
      scatterlines!(ax1, psol_param_max_profile_data, label="max x1 profile nsbranchII[1]")
      ax1.xlabel = "c"
      ax1.ylabel = L"\max x_1"
      lines!(ax2, cos.(2 * pi * range(0, 1, 100)), sin.(2 * pi * range(0, 1, 100)))
      # plot multipliers
      μ = nsbranchII[1].points[1].stability         # hide
      μ_psol = Observable(Point2f.(vec(real(μ)), vec(imag(μ))))         # hide
      scatter!(ax2, μ_psol)         # hide
      # active point
      idx = 1
      abscissa = nsbranchII[1].points[idx].parameters[2]
      ordinate = maximum(reduce(hcat, nsbranchII[1].points[idx].profile)[1, :])
      points = Observable(Point2f([abscissa; ordinate]))         # hide
      scatter!(ax1, points, color = :black)         # hide
      # xlims!(ax2, (-1.1, -0.9))
      abscissa_profile = nsbranchII[1].points[idx].period*nsbranchII[1].points[idx].mesh
      ordinate_profile = hcat(nsbranchII[1].points[idx].profile...)[1,:]        # hide
      points_profile = Observable(Point2f.(
          [r for r in eachrow([abscissa_profile ordinate_profile])]))         # hide
      scatter!(ax3, points_profile)         # hide
      on(index_slider) do idx         # hide
        abscissa = nsbranchII[1].points[idx].parameters[2]
        ordinate = maximum(reduce(hcat, nsbranchII[1].points[idx].profile)[1, :])
        points[] = Point2f([abscissa; ordinate])         # hide
        notify(points)         # hide
        abscissa_profile = nsbranchII[1].points[idx].period*nsbranchII[1].points[idx].mesh
        ordinate_profile = hcat(nsbranchII[1].points[idx].profile...)[1,:]        # hide
        points_profile[] = Point2f.(
          [r for r in eachrow([abscissa_profile ordinate_profile])])         # hide
        notify(points_profile)         # hide
        μ = nsbranchII[1].points[idx].stability         # hide
        μ_psol[] = Point2f.(vec(real(μ)), vec(imag(μ)))         # hide
        notify(μ_psol)         # hide
      end         # hide
      slider = DOM.div("step: ", index_slider, index_slider.value)         # hide
      return Bonito.record_states(session, DOM.div(slider, fig))         # hide
  end
end
##!CODEBOUNDARY

#
# ## In pursuit of period-doubling orbits.
#
# In order to find a branch of period-doubling orbits we select
# a Hopf point near c ≈ 0.7 hopf_branchI and continue the emenating
# branch of periodic orbits in the parameter $c$.
#

##!CODEBOUNDARY @example NeuralMassModel
indx = findlast(p -> 0.69 < p.parameters[2] < 0.7, hopf_branches[1].points)
hopf_point = hopf_branches[1].points[indx]

# continu psol branch emanating from hopf point
psol_branch = SetupPsolBranch(jet, con_par, hopf_point, τs; 
                              parameterbounds=parameterbounds, 
                              MaxNumberofSteps=200, ntst=30, ncol=3, δ=0.001, δmax=0.1)
reverse_branch!(psol_branch)
continue!(psol_branch)
##!CODEBOUNDARY

#
# ## Plot the periodic orbit branch
#
# Here we display the periodic orbit branch emanating from the selected Hopf point. We also show the Hopf branches for context.

##!CODEBOUNDARY @example NeuralMassModel
fig = Figure();
ax1 = Axis(fig[1,1])
for hopf_branch in hopf_branches
  lines!(ax1, extract_pars(hopf_branch.points), color=:lightgray, label="Hopf branch")
end
lines!(ax1, extract_pars(psol_branch.points), label="Periodic orbit branch", color=:green)
ax1.xlabel = "c"
ax1.ylabel = "a"
axislegend(ax1, merge=true)
fig
##!CODEBOUNDARY

#
# ## Detect bifurcations along the periodic orbit branch
#
# First we calculate the multipliers along the branch with the function `multiples!`. Then we detect the special points with the function `detectSpecialPoints`.
#

##!CODEBOUNDARY @example NeuralMassModel
multipliers!(jet, psol_branch, τs)
psol_branch = detectSpecialPoints(psol_branch)
##!CODEBOUNDARY

#
# ## Plot the periodic orbit branch
#
# We plot the again the periodic orbit branch, but now with the special points marked. We see that the branch contains fold, period-doubling and Neimark-Sacker bifurcations. Additionally, we plot one of the Neimark-Sacker branches emanating from the first double Hopf point and see that it goes through a Neimark-Sacker bifurcation on the periodc orbit branch.

##!CODEBOUNDARY @example NeuralMassModel
fig = Figure();
ax1 = Axis(fig[1,1], title="Periodic orbit branch with bifurcation points")
for hopf_branch in hopf_branches
  lines!(ax1, extract_pars(hopf_branch.points), color=:lightgray, label="Hopf branch")
end
lines!(ax1, extract_pars(psol_branch.points), label="Periodic orbit branch", color=:green)
lines!(ax1, extract_pars(nsbranchI[2].points), color=:gray, label="NS branch")
# add special points
scatter!(ax1, extract_pars(psol_branch.points[psol_branch.specialpoints.indx_pd]),
         color=:magenta, marker=:xcross, markersize=5, label="Period doubling")
scatter!(ax1, extract_pars(psol_branch.points[psol_branch.specialpoints.indx_ns]),
         color=:purple, marker=:diamond, markersize=5, label="Neimark-Sacker")
scatter!(ax1, extract_pars(psol_branch.points[psol_branch.specialpoints.indx_fold]),
         color=:orange, marker=:circle, markersize=5, label="Fold")
axislegend(ax1, position=:rt, labelsize=8)
ax1.xlabel = "c"
ax1.ylabel = "a"
fig
##!CODEBOUNDARY

#
# ## Inspect multipliers on psol branch
#
# Here we persent another visalization of the periodic solution branch along with its special points. In the first subplot, we plot the relationship between the parameter \( c \) and the maximum value of the first state variable, \( \max x_1 \), along the branch. The detected special points folds, period-doubling bifurcations, and Neimark-Sacker bifurcations are marked.
#
# In the second subplot, we show the multipliers at the Neimark-Sacker point, in the complex plane. The closest pair of multiplier to the unit circle is highlighted in red.
#
# The last two subplots show the time profiles of the first and second state variables, \( x_1 \) and \( x_2 \), at the Neimark-Sacker point.
#

##!CODEBOUNDARY @example NeuralMassModel
fig = Figure();
ax1 = Axis(fig[1, 1], title="Psol branch")
ax2 = Axis(fig[1, 2], title="Multipliers at Neimark-Sacker point")
ax3 = Axis(fig[2, 1], title="Profile")
ax4 = Axis(fig[2, 2], title="Profile")
psol_param_max_profile_data = hcat(
  [p.parameters[2] for p in psol_branch.points],
  [maximum(reduce(hcat, p.profile)[1, :]) for p in psol_branch.points])
lines!(ax1, psol_param_max_profile_data)
ax1.xlabel = L"c"
ax1.ylabel = L"\max x_1"
lines!(ax2, cos.(2 * pi * range(0, 1, 100)), sin.(2 * pi * range(0, 1, 100)))
# plot multipliers at Neimark-Sacker point
ns_point = psol_branch.points[psol_branch.specialpoints.indx_ns[1]]
μ = ns_point.stability
μ_psol = Point2f.(vec(real(μ)), vec(imag(μ)))
scatter!(ax2, μ_psol, markersize=5)
# High light the closest multiplier to the unit circle
# Calculate distances from the unit circle
distances = abs.(abs.(μ) .- 1)
# Find the index of the closest multiplier
closest_index = argmin(distances)
# Extract the closest multiplier
closest_multiplier = Point2f(real(μ[closest_index]), imag(μ[closest_index]))
closest_multiplier_conj = Point2f(real(μ[closest_index]), -imag(μ[closest_index]))
scatter!(ax2, closest_multiplier, color=:red, markersize=5)
scatter!(ax2, closest_multiplier_conj, color=:red, markersize=5)
ax2.xlabel = "Re(μ)"
ax2.ylabel = "Im(μ)"
# add special points
scatter!(ax1, [p.parameters[2] for p in psol_branch.points[psol_branch.specialpoints.indx_fold]],
  [maximum(reduce(hcat, p.profile)[1, :]) for p in psol_branch.points[psol_branch.specialpoints.indx_fold]],
  color=:orange, marker=:circle, markersize=8, label="Fold")
scatter!(ax1, [p.parameters[2] for p in psol_branch.points[psol_branch.specialpoints.indx_pd]],
  [maximum(reduce(hcat, p.profile)[1, :]) for p in psol_branch.points[psol_branch.specialpoints.indx_pd]],
  color=:magenta, marker=:xcross, markersize=8, label="Period doubling")
scatter!(ax1, [p.parameters[2] for p in psol_branch.points[psol_branch.specialpoints.indx_ns]],
  [maximum(reduce(hcat, p.profile)[1, :]) for p in psol_branch.points[psol_branch.specialpoints.indx_ns]],
    color=:purple, marker=:diamond, markersize=8, label="Neimark-Sacker")
axislegend(ax1, position=:rt, labelsize=8)
# plot x_1 profile at Neimark-Sacker point
lines!(ax3, ns_point.period*ns_point.mesh, hcat(ns_point.profile...)[1,:])
lines!(ax4, ns_point.period*ns_point.mesh, hcat(ns_point.profile...)[2,:])
ax3.xlabel = L"t"
ax3.ylabel = L"x_1"
ax4.xlabel = L"t"
ax4.ylabel = L"x_2"
fig
##!CODEBOUNDARY

#
# ## Continuation of Limit Point of Cycles (LPC)
#

##!CODEBOUNDARY @setup NeuralMassModel
if ENV["INTERACTIVE"] == "true"
  lpc_guess = psol_branch.points[psol_branch.specialpoints.indx_fold[1]]
  lpc_branchI = SetupLPCBranch(jet, lpc_guess, τs; 
                  parameterbounds=parameterbounds, δ=0.1, δmax=0.1,
                              con_par = [par_indx.a; par_indx.c],
                              MaxNumberofSteps=1000)
  continue!(lpc_branchI)
  reverse_branch!(lpc_branchI)
  continue!(lpc_branchI)

  # plot parameters
  fig = Figure();
  ax1 = Axis(fig[1, 1], xlabel="c", ylabel="a", title="LPC branch")
  lines!(ax1, extract_pars(lpc_branchI.points))
  # save("lpc_branch.png", fig)
  fig
end
##!CODEBOUNDARY

#
# ## Computation of normal form coefficients for LPC branch
#

##!CODEBOUNDARY @setup NeuralMassModel
if ENV["INTERACTIVE"] == "true"
normal_form_coefficients!(jet, lpc_branchI, τs)
end
##!CODEBOUNDARY

#
# # Compuation of multipliers for the LPC branch
#

##!CODEBOUNDARY @setup NeuralMassModel
if ENV["INTERACTIVE"] == "true"
# compute multipliers
multipliers!(jet, lpc_branchI, τs)
end
##!CODEBOUNDARY

#
# ## Figure showing normal form coefficient along the parameter bifurcation curve
#
# Below is a figure showing the sign of the normal form coefficients along the LPC branch. The orange points are where the normal form coefficients are negative, and the blue points are where the normal form coefficients are positive. We see that the normal form coefficients change sign at the cusp bifurcation point in the top left corner of the plot. However, see also see a change of sign between the two cusp points and a missing a change of sign at the cusp point in the lower left corner. To see why this happens we will visualize the normal form coefficients in 3D in the next figure.
#

##!CODEBOUNDARY @setup NeuralMassModel
if ENV["INTERACTIVE"] == "true"
# compute multipliers
# using GLMakie
# CairoMakie.activate!()
# GLMakie.activate!()
# 2D Plot
# CairoMakie.activate!()
fig2D = Figure()
ax2 = Axis(fig2D[1,1], xlabel="c", ylabel="a", title="LPC branch (2D)")
lines!(ax2, [extract_pars(lpc_branchI.points); get_nmfm_coefficients(lpc_branchI)'])
scatter!(ax2, extract_pars(lpc_branchI.points[get_nmfm_coefficients(lpc_branchI).<0.0]), label="Negative normal form coefficients", markersize=5, color=:orangered)
scatter!(ax2, extract_pars(lpc_branchI.points[get_nmfm_coefficients(lpc_branchI).>0.0]), label="Positive normal form coefficients", markersize=5, color=:royalblue)
axislegend(ax2)
save("lpc_branch_2D.png", fig2D)
fig2D
end
##!CODEBOUNDARY

#
# ## Figure showing normal form coefficient along the parameter bifurcation curve
#
# Here we will show the normal form coefficients in 3D. Again, we see, as expected, that the normal form coefficients change sign at the cusp bifurcation point in the top left corner of the plot. We also see that the sign change in the middle of the two bifurcations points is due to the normal form coefficient goes to infinity at that point. Indicating a degenerated point along the fold curve. Secondly, the cusp in the lower left corner is degenerated as well. Inspection of the multipliers will demonstrate that the point have additional point on the unit circle at -1.
#

##!CODEBOUNDARY @setup NeuralMassModel
if ENV["INTERACTIVE"] == "true"
# compute multipliers
# 3D Plot
fig3D = Figure(;fontsize=8)
ax1 = Axis3(fig3D[1, 1], xlabel="c", ylabel="a",
      zlabel="normal form coefficient c", title="LPC branch (3D)")
lines!([extract_pars(lpc_branchI.points); get_nmfm_coefficients(lpc_branchI)'])
scatter!(extract_pars(lpc_branchI.points[get_nmfm_coefficients(lpc_branchI).<0.0]), markersize=5)
scatter!(extract_pars(lpc_branchI.points[get_nmfm_coefficients(lpc_branchI).>0.0]), markersize=5)
ax1.azimuth = 4.040000000000005
ax1.elevation = 0.26020352061721297
zlims!(ax1, (-1, 1))
save("lpc_branch_3D.png", fig3D)
fig3D
end
##!CODEBOUNDARY

#
# ## Plot orbits along the LPC branch
#

##!CODEBOUNDARY @setup NeuralMassModel
if ENV["INTERACTIVE"] == "true"
figHopf = Figure();
ax1 = Axis(figHopf[1, 1], title = "Bifurcation diagram")
lines!(ax1, extract_pars(stst_branch.points), label="Steady-state branch", color=:black)
for (i,hopf_branch) in enumerate(hopf_branches)
  lines!(ax1, extract_pars(hopf_branch.points), label="Hopf branch $i", color=:royalblue)
end
for (i, p) in enumerate(double_hopf_points)
  scatter!(ax1, [p.parameters[2]], [p.parameters[1]], color=:orangered, label="HH")
  text!(ax1, p.parameters[2], p.parameters[1], text="HH$i")
end
for (i,p) in enumerate(genh_points)
  scatter!(ax1, [p.parameters[2]], [p.parameters[1]], color=:purple, label="GH")
  text!(ax1, p.parameters[2], p.parameters[1], text="GH$i")
end
# xlims!(ax1, (0.0, 1.2))
# ylims!(ax1, (-1.0, 1.0))
ax1.xlabel = "c"
ax1.ylabel = "a"
axislegend(ax1, merge=true)
lines!(ax1, extract_pars(nsbranch[1].points), color=:forestgreen)
lines!(ax1, extract_pars(nsbranch[2].points), color=:forestgreen)
figHopf
end
##!CODEBOUNDARY

#

##!CODEBOUNDARY @setup NeuralMassModel
if ENV["INTERACTIVE"] == "true"
  fig = Figure()
  ax1 = Axis3(fig[1, 1], title="LPC branch")
  for p in lpc_branchI.points
    lines!(ax1, p.parameters[1]*ones(length(p.mesh)), p.period*p.mesh, hcat(p.profile...)[1,:],
          color=:royalblue, linewidth=0.5)
  end
end
##!CODEBOUNDARY

#

##!CODEBOUNDARY @setup NeuralMassModel
if ENV["INTERACTIVE"] == "true"
  using WGLMakie
  using Bonito, Markdown
  using Bonito:Slider
  App() do session::Session
      index_slider = Slider(1:length(lpc_branchI.points))
      fig = Figure(; size=(800, 600))
      ax1 = Axis(fig[1, 1], title="LPC branch")
      ax2 = Axis(fig[1, 2], title="Multipliers")
      ax3 = Axis(fig[2, 1], title="Profile")
      ax4 = Axis(fig[2, 2], title="Profile")
      psol_param_max_profile_data = hcat(
        [p.parameters[2] for p in lpc_branchI.points],
        [maximum(reduce(hcat, p.profile)[1, :]) for p in lpc_branchI.points])
      scatterlines!(ax1, psol_param_max_profile_data, label="max x1 profile lpc_branchI")
      ax1.xlabel = "c"
      ax1.ylabel = L"\max x_1"
      lines!(ax2, cos.(2 * pi * range(0, 1, 100)), sin.(2 * pi * range(0, 1, 100)))
      # plot multipliers
      μ = lpc_branchI.points[1].stability         # hide
      μ_psol = Observable(Point2f.(vec(real(μ)), vec(imag(μ))))         # hide
      scatter!(ax2, μ_psol)         # hide
      # active point
      idx = 1
      abscissa = lpc_branchI.points[idx].parameters[2]
      ordinate = maximum(reduce(hcat, lpc_branchI.points[idx].profile)[1, :])
      points = Observable(Point2f([abscissa; ordinate]))         # hide
      scatter!(ax1, points, color = :black)         # hide
      xlims!(ax2, (-1.1, -0.9))
      abscissa_profile = lpc_branchI.points[idx].period*lpc_branchI.points[idx].mesh
      ordinate_profile = hcat(lpc_branchI.points[idx].profile...)[1,:]        # hide
      points_profile = Observable(Point2f.(
          [r for r in eachrow([abscissa_profile ordinate_profile])]))         # hide
      scatter!(ax3, points_profile)         # hide
      on(index_slider) do idx         # hide
        abscissa = lpc_branchI.points[idx].parameters[2]
        ordinate = maximum(reduce(hcat, lpc_branchI.points[idx].profile)[1, :])
        points[] = Point2f([abscissa; ordinate])         # hide
        notify(points)         # hide
        abscissa_profile = lpc_branchI.points[idx].period*lpc_branchI.points[idx].mesh
        ordinate_profile = hcat(lpc_branchI.points[idx].profile...)[1,:]        # hide
        points_profile[] = Point2f.(
          [r for r in eachrow([abscissa_profile ordinate_profile])])         # hide
        notify(points_profile)         # hide
        μ = lpc_branchI.points[idx].stability         # hide
        μ_psol[] = Point2f.(vec(real(μ)), vec(imag(μ)))         # hide
        notify(μ_psol)         # hide
      end         # hide
      slider = DOM.div("step: ", index_slider, index_slider.value)         # hide
      return Bonito.record_states(session, DOM.div(slider, fig))         # hide
  end
end
##!CODEBOUNDARY

#
# ## Continuation of Period Doubling (PD)
#
# Next we continue the period doubling branch emanating from the period doubling point on the periodic orbit branch. We will continue the branch in the parameters \( a \) and \( c \).

##!CODEBOUNDARY @example NeuralMassModel
pd_guess = psol_branch.points[psol_branch.specialpoints.indx_pd][1]
pd_branchI = SetupPDBranch(
  jet,
  pd_guess, 
  τs;
  parameterbounds=parameterbounds,
  δ=0.001, δmax=0.1,
  MaxNumberofSteps=500,
  con_par = [par_indx.a; par_indx.c])
continue!(pd_branchI)
reverse_branch!(pd_branchI)
continue!(pd_branchI)
##!CODEBOUNDARY

#
# ## Visaualize the period doubling branch
#
# Here we show the continued period doubling branch emanating from the period doubling point on the periodic orbit branch. We also show the periodic orbit branch for context. We see that the two detected period doubling points on the period orbit branch are connected by the period doubling branch.

##!CODEBOUNDARY @example NeuralMassModel
fig2D = Figure();
ax1 = Axis(fig2D[1,1], xlabel="c", ylabel="a", title="Period doubling branch (2D)")
lines!(ax1, extract_pars(psol_branch.points), label="Periodic orbit branch", color=:green)
# add special points
axislegend(ax1, position=:rt, labelsize=8)
lines!(ax1, extract_pars(pd_branchI.points))
scatter!(ax1, extract_pars(psol_branch.points[psol_branch.specialpoints.indx_pd]),
         color=:magenta, marker=:xcross, markersize=5, label="Period doubling")
fig2D
##!CODEBOUNDARY

#
# ## Computation of normal form coefficients for period doubling branch
#
# Next, we compute the normal form coefficients for the period doubling branch.

##!CODEBOUNDARY @example NeuralMassModel
normal_form_coefficients!(jet, pd_branchI, τs)
##!CODEBOUNDARY

#
# ## Computaion of multipliers for period doubling branch
#

##!CODEBOUNDARY @example NeuralMassModel
# compute multipliers
multipliers!(jet, pd_branchI, τs)
##!CODEBOUNDARY

#
# ## Figure showing normal form coefficient along the parameter bifurcation curve
#
# Here we show the normal form coefficients in 2D. The orangered points are where the normal form coefficients are negative, and the royalblue points are where the normal form coefficients are positive. Additionally, we plot the Neimark-Sacker branch emanating from the first period doubling point. We see that the Neimark-Sacker branch touches the period doubling branch at a sign change in the normal form coefficient.

##!CODEBOUNDARY @example NeuralMassModel
fig2D = Figure();
ax1 = Axis(fig2D[1,1], xlabel="c", ylabel="a", title="Period doubling branch (2D)")
lines!(ax1, extract_pars(pd_branchI.points))
scatter!(ax1, extract_pars(pd_branchI.points[real(get_nmfm_coefficients(pd_branchI)).<0.0]), label="Negative normal form coefficients", markersize=8, color=:orangered)
scatter!(ax1, extract_pars(pd_branchI.points[real(get_nmfm_coefficients(pd_branchI)).>0.0]), label="Positive normal form coefficients", markersize=8, color=:royalblue)
lines!(ax1, extract_pars(nsbranchI[2].points), color=:gray, label="NS branch")
axislegend(ax1)
fig2D
##!CODEBOUNDARY

#
# ## Figure showing normal form coefficient along the parameter bifurcation curve in 3d
#
# By creating a three dimensional visualization of the normal form coefficients along the period doubling branch, we can better understand the bifurcation structure. We see that there are period doubling points on the branch where the sign of the normal form coefficient vanishes. Although, we see a sign changes on the period-doubling branch near \(a \approx 0 \) these are due to round off errors. We also see that the normal form coefficient goes to infinity at the point where the Neimark-Sacker branch touches the period doubling branch.

##!CODEBOUNDARY @example NeuralMassModel
fig3D = Figure(;fontsize=8);
ax1 = Axis3(fig3D[1, 1], xlabel="c", ylabel="a", 
      zlabel="normal form coefficient c", title="Period doubling branch (3D)")
lines!([extract_pars(pd_branchI.points); real(get_nmfm_coefficients(pd_branchI))'])
scatter!(extract_pars(pd_branchI.points[real(get_nmfm_coefficients(pd_branchI)).<0.0]), markersize=5)
scatter!(extract_pars(pd_branchI.points[real(get_nmfm_coefficients(pd_branchI)).>0.0]), markersize=5)
ax1.azimuth = 4.040000000000005
ax1.elevation = 0.26020352061721297
zlims!(ax1, (-0.006, 0.006))
fig3D
##!CODEBOUNDARY

#
# ## Generalized Period Doubling (GPD)
#
# We locate the two points on the period doubling branch where the sign of the critical normal form coefficient at the period-doubling bifurcation vanishes. These points are potential candidates for generalized period-doubling points. We then filter the points to only include those with a critical normal form coefficient where the parameter \( c \) is greater than 0.55.

##!CODEBOUNDARY @example NeuralMassModel
pd_branchI = LocateGPDPoints(pd_branchI)
# filter GPD with parameter c > 0.5
gpd_indx = pd_branchI.specialpoints.gpd_indx
gpd_indx = [i for (i, p) in enumerate(pd_branchI.points[gpd_indx]) if p.parameters[end] > 0.55]
gpd_indx = pd_branchI.specialpoints.gpd_indx[gpd_indx]
gpd_points = pd_branchI.points[gpd_indx.+1]
##!CODEBOUNDARY

#
# ## In pursuit of a limit point a cycle (LPC) branch eamanating from a GPD point
#
# Start from an existing period-doubling (PD) branch (`pd_branchI`) that contains a known period-doubling bifurcation, and extract a particular periodic solution (`psol1`). This periodic solution will serve as the initial point for further continuation steps.
#
# **Explanation:**
#
# - Access a solution near a generalized period-doubling point using `gpd_indx`.
# - After extracting the solution, plot it alongside known points (`gpd_points`) and the PD branch to verify its location in parameter space.
#

##!CODEBOUNDARY @example NeuralMassModel
# select psol1 from the PD branch in the middle of the two GPD points
psol_indx = if gpd_indx[2] > gpd_indx[1]
    (gpd_indx[2] + gpd_indx[1]) ÷ 2
else
    (gpd_indx[1] + gpd_indx[2]) ÷ 2
end
psol1 = [pd_branchI.points[psol_indx]]

# # psol with double period and
psol1_double_period = psol1[1]
psol1_double_period = @set psol1_double_period.period *= 2
# # lines(hcat(psol1_double_period.profile ...)[1,:])
psol1_double_period = @set psol1_double_period.profile = [psol1_double_period.profile; psol1_double_period.profile[2:end]]
psol1_double_period = @set psol1_double_period.profile = psol1_double_period.profile[1:2:end]

ax1 = fig2D.content[1]
scatter!(ax1, extract_pars(gpd_points), color=:red)
scatter!(ax1, extract_pars(psol1), color=:magenta)
# save("pd_branch_2D_2.png", fig2D)
fig2D

##!CODEBOUNDARY

#
# ## Setting Up and Continuing a New Periodic Solution Branch
#
# **Objective:**
# With `psol1` selected, set up and continue a new periodic solution branch (`psol_branchII`). This allows tracing out how solutions evolve as parameters vary.
#
# **Explanation:**
#
# - `SetupPsolBranch` initializes a branch starting from `psol1`.
# - After specifying parameter bounds and step sizes, `continue!` follows the branch through parameter space.
# - Visualizing the solution amplitude or parameter variation helps understand how solutions change along the branch.
#

##!CODEBOUNDARY @example NeuralMassModel
psol_branchII = SetupPsolBranch(jet, con_par, psol1[1], τs; parameterbounds=parameterbounds, MaxNumberofSteps=200, δ=0.001, δmax=0.01)
# psol_branchII = SetupPsolBranch(jet, con_par, psol1_double_period, τs;
# parameterbounds=parameterbounds, MaxNumberofSteps=100, δ=0.01, δmax=0.1)
# continue!(psol_branchII)
reverse_branch!(psol_branchII)
continue!(psol_branchII)
##!CODEBOUNDARY

#
# ## Detecting Special Points
#

##!CODEBOUNDARY @example NeuralMassModel
multipliers!(jet, psol_branchII, τs)
psol_branchII = detectSpecialPoints(psol_branchII)
##!CODEBOUNDARY

#
# ## Plot
#

##!CODEBOUNDARY @example NeuralMassModel
fig2D = Figure();
ax2 = Axis(fig2D[1,1], xlabel="c", ylabel="a", title="Period doubling branch (2D)")
lines!(ax2, extract_pars(psol_branchII.points), color=:royalblue)
scatter!(ax2, extract_pars(psol_branchII.points[psol_branchII.specialpoints.indx_fold]), color=:purple)
lines!(ax2, extract_pars(pd_branchI.points))
fig2D
##!CODEBOUNDARY

##!CODEBOUNDARY @example NeuralMassModel
fig = Figure();
ax1 = Axis(fig[1, 1], title="Psol branch II")
lines!(ax1,
    [p.parameters[2] for p in psol_branchII.points],
    [maximum(hcat(p.profile...)[1, :]) for p in psol_branchII.points]
)
scatter!(ax1, [p.parameters[2] for p in psol_branchII.points[psol_branchII.specialpoints.indx_fold]],
  [maximum(reduce(hcat, p.profile)[1, :]) for p in psol_branchII.points[psol_branchII.specialpoints.indx_fold]],
  color=:orange, marker=:circle, markersize=8, label="Fold")
axislegend(ax1, position=:rt, labelsize=8)
ax1.xlabel = "c"
ax1.ylabel = "a"
# xlims!(ax1, (0.68, 0.70))
# ylims!(ax1, (0.25, 0.4))
# save("psol_branchII.png", fig)
fig
##!CODEBOUNDARY

#
# ## Detecting Special Points and Preparing for LPC Branch Computation
#
# **Objective:**
# Once `psol_branchII` is obtained, detect special bifurcation points. Identifying these points is crucial for further analysis, such as setting up a Limit Point of Cycle (LPC) branch from a fold point.
#
# **Explanation:**
#
# - `multipliers!` computes Floquet multipliers for stability analysis.
# - `detectSpecialPoints` identifies bifurcation points (e.g., folds).
# - Using a detected fold, we set up an LPC branch to investigate deeper bifurcation structures.
#

##!CODEBOUNDARY @example NeuralMassModel
fig = Figure();
ax = Axis(fig[1, 1], title="Psol branch II")
lines!(extract_pars(pd_branchI.points), color=:royalblue)
lines!(extract_pars(psol_branchII.points), color=:royalblue)
scatter!(extract_pars(gpd_points), color=:red)
scatter!(extract_pars(psol1), color=:lightgreen)
fig
##!CODEBOUNDARY

#
# ## Continue lpc branch with double period
#
# Here we continue the lpc branch in two parameters, \( a \) and \( c \), starting from the initial point `lpc_guess`. We set the maximum number of steps to 10 and the continuation parameter to \( a \) and \( c \). We also set a callback function to stop the continuation when the distance between the generalized period doubling points and solution is less than \(10^{-4}\). We see that the two generalized period doubling points are quite beautifully connected by the lpc branch.

##!CODEBOUNDARY @example NeuralMassModel
lpc_guess = psol_branchII.points[psol_branchII.specialpoints.indx_fold[end]]
lpc_branchII = SetupLPCBranch(jet, lpc_guess, τs;
    # parameterbounds=(min = [0.0, 0.0], max = [1.2, gpd_points[1].parameters[2]]),
    δ=0.01, δmax=0.1, MaxNumberofSteps=600,
    con_par=[par_indx.a; par_indx.c])

function callback_lpc_branchII(x, gpd_points)
    distance1 = norm(gpd_points[1].parameters - x.parameters)
    distance2 = norm(gpd_points[2].parameters - x.parameters)
    if distance1 < 10^-4 || distance2 < 10^-4
        println("Callback: distance between GPD and solution is less than 1e-4")
        return true
    end
    return false
end

continue!(lpc_branchII,callback = x -> callback_lpc_branchII(x, gpd_points))
reverse_branch!(lpc_branchII)
continue!(lpc_branchII,callback = x -> callback_lpc_branchII(x, gpd_points))

fig = Figure();
ax1 = Axis(fig[1, 1], title="LPC branch emanating from, and connecting, the two GPD points")
lines!(ax1, extract_pars(pd_branchI.points), linewidth=2)
lines!(ax1, extract_pars(lpc_branchII.points), linewidth=2)
scatter!(extract_pars(gpd_points), color=:red)
fig

##!CODEBOUNDARY

#
# ## Normal Form Coefficients on LPC Branch
#
# We compute the normal form coefficients for the LPC branch.
#

##!CODEBOUNDARY @example NeuralMassModel
# Optional filtering:
normal_form_coefficients!(jet, lpc_branchII, τs)
##!CODEBOUNDARY

#
# ## Visualizing Results in 3D Parameter Space
#
# **Objective:**
# Visualize the bifurcation structures and normal form coefficients in three-dimensional parameter space. Such visualization provides a comprehensive understanding of how parameters, bifurcation points, and normal form coefficients interact.
#
# **Explanation:**
#
# - Plotting in 3D shows how normal form coefficients vary over a parameter surface.
# - By overlaying branches, fold points, and special points, we gain a clearer picture of the overall bifurcation landscape.
#

##!CODEBOUNDARY @example NeuralMassModel
fig = Figure();
ax1 = Axis(fig[1, 1], title="LPC branch emanating from, and connecting, the two GPD points")
lines!(ax1, extract_pars(pd_branchI.points), linewidth=2, color=:lightgray)
scatter!(ax1,
        extract_pars(lpc_branchII.points[get_nmfm_coefficients(lpc_branchII).<0.0]),      
        markersize=5)
scatter!(ax1, 
         extract_pars(lpc_branchII.points[get_nmfm_coefficients(lpc_branchII).>0.0]),
          markersize=5)
scatter!(ax1, extract_pars(gpd_points), color=:red, markersize=8)
ax1.xlabel = "c"
ax1.ylabel = "a"
# save("second_lpc_branch.png", fig) # hide
fig
##!CODEBOUNDARY

#
# **Conclusion:**
# This unified workflow demonstrates how to extract a periodic solution from a period-doubling branch, continue a new branch of solutions, detect special points, initiate an LPC branch, optionally filter solutions, compute normal form coefficients, and visualize the results. This process reveals intricate details about the bifurcation structure and stability characteristics of the analyzed dynamical system.
