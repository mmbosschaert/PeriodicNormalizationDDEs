# Neural Mass Model

This tutorial demonstrates the correctness of our derivations for critical periodic normal form coefficients associated with the fold, period-doubling, and Neimark-Sacker bifurcations. For the periodic fold bifurcation, we follow a branch of limit points of cycles, showing that the critical normal form coefficient's sign correctly vanishes when passing a cusp point. Similarly, for the period-doubling bifurcation, we continue a curve of period-doubling bifurcations and identify points where the critical normal form coefficient vanishes, indicating a generalized period-doubling point. Normal form analysis predicts a curve of limit points of cycles with double the period emanating from these points, and we will demonstrate this prediction. Lastly, for the Neimark-Sacker bifurcation, we confirm that the normal form coefficients of the Neimark-Sacker curves emanating from a double Hopf point align with the predictions derived from the normal form coefficients at the double Hopf point.

Using a neural mass model as a case study, validate the aforementioned predictions. The model is defined by a system of coupled delay differential equations with two time delays, representing the interaction of neural populations in a simplified yet biologically plausible framework:


```math
\dot{x}_1 = - x_1 - a \cdot g(b \cdot x_{1, \tau_1}) + c \cdot g(d \cdot x_{2, \tau_2})
\dot{x}_2 = - x_2 - a \cdot g(b \cdot x_{2, \tau_1}) + c \cdot g(d \cdot x_{1, \tau_2})
```


This tutorial guides readers through the process of verifying these derivations by integrating symbolic calculations with numerical continuation implemented in Julia. Along the way, we demonstrate how to set up the model, compute periodic solutions, and analyze bifurcation diagrams. By providing a detailed, step-by-step framework, this tutorial not only validates the correctness of our theoretical derivations but also equips researchers with the tools to investigate other models that incorporate delays.

## Load package

We start by loading the necessary packages.


```@example NeuralMassModel
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
```


## Define constants

Here, we define the constants used throughout the model. The `dims` variable represents the number of components. The delays, `τ₁` and `τ₂`, are provided in a vector `τs`, which includes a zero delay for convenience in calculations.


```@example NeuralMassModel
const gᵤ = 0.1
const gᵥ = 0.52
const β = 0.1
const ap = [1, 2]
const dims = 2
const τs = [0.0; 1.0]
```


## Model definition

The neural mass model is defined as a system of delay differential equations. The function `g(x)` acts as a nonlinear transfer function, using a hyperbolic tangent for sigmoid-like behavior. The main function, `neuralMassModel`, defines the system's dynamics in terms of delayed variables, parameters `a` and `c`, and the nonlinear interactions between two coupled units.


```@example NeuralMassModel
function acs(u, p)
  ζ, τ = p
  x, y = u[:, 1]
  xτ, yτ = u[:, 2]
  du = similar(u[:, 1])
  du[1] = τ * y
  du[2] = τ * (-x - gᵤ * xτ - 2 * ζ * y - gᵥ * yτ + β * xτ^3)
  du
end
```


## Model definition for simulation

For simulation purposes, the model is adapted to use a standard form compatible with solvers for delay differential equations. The function `neuralMassModel!` computes the derivatives in place, given the current state `u`, the delayed state `h`, parameters `p`, and the current time `t`.


```@example NeuralMassModel
const out = zeros(dims)
function acs!(du, u, h, p, t)
  ζ, τ = p
  x, y = u
  h(out, p, t - 1.0)
  du[1] = τ * y
  du[2] = τ * (-x - gᵤ * out[1] - 2.0 * ζ * y - gᵥ * out[2] + β * out[1]^3)
end
const out2 = zeros(3)
function acs_cross_section!(du, u, h, p, t)
  ζ, τ = p
  x, y, _ = u
  h(out2, p, t - 1.0)
  du[1] = τ * y
  du[2] = τ * (-x - gᵤ * out2[1] - 2.0 * ζ * y - gᵥ * out2[2] + β * out2[1]^3)
  du[3] = out2[2]
end
```


## Generate multi-linear forms

We use the `getJet` function to create a structure containing the multi-linear forms of the neural mass model. These forms are used in calculations of normal form coefficients, continuation of equilibria and periodic orbits, and their stability analysis. By calculating the multi-linear forms symbolically, we ensure precise and efficient computation for subsequent analysis steps.


```@example NeuralMassModel
# Calculate multi-linear forms
jet = getJet(acs, dims, τs)
```


## Equilibria continuation

Continuing the tutorial, we now perform a continuation of equilibria for the
neural mass model. This involves tracking steady states as a chosen parameter
(`c`) varies. First, we define a named tuple `par_indx` for easy access to
parameters and set `c` as the continuation parameter. Initial values for `a`
and `c` are stored in `params`, and parameter bounds are defined with
`parameterbounds` to constrain the values during continuation. The steady-state
branch is initialized using `SetupStstBranch`, which prepares the model,
parameter settings, and state variables for continuation. The process is then
executed with `continue!`, calculating equilibria over the specified range of
the continuation parameter. This step allows us to explore how the system's
steady states respond to changes in the parameter `c`.


```@example NeuralMassModel
# Create a named tuple for accessing parameters
par_indx = (ζ=1, τ=2)

# set parameters for double Hopf point derived analytically
ζ = -0.016225
τ = 5.89802
params = [ζ; τ]

# calculate stability of equilibrium
stst1 = stst(zeros(dims), params)
stst1 = stability(jet, stst1, τs)
sort(stst1.stability; by=x -> (floor(real(x); digits=6), -floor(imag(x); digits=6)), rev=true)

parameterbounds = (min=[-0.5; 5.0], max=[0.5; 15.0])

# calculate normal form coefficients
hoho1 = point_to_hoho(jet, stst1, τs)
hoho1 = PN.LocateDoubleHopf(jet, hoho1, τs; MaxIter=100, tol=1e-10)
hoho1 = normalform(jet, hoho1, τs);

# show normal form coefficients
hoho1.nmfm.δ
hoho1.nmfm.θ

PeriodicNormalizationDDEs.determine_unfolding_case([hoho1])
```


## Continuation of Hopf branches

We continue the Hopf branches crossing the double Hopf point in both directions. We set up the branches using the function `SetupHopfBranch` and continue them using the function `continue!`.

```@example NeuralMassModel
hopf1 = hopf_from_hoho(jet, hoho1, τs)
hopf_brI = SetupHopfBranch(jet, hopf1, τs; parameterbounds=parameterbounds, δ=0.001, δmax=0.01, MaxNumberofSteps=20000);
continue!(hopf_brI)
reverse_branch!(hopf_brI)
continue!(hopf_brI)
```


## Parameter Manipulation

In this section, we define a convenience function to extract parameters from branch points, enabling easy access and manipulation of the parameter values. To align with the continuation setup, the function reverses the order of the parameters, placing `c` as the first coordinate.


```@example NeuralMassModel
extract_pars(p) = reduce(hcat, [p.parameters[1:2] for p in p])
```


## Plot Hopf branches

Here we plot the continued Hopf branch. We see that it loops back two the double Hopf points. Additionally, we see that there is a second double Hopf point on the Hopf branch.

```@example NeuralMassModel
fig = Figure();
ax = Axis(fig[1, 1], title="Hopf branch with double Hopf points")
lines!(ax, extract_pars(hopf_brI.points), color=:royalblue)
ax.xlabel = "ζ"
ax.ylabel = "τ"
fig
```


## Calculate stability along the Hopf curves

Here we calculate the stability of the Hopf points on the Hopf branches.


```@example NeuralMassModel
## calculate stability along the Hopf curves
stability(jet, hopf_brI, τs)
```

## Detect special points

Using the spectral information derived in the previous step we use the function
`find_corrected_double_hopf_points` to locates double Hopf points present on the
Hopf branches. We set the maximum number of iterations to 100 and the tolerance to 1e-15.
After this we order the double hopf points such the order remains
the same for each execution of the code.
We see that the first double Hopf point is the one we found analytically.

```@example NeuralMassModel
double_hopf_points = PN.find_corrected_double_hopf_points(
  [hopf_brI], jet, τs; 
  max_iter=100, tol=1e-15)
double_hopf_points = sort(double_hopf_points, by = x -> x.parameters[2])
```


## Detect generalized Hopf points

Similar we use the function `find_genh_points` to locate generalized Hopf points present on the Hopf branches. We set the maximum number of iterations to 100 and the tolerance to 1e-10.

```@example NeuralMassModel
genh_points = PN.find_genh_points(
  [hopf_brI], jet, τs, par_indx.ζ;
  max_iter=100, tol=1e-10)
genh_points = sort(genh_points, by = x -> x.parameters[2])
```


## Visualizing Hopf branches with bifurcation points

We creat a figure with the continued Hopf branches together with the several double Hopf and generalized Hopf pionts

```@example NeuralMassModel
figHopf = Figure();
ax1 = Axis(figHopf[1, 1], title = "Bifurcation diagram")
lines!(ax1, extract_pars(hopf_brI.points), label="Hopf branch", color=:royalblue)
for (i, p) in enumerate(double_hopf_points)
  scatter!(ax1, [p.parameters[1]], [p.parameters[2]], color=:orangered)
  text!(ax1, p.parameters[1], p.parameters[2], text="HH$i")
end
for (i,p) in enumerate(genh_points)
  scatter!(ax1, [p.parameters[1]], [p.parameters[2]], color=:purple)
  text!(ax1, p.parameters[1], p.parameters[2], text="GH$i")
end
ax1.xlabel = "ζ"
ax1.ylabel = "τ"
axislegend(ax1, merge=true)
figHopf
```



## Compute normal form coefficients at the double Hopf point

Here we use the function `compute_nmfm_coefficients` to compute the normal form coefficients at the double Hopf points.

```@example NeuralMassModel
double_hopf_points = PN.compute_nmfm_coefficients(double_hopf_points, jet, τs)
```



## Continuation of Neimark-Sacker branches

Here we start with continuation of the Neimark-Sacker branches emanating from the double Hopf points. We set up the branches using the function `SetupNSBranch` and continue them using the function `continue!`. Additionally, we define a callback function to stop the continuation when the distance between the second double Hopf point and the solution is less than 1e-4. Indeed as we will see in the next figure the two double Hopf points are connected via a Neimark-Sacker curve. We finish this code block with continuation of the second Neimark-Sacker curve emanating from the second double Hopf point.


```@example NeuralMassModel
ϵ₁ = 0.08
ϵ₂ = 0.1
ntst = 20
ncol = 3
nsbranchI = SetupNSBranch(
  jet, hoho1, ϵ₁, ϵ₂, ntst, ncol, τs;
  parameterbounds=parameterbounds, δ=0.01, δmin=1e-06, δmax=[0.1; 0.1],
  MaxNumberofSteps=[300; 600], con_par = [par_indx.ζ; par_indx.τ]
);
continue!(nsbranchI[1])

# create callback function to stop continuation when distance between the second double Hopf point and solution is less than 1e-4
function callback_ns_branchI(x, double_hopf_points)
    distance = norm(double_hopf_points[2].parameters - x.parameters)
    if distance < 10^-4
        println("Callback: distance between double Hopf piont and solution is less than 1e-4")
        return true
    end
    return false
end

# start continuation of the second Neimark-Sacker branch emanating from the first double Hopf point
continue!(nsbranchI[2], callback = x -> callback_ns_branchI(x, double_hopf_points))
# start continuation of the second Neimark-Sacker branch emanating from the second double Hopf point
psol_guessII = doubleHopfToPsol(jet, double_hopf_points[2], ϵ₁, ϵ₂, ntst, ncol, τs);
nsbranchII = SetupNSBranch(jet, psol_guessII[2], τs; parameterbounds=parameterbounds, δ=0.02, δmax=0.1, MaxNumberofSteps=300,
                        con_par = [par_indx.ζ; par_indx.τ]);
continue!(nsbranchII)
```




## Next we compute the

compute normalform coefficients and multipliers


```@example NeuralMassModel
for nsbranch in [nsbranchI; nsbranchII]
  normal_form_coefficients!(jet, nsbranch, τs)
  multipliers!(jet, nsbranch, τs)
end
```



Plot Neimark-Sacker curvers along with normal form coefficients


```@example NeuralMassModel
fig = Figure();
ax1 = Axis(fig[1, 1], title="Neimark-Sacker curves")
lines!(ax1, extract_pars(hopf_brI.points), label="Hopf branch", color=:lightgray)
for branch in [nsbranchI; nsbranchII]
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
  scatter!(ax1, [p.parameters[1]], [p.parameters[2]], 
           color=:purple,
           markersize=5,
           label="double Hopf point",)
end
ax1.xlabel = "ζ"
ax1.ylabel = "τ"
axislegend(ax1, merge=:true, position=:lt)
fig
```


## Continue limit point of cycles emanating from the generalized Hopf points

We continue the limit point of cycles emanating from the generalized Hopf points.


```@example NeuralMassModel
ϵ = 0.01
ntst = 20
ncol = 3
MaxNumberofSteps = 10^6
δmax = 0.1

# Create the initial guess from the 1-th generalized Hopf point
lpc_guess = generalizedHopfToPsol(jet, genh_points[1], ϵ, ntst, ncol, τs)

# Set up the LPC branch
lpcbranchI = SetupLPCBranch(
    jet, lpc_guess, τs;
    parameterbounds = parameterbounds,
    δ = 0.01,
    δmax = δmax,
    MaxNumberofSteps = MaxNumberofSteps,
    NumberOfFails = 2,
    con_par = [par_indx.ζ; par_indx.τ]
)
    
function callback_LPC_branchI(x, genh_points)
    distance = norm(genh_points[2].parameters - x.parameters)
    if distance < 10^-4
        println("Callback: distance between the generalized Hopf piont and solution is less than 1e-4")
        return true
    end
    return false
end
    
continue!(lpcbranchI, callback = x -> callback_LPC_branchI(x, genh_points))

# Create the initial guess from the 3-th generalized Hopf point
lpc_guess = generalizedHopfToPsol(jet, genh_points[3], ϵ, ntst, ncol, τs)

# Set up the LPC branch
lpcbranchII = SetupLPCBranch(
    jet, lpc_guess, τs;
    parameterbounds = parameterbounds,
    δ = 0.01,
    δmax = δmax,
    MaxNumberofSteps = MaxNumberofSteps,
    NumberOfFails = 2,
    con_par = [par_indx.ζ; par_indx.τ]
)

continue!(lpcbranchII)
lpc_branches = [lpcbranchI, lpcbranchII]
```


## Plot limit point of cycles branches

Here we plot the three continued limit point of cycles branches. We see that there is a cusp point on the first lpc branch where we should expect the normal form coefficients to change sign.


```@example NeuralMassModel
fig = Figure();
set_theme!(theme_light())
ax1 = Axis(fig[1, 1])
lines!(ax1, extract_pars(hopf_brI.points), label="Hopf branch", color=:lightgray)
for lpc_branch in lpc_branches
  lines!(ax1, extract_pars(lpc_branch.points), color=:green, label="LPC branch")
end
for (i,p) in enumerate(genh_points)
  scatter!(ax1, [p.parameters[1]], [p.parameters[2]], color=:purple, label="GH")
end
axislegend(ax1, merge=:true)
fig
```


## Calculate normal form coefficients

We calculate the normal form coefficients for each limit point of cycles branch.


calculate normal form coefficients
```
for branch in lpc_branches
  normal_form_coefficients!(jet, branch, τs)
end
```


## Plot sign of normal form coefficients along the branches

In this visualization, we compare the signs of the normal form coefficients along each limit-point cycle (LPC) branch, coloring positive coefficients in blue and negative coefficients in orange-red to draw attention to sign changes. We see that the there is a indeed a cusp point on the first LPC branch where the normal form coefficients change sign.

```@example NeuralMassModel
fig = Figure();
set_theme!(theme_light())
ax1 = Axis(fig[1, 1], title="Normal form coefficients along the LPC branches")
lines!(ax1, extract_pars(hopf_brI.points), label="Hopf branch", color=:lightgray)
for branch in lpc_branches
  if ~isempty(branch.points[get_nmfm_coefficients(branch).<0.0])
    scatter!(
      ax1,
      extract_pars(branch.points[get_nmfm_coefficients(branch).<0.0]),
      markersize=5,
      color=:orangered,
      label="LPC branch negative sign"
    )
  end
  if ~isempty(branch.points[get_nmfm_coefficients(branch).>0.0])
    scatter!(
      ax1,
      extract_pars(branch.points[get_nmfm_coefficients(branch).>0.0]),
      markersize=5,
      color=:royalblue,
      label="LPC branch positive sign"
    )
  end
end
for (i,p) in enumerate(genh_points)
  scatter!(ax1, [p.parameters[1]], [p.parameters[2]], 
           color=:purple,
           markersize=5,
           label="Generalized Hopf point",)
end
axislegend(ax1, merge=:true, position=:rb)
fig
```
