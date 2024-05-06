struct stst
  coords::Vector{Float64}
  parameters::Vector{Float64}
  stability::Union{Vector{ComplexF64},Nothing}
end
stst(coords, parameters) = stst(coords, parameters, nothing)

# Define custom show function for DoubleHopf
Base.show(io::IO, stst1::stst) = begin
  println(io, "Coordinates: $(stst1.coords)")
  println(io, "Parameters: $(stst1.parameters)")
  if stst1.stability === nothing
    println(io, "Stability: not calculated")
  else
    eigenvalues = sort(stst1.stability; by=x -> (floor(real(x); digits=6), -floor(imag(x); digits=6)), rev=true)
    println(io, "Stability:")
    display(eigenvalues)
  end
end

function vec(stst1::stst, con_par)
  vcat([stst1.coords, stst1.parameters[con_par]]...)
end

function vec_to_point(vec, prev_stst::stst, con_par)
  dims = length(vec) - 1
  params = prev_stst.parameters
  stst(vec[1:dims], [params[1:con_par-1]; vec[dims+1]; params[con_par+1:end]])
end

function SetupStstBranch(model, con_par, eq, params, τs, dims; parameterbounds=nothing, δ=0.001, δmin=1e-06, δmax=0.1, MaxNumberofSteps=250, NumberOfFails=4)
  # generate code for continuation of steady points
  f = (x, _) -> model(repeat(x[1:dims], 1, length(τs)), [params[1:con_par-1]; x[dims+1]; params[con_par+1:end]])
  @variables xx[1:dims+1] xx_prev[1:dims+1]
  xx = collect(xx)
  xx_prev = collect(xx_prev)
  f = build_function(f(xx, xx_prev), xx, xx_prev, expression=Val{false})[1]
  df = build_function(Symbolics.jacobian(f(xx, xx_prev), xx), xx, xx_prev, expression=Val{false})[1]

  stst1 = stst(eq, params)
  point_to_vec(p) = [p.coords; p.parameters[con_par]]
  vec_to_point(v, _) = stst(v[1:dims], [params[1:con_par-1]; v[dims+1]; params[con_par+1:end]])
  x₀ = point_to_vec(stst1)

  v₀ = [df(x₀, x₀); rand(dims + 1)'] \ [zeros(dims); 1.0]
  v₀ /= norm(v₀)

  stst_branch = (points=stst[],
    tangents=[],
    stepsizes=[],
    f=f,
    df=df,
    parameterbounds=parameterbounds,
    δ=δ,
    δmin=δmin,
    δmax=δmax,
    MaxNumberofSteps=MaxNumberofSteps,
    con_par=con_par,
    NumberOfFails=NumberOfFails)
  push!(stst_branch.points, stst1)
  push!(stst_branch.tangents, v₀)
  push!(stst_branch.stepsizes, 0.0)
  stst_branch
end

function LocateHopfPoints(jet, branch, τs)
  # search for Hopf bifucation points along the steady-state branch stst
  num_unstable = [sum(real(p.stability) .> 0) for p in branch.points]

  # search for changes in number of unstable eigenvalues
  ind_hopf = findall(p -> abs(p) == 2, num_unstable[2:end] - num_unstable[1:end-1])

  # correct potential Hopf points
  n = length(branch.points[1].coords)
  @variables xx[1:3n+2] xx_prev[1:3n+2]
  xx = vcat(xx...)
  xx_prev = vcat(xx_prev...)
  res = similar(xx)
  Δre, Δim = characteristic_matrices_unevaluated_re_im(jet.system, 2, τs)

  params = branch.points[1].parameters
  con_par = branch.con_par
  params[1:con_par-1]
  xx[n+1]
  params[con_par+1:end]
  Hopf_res!(res, jet.system, τs, Δre, Δim, [xx[1:n]; params[1:con_par-1]; xx[n+1]; params[con_par+1:end]; xx[n+2:end]], [xx_prev[1:n]; params[1:con_par-1]; xx_prev[n+1]; params[con_par+1:end]; xx_prev[n+2:end]], n)
  hopfJac = Symbolics.jacobian(res, xx)

  f = build_function(res, xx, xx_prev, expression=Val{false})[1]
  Df = build_function(hopfJac, xx, xx_prev, expression=Val{false})[1]

  hopf_points = []
  for point in branch.points[ind_hopf]
    ind_omega = sortperm(abs.(real(point.stability)))
    ω₀ = abs(imag(point.stability[ind_omega[1]]))
    F = eigen(jet.Δ(repeat(point.coords, 1, length(τs)), [params[1:con_par-1]; point.parameters[con_par]; params[con_par+1:end]], ω₀ * im))
    hopf = Hopf(point.coords, point.parameters, F.vectors[:, 1], ω₀, nothing, nothing)

    x₀ = vcat([hopf.coords, hopf.parameters[con_par], real(hopf.v), imag(hopf.v), hopf.ω]...)

    x₁, _ = DDEBifTool.newton(f, Df, x₀, x₀, tol=1e-14)
    hopf = Hopf(x₁[1:n], [params[1:con_par-1]; x₁[n+1]; params[con_par+1:end]], x₁[n+2:2n+1] + x₁[2n+2:3n+1] * im, x₁[end], nothing, nothing)
    push!(hopf_points, hopf)
  end

  hopf_points
end
