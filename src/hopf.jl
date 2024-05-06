mutable struct Hopf
  coords::Vector{Float64}
  parameters::Vector{Float64}
  v::Union{Vector{ComplexF64},Vector{Complex{Num}}}
  ω::Float64
  stability::Union{Vector{ComplexF64},Nothing}
  nmfm::Union{ComplexF64,Nothing}
end
Hopf(coords, parameters, v, ω) = Hopf(coords, parameters, v, ω, nothing, nothing)

function vec(hopf::Hopf, _)::Vector{Float64}
  vcat(hopf.coords, hopf.parameters, real(hopf.v), imag(hopf.v), hopf.ω)
end

function vec_to_point(v::Vector{Float64}, ::Hopf, _)
  dims = div(length(v) - 1, 4)
  Hopf(v[1:dims], v[dims+1:2dims], v[2dims+1:3dims] + v[3dims+1:4dims] * im, v[4dims+1])
end

function Hopf_res!(res, model, τs, Δre, Δim, xx, xx_prev, n)
  x, α, vre, vim, ω = xx[1:n], xx[n+1:n+2], xx[n+3:2n+2], xx[2n+3:3n+2], xx[3n+3]
  vre_prev, vim_prev = xx_prev[n+3:2n+2], xx_prev[2n+3:3n+2]

  res[1:n] = model(repeat(x, 1, length(τs)), α)
  res[n+1:2n] = Δre(repeat(x, 1, length(τs)), α, ω) * vre - Δim(repeat(x, 1, length(τs)), α, ω) * vim
  res[2n+1:3n] = Δre(repeat(x, 1, length(τs)), α, ω) * vim + Δim(repeat(x, 1, length(τs)), α, ω) * vre
  res[3n+1] = dot(vre, vre_prev) - dot(vim, vim_prev) - 1
  res[3n+2] = dot(vim, vre_prev) + dot(vre, vim_prev)
end

function Hopf_res!(res, model, τs, Δre, Δim, xx, hopf_prev::Hopf, n)
  x, α, vre, vim, ω = xx[1:n], xx[n+1:n+2], xx[n+3:2n+2], xx[2n+3:3n+2], xx[3n+3]
  vre_prev, vim_prev = real(hopf_prev.v), imag(hopf_prev.v)

  res[1:n] = model(repeat(x, 1, length(τs)), α)
  res[n+1:2n] = Δre(repeat(x, 1, length(τs)), α, ω) * vre - Δim(repeat(x, 1, length(τs)), α, ω) * vim
  res[2n+1:3n] = Δre(repeat(x, 1, length(τs)), α, ω) * vim + Δim(repeat(x, 1, length(τs)), α, ω) * vre
  res[3n+1] = dot(vre, vre) + dot(vim, vim) - 1
  res[3n+2] = dot(vim, vre_prev) + dot(vre, vim_prev)
end

function Hopf_res(model, τs, Δre, Δim, xx, hopf_prev::Hopf, n)
  res = zeros(3n + 3)
  Hopf_res!(res, model, τs, Δre, Δim, xx, hopf_prev::Hopf, n)
  res
end

# function Hopf_res!(res, model, τs, Δ, xx, xx_prev, n)
# 
#     x, α, v, ω = xx[1:n], xx[n+1:n+2], xx[n+3:2n+2] + xx[2n+3:3n+2]*im, xx[3n+3]
#     v_prev = xx_prev[n+3:2n+2] + xx_prev[2n+3:3n+2]*im
# 
#     res[1:n] = model(repeat(x,1,length(τs)),α)
#     res[n+1:2n]  = real(Δ(repeat(x,1,length(τs)),α,ω*im)*v)
#     res[2n+1:3n] = imag(Δ(repeat(x,1,length(τs)),α,ω*im)*v)
#     res[3n+1] = real(dot(v,v_prev)) - 1
#     res[3n+2] = imag(dot(v,v_prev))
# end

function normalform(jet, hopf::Hopf, τs)

  m = length(τs)
  φ = repeat(hopf.coords, 1, m)
  α = hopf.parameters

  λ = hopf.ω * im
  Δ(λ) = jet.Δ(φ, α, λ)
  Δ′(λ) = jet.Δ′(φ, α, λ)

  q = qr(Δ(λ)').Q[:, end]
  p = qr(Δ(λ)).Q[:, end]

  _, s, V = svd(Δ(λ))
  indxmin = last(findmin(s))
  q = V[:, indxmin]

  _, s, V = svd(transpose(Δ(λ)))
  indxmin = last(findmin(s))
  p = V[:, indxmin]

  # normalize
  p /= transpose(p) * (Δ′(λ) * q)

  # multi-linear forms at bt point
  Ξ(h) = vcat([h(-τ) for τ ∈ τs]...)
  B(v₁, v₂) = jet.D2(φ, α, Ξ(v₁), Ξ(v₂))
  C(v₁, v₂, v₃) = jet.D3(φ, α, Ξ(v₁), Ξ(v₂), Ξ(v₃))

  ϕ(θ) = exp(λ * θ) * q
  h2000(θ) = exp(2 * λ * θ) * (Δ(2 * λ) \ B(ϕ, ϕ))
  h1100(_) = Δ(0) \ B(ϕ, conj ∘ ϕ)

  c₁ = first((1 / 2) * transpose(p) * (B(conj ∘ ϕ, h2000) + 2 * B(ϕ, h1100) + C(ϕ, ϕ, conj ∘ ϕ)))

  ℓ₁ = real(c₁) / hopf.ω
  ℓ₁
end

function detect_genh(jet, branch, τs)
  for hopf in branch
    hopf.nmfm = normalform(jet, hopf, τs)
  end
  nmfm = [hopf.nmfm for hopf in branch]
  findall(!iszero, sign.(nmfm)[2:end] - sign.(nmfm)[1:end-1])
end

# add two Hopf points
+(hopf1::Hopf, hopf2::Hopf) = Hopf(hopf1.coords + hopf2.coords, hopf1.parameters + hopf2.parameters, hopf1.v + hopf2.v, hopf1.ω + hopf2.ω)

# multiply Hopf point with scalar
*(a::Float64, hopf::Hopf) = Hopf(a * hopf.coords, a * hopf.parameters, a * hopf.v, a * hopf.ω)

# function to locate generalized Hopf points
# input are two Hopf points where the first Lyapunov coefficient changes sign
# output is the Hopf point in between with Lyaupnov coefficient close to zero
function locate_genh(jet, hopf_br, indx, par_dir, τs; MaxIter=100, tol=1e-10)

  hopf1 = hopf_br.points[indx]
  hopf2 = hopf_br.points[indx+1]

  parameter_index = length(hopf1.coords) + par_dir

  # create function and Jacobian for correction in one parameter
  fpar(p, hopf_ref) = hopf_br.f(vcat(p[1:parameter_index-1], hopf_ref.parameters[par_dir], p[parameter_index:end]), hopf_ref)[1:end-1]
  dfpar(p, hopf_ref) = hopf_br.df(vcat(p[1:parameter_index-1], hopf_ref.parameters[par_dir], p[parameter_index:end]), hopf1)[1:end-1, 1:end.!=parameter_index]

  for i = 1:MaxIter
    hopf_new = middle_point(fpar, dfpar, hopf1, hopf2, parameter_index, par_dir, jet, τs)
    if real(sign(hopf1.nmfm * hopf_new.nmfm)) < 0.0
      hopf2 = hopf_new
    else
      hopf1 = hopf_new
    end

    if abs(real(hopf_new.nmfm)) < tol
      break
    end
  end

  hopf_new = stability(jet, hopf_new, τs)
  genh_new = point_to_genhopf(hopf_new)
  genh_new = normalform(jet, genh_new, τs)

  # return generalized Hopf point
  genh_new
end

# This function calculates the "middle point" between two Hopf points
function middle_point(fpar, dfpar, hopf1::Hopf, hopf2::Hopf, parameter_index, par_dir, jet, τs)
  # construct Hopf point in between
  hopfmiddle = 0.5 * (hopf1 + hopf2)

  # convert to vector and remove fixed parameter
  pmiddle = vec(hopfmiddle, nothing)[1:end.!=parameter_index]

  # correct new point
  x_new, iterations, convergenced = newton(fpar, dfpar, pmiddle, hopfmiddle; tol=1e-10, maxIter=100)

  # insert fixed parameter
  insert!(x_new, parameter_index, hopfmiddle.parameters[par_dir])

  # convert to Hopf point
  hopf_new = DDEBifTool.vec_to_point(x_new, hopfmiddle, nothing)

  # calculate l1 norm of normal form coefficients
  hopf_new.nmfm = normalform(jet, hopf_new, τs)
  hopf_new
end

function LocateDoubleHopf(jet, hoho_guess, τs; MaxIter=100, tol=1e-10)
  dims = length(hoho_guess.coords)

  # convert point to DoubleHopf if not already
  if typeof(hoho_guess) !== DoubleHopf
    @show hoho_guess = point_to_hoho(jet, hoho_guess, τs)
  end

  f, df = DDEBifTool.defining_system_DoubleHopf(jet, τs, dims)
  x₀, _, converged = newton(f, df, vec(hoho_guess, nothing), hoho_guess; tol=1e-10, maxIter=100)
  if !converged
    println("Newton did not converge")
  end

  # convert corrected double Hopf point
  hoho = DDEBifTool.vec_to_point(x₀, hoho_guess, nothing)

  # add stability information
  hoho = stability(jet, hoho, τs)

  # return double Hopf point
  hoho
end

function defining_system_Hopf(jet, τs, dims)
  # create f and jac for continuation of Hopf points
  @variables xx[1:3dims+3] v_prev[1:2dims]
  xx = Symbolics.scalarize(xx)
  v_prev = Symbolics.scalarize(v_prev)
  res = similar(xx)
  hopf_prev = Hopf(zeros(dims), zeros(2), v_prev[1:dims] + v_prev[dims+1:2dims] * im, 0.0, nothing, nothing)
  Δre, Δim = characteristic_matrices_unevaluated_re_im(jet.system, dims, τs)
  Hopf_res!(res, jet.system, τs, Δre, Δim, xx, hopf_prev, dims)
  res[end] = 0.0

  hopfJac = Symbolics.jacobian(res, xx)

  f = build_function(res, xx, v_prev, expression=Val{false})[1]
  df = build_function(hopfJac, xx, v_prev, expression=Val{false})[1]

  f, df, Δre, Δim
end

function SetupHopfBranch(jet, hopf_point, τs; parameterbounds=nothing, δ=0.001, δmin=1e-06, δmax=0.01, MaxNumberofSteps=250, NumberOfFails=4)
  dims = length(hopf_point.coords)

  # create f, df, Δre, Δim for continuation of Hopf points
  f, df, Δre, Δim = defining_system_Hopf(jet, τs, dims)

  # hopf branch I
  x₀ = vec(hopf_point, nothing)
  x₀prev = x₀
  jac = df(x₀, x₀prev)
  # tangent vector
  v₀ = qr(jac').Q[:, end]

  hopf_branch = (points=Hopf[], tangents=[], stepsizes=[],
    f=(x, hopf) -> Hopf_res(jet.system, τs, Δre, Δim, x, hopf, dims),
    df=(x, hopf) -> df(x, [real(hopf.v); imag(hopf.v)]),
    parameterbounds=parameterbounds,
    δ=δ,
    δmin=δmin,
    δmax=δmax,
    MaxNumberofSteps=MaxNumberofSteps,
    con_par=nothing,
    NumberOfFails=NumberOfFails
  )

  push!(hopf_branch.points, hopf_point)
  push!(hopf_branch.tangents, v₀)
  push!(hopf_branch.stepsizes, 0.0)
  hopf_branch
end

function hopf_from_hoho(jet, hoho1, τs; ϵ=1e-03)
  # TODO: add second curve
  β₁ = -ϵ
  # coords = hoho1.coords + [sqrt(-β₁/real(hoho1.nmfm.g2100)); 0.0]
  # β₂ = 0.00001
  # coords = [0.0; sqrt(-β₂/real(hoho1.nmfm.g0021))]
  coords = [0.0; 0.0]

  # hopf1 = Hopf( coords, hoho1.parameters + hoho1.nmfm.K*[0.0; β₂], hoho1.nmfm.q[2], hoho1.ω₂)
  hopf1 = Hopf(coords, hoho1.parameters + hoho1.nmfm.K * [β₁; 0.0], hoho1.nmfm.q[1], hoho1.ω₁)

  m = length(τs)
  φ = repeat(hopf1.coords, 1, m)
  α = hopf1.parameters
  λ₁ = hoho1.ω₁ * im
  Δ(λ) = jet.Δ(φ, α, λ)
  _, s, V = svd(Δ(λ₁))
  indxmin = last(findmin(s))
  q1 = V[:, indxmin]

  hopf_prev = hopf1
  c = sqrt(1 / (dot(real(hopf1.v), real(hopf_prev.v)) - dot(imag(hopf1.v), imag(hopf_prev.v))))
  c = sqrt(1 / (dot(real(hopf1.v), real(hopf1.v)) - dot(imag(hopf1.v), imag(hopf1.v))))
  hopf1 = @set hopf1.v = c * q1
end
