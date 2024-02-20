mutable struct DoubleHopfNormalform
  q
  g2100
  g1011
  g1110
  g0021
  θ
  δ
  b
  K
  h0000
  h0011
  h0020
  h2000
  h1100
end

mutable struct DoubleHopf
  coords::Vector{Float64}
  parameters::Vector{Float64}
  v₁::Union{Vector{ComplexF64},Vector{Complex{Num}}}
  v₂::Union{Vector{ComplexF64},Vector{Complex{Num}}}
  ω₁::Float64
  ω₂::Float64
  stability::Union{Vector{ComplexF64},Nothing}
  nmfm::Union{DoubleHopfNormalform,Nothing}
end
DoubleHopf(coords, parameters, v₁, v₂, ω₁, ω₂) = DoubleHopf(coords, parameters, v₁, v₂, ω₁, ω₂, nothing, nothing)

function point_to_hoho(jet, p, τs)
  if p.stability === nothing
    println("Need to calculate stability")
    return nothing
  else
    ω_indx = findall(x -> abs(real(x)) < 1e-02, p.stability)
    if length(ω_indx) !== 4
      println("Unable to extract frequencies")
    end
    ω₁, ω₂ = sort(abs.(imag(p.stability[ω_indx])))[[1, 3]]

    m = length(τs)
    φ = repeat(p.coords, 1, m)
    Δ(λ) = jet.Δ(φ, p.parameters, λ)

    λ₁ = ω₁ * im
    λ₂ = ω₂ * im

    _, s, V = svd(Δ(λ₁))
    indxmin = last(findmin(s))
    q1 = V[:, indxmin]

    _, s, V = svd(Δ(λ₂))
    indxmin = last(findmin(s))
    q2 = V[:, indxmin]

    DoubleHopf(p.coords, p.parameters, q1, q2, ω₁, ω₂)
  end
end

function locate_double_hopf(branch)
  num_unstable = zeros(Int, length(branch))
  for (i, p) in enumerate(branch)
    p.stability
    # remove ω from stability field
    ind1 = findfirst(λ -> imag(λ) ≈ p.ω, p.stability)
    ind2 = findfirst(λ -> imag(λ) ≈ -p.ω, p.stability)
    eigenvalues = p.stability[1:end.!=ind1.&&1:end.!=ind2]
    # count unstable eigenvalues
    num_unstable[i] = sum(real(eigenvalues) .> 0.0)
  end
  ind_double_hopf = findall(p -> abs(p) == 2, num_unstable[2:end] - num_unstable[1:end-1])
  ind_double_hopf
end

function vec(hoho::DoubleHopf, _)::Vector{Float64}
  vcat(hoho.coords, hoho.parameters, real(hoho.v₁), imag(hoho.v₁), real(hoho.v₂), imag(hoho.v₂), hoho.ω₁, hoho.ω₂)
end

function vec_to_point(v::Vector{Float64}, ::DoubleHopf, _)
  dims = div(length(v) - 4, 5)
  DoubleHopf(v[1:dims], v[dims+1:dims+2], v[dims+3:2dims+2] + v[2dims+3:3dims+2] * im, v[3dims+3:4dims+2] + v[4dims+3:5dims+2] * im, v[5dims+3], v[5dims+4])
end

function doubleHopf_res!(res, model, τs, Δre, Δim, xx, hopf_prev::DoubleHopf, n)
  x, α, v1re, v1im, v2re, v2im, ω₁, ω₂ = xx[1:n], xx[n+1:n+2], xx[n+3:2n+2], xx[2n+3:3n+2], xx[3n+3:4n+2], xx[4n+3:5n+2], xx[5n+3], xx[5n+4]
  v1re_prev, v1im_prev = real(hopf_prev.v₁), imag(hopf_prev.v₁)
  v2re_prev, v2im_prev = real(hopf_prev.v₂), imag(hopf_prev.v₂)

  res[1:n] = model(repeat(x, 1, length(τs)), α)
  res[n+1:2n] = Δre(repeat(x, 1, length(τs)), α, ω₁) * v1re - Δim(repeat(x, 1, length(τs)), α, ω₁) * v1im
  res[2n+1:3n] = Δre(repeat(x, 1, length(τs)), α, ω₁) * v1im + Δim(repeat(x, 1, length(τs)), α, ω₁) * v1re
  res[3n+1:4n] = Δre(repeat(x, 1, length(τs)), α, ω₂) * v2re - Δim(repeat(x, 1, length(τs)), α, ω₂) * v2im
  res[4n+1:5n] = Δre(repeat(x, 1, length(τs)), α, ω₂) * v2im + Δim(repeat(x, 1, length(τs)), α, ω₂) * v2re
  res[5n+1] = dot(v1re, v1re) + dot(v1im, v1im) - 1
  res[5n+2] = dot(v1im, v1re_prev) + dot(v1re, v1im_prev)
  res[5n+3] = dot(v2re, v2re) + dot(v2im, v2im) - 1
  res[5n+4] = dot(v2im, v2re_prev) + dot(v2re, v2im_prev)
end

function doubleHopf_res(model, τs, Δre, Δim, xx, hopf_prev::DoubleHopf, n)
  res = zeros(5n + 4)
  doubleHopf_res!(res, model, τs, Δre, Δim, xx, hopf_prev::DoubleHopf, n)
  res
end

function defining_system_DoubleHopf(jet, τs, dims)
  # create f and jac for continuation of Hopf points
  @variables xx[1:5dims+4] v1_prev[1:2dims] v2_prev[1:2dims]
  xx = Symbolics.scalarize(xx)
  v1_prev = Symbolics.scalarize(v1_prev)
  v2_prev = Symbolics.scalarize(v2_prev)
  res = similar(xx)
  hoho_prev = DoubleHopf(zeros(dims), zeros(2), v1_prev[1:dims] + v1_prev[dims+1:2dims] * im, v2_prev[1:dims] + v2_prev[dims+1:2dims] * im, 0.0, 0.0)
  Δre, Δim = characteristic_matrices_unevaluated_re_im(jet.system, dims, τs)
  doubleHopf_res!(res, jet.system, τs, Δre, Δim, xx, hoho_prev, dims)

  # calculate jacobian symbolically
  jac = Symbolics.jacobian(res, xx)

  df = build_function(jac, xx, v1_prev, v2_prev, expression=Val{false})[1]

  f = (x, hoho) -> doubleHopf_res(jet.system, τs, Δre, Δim, x, hoho, dims)
  Df = (x, hoho) -> df(x, [real(hoho.v₁); imag(hoho.v₁)], [real(hoho.v₂); imag(hoho.v₂)])

  f, Df
end

function normalform(jet, hoho::DoubleHopf, τs)

  m = length(τs)
  φ = repeat(hoho.coords, 1, m)
  α = hoho.parameters

  λ₁ = hoho.ω₁ * im
  λ₂ = hoho.ω₂ * im

  Δ(λ) = jet.Δ(φ, α, λ)
  Δ′(λ) = jet.Δ′(φ, α, λ)
  # q = nullspace(Δ(λ); atol=1e-14)
  # p = nullspace(Δ(λ)';atol=1e-14)

  _, s, V = svd(Δ(λ₁))
  indxmin = last(findmin(s))
  q1 = V[:, indxmin]

  _, s, V = svd(transpose(Δ(λ₁)))
  indxmin = last(findmin(s))
  p1 = V[:, indxmin]

  _, s, V = svd(Δ(λ₂))
  indxmin = last(findmin(s))
  q2 = V[:, indxmin]

  _, s, V = svd(transpose(Δ(λ₂)))
  indxmin = last(findmin(s))
  p2 = V[:, indxmin]

  # normalize

  p1 /= transpose(p1) * (Δ′(λ₁) * q1)
  p2 /= transpose(p2) * (Δ′(λ₂) * q2)

  q = [q1, q2]
  p = [p1, p2]
  # first(transpose(p1)*(Δ′(λ₁)*q1)) ≈ 1.0
  # first(transpose(p2)*(Δ′(λ₂)*q2)) ≈ 1.0

  # multi-linear forms at bt point
  Ξ(h) = vcat([h(-τ) for τ ∈ τs]...)
  B(v₁, v₂) = jet.D2(φ, α, Ξ(v₁), Ξ(v₂))
  C(v₁, v₂, v₃) = jet.D3(φ, α, Ξ(v₁), Ξ(v₂), Ξ(v₃))
  A1(v₁, p₁) = jet.D11(φ, α, Ξ(v₁), p₁)
  J1 = jet.D01(φ, α)

  ϕ1(θ) = exp(λ₁ * θ) * q1
  ϕ2(θ) = exp(λ₂ * θ) * q2
  ϕs = [ϕ1, ϕ2]

  # Quadratic center manifold
  h1100(_) = Δ(0) \ B(ϕ1, conj ∘ ϕ1)
  h2000(θ) = exp(2 * λ₁ * θ) * (Δ(2 * λ₁) \ B(ϕ1, ϕ1))
  h1010(θ) = exp((λ₁ + λ₂) * θ) * (Δ(λ₁ + λ₂) \ B(ϕ1, ϕ2))
  h1001(θ) = exp((λ₁ - λ₂) * θ) * (Δ(λ₁ - λ₂) \ B(ϕ1, conj ∘ ϕ2))
  h0020(θ) = exp(2λ₂ * θ) * (Δ(2λ₂) \ B(ϕ2, ϕ2))
  h0011(_) = Δ(0) \ B(ϕ2, conj ∘ ϕ2)

  g2100 = (1 / 2) * transpose(p1) * (2 * B(h1100, ϕ1) + B(h2000, conj ∘ ϕ1) + C(ϕ1, ϕ1, conj ∘ ϕ1))
  g1011 = transpose(p1) * (B(h0011, ϕ1) + B(h1001, ϕ2) + B(h1010, conj ∘ ϕ2) + C(ϕ1, ϕ2, conj ∘ ϕ2))
  g1110 = transpose(p2) * (B(conj ∘ h1001, ϕ1) + B(h1010, conj ∘ ϕ1) + B(h1100, ϕ2) + C(ϕ1, conj ∘ ϕ1, ϕ2))
  g0021 = (1 / 2) * transpose(p2) * (2 * B(h0011, ϕ2) + B(h0020, conj ∘ ϕ2) + C(ϕ2, ϕ2, conj ∘ ϕ2))

  θ = real(g1011) / real(g0021)
  δ = real(g1110) / real(g2100)

  Id = diagm(ones(2))
  Γ = [transpose(p[i]) * (A1(ϕs[i], Id[:, j]) + B(ϕs[i], _ -> Δ(0) \ (J1 * Id[:, j]))) for i = 1:2, j = 1:2]
  K = inv(real(Γ))

  h000001 = Δ(0) \ (J1 * K[:, 1])
  h000010 = Δ(0) \ (J1 * K[:, 2])
  h0000 = [_ -> h000001, _ -> h000010]
  b = [imag(transpose(p[i]) * (A1(ϕs[i], K[:, j]) + B(ϕs[i], h0000[j]))) for i = 1:2, j = 1:2]

  nmfm = DoubleHopfNormalform(
    q,
    g2100,
    g1011,
    g1110,
    g0021,
    θ,
    δ,
    b,
    K,
    h0000,
    h0011,
    h0020,
    h2000,
    h1100
  )

  hoho = @set hoho.nmfm = nmfm
end


