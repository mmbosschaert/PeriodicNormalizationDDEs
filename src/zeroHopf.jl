function locate_zero_hopf(branch)
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
  ind_zero_hopf = findall(p -> abs(p) == 1, num_unstable[2:end] - num_unstable[1:end-1])
  ind_zero_hopf
end

mutable struct FoldHopfNormalform
  g200
  g110
  g011
  g300
  g111
  g210
  g021
  b
  c
  d
  e
  s
  θ
  q0
  q1
  h200
  h011
  h020
  h110
  h000μ
  K
  ω₀
  ω₁
  ω₂
end

mutable struct FoldHopf
  coords::Vector{Float64}
  parameters::Vector{Float64}
  v₁::Union{Vector{Float64},Vector{Num}}
  v₂::Union{Vector{ComplexF64},Vector{Complex{Num}}}
  ω₀::Float64
  stability::Union{Vector{ComplexF64},Nothing}
  nmfm::Union{FoldHopfNormalform,Nothing}
end
FoldHopf(coords, parameters, v₁, v₂, ω₀) = FoldHopf(coords, parameters, v₁, v₂, ω₀, nothing, nothing)

# Define custom show function for FoldHopfNormalform
Base.show(io::IO, nf::FoldHopfNormalform) = begin
  println(io, "g200: $(nf.g200)")
  println(io, "g110: $(nf.g110)")
  println(io, "g011: $(nf.g011)")
  println(io, "g300: $(nf.g300)")
  println(io, "g111: $(nf.g111)")
  println(io, "g210: $(nf.g210)")
  println(io, "g021: $(nf.g021)")
  println(io, "b: $(nf.b)")
  println(io, "c: $(nf.c)")
  println(io, "d: $(nf.d)")
  println(io, "e: $(nf.e)")
  println(io, "s: $(nf.s)")
  println(io, "θ: $(nf.θ)")
  println(io, "q0: $(nf.q0)")
  println(io, "q1: $(nf.q1)")
  println(io, "h200: $(nf.h200)")
  println(io, "h011: $(nf.h011)")
  println(io, "h020: $(nf.h020)")
  println(io, "h110: $(nf.h110)")
  #println(io, "h000μ: $(nf.h000μ)")
  println(io, "K: $(nf.K)")
  println(io, "ω₀: $(nf.ω₀)")
  println(io, "ω₁: $(nf.ω₁)")
  println(io, "ω₂: $(nf.ω₂)")
end

# Define custom show function for DoubleHopf
Base.show(io::IO, hopf::FoldHopf) = begin
  println(io, "Coordinates: $(hopf.coords)")
  println(io, "Parameters: $(hopf.parameters)")
  println(io, "V₁: $(hopf.v₁)")
  println(io, "V₂: $(hopf.v₂)")
  println(io, "ω₀: $(hopf.ω₀)")
  println(io, "Stability: $(hopf.stability)")
  println(io, "Normal form: $(hopf.nmfm)")
end

function point_to_zeho(jet, p, τs)
  if p.stability === nothing
    println("Need to calculate stability")
    return nothing
  else
    # Find the eigenvalue closest to zero in the complex plane
    ~, closest_to_zero_index = findmin(abs, p.stability)

    # Remove the eigenvalue closest to zero
    filtered_eigenvalues = deleteat!(copy(p.stability), closest_to_zero_index)

    freqs = closest_eigenvalues_to_imaginary_axis(filtered_eigenvalues; num_closest=2)
    ω₀ = abs.(imag(freqs[1]))

    m = length(τs)
    φ = repeat(p.coords, 1, m)
    Δ(λ) = jet.Δ(φ, p.parameters, λ)

    λ₁ = 0.0
    λ₂ = ω₀ * im

    _, s, V = svd(Δ(λ₁))
    indxmin = last(findmin(s))
    q1 = V[:, indxmin]

    _, s, V = svd(Δ(λ₂))
    indxmin = last(findmin(s))
    q2 = V[:, indxmin]

    FoldHopf(p.coords, p.parameters, q1, q2, ω₀)
  end
end

# define border inverse of characteristic matrix
function Δᴵᴺⱽ(λ, q, p, y, Δ) 
  ([Δ(λ) q; [transpose(p) 0]]\[y; 0])[1:end-1]
end
function Aᴵᴺⱽ(λ, q, p, η, κ, Δ, Δ′, Δ′′)
  ξ = Δᴵᴺⱽ(λ, q, p, η + κ * Δ′(λ) * q, Δ)
  γ = first(transpose(p) * (-Δ′(λ) * ξ + 0.5 * κ * Δ′′(λ) * q))
  θ -> exp(λ * θ) * (ξ + γ * q - κ * θ * q)
end

function normalform(jet, zeho::FoldHopf, τs)

  m = length(τs)
  φ = repeat(zeho.coords, 1, m)
  α = zeho.parameters

  λ₀ = 0.0
  λ₁ = zeho.ω₀ * im

  Δ(λ) = jet.Δ(φ, α, λ)
  Δ′(λ) = jet.Δ′(φ, α, λ)
  Δ′′(λ) = jet.Δ′′(φ, α, λ)
  # q = nullspace(Δ(λ); atol=1e-14)
  # p = nullspace(Δ(λ)';atol=1e-14)

  _, s, V = svd(Δ(λ₀))
  indxmin = last(findmin(s))
  q0 = V[:, indxmin]

  _, s, V = svd(transpose(Δ(λ₀)))
  indxmin = last(findmin(s))
  p0 = V[:, indxmin]

  _, s, V = svd(Δ(λ₁))
  indxmin = last(findmin(s))
  q1 = V[:, indxmin]

  _, s, V = svd(transpose(Δ(λ₁)))
  indxmin = last(findmin(s))
  p1 = V[:, indxmin]

  # normalize
  p0 /= transpose(p0) * (Δ′(λ₀) * q0)
  p1 /= transpose(p1) * (Δ′(λ₁) * q1)

  @assert first(transpose(p0) * (Δ′(λ₀) * q0)) ≈ 1.0
  @assert first(transpose(p1) * (Δ′(λ₁) * q1)) ≈ 1.0

  # multi-linear forms at bt point
  Ξ(h) = vcat([h(-τ) for τ ∈ τs]...)
  B(v₁, v₂) = jet.D2(φ, α, Ξ(v₁), Ξ(v₂))
  C(v₁, v₂, v₃) = jet.D3(φ, α, Ξ(v₁), Ξ(v₂), Ξ(v₃))
  A1(v₁, p₁) = jet.D11(φ, α, Ξ(v₁), p₁)
  J1 = jet.D01(φ, α)

  ϕ0(_) = q0
  ϕ1(θ) = exp(λ₁ * θ) * q1

  # Normal form coefficients
  # Quadratic
  g200 = 1 / 2 * transpose(p0) * B(ϕ0, ϕ0)
  g110 = transpose(p1) * B(ϕ0, ϕ1)
  g011 = transpose(p0) * B(ϕ1, conj ∘ ϕ1)

  # Quadratic center manifold
  h200 = Aᴵᴺⱽ(0.0, q0, p0, B(ϕ0, ϕ0), -2g200, Δ, Δ′, Δ′′)
  h020(θ) = exp(2 * λ₁ * θ) * (Δ(2 * λ₁) \ B(ϕ1, ϕ1))
  h110 = Aᴵᴺⱽ(λ₁, q1, p1, B(ϕ0, ϕ1), -g110, Δ, Δ′, Δ′′)
  h011 = Aᴵᴺⱽ(0.0, q0, p0, B(ϕ1, conj ∘ ϕ1), -g011, Δ, Δ′, Δ′′)

  # Cubic normal form coefficients
  g300 = (1 / 6) * transpose(p0) * (3 * B(ϕ0, h200) + C(ϕ0, ϕ0, ϕ0))
  g111 = transpose(p0) * (B(ϕ0, h011) + B(conj ∘ ϕ1, h110) + B(ϕ1, conj ∘ h110) + C(ϕ0, ϕ1, conj ∘ ϕ1))
  g210 = (1 / 2) * transpose(p1) * (B(ϕ1, h200) + 2 * B(ϕ0, h110) + C(ϕ0, ϕ0, ϕ1))
  g021 = (1 / 2) * transpose(p1) * (B(conj ∘ ϕ1, h020) + 2 * B(ϕ1, h011) + C(ϕ1, ϕ1, conj ∘ ϕ1))

  # Gavrilov normalform coefficients
  b = g200
  c = g011
  d = g110 - λ₁ * g300 / g200
  e = real(g210 + g110 * (real(g021) / g011 - 3 * g300 / (2 * g200) + g111 / (2 * g011)) - g021 * g200 / g011)

  s = b * c
  θ = real(g110) / g200
  println("s: $s, θ: $θ, e: $e")

  # parameter-dependent normal form coefficients
  γ = (transpose(p0) * J1)'
  s1 = γ / dot(γ, γ)
  s2 = [γ[2]; -γ[1]]

  r1(_) = Δᴵᴺⱽ(0.0, q0, p0, J1 * s1, Δ)
  r2(_) = Δᴵᴺⱽ(0.0, q0, p0, J1 * s2, Δ)
  r3(θ) = Δᴵᴺⱽ(0.0, q0, p0, (Δ′(0.0) * q0), Δ) - θ * q0

  LL = [dot(p0, A1(ϕ0, s2)+B(ϕ0, r2)) transpose(p0)*B(ϕ0, ϕ0)
        transpose(p1)*(A1(ϕ1, s2)+B(ϕ1, r2)) transpose(p1)*B(ϕ1, ϕ0)]

  RR = -[transpose(p0) * (A1(ϕ0, s1) + B(ϕ0, r1) - B(ϕ0, r3))
    transpose(p1) * (A1(ϕ1, s1) + B(ϕ1, r1) - B(ϕ1, r3))]

  δ = real(LL) \ [real(RR) [0; 1]]

  K10 = s1 + δ[1] * s2
  K01 = δ[3] * s2
  K = [K10 K01]

  h00010(θ) = r1(θ) + δ[1] * r2(θ) + δ[2] * q0 - r3(θ)
  h00001(θ) = δ[3] * r2(θ) + δ[4] * q0

  h000μ = [h00010, h00001]

  ω₁ = imag(transpose(p1) * B(ϕ1, h00010) + transpose(p1) * A1(ϕ1, K10))
  ω₂ = imag(transpose(p1) * B(ϕ1, h00001) + transpose(p1) * A1(ϕ1, K01))

  nmfm = FoldHopfNormalform(
    g200,
    g110,
    g011,
    g300,
    g111,
    g210,
    g021,
    b,
    c,
    d,
    e,
    s,
    θ,
    q0,
    q1,
    h200,
    h011,
    h020,
    h110,
    h000μ,
    K,
    zeho.ω₀,
    ω₁,
    ω₂
  )

  zeho = @set zeho.nmfm = nmfm
end
