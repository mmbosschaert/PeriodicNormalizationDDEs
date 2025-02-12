mutable struct GenHopfNormalform
  γ110
  γ101
  γ210
  γ201
  c₁
  c₂
  ℓ₁
  ℓ₂
  q
  h2000
  h1100
  h0010
  h0001
  ω₂
end

# TODO: implement this
mutable struct GenHopfNormalformHigherOrder
  K10
  K01
  K02
  K11
  K03
  c₁
  c₂
  c₃
  ℓ₁
  ℓ₂
  ℓ₃
  a3201
  q
  h2000
  h1100
  h3000
  h2100
  h4000
  h3100
  h2200
  h5000
  h4100
  h3200
  h6000
  h5100
  h4200
  h3300
  h4300
  h7000
  h6100
  h5200
  h0010
  h0001
  h0002
  h0011
  h0003
  h1002
  h1011
  h1003
  h1001
  h1010
  h2001
  h2010
  h2002
  h1101
  h1110
  h1102
  h3001
  h3010
  h2101
  h2110
  h4001
  h2201
  h3101
  h4101
  h5001
  h3201
  h2102
  h3002
  b101
  b110
  b201
  b102
end

mutable struct GenHopf
  coords::Vector{Float64}
  parameters::Vector{Float64}
  v::Union{Vector{ComplexF64},Vector{Complex{Num}}}
  ω::Float64
  stability::Union{Vector{ComplexF64},Nothing}
  nmfm::Union{GenHopfNormalform,GenHopfNormalformHigherOrder,Nothing}
end
GenHopf(coords, parameters, v, ω) = GenHopf(coords, parameters, v, ω, nothing, nothing)

# Define custom show function for GenHopfNormalform
Base.show(io::IO, nf::GenHopfNormalform) = begin
  println(io, "γ110: $(nf.γ110)")
  println(io, "γ101: $(nf.γ101)")
  println(io, "γ210: $(nf.γ210)")
  println(io, "γ201: $(nf.γ201)")
  println(io, "c₁: $(nf.c₁)")
  println(io, "c₂: $(nf.c₂)")
  println(io, "ℓ₁: $(nf.ℓ₁)")
  println(io, "ℓ₂: $(nf.ℓ₂)")
  println(io, "q: $(nf.q)")
  println(io, "h2000: $(nf.h2000)")
  println(io, "h1100: $(nf.h1100)")
  println(io, "h0010: $(nf.h0010)")
  println(io, "h0001: $(nf.h0001)")
  println(io, "ω₂: $(nf.ω₂)")
end

# TODO: implement this
# Define custom show function for GenHopfNormalformHigherOrder
Base.show(io::IO, nf::GenHopfNormalformHigherOrder) = begin
  println(io, "K10: $(nf.K10)")
  println(io, "K01: $(nf.K01)")
  println(io, "K02: $(nf.K02)")
  println(io, "K11: $(nf.K11)")
  println(io, "c₁: $(nf.c₁)")
  println(io, "c₂: $(nf.c₂)")
  println(io, "ℓ₁: $(nf.ℓ₁)")
  println(io, "ℓ₂: $(nf.ℓ₂)")
  println(io, "a3201: $(nf.a3201)")
  println(io, "q: $(nf.q)")
  println(io, "h2000: $(nf.h2000)")
  println(io, "h1100: $(nf.h1100)")
  println(io, "h3000: $(nf.h3000)")
  println(io, "h2100: $(nf.h2100)")
  println(io, "h4000: $(nf.h4000)")
  println(io, "h3100: $(nf.h3100)")
  println(io, "h2200: $(nf.h2200)")
  println(io, "h5000: $(nf.h5000)")
  println(io, "h4100: $(nf.h4100)")
  println(io, "h3200: $(nf.h3200)")
  println(io, "h0010: $(nf.h0010)")
  println(io, "h0001: $(nf.h0001)")
  println(io, "h0002: $(nf.h0002)")
  println(io, "h0011: $(nf.h0011)")
  println(io, "h1002: $(nf.h1002)")
  println(io, "h1011: $(nf.h1011)")
  println(io, "h1001: $(nf.h1001)")
  println(io, "h1010: $(nf.h1010)")
  println(io, "h2001: $(nf.h2001)")
  println(io, "h1101: $(nf.h1101)")
  println(io, "h3001: $(nf.h3001)")
  println(io, "h3010: $(nf.h3010)")
  println(io, "h2101: $(nf.h2101)")
  println(io, "h2110: $(nf.h2110)")
  println(io, "h4001: $(nf.h4001)")
  println(io, "h2201: $(nf.h2201)")
  println(io, "h3101: $(nf.h3101)")
  println(io, "h4101: $(nf.h4101)")
  println(io, "h5001: $(nf.h5001)")
  println(io, "h3201: $(nf.h3201)")
  println(io, "h2102: $(nf.h2102)")
  println(io, "h3002: $(nf.h3002)")
  println(io, "b101: $(nf.b101)")
  println(io, "b201: $(nf.b201)")
  println(io, "b102: $(nf.b102)")
end

# Define custom show function for GenHopf
Base.show(io::IO, hopf::GenHopf) = begin
  println(io, "Coordinates: $(hopf.coords)")
  println(io, "Parameters: $(hopf.parameters)")
  println(io, "V:")
  display(hopf.v)
  println(io, "ω: $(hopf.ω)")
  if hopf.stability === nothing
    println(io, "Stability: not calculated")
  else
    eigenvalues = sort(hopf.stability; by=x -> (floor(real(x); digits=6), -floor(imag(x); digits=6)), rev=true)
    println(io, "Stability:")
    display(eigenvalues)
  end
  if hopf.nmfm === nothing
    println(io, "Normal form: not calculated")
  else
    println(io, "Normal form:")
    println(io, hopf.nmfm)
  end
end

function point_to_genhopf(p::Hopf)
  if p.stability === nothing
    println("Need to calculate stability")
    return nothing
  else
    GenHopf(p.coords, p.parameters, p.stability, p.ω)
  end
end

function point_to_genhopf(jet, p::stst, τs)
  if p.stability === nothing
    println("Need to calculate stability")
    return nothing
  else
    freqs = DDEBifTool.closest_eigenvalues_to_imaginary_axis(p.stability)
    ω = first(sort(abs.(imag(freqs))))

    m = length(τs)
    φ = repeat(p.coords, 1, m)
    Δ(λ) = jet.Δ(φ, p.parameters, λ)

    λ = ω * im

    _, s, V = svd(Δ(λ))
    indxmin = last(findmin(s))
    V = V[:, indxmin]
    GenHopf(p.coords, p.parameters, V, ω, p.stability, nothing)
  end
end

function normalform(jet, hopf::GenHopf, τs)

  m = length(τs)
  φ = repeat(hopf.coords, 1, m)
  α = hopf.parameters

  λ = hopf.ω * im
  Δ(λ) = jet.Δ(φ, α, λ)
  Δ′(λ) = jet.Δ′(φ, α, λ)
  Δ′′(λ) = jet.Δ′′(φ, α, λ)

  _, s, V = svd(Δ(λ))
  indxmin = last(findmin(s)) # should this be abs(s) instead?
  q = V[:, indxmin]

  _, s, V = svd(transpose(Δ(λ)))
  indxmin = last(findmin(s))
  p = V[:, indxmin]

  # normalize
  p /= transpose(p) * (Δ′(λ) * q)

  # define border inverse of characteristic matrix
  # Δᴵᴺⱽ(λ, y) = ([Δ(λ) p; [transpose(q) 0]]\[y; 0])[1:end-1]
  function Δᴵᴺⱽ(λ, y)
    sol = [Δ(λ) p; [transpose(q) 0]] \ [y; 0]
    if abs(sol[end]) > 1e-10
      error("Border inverse failed")
    end
    sol[1:end-1]
  end

  function Aᴵᴺⱽ(λ, η, κ)
    ξ = Δᴵᴺⱽ(λ, η + κ * Δ′(λ) * q)
    γ = first(transpose(p) * (-Δ′(λ) * ξ + 0.5 * κ * Δ′′(λ) * q))
    θ -> exp(λ * θ) * (ξ + γ * q - κ * θ * q)
  end

  # multi-linear forms at bt point
  Ξ(h) = vcat([h(-τ) for τ ∈ τs]...)
  B(v₁, v₂) = jet.D2(φ, α, Ξ(v₁), Ξ(v₂))
  C(v₁, v₂, v₃) = jet.D3(φ, α, Ξ(v₁), Ξ(v₂), Ξ(v₃))
  D(v₁, v₂, v₃, v₄) = jet.D4(φ, α, Ξ(v₁), Ξ(v₂), Ξ(v₃), Ξ(v₄))
  E(v₁, v₂, v₃, v₄, v₅) = jet.D5(φ, α, Ξ(v₁), Ξ(v₂), Ξ(v₃), Ξ(v₄), Ξ(v₅))
  J₁ = jet.D01(φ, α)
  A₁(v₁, p₁) = jet.D11(φ, α, Ξ(v₁), p₁)
  B₁(v₁, v₂, p₁) = jet.D21(φ, α, Ξ(v₁), Ξ(v₂), p₁)
  C₁(v₁, v₂, v₃, p₁) = jet.D31(φ, α, Ξ(v₁), Ξ(v₂), Ξ(v₃), p₁)

  ϕ(θ) = exp(λ * θ) * q
  h2000(θ) = exp(2.0 * λ * θ) * (Δ(2.0 * λ) \ B(ϕ, ϕ))
  h1100(_) = Δ(0.0) \ B(ϕ, conj ∘ ϕ)

  c₁ = first(0.5 * transpose(p) * (B(conj ∘ ϕ, h2000) + 2 * B(ϕ, h1100) + C(ϕ, ϕ, conj ∘ ϕ)))
  ℓ₁ = real(c₁) / hopf.ω

  h3000(θ) = exp(3 * λ * θ) * (Δ(3 * λ) \ (3 * B(ϕ, h2000) + C(ϕ, ϕ, ϕ)))
  h2100 = Aᴵᴺⱽ(λ, B(conj ∘ ϕ, h2000) + 2 * B(ϕ, h1100) + C(ϕ, ϕ, conj ∘ ϕ), -2c₁)

  h2200(_) = Δ(0.0) \ (2B(conj ∘ ϕ, h2100) + 2B(ϕ, conj ∘ h2100) + B(conj ∘ h2000, h2000) + 2B(h1100, h1100) + C(ϕ, ϕ, conj ∘ h2000) + 4C(ϕ, conj ∘ ϕ, h1100) + C(conj ∘ ϕ, conj ∘ ϕ, h2000) + D(ϕ, ϕ, conj ∘ ϕ, conj ∘ ϕ))

  h3100(θ) = exp(2 * λ * θ) * (Δ(2 * λ) \ (B(conj ∘ ϕ, h3000) + 3B(ϕ, h2100) + 3B(h1100, h2000) + 3C(ϕ, conj ∘ ϕ, h2000) + 3C(ϕ, ϕ, h1100) + D(ϕ, ϕ, ϕ, conj ∘ ϕ))) - 6c₁ * (Δ(2λ) \ (Δ′(2λ) - θ * Δ(2λ))) * h2000(θ)

  c₂ = 1 / 12 * transpose(p) * (2B(conj ∘ ϕ, h3100) + 3B(ϕ, h2200) + B(conj ∘ h2000, h3000) + 6B(h1100, h2100) + 3B(conj ∘ h2100, h2000) +
                                6C(conj ∘ ϕ, h2000, h1100) + 6C(ϕ, conj ∘ ϕ, h2100) + C(conj ∘ ϕ, conj ∘ ϕ, h3000) + 3C(ϕ, ϕ, conj ∘ h2100) + 3C(ϕ, h2000, conj ∘ h2000) +
                                6C(ϕ, h1100, h1100) + 6D(ϕ, ϕ, conj ∘ ϕ, h1100) + 3D(ϕ, conj ∘ ϕ, conj ∘ ϕ, h2000) +
                                D(ϕ, ϕ, ϕ, conj ∘ h2000) + E(ϕ, ϕ, ϕ, conj ∘ ϕ, conj ∘ ϕ))
  ℓ₂ = real(c₂) / hopf.ω

  # parameter-related normal form coefficients
  v10 = [1.0; 0.0]
  v01 = [0.0; 1.0]
  h0010(_) = Δ(0.0) \ (J₁ * v10)
  h0001(_) = Δ(0.0) \ (J₁ * v01)

  γ110 = first(transpose(p) * (A₁(ϕ, v10) + B(ϕ, h0010)))
  γ101 = first(transpose(p) * (A₁(ϕ, v01) + B(ϕ, h0001)))

  h1010 = Aᴵᴺⱽ(λ, A₁(ϕ, v10) + B(ϕ, h0010), -γ110)
  h1001 = Aᴵᴺⱽ(λ, A₁(ϕ, v01) + B(ϕ, h0001), -γ101)

  h2010(θ) = exp(2λ * θ) * Δ(2λ) \ (A₁(h2000, v10) + 2B(ϕ, h1010) + B(h2000, h0010) + B₁(ϕ, ϕ, v10) + C(ϕ, ϕ, h0010)) - 2γ110 * (Δ(2λ) \ (Δ′(2λ) - θ * Δ(2λ))) * h2000(θ)
  h2001(θ) = exp(2λ * θ) * Δ(2λ) \ (A₁(h2000, v01) + 2B(ϕ, h1001) + B(h2000, h0001) + B₁(ϕ, ϕ, v01) + C(ϕ, ϕ, h0001)) - 2γ101 * (Δ(2λ) \ (Δ′(2λ) - θ * Δ(2λ)) * h2000(θ))


  h1110(θ) = Δ(0.0) \ (A₁(h1100, v10) + 2real(B(conj ∘ ϕ, h1010)) + B(h1100, h0010) + B₁(ϕ, conj ∘ ϕ, v10) + C(ϕ, conj ∘ ϕ, h0010)) - 2 * real(γ110) * (Δ(0.0) \ (Δ′(0.0) - θ * Δ(0.0))) * h1100(θ)
  h1101(θ) = Δ(0.0) \ (A₁(h1100, v01) + 2real(B(conj ∘ ϕ, h1001)) + B(h1100, h0001) + B₁(ϕ, conj ∘ ϕ, v01) + C(ϕ, conj ∘ ϕ, h0001)) - 2 * real(γ101) * (Δ(0.0) \ (Δ′(0.0) - θ * Δ(0.0))) * h1100(θ)

  γ210 = 0.5 * transpose(p) * (A₁(h2100, v10) + B(conj ∘ ϕ, h2010) + 2B(ϕ, h1110) + B(h2100, h0010) + B(h2000, conj ∘ h1010) + 2B(h1100, h1010) + B₁(h2000, conj ∘ ϕ, v10) + 2B₁(ϕ, h1100, v10) + 2C(ϕ, conj ∘ ϕ, h1010) + C(h2000, conj ∘ ϕ, h0010) + C(ϕ, ϕ, conj ∘ h1010) + 2C(ϕ, h1100, h0010) + C₁(ϕ, ϕ, conj ∘ ϕ, v10) + D(ϕ, ϕ, conj ∘ ϕ, h0010)) |> first
  γ201 = 0.5 * transpose(p) * (A₁(h2100, v01) + B(conj ∘ ϕ, h2001) + 2B(ϕ, h1101) + B(h2100, h0001) + B(h2000, conj ∘ h1001) + 2B(h1100, h1001) + B₁(h2000, conj ∘ ϕ, v01) + 2B₁(ϕ, h1100, v01) + 2C(ϕ, conj ∘ ϕ, h1001) + C(h2000, conj ∘ ϕ, h0001) + C(ϕ, ϕ, conj ∘ h1001) + 2C(ϕ, h1100, h0001) + C₁(ϕ, ϕ, conj ∘ ϕ, v01) + D(ϕ, ϕ, conj ∘ ϕ, h0001)) |> first

  # (;γ110,γ101,γ210,γ201,c₁,c₂,ℓ₁,ℓ₂,q,h2000,h1100,h2010,h2001,h1101,h1110,h2100,h2200,h3100,h1010,h1001,h0010,h0001)

  ω₂ = imag([γ110 γ101] * (real([γ110 γ101; γ210 γ201]) \ [0.0; 1.0])) |> first

  nmfm = GenHopfNormalform(
    γ110,
    γ101,
    γ210,
    γ201,
    c₁,
    c₂,
    ℓ₁,
    ℓ₂,
    q,
    h2000,
    h1100,
    h0010,
    h0001,
    ω₂)

  hopf = @set hopf.nmfm = nmfm
end

function normalform_beta(jet, hopf::GenHopf, τs)
  m = length(τs)
  φ = repeat(hopf.coords, 1, m)
  α = hopf.parameters

  λ = hopf.ω * im
  Δ(λ) = jet.Δ(φ, α, λ)
  Δ′(λ) = jet.Δ′(φ, α, λ)
  Δ′′(λ) = jet.Δ′′(φ, α, λ)
  Δ′′′(λ) = jet.Δ′′′(φ, α, λ)
  Δ⁴(λ) = jet.Δ⁴(φ, α, λ)

  _, s, V = svd(Δ(λ))
  indxmin = last(findmin(s))
  q = V[:, indxmin]

  _, s, V = svd(transpose(Δ(λ)))
  indxmin = last(findmin(s))
  p = V[:, indxmin]

  # normalize
  p /= transpose(p) * (Δ′(λ) * q)

  # define border inverse of characteristic matrix
  Δᴵᴺⱽ(λ, y) = ([Δ(λ) p; [q' 0]]\[y; 0])[1:end-1]

  function Aᴵᴺⱽ(λ, η, κ)
    ξ = Δᴵᴺⱽ(λ, η + κ * Δ′(λ) * q)
    γ = first(transpose(p) * (-Δ′(λ) * ξ + 0.5 * κ * Δ′′(λ) * q))
    θ -> exp(λ * θ) * (ξ + γ * q - κ * θ * q)
  end

  function A2ᴵᴺⱽ(λ, w, κ)
    ξ = Δᴵᴺⱽ(λ, Δ′(λ) * w - 0.5 * κ * Δ′′(λ) * q)
    γ = first(transpose(p) * (-Δ′(λ) * ξ + 0.5 * Δ′′(λ) * w - (1 / 6) * κ * Δ′′′(λ) * q))
    θ -> exp(λ * θ) * (ξ + γ * q - θ * w + 0.5 * κ * θ^2 * q)
  end

  function A3ᴵᴺⱽ(λ, w, η, κ)
    ξ = Δᴵᴺⱽ(λ, Δ′(λ) * w - 0.5 * Δ′′(λ) * η + (1 / 6) * Δ′′′(λ) * (κ * q))
    γ = first(transpose(p) * (-Δ′(λ) * ξ + 0.5 * Δ′′(λ) * w - (1 / 6) * Δ′′′(λ) * η + (1 / 24) * Δ⁴(λ) * (κ * q)))
    θ -> exp(λ * θ) * (ξ + γ * q - θ * w + 0.5 * θ^2 * η - (1 / 6) * θ^3 * (κ * q))
  end

  # multi-linear forms at bt point
  Ξ(h) = vcat([h(-τ) for τ ∈ τs]...)
  B(v₁, v₂) = jet.D2(φ, α, Ξ(v₁), Ξ(v₂))
  C(v₁, v₂, v₃) = jet.D3(φ, α, Ξ(v₁), Ξ(v₂), Ξ(v₃))
  D(v₁, v₂, v₃, v₄) = jet.D4(φ, α, Ξ(v₁), Ξ(v₂), Ξ(v₃), Ξ(v₄))
  E(v₁, v₂, v₃, v₄, v₅) = jet.D5(φ, α, Ξ(v₁), Ξ(v₂), Ξ(v₃), Ξ(v₄), Ξ(v₅))
  K(v₁, v₂, v₃, v₄, v₅, v₆) = jet.D6(φ, α, Ξ(v₁), Ξ(v₂), Ξ(v₃), Ξ(v₄), Ξ(v₅), Ξ(v₆))
  L(v₁, v₂, v₃, v₄, v₅, v₆, v₇) = jet.D7(φ, α, Ξ(v₁), Ξ(v₂), Ξ(v₃), Ξ(v₄), Ξ(v₅), Ξ(v₆), Ξ(v₇))

  J₁ = jet.D01(φ, α)
  A₁(v₁, p₁) = jet.D11(φ, α, Ξ(v₁), p₁)
  B₁(v₁, v₂, p₁) = jet.D21(φ, α, Ξ(v₁), Ξ(v₂), p₁)
  C₁(v₁, v₂, v₃, p₁) = jet.D31(φ, α, Ξ(v₁), Ξ(v₂), Ξ(v₃), p₁)
  D₁(v₁, v₂, v₃, v₄, p₁) = jet.D41(φ, α, Ξ(v₁), Ξ(v₂), Ξ(v₃), Ξ(v₄), p₁)
  E₁(v₁, v₂, v₃, v₄, v₅, p₁) = jet.D51(φ, α, Ξ(v₁), Ξ(v₂), Ξ(v₃), Ξ(v₄), Ξ(v₅), p₁)

  J₂(p₁, p₂) = jet.D02(φ, α, p₁, p₂)
  A₂(v₁, p₁, p₂) = jet.D12(φ, α, Ξ(v₁), p₁, p₂)
  B₂(v₁, v₂, p₁, p₂) = jet.D22(φ, α, Ξ(v₁), Ξ(v₂), p₁, p₂)
  C₂(v₁, v₂, v₃, p₁, p₂) = jet.D32(φ, α, Ξ(v₁), Ξ(v₂), Ξ(v₃), p₁, p₂)

  J₃(p₁, p₂, p₃) = jet.D03(φ, α, p₁, p₂, p₃)
  A₃(v₁, p₁, p₂, p₃) = jet.D13(φ, α, Ξ(v₁), p₁, p₂, p₃)
  B₃(v₁, v₂, p₁, p₂, p₃) = jet.D32(φ, α, Ξ(v₁), Ξ(v₂), p₁, p₂, p₃)
  C₃(v₁, v₂, v₃, p₁, p₂, p₃) = jet.D33(φ, α, Ξ(v₁), Ξ(v₂), Ξ(v₃), p₁, p₂, p₃)


  ϕ(θ) = exp(λ * θ) * q
  h2000(θ) = exp(2.0 * λ * θ) * (Δ(2.0 * λ) \ B(ϕ, ϕ))
  h1100(_) = Δ(0.0) \ B(ϕ, conj ∘ ϕ)

  c₁ = first(0.5 * transpose(p) * (B(conj ∘ ϕ, h2000) + 2 * B(ϕ, h1100) + C(ϕ, ϕ, conj ∘ ϕ)))
  ℓ₁ = real(c₁) / hopf.ω

  h3000(θ) = exp(3 * λ * θ) * (Δ(3 * λ) \ (3 * B(ϕ, h2000) + C(ϕ, ϕ, ϕ)))
  h2100 = Aᴵᴺⱽ(λ, B(conj ∘ ϕ, h2000) + 2 * B(ϕ, h1100) + C(ϕ, ϕ, conj ∘ ϕ), -2c₁)

  h2200(_) = Δ(0.0) \ (2B(conj ∘ ϕ, h2100) + 2B(ϕ, conj ∘ h2100) + B(conj ∘ h2000, h2000) + 2B(h1100, h1100) + C(ϕ, ϕ, conj ∘ h2000) + 4C(ϕ, conj ∘ ϕ, h1100) + C(conj ∘ ϕ, conj ∘ ϕ, h2000) + D(ϕ, ϕ, conj ∘ ϕ, conj ∘ ϕ))

  h3100(θ) = exp(2 * λ * θ) * (Δ(2 * λ) \ (B(conj ∘ ϕ, h3000) + 3B(ϕ, h2100) + 3B(h1100, h2000) + 3C(ϕ, conj ∘ ϕ, h2000) + 3C(ϕ, ϕ, h1100) + D(ϕ, ϕ, ϕ, conj ∘ ϕ))) - 6c₁ * (Δ(2λ) \ (Δ′(2λ) - θ * Δ(2λ))) * h2000(θ)
  h4000(θ) = exp(4 * λ * θ) * Δ(4λ) \ (4B(ϕ, h3000) + 3B(h2000, h2000) + 6C(ϕ, ϕ, h2000) + D(ϕ, ϕ, ϕ, ϕ))

  M3200 = (2B(conj ∘ ϕ, h3100) + 3B(ϕ, h2200) + B(conj ∘ h2000, h3000) + 6B(h1100, h2100) + 3B(conj ∘ h2100, h2000) +
           6C(conj ∘ ϕ, h2000, h1100) + 6C(ϕ, conj ∘ ϕ, h2100) + C(conj ∘ ϕ, conj ∘ ϕ, h3000) + 3C(ϕ, ϕ, conj ∘ h2100) + 3C(ϕ, h2000, conj ∘ h2000) +
           6C(ϕ, h1100, h1100) + 6D(ϕ, ϕ, conj ∘ ϕ, h1100) + 3D(ϕ, conj ∘ ϕ, conj ∘ ϕ, h2000) +
           D(ϕ, ϕ, ϕ, conj ∘ h2000) + E(ϕ, ϕ, ϕ, conj ∘ ϕ, conj ∘ ϕ))

  c₂ = 1 / 12 * transpose(p) * M3200

  ℓ₂ = real(c₂) / hopf.ω

  h3200(θ) = Aᴵᴺⱽ(λ, M3200, -12c₂)(θ) - 6 * im * imag(c₁) * A2ᴵᴺⱽ(λ, h2100(0), -2c₁)(θ)
  h4100(θ) = (exp(3 * λ * θ) * Δ(3λ) \ (4B(ϕ, h3100) + B(conj ∘ ϕ, h4000) + 4B(h1100, h3000) + 6B(h2000, h2100) + 6C(ϕ, ϕ, h2100)
                                        + 4C(ϕ, conj ∘ ϕ, h3000) + 12C(ϕ, h1100, h2000) + 3C(conj ∘ ϕ, h2000, h2000) + 4D(ϕ, ϕ, ϕ, h1100)
                                        + 6D(ϕ, ϕ, conj ∘ ϕ, h2000) + E(ϕ, ϕ, ϕ, ϕ, conj ∘ ϕ)) - 12 * c₁ * (Δ(3λ) \ (Δ′(3λ) - θ * Δ(3λ))) * h3000(θ))

  h5000(θ) = (exp(5 * λ * θ) * Δ(5λ) \ (5B(ϕ, h4000) + 10B(h2000, h3000) + 10C(ϕ, ϕ, h3000)
                                        + 15C(ϕ, h2000, h2000) + 10D(ϕ, ϕ, ϕ, h2000) + E(ϕ, ϕ, ϕ, ϕ, ϕ)))

  h6000(θ) = (exp(6 * λ * θ) * Δ(6λ) \ (6B(ϕ, h5000) + 15B(h2000, h4000) + 10B(h3000, h3000)
                                        + 15C(ϕ, ϕ, h4000) + 60C(ϕ, h2000, h3000) + 15C(h2000, h2000, h2000)
                                        + 20D(ϕ, ϕ, ϕ, h3000) + 45D(ϕ, ϕ, h2000, h2000) + 15E(ϕ, ϕ, ϕ, ϕ, h2000)
                                        + K(ϕ, ϕ, ϕ, ϕ, ϕ, ϕ)))

  h5100(θ) = (exp(4 * λ * θ) * Δ(4λ) \ (5B(ϕ, h4100) + B(conj ∘ (ϕ), h5000) + 5B(h1100, h4000)
                                        + 10B(h2000, h3100) + 10B(h2100, h3000) + 10C(ϕ, ϕ, h3100) + 5C(ϕ, conj ∘ (ϕ), h4000)
                                        + 20C(ϕ, h1100, h3000) + 30C(ϕ, h2000, h2100) + 10C(conj ∘ (ϕ), h2000, h3000)
                                        + 15C(h1100, h2000, h2000) + 10D(ϕ, ϕ, ϕ, h2100) + 10D(ϕ, ϕ, conj ∘ (ϕ), h3000)
                                        + 30D(ϕ, ϕ, h1100, h2000) + 15D(ϕ, conj ∘ (ϕ), h2000, h2000) + 5E(ϕ, ϕ, ϕ, ϕ, h1100)
                                        + 10E(ϕ, ϕ, ϕ, conj ∘ (ϕ), h2000) + K(ϕ, ϕ, ϕ, ϕ, ϕ, conj ∘ (ϕ)))
              -
              20 * c₁ * (Δ(4λ) \ (Δ′(4λ) - θ * Δ(4λ))) * h4000(θ))

  h4200(θ) = (exp(2 * λ * θ) * Δ(2λ) \ (4B(ϕ, h3200) + 2B(conj ∘ (ϕ), h4100) + B(conj ∘ (h2000), h4000) + 8B(h1100, h3100) + 4B(conj ∘ (h2100), h3000)
                                        + 6B(h2000, h2200) + 6B(h2100, h2100) + 6C(ϕ, ϕ, h2200) + 8C(ϕ, conj ∘ (ϕ), h3100) + 4C(ϕ, conj ∘ (h2000), h3000)
                                        + 24C(ϕ, h1100, h2100) + 12C(ϕ, conj ∘ (h2100), h2000) + C(conj ∘ (ϕ), conj ∘ (ϕ), h4000) + 8C(conj ∘ (ϕ), h1100, h3000)
                                        + 12C(conj ∘ (ϕ), h2000, h2100) + 3C(conj ∘ (h2000), h2000, h2000) + 12C(h1100, h1100, h2000) + 4D(ϕ, ϕ, ϕ, conj ∘ (h2100))
                                        + 12D(ϕ, ϕ, conj ∘ (ϕ), h2100) + 6D(ϕ, ϕ, conj ∘ (h2000), h2000) + 12D(ϕ, ϕ, h1100, h1100) + 4D(ϕ, conj ∘ (ϕ), conj ∘ (ϕ), h3000)
                                        + 24D(ϕ, conj ∘ (ϕ), h1100, h2000) + 3D(conj ∘ (ϕ), conj ∘ (ϕ), h2000, h2000) + E(ϕ, ϕ, ϕ, ϕ, conj ∘ (h2000))
                                        + 8E(ϕ, ϕ, ϕ, conj ∘ (ϕ), h1100) + 6E(ϕ, ϕ, conj ∘ (ϕ), conj ∘ (ϕ), h2000) + K(ϕ, ϕ, ϕ, ϕ, conj ∘ (ϕ), conj ∘ (ϕ)))
              -
              48 * c₂ * (Δ(2λ) \ (Δ′(2λ) - θ * Δ(2λ))) * h2000(θ)
              -
              8 * (3 * c₁ + conj(c₁)) * exp(2λ * θ) * (Δ(2λ) \ ((Δ′(2λ) - θ * Δ(2λ)) * h3100(0) + 3 * c₁ * (Δ′′(2λ) - θ^2 * Δ(2λ)) * h2000(0))))

  h3300(θ) = (Δ(0) \ (3B(ϕ, conj ∘ (h3200)) + 3B(conj ∘ (ϕ), h3200) + 3B(conj ∘ (h2000), h3100) + B(conj ∘ (h3000), h3000)
                      + 9B(h1100, h2200) + 9B(h2100, conj ∘ (h2100)) + 3B(conj ∘ (h3100), h2000) + 3C(ϕ, ϕ, conj ∘ (h3100))
                      + 9C(ϕ, conj ∘ (ϕ), h2200) + 9C(ϕ, conj ∘ (h2000), h2100) + 3C(ϕ, conj ∘ (h3000), h2000) + 18C(ϕ, h1100, conj ∘ (h2100))
                      + 3C(conj ∘ (ϕ), conj ∘ (ϕ), h3100) + 3C(conj ∘ (ϕ), conj ∘ (h2000), h3000) + 18C(conj ∘ (ϕ), h1100, h2100)
                      + 9C(conj ∘ (ϕ), conj ∘ (h2100), h2000) + 9C(conj ∘ (h2000), h1100, h2000) + 6C(h1100, h1100, h1100)
                      + D(ϕ, ϕ, ϕ, conj ∘ (h3000)) + 9D(ϕ, ϕ, conj ∘ (ϕ), conj ∘ (h2100)) + 9D(ϕ, ϕ, conj ∘ (h2000), h1100)
                      + 9D(ϕ, conj ∘ (ϕ), conj ∘ (ϕ), h2100) + 9D(ϕ, conj ∘ (ϕ), conj ∘ (h2000), h2000) + 18D(ϕ, conj ∘ (ϕ), h1100, h1100)
                      + D(conj ∘ (ϕ), conj ∘ (ϕ), conj ∘ (ϕ), h3000) + 9D(conj ∘ (ϕ), conj ∘ (ϕ), h1100, h2000) + 3E(ϕ, ϕ, ϕ, conj ∘ (ϕ), conj ∘ (h2000))
                      + 9E(ϕ, ϕ, conj ∘ (ϕ), conj ∘ (ϕ), h1100) + 3E(ϕ, conj ∘ (ϕ), conj ∘ (ϕ), conj ∘ (ϕ), h2000)
                      + K(ϕ, ϕ, ϕ, conj ∘ (ϕ), conj ∘ (ϕ), conj ∘ (ϕ)))
              -
              72 * real(c₂) * (Δ(0) \ (Δ′(0) - θ * Δ(0))) * h1100(θ))

  M4300 = (4B(ϕ, h3300) + 3B(conj ∘ (ϕ), h4200) + 3B(conj ∘ (h2000), h4100) + B(conj ∘ (h3000), h4000) + 12B(h1100, h3200)
           + 12B(conj ∘ (h2100), h3100) + 4B(conj ∘ (h3100), h3000) + 6B(h2000, conj ∘ (h3200)) + 18B(h2100, h2200) + 6C(ϕ, ϕ, conj ∘ (h3200))
           + 12C(ϕ, conj ∘ (ϕ), h3200) + 12C(ϕ, conj ∘ (h2000), h3100) + 4C(ϕ, conj ∘ (h3000), h3000) + 36C(ϕ, h1100, h2200)
           + 36C(ϕ, conj ∘ (h2100), h2100) + 12C(ϕ, conj ∘ (h3100), h2000) + 3C(conj ∘ (ϕ), conj ∘ (ϕ), h4100) + 3C(conj ∘ (ϕ), conj ∘ (h2000), h4000)
           + 24C(conj ∘ (ϕ), h1100, h3100) + 12C(conj ∘ (ϕ), conj ∘ (h2100), h3000) + 18C(conj ∘ (ϕ), h2000, h2200) + 18C(conj ∘ (ϕ), h2100, h2100)
           + 12C(conj ∘ (h2000), h1100, h3000) + 18C(conj ∘ (h2000), h2000, h2100) + 3C(conj ∘ (h3000), h2000, h2000) + 36C(h1100, h1100, h2100)
           + 36C(h1100, conj ∘ (h2100), h2000) + 4D(ϕ, ϕ, ϕ, conj ∘ (h3100)) + 18D(ϕ, ϕ, conj ∘ (ϕ), h2200) + 18D(ϕ, ϕ, conj ∘ (h2000), h2100)
           + 6D(ϕ, ϕ, conj ∘ (h3000), h2000) + 36D(ϕ, ϕ, h1100, conj ∘ (h2100)) + 12D(ϕ, conj ∘ (ϕ), conj ∘ (ϕ), h3100)
           + 12D(ϕ, conj ∘ (ϕ), conj ∘ (h2000), h3000) + 72D(ϕ, conj ∘ (ϕ), h1100, h2100) + 36D(ϕ, conj ∘ (ϕ), conj ∘ (h2100), h2000)
           + 36D(ϕ, conj ∘ (h2000), h1100, h2000) + 24D(ϕ, h1100, h1100, h1100) + D(conj ∘ (ϕ), conj ∘ (ϕ), conj ∘ (ϕ), h4000)
           + 12D(conj ∘ (ϕ), conj ∘ (ϕ), h1100, h3000) + 18D(conj ∘ (ϕ), conj ∘ (ϕ), h2000, h2100) + 9D(conj ∘ (ϕ), conj ∘ (h2000), h2000, h2000)
           + 36D(conj ∘ (ϕ), h1100, h1100, h2000) + E(ϕ, ϕ, ϕ, ϕ, conj ∘ (h3000)) + 12E(ϕ, ϕ, ϕ, conj ∘ (ϕ), conj ∘ (h2100))
           + 12E(ϕ, ϕ, ϕ, conj ∘ (h2000), h1100) + 18E(ϕ, ϕ, conj ∘ (ϕ), conj ∘ (ϕ), h2100) + 18E(ϕ, ϕ, conj ∘ (ϕ), conj ∘ (h2000), h2000)
           + 36E(ϕ, ϕ, conj ∘ (ϕ), h1100, h1100) + 4E(ϕ, conj ∘ (ϕ), conj ∘ (ϕ), conj ∘ (ϕ), h3000) + 36E(ϕ, conj ∘ (ϕ), conj ∘ (ϕ), h1100, h2000)
           + 3E(conj ∘ (ϕ), conj ∘ (ϕ), conj ∘ (ϕ), h2000, h2000) + 3K(ϕ, ϕ, ϕ, ϕ, conj ∘ (ϕ), conj ∘ (h2000)) + 12K(ϕ, ϕ, ϕ, conj ∘ (ϕ), conj ∘ (ϕ), h1100)
           + 6K(ϕ, ϕ, conj ∘ (ϕ), conj ∘ (ϕ), conj ∘ (ϕ), h2000) + L(ϕ, ϕ, ϕ, ϕ, conj ∘ (ϕ), conj ∘ (ϕ), conj ∘ (ϕ)))

  c₃ = (1 / 144) * transpose(p) * M4300
  ℓ₃ = real(c₃) / hopf.ω

  h4300(θ) = (Aᴵᴺⱽ(λ, M4300, -144c₃)(θ) - 72 * (2c₂ + conj(c₂)) * A2ᴵᴺⱽ(λ, h2100(0), -2c₁)(θ)
              -
              12 * im * imag(c₁) * A3ᴵᴺⱽ(λ, h3200(0), -(12c₂ * q + 6 * im * imag(c₁) * h2100(0)), 12 * im * imag(c₁) * c₁)(θ))

  h7000(θ) = (exp(7 * λ * θ) * Δ(7λ) \ (7B(ϕ, h6000) + 21B(h2000, h5000) + 35B(h3000, h4000)
                                        + 21C(ϕ, ϕ, h5000) + 105C(ϕ, h2000, h4000) + 70C(ϕ, h3000, h3000)
                                        + 105C(h2000, h2000, h3000) + 35D(ϕ, ϕ, ϕ, h4000)
                                        + 210D(ϕ, ϕ, h2000, h3000) + 105D(ϕ, h2000, h2000, h2000)
                                        + 35E(ϕ, ϕ, ϕ, ϕ, h3000) + 105E(ϕ, ϕ, ϕ, h2000, h2000)
                                        + 21K(ϕ, ϕ, ϕ, ϕ, ϕ, h2000) + L(ϕ, ϕ, ϕ, ϕ, ϕ, ϕ, ϕ)))


  h6100(θ) = (exp(5 * λ * θ) * Δ(5λ) \ (6B(ϕ, h5100) + B(conj ∘ (ϕ), h6000) + 6B(h1100, h5000)
                                        + 15B(h2000, h4100) + 15B(h2100, h4000) + 20B(h3000, h3100)
                                        + 15C(ϕ, ϕ, h4100) + 6C(ϕ, conj ∘ (ϕ), h5000) + 30C(ϕ, h1100, h4000)
                                        + 60C(ϕ, h2000, h3100) + 60C(ϕ, h2100, h3000) + 15C(conj ∘ (ϕ), h2000, h4000)
                                        + 10C(conj ∘ (ϕ), h3000, h3000) + 60C(h1100, h2000, h3000)
                                        + 45C(h2000, h2000, h2100) + 20D(ϕ, ϕ, ϕ, h3100) + 15D(ϕ, ϕ, conj ∘ (ϕ), h4000)
                                        + 60D(ϕ, ϕ, h1100, h3000) + 90D(ϕ, ϕ, h2000, h2100)
                                        + 60D(ϕ, conj ∘ (ϕ), h2000, h3000) + 90D(ϕ, h1100, h2000, h2000)
                                        + 15D(conj ∘ (ϕ), h2000, h2000, h2000) + 15E(ϕ, ϕ, ϕ, ϕ, h2100)
                                        + 20E(ϕ, ϕ, ϕ, conj ∘ (ϕ), h3000) + 60E(ϕ, ϕ, ϕ, h1100, h2000)
                                        + 45E(ϕ, ϕ, conj ∘ (ϕ), h2000, h2000) + 6K(ϕ, ϕ, ϕ, ϕ, ϕ, h1100)
                                        + 15K(ϕ, ϕ, ϕ, ϕ, conj ∘ (ϕ), h2000) + L(ϕ, ϕ, ϕ, ϕ, ϕ, ϕ, conj ∘ (ϕ)))
              -
              30 * c₁ * (Δ(5λ) \ (Δ′(5λ) - θ * Δ(5λ))) * h5000(θ))


  h5200(θ) = (exp(3 * λ * θ) * Δ(3λ) \ (5B(ϕ, h4200) + 2B(conj ∘ (ϕ), h5100)
                                        + B(conj ∘ (h2000), h5000) + 10B(h1100, h4100) + 5B(conj ∘ (h2100), h4000)
                                        + 10B(h2000, h3200) + 20B(h2100, h3100) + 10B(h2200, h3000)
                                        + 10C(ϕ, ϕ, h3200) + 10C(ϕ, conj ∘ (ϕ), h4100) + 5C(ϕ, conj ∘ (h2000), h4000)
                                        + 40C(ϕ, h1100, h3100) + 20C(ϕ, conj ∘ (h2100), h3000)
                                        + 30C(ϕ, h2000, h2200) + 30C(ϕ, h2100, h2100) + C(conj ∘ (ϕ), conj ∘ (ϕ), h5000)
                                        + 10C(conj ∘ (ϕ), h1100, h4000) + 20C(conj ∘ (ϕ), h2000, h3100) + 20C(conj ∘ (ϕ), h2100, h3000)
                                        + 10C(conj ∘ (h2000), h2000, h3000) + 20C(h1100, h1100, h3000)
                                        + 60C(h1100, h2000, h2100) + 15C(conj ∘ (h2100), h2000, h2000)
                                        + 10D(ϕ, ϕ, ϕ, h2200) + 20D(ϕ, ϕ, conj ∘ (ϕ), h3100) + 10D(ϕ, ϕ, conj ∘ (h2000), h3000)
                                        + 60D(ϕ, ϕ, h1100, h2100) + 30D(ϕ, ϕ, conj ∘ (h2100), h2000)
                                        + 5D(ϕ, conj ∘ (ϕ), conj ∘ (ϕ), h4000) + 40D(ϕ, conj ∘ (ϕ), h1100, h3000) + 60D(ϕ, conj ∘ (ϕ), h2000, h2100)
                                        + 15D(ϕ, conj ∘ (h2000), h2000, h2000) + 60D(ϕ, h1100, h1100, h2000)
                                        + 10D(conj ∘ (ϕ), conj ∘ (ϕ), h2000, h3000) + 30D(conj ∘ (ϕ), h1100, h2000, h2000)
                                        + 5E(ϕ, ϕ, ϕ, ϕ, conj ∘ (h2100)) + 20E(ϕ, ϕ, ϕ, conj ∘ (ϕ), h2100) + 10E(ϕ, ϕ, ϕ, conj ∘ (h2000), h2000)
                                        + 20E(ϕ, ϕ, ϕ, h1100, h1100) + 10E(ϕ, ϕ, conj ∘ (ϕ), conj ∘ (ϕ), h3000)
                                        + 60E(ϕ, ϕ, conj ∘ (ϕ), h1100, h2000) + 15E(ϕ, conj ∘ (ϕ), conj ∘ (ϕ), h2000, h2000)
                                        + K(ϕ, ϕ, ϕ, ϕ, ϕ, conj ∘ (h2000)) + 10K(ϕ, ϕ, ϕ, ϕ, conj ∘ (ϕ), h1100) + 10K(ϕ, ϕ, ϕ, conj ∘ (ϕ), conj ∘ (ϕ), h2000)
                                        + L(ϕ, ϕ, ϕ, ϕ, ϕ, conj ∘ (ϕ), conj ∘ (ϕ)))
              -
              120 * c₂ * (Δ(3λ) \ (Δ′(3λ) - θ * Δ(3λ))) * h3000(θ)
              -
              (40 * c₁ + 10 * conj(c₁)) * exp(3λ * θ) * (Δ(3λ) \ ((Δ′(3λ) - θ * Δ(3λ)) * h4100(0) + 6 * c₁ * (Δ′′(3λ) - θ^2 * Δ(3λ)) * h3000(0))))



  # parameter-related normal form coefficients
  e₁ = [1.0; 0.0]
  e₂ = [0.0; 1.0]

  Γ₁(u) = A₁(u, e₁) + B(u, (θ -> Δ(0) \ (J₁ * e₁)))
  Γ₂(u) = A₁(u, e₂) + B(u, (θ -> Δ(0) \ (J₁ * e₂)))

  Λ₁(u, v, w) = Γ₁(u) + 2B(v, Aᴵᴺⱽ(λ, Γ₁(w), 0)) + B₁(v, w, e₁) + C(v, w, (θ -> Δ(0) \ (J₁ * e₁)))
  Λ₂(u, v, w) = Γ₂(u) + 2B(v, Aᴵᴺⱽ(λ, Γ₂(w), 0)) + B₁(v, w, e₂) + C(v, w, (θ -> Δ(0) \ (J₁ * e₂)))
  Π₁(u, v, w) = Γ₁(u) + 2real(B(v, Aᴵᴺⱽ(λ, Γ₁(w), 0))) + B₁(v, w, e₁) + C(v, w, (θ -> Δ(0) \ (J₁ * e₁)))
  Π₂(u, v, w) = Γ₂(u) + 2real(B(v, Aᴵᴺⱽ(λ, Γ₂(w), 0))) + B₁(v, w, e₂) + C(v, w, (θ -> Δ(0) \ (J₁ * e₂)))

  P11 = real(first(transpose(p) * Γ₁(ϕ)))
  P12 = real(first(transpose(p) * Γ₂(ϕ)))

  rP2 = real(first(transpose(p) * (-4B(ϕ, (θ -> Δ(0) \ real(im * B(conj ∘ ϕ, Aᴵᴺⱽ(λ, zeros(length(q)), 1)))))
                                   -
                                   2 * im * B(conj ∘ ϕ, (θ -> ((Δ(2λ) \ (Δ′(2λ) - θ * Δ(2λ))) * h2000(θ) + exp(2.0 * λ * θ) * (Δ(2.0 * λ) \ B(ϕ, Aᴵᴺⱽ(λ, zeros(length(q)), 1))))))
                                   +
                                   im * B(h2000, conj ∘ Aᴵᴺⱽ(λ, zeros(length(q)), 1)) - 2 * im * B(h1100, Aᴵᴺⱽ(λ, zeros(length(q)), 1)) + im * C(ϕ, ϕ, conj ∘ Aᴵᴺⱽ(λ, zeros(length(q)), 1))
                                   -
                                   2 * im * C(ϕ, conj ∘ ϕ, Aᴵᴺⱽ(λ, zeros(length(q)), 1)))))

  P21 = (0.5 * real(first(transpose(p) * (Γ₁(h2100) + 2B(ϕ, (θ -> (Δ(0) \ Π₁(h1100, conj ∘ ϕ, ϕ)))) + B(conj ∘ ϕ, (θ -> (exp(2.0 * λ * θ) * (Δ(2λ) \ Λ₁(h2000, ϕ, ϕ)))))
                                          + B(h2000, conj ∘ Aᴵᴺⱽ(λ, Γ₁(ϕ), 0)) + 2B(h1100, Aᴵᴺⱽ(λ, Γ₁(ϕ), 0)) + 2B₁(ϕ, h1100, e₁) + B₁(conj ∘ ϕ, h2000, e₁)
                                          + C(ϕ, ϕ, conj ∘ Aᴵᴺⱽ(λ, Γ₁(ϕ), 0)) + 2C(ϕ, conj ∘ ϕ, Aᴵᴺⱽ(λ, Γ₁(ϕ), 0)) + 2C(ϕ, h1100, (θ -> Δ(0) \ (J₁ * e₁)))
                                          + C(conj ∘ ϕ, h2000, (θ -> Δ(0) \ (J₁ * e₁))) + C₁(ϕ, ϕ, conj ∘ ϕ, e₁) + D(ϕ, ϕ, conj ∘ ϕ, (θ -> Δ(0) \ (J₁ * e₁)))
  ))) + imag(first(transpose(p) * Γ₁(ϕ))) * 0.5 * rP2)


  P22 = (0.5 * real(first(transpose(p) * (Γ₂(h2100) + 2B(ϕ, (θ -> Δ(0) \ Π₂(h1100, conj ∘ ϕ, ϕ))) + B(conj ∘ ϕ, (θ -> exp(2.0 * λ * θ) * (Δ(2λ) \ Λ₂(h2000, ϕ, ϕ))))
                                          + B(h2000, conj ∘ Aᴵᴺⱽ(λ, Γ₂(ϕ), 0)) + 2B(h1100, Aᴵᴺⱽ(λ, Γ₂(ϕ), 0)) + 2B₁(ϕ, h1100, e₂) + B₁(conj ∘ ϕ, h2000, e₂)
                                          + C(ϕ, ϕ, conj ∘ Aᴵᴺⱽ(λ, Γ₂(ϕ), 0)) + 2C(ϕ, conj ∘ ϕ, Aᴵᴺⱽ(λ, Γ₂(ϕ), 0)) + 2C(ϕ, h1100, (θ -> Δ(0) \ (J₁ * e₂)))
                                          + C(conj ∘ ϕ, h2000, (θ -> Δ(0) \ (J₁ * e₂))) + C₁(ϕ, ϕ, conj ∘ ϕ, e₂) + D(ϕ, ϕ, conj ∘ ϕ, (θ -> Δ(0) \ (J₁ * e₂))))
  ))) + imag(first(transpose(p) * Γ₂(ϕ))) * 0.5 * rP2

  P = [P11 P12; P21 P22]

  Q210 = 0.5 * real(first(transpose(p) * (4B(ϕ, (θ -> ((Δ(0) \ (Δ′(0) - θ * Δ(0))) * h1100(θ) - (Δ(0) \ real(B(conj ∘ ϕ, Aᴵᴺⱽ(λ, zeros(length(q)), 1)))))))
                                          + 2 * B(conj ∘ ϕ, (θ -> ((Δ(2λ) \ (Δ′(2λ) - θ * Δ(2λ))) * h2000(θ) + (exp(2.0 * λ * θ) * (Δ(2.0 * λ) \ B(ϕ, Aᴵᴺⱽ(λ, zeros(length(q)), 1)))))))
                                          + B(h2000, conj ∘ Aᴵᴺⱽ(λ, zeros(length(q)), 1)) + 2 * B(h1100, Aᴵᴺⱽ(λ, zeros(length(q)), 1)) + C(ϕ, ϕ, conj ∘ Aᴵᴺⱽ(λ, zeros(length(q)), 1))
                                          + 2 * C(ϕ, conj ∘ ϕ, Aᴵᴺⱽ(λ, zeros(length(q)), 1))
  )))


  K10 = P \ [1; Q210]
  K01 = P \ [0; 1]


  h0010(_) = Δ(0.0) \ (J₁ * K10)
  h0001(_) = Δ(0.0) \ (J₁ * K01)

  b110 = imag(first(transpose(p) * (A₁(ϕ, K10) + B(ϕ, h0010))))
  b101 = imag(first(transpose(p) * (A₁(ϕ, K01) + B(ϕ, h0001))))

  h1010 = Aᴵᴺⱽ(λ, A₁(ϕ, K10) + B(ϕ, h0010), -(1 + im * b110))
  h1001 = Aᴵᴺⱽ(λ, A₁(ϕ, K01) + B(ϕ, h0001), -im * b101)

  h2010(θ) = exp(2λ * θ) * Δ(2λ) \ (A₁(h2000, K10) + 2B(ϕ, h1010) + B(h2000, h0010) + B₁(ϕ, ϕ, K10) + C(ϕ, ϕ, h0010)) - 2(1 + im * b110) * (Δ(2λ) \ (Δ′(2λ) - θ * Δ(2λ))) * h2000(θ)
  h2001(θ) = exp(2λ * θ) * Δ(2λ) \ (A₁(h2000, K01) + 2B(ϕ, h1001) + B(h2000, h0001) + B₁(ϕ, ϕ, K01) + C(ϕ, ϕ, h0001)) - 2 * im * b101 * (Δ(2λ) \ (Δ′(2λ) - θ * Δ(2λ))) * h2000(θ)

  h1110(θ) = Δ(0.0) \ (A₁(h1100, K10) + 2real(B(conj ∘ ϕ, h1010)) + B(h1100, h0010) + B₁(ϕ, conj ∘ ϕ, K10) + C(ϕ, conj ∘ ϕ, h0010)) - 2 * (Δ(0.0) \ (Δ′(0.0) - θ * Δ(0.0))) * h1100(θ)
  h1101(_) = Δ(0.0) \ (A₁(h1100, K01) + 2real(B(conj ∘ ϕ, h1001)) + B(h1100, h0001) + B₁(ϕ, conj ∘ ϕ, K01) + C(ϕ, conj ∘ ϕ, h0001))

  M2101 = (A₁(h2100, K10) + B(conj ∘ ϕ, h2010) + 2B(ϕ, h1110) + B(h2100, h0010) + B(h2000, conj ∘ h1010)
           + 2B(h1100, h1010) + B₁(h2000, conj ∘ ϕ, K10) + 2B₁(ϕ, h1100, K10) + 2C(ϕ, conj ∘ ϕ, h1010)
           + C(h2000, conj ∘ ϕ, h0010) + C(ϕ, ϕ, conj ∘ h1010) + 2C(ϕ, h1100, h0010) + C₁(ϕ, ϕ, conj ∘ ϕ, K10) + D(ϕ, ϕ, conj ∘ ϕ, h0010))

  M2110 = (A₁(h2100, K01) + B(conj ∘ ϕ, h2001) + 2B(ϕ, h1101) + B(h2100, h0001) + B(h2000, conj ∘ h1001)
           + 2B(h1100, h1001) + B₁(h2000, conj ∘ ϕ, K01) + 2B₁(ϕ, h1100, K01) + 2C(ϕ, conj ∘ ϕ, h1001)
           + C(h2000, conj ∘ ϕ, h0001) + C(ϕ, ϕ, conj ∘ h1001) + 2C(ϕ, h1100, h0001) + C₁(ϕ, ϕ, conj ∘ ϕ, K01) + D(ϕ, ϕ, conj ∘ ϕ, h0001))

  b201 = 0.5 * imag(first(transpose(p) * (M2101)))
  b210 = 0.5 * imag(first(transpose(p) * (M2110)))

  h2101(θ) = Aᴵᴺⱽ(λ, M2101, -2 * (1 + im * b201))(θ) - im * b101 * A2ᴵᴺⱽ(λ, h2100(0), -2 * c₁)(θ) - 2 * c₁ * A2ᴵᴺⱽ(λ, h1001(0), -im * b101)(θ)
  h2110(θ) = Aᴵᴺⱽ(λ, M2110, -2 * im * b210)(θ) - (3 + im * b110) * A2ᴵᴺⱽ(λ, h2100(0), -2 * c₁)(θ) - 2 * c₁ * A2ᴵᴺⱽ(λ, h1010(0), -(1 + im * b110))(θ)

  h3001(θ) = (exp(3λ * θ) * Δ(3λ) \ (A₁(h3000, K01) + 3B(ϕ, h2001) + B(h0001, h3000) + 3B(h1001, h2000) + 3B₁(ϕ, h2000, K01) + 3C(ϕ, ϕ, h1001)
                                     + 3C(ϕ, h0001, h2000) + C₁(ϕ, ϕ, ϕ, K01) + D(ϕ, ϕ, ϕ, h0001)) - 3 * im * b101 * (Δ(3λ) \ (Δ′(3λ) - θ * Δ(3λ))) * h3000(θ))
  h3010(θ) = (exp(3λ * θ) * Δ(3λ) \ (A₁(h3000, K10) + 3B(ϕ, h2010) + B(h0010, h3000) + 3B(h1010, h2000) + 3B₁(ϕ, h2000, K10) + 3C(ϕ, ϕ, h1010)
                                     + 3C(ϕ, h0010, h2000) + C₁(ϕ, ϕ, ϕ, K10) + D(ϕ, ϕ, ϕ, h0010)) - 3 * (1 + im * b110) * (Δ(3λ) \ (Δ′(3λ) - θ * Δ(3λ))) * h3000(θ))
  h4001(θ) = (exp(4λ * θ) * Δ(4λ) \ (A₁(h4000, K01) + 4B(ϕ, h3001) + B(h0001, h4000)
                                     + 4B(h1001, h3000) + 6B(h2000, h2001) + 4B₁(ϕ, h3000, K01)
                                     + 3B₁(h2000, h2000, K01) + 6C(ϕ, ϕ, h2001) + 4C(ϕ, h0001, h3000)
                                     + 12C(ϕ, h1001, h2000) + 3C(h0001, h2000, h2000)
                                     + 6C₁(ϕ, ϕ, h2000, K01) + 4D(ϕ, ϕ, ϕ, h1001) + 6D(ϕ, ϕ, h0001, h2000)
                                     + D₁(ϕ, ϕ, ϕ, ϕ, K01) + E(ϕ, ϕ, ϕ, ϕ, h0001))
              -
              4 * im * b101 * (Δ(4λ) \ (Δ′(4λ) - θ * Δ(4λ))) * h4000(θ))

  h3101(θ) = (exp(2λ * θ) * Δ(2λ) \ (A₁(h3100, K01) + 3B(ϕ, h2101) + B(conj ∘ ϕ, h3001) + B(h0001, h3100) + B(conj ∘ h1001, h3000) + 3B(h1001, h2100)
                                     + 3B(h1100, h2001) + 3B(h1101, h2000) + 3B₁(ϕ, h2100, K01) + B₁(conj ∘ ϕ, h3000, K01) + 3B₁(h1100, h2000, K01)
                                     + 3C(ϕ, ϕ, h1101) + 3C(ϕ, conj ∘ ϕ, h2001) + 3C(ϕ, h0001, h2100) + 3C(ϕ, conj ∘ h1001, h2000) + 6C(ϕ, h1001, h1100)
                                     + C(conj ∘ ϕ, h0001, h3000) + 3C(conj ∘ ϕ, h1001, h2000) + 3C(h0001, h1100, h2000) + 3C₁(ϕ, ϕ, h1100, K01)
                                     + 3C₁(ϕ, conj ∘ ϕ, h2000, K01) + D(ϕ, ϕ, ϕ, conj ∘ h1001) + 3D(ϕ, ϕ, conj ∘ ϕ, h1001) + 3D(ϕ, ϕ, h0001, h1100)
                                     + 3D(ϕ, conj ∘ ϕ, h0001, h2000) + D₁(ϕ, ϕ, ϕ, conj ∘ ϕ, K01) + E(ϕ, ϕ, ϕ, conj ∘ ϕ, h0001))
              -
              (6 * (1 + im * b201) * (Δ(2λ) \ (Δ′(2λ) - θ * Δ(2λ))) * h2000(θ)
               -
               6 * c₁ * exp(2λ * θ) * (Δ(2λ) \ ((Δ′(2λ) - θ * Δ(2λ)) * h2001(0) + im * b101 * (Δ′′(2λ) - θ^2 * Δ(2λ)) * h2000(0)))
               -
               2 * im * b101 * exp(2λ * θ) * (Δ(2λ) \ ((Δ′(2λ) - θ * Δ(2λ)) * h3100(0) + 3c₁ * (Δ′′(2λ) - θ^2 * Δ(2λ)) * h2000(0)))))

  h2201(θ) = (Δ(0) \ (A₁(h2200, K01) + 2B(ϕ, conj ∘ h2101) + 2B(conj ∘ ϕ, h2101) + B(h0001, h2200)
                      + 2B(conj ∘ h1001, h2100) + B(conj ∘ h2000, h2001) + B(conj ∘ h2001, h2000)
                      + 2B(h1001, conj ∘ h2100) + 4B(h1100, h1101) + 2B₁(ϕ, conj ∘ h2100, K01)
                      + B₁(conj ∘ h2000, h2000, K01) + 2B₁(h1100, h1100, K01) + 2B₁(conj ∘ ϕ, h2100, K01)
                      + C(ϕ, ϕ, conj ∘ h2001) + 4C(ϕ, conj ∘ ϕ, h1101) + 2C(ϕ, h0001, conj ∘ h2100)
                      + 4C(ϕ, conj ∘ h1001, h1100) + 2C(ϕ, conj ∘ h2000, h1001) + C(conj ∘ ϕ, conj ∘ ϕ, h2001)
                      + 2C(conj ∘ ϕ, h0001, h2100) + 2C(conj ∘ ϕ, conj ∘ h1001, h2000) + 4C(conj ∘ ϕ, h1001, h1100)
                      + C(h0001, conj ∘ h2000, h2000) + 2C(h0001, h1100, h1100) + C₁(ϕ, ϕ, conj ∘ h2000, K01)
                      + 4C₁(ϕ, conj ∘ ϕ, h1100, K01) + C₁(conj ∘ ϕ, conj ∘ ϕ, h2000, K01)
                      + 2D(ϕ, ϕ, conj ∘ ϕ, conj ∘ h1001) + D(ϕ, ϕ, h0001, conj ∘ h2000) + 2D(ϕ, conj ∘ ϕ, conj ∘ ϕ, h1001)
                      + 4D(ϕ, conj ∘ ϕ, h0001, h1100) + D(conj ∘ ϕ, conj ∘ ϕ, h0001, h2000) + D₁(ϕ, ϕ, conj ∘ ϕ, conj ∘ ϕ, K01)
                      + E(ϕ, ϕ, conj ∘ ϕ, conj ∘ ϕ, h0001))
              -
              8 * (Δ(0) \ (Δ′(0) - θ * Δ(0))) * h1100(θ))

  M3201 = (A₁(h3200, K01) + 3B(ϕ, h2201) + 2B(conj ∘ ϕ, h3101) + B(h0001, h3200)
           + 2B(conj ∘ h1001, h3100) + B(conj ∘ h2000, h3001) + B(conj ∘ h2001, h3000)
           + 3B(h1001, h2200) + 6B(h1100, h2101) + 6B(h1101, h2100)
           + 3B(conj ∘ h2100, h2001) + 3B(conj ∘ h2101, h2000) + 3B₁(ϕ, h2200, K01)
           + 2B₁(conj ∘ ϕ, h3100, K01) + B₁(conj ∘ h2000, h3000, K01) + 6B₁(h1100, h2100, K01)
           + 3B₁(conj ∘ h2100, h2000, K01) + 3C(ϕ, ϕ, conj ∘ h2101) + 6C(ϕ, conj ∘ ϕ, h2101)
           + 3C(ϕ, h0001, h2200) + 6C(ϕ, conj ∘ h1001, h2100) + 3C(ϕ, conj ∘ h2000, h2001)
           + 3C(ϕ, conj ∘ h2001, h2000) + 6C(ϕ, h1001, conj ∘ h2100) + 12C(ϕ, h1100, h1101)
           + C(conj ∘ ϕ, conj ∘ ϕ, h3001) + 2C(conj ∘ ϕ, h0001, h3100) + 2C(conj ∘ ϕ, conj ∘ h1001, h3000)
           + 6C(conj ∘ ϕ, h1001, h2100) + 6C(conj ∘ ϕ, h1100, h2001) + 6C(conj ∘ ϕ, h1101, h2000)
           + C(h0001, conj ∘ h2000, h3000) + 6C(h0001, h1100, h2100) + 3C(h0001, conj ∘ h2100, h2000)
           + 6C(conj ∘ h1001, h1100, h2000) + 3C(conj ∘ h2000, h1001, h2000)
           + 6C(h1001, h1100, h1100) + 3C₁(ϕ, ϕ, conj ∘ h2100, K01) + 6C₁(ϕ, conj ∘ ϕ, h2100, K01)
           + 3C₁(ϕ, conj ∘ h2000, h2000, K01) + 6C₁(ϕ, h1100, h1100, K01) + C₁(conj ∘ ϕ, conj ∘ ϕ, h3000, K01)
           + 6C₁(conj ∘ ϕ, h1100, h2000, K01) + D(ϕ, ϕ, ϕ, conj ∘ h2001) + 6D(ϕ, ϕ, conj ∘ ϕ, h1101)
           + 3D(ϕ, ϕ, h0001, conj ∘ h2100) + 6D(ϕ, ϕ, conj ∘ h1001, h1100) + 3D(ϕ, ϕ, conj ∘ h2000, h1001)
           + 3D(ϕ, conj ∘ ϕ, conj ∘ ϕ, h2001) + 6D(ϕ, conj ∘ ϕ, h0001, h2100) + 6D(ϕ, conj ∘ ϕ, conj ∘ h1001, h2000)
           + 12D(ϕ, conj ∘ ϕ, h1001, h1100) + 3D(ϕ, h0001, conj ∘ h2000, h2000)
           + 6D(ϕ, h0001, h1100, h1100) + D(conj ∘ ϕ, conj ∘ ϕ, h0001, h3000) + 3D(conj ∘ ϕ, conj ∘ ϕ, h1001, h2000)
           + 6D(conj ∘ ϕ, h0001, h1100, h2000) + D₁(ϕ, ϕ, ϕ, conj ∘ h2000, K01) + 6D₁(ϕ, ϕ, conj ∘ ϕ, h1100, K01)
           + 3D₁(ϕ, conj ∘ ϕ, conj ∘ ϕ, h2000, K01) + 2E(ϕ, ϕ, ϕ, conj ∘ ϕ, conj ∘ h1001) + E(ϕ, ϕ, ϕ, h0001, conj ∘ h2000)
           + 3E(ϕ, ϕ, conj ∘ ϕ, conj ∘ ϕ, h1001) + 6E(ϕ, ϕ, conj ∘ ϕ, h0001, h1100) + 3E(ϕ, conj ∘ ϕ, conj ∘ ϕ, h0001, h2000)
           + E₁(ϕ, ϕ, ϕ, conj ∘ ϕ, conj ∘ ϕ, K01) + K(ϕ, ϕ, ϕ, conj ∘ ϕ, conj ∘ ϕ, h0001))

  g3201 = (1 / 12) * first(transpose(p) * M3201)
  a3201 = real(g3201)

  h3201(θ) = (Aᴵᴺⱽ(λ, M3201, -12g3201)(θ) - 12c₂ * A2ᴵᴺⱽ(λ, h1001(0), -im * b101)(θ) - (18 + 6 * im * b201) * A2ᴵᴺⱽ(λ, h2100(0), -2c₁)(θ)
              -
              6 * im * imag(c₁) * A3ᴵᴺⱽ(λ, h2101(0), -(2(1 + im * b201) * q + im * b101 * h2100(0) + 2c₁ * h1001(0)), 4 * im * c₁ * b101)(θ)
              -
              im * b101 * A3ᴵᴺⱽ(λ, h3200(0), -(12c₂ * q + 6 * im * imag(c₁) * h2100(0)), 12 * im * imag(c₁) * c₁)(θ))


  h5001(θ) = (exp(5λ * θ) * Δ(5λ) \ (A₁(h5000, K01) + 5B(ϕ, h4001) + B(h0001, h5000)
                                     + 5B(h1001, h4000) + 10B(h2000, h3001) + 10B(h2001, h3000)
                                     + 5B₁(ϕ, h4000, K01) + 10B₁(h2000, h3000, K01) + 10C(ϕ, ϕ, h3001)
                                     + 5C(ϕ, h0001, h4000) + 20C(ϕ, h1001, h3000) + 30C(ϕ, h2000, h2001)
                                     + 10C(h0001, h2000, h3000) + 15C(h1001, h2000, h2000)
                                     + 10C₁(ϕ, ϕ, h3000, K01) + 15C₁(ϕ, h2000, h2000, K01)
                                     + 10D(ϕ, ϕ, ϕ, h2001) + 10D(ϕ, ϕ, h0001, h3000) + 30D(ϕ, ϕ, h1001, h2000)
                                     + 15D(ϕ, h0001, h2000, h2000) + 10D₁(ϕ, ϕ, ϕ, h2000, K01)
                                     + 5E(ϕ, ϕ, ϕ, ϕ, h1001) + 10E(ϕ, ϕ, ϕ, h0001, h2000)
                                     + E₁(ϕ, ϕ, ϕ, ϕ, ϕ, K01) + K(ϕ, ϕ, ϕ, ϕ, ϕ, h0001))
              -
              5 * im * b101 * (Δ(5λ) \ (Δ′(5λ) - θ * Δ(5λ))) * h5000(θ))

  h4101(θ) = (exp(3λ * θ) * Δ(3λ) \ (A₁(h4100, K01) + 4B(ϕ, h3101) + B(conj ∘ (ϕ), h4001)
                                     + B(h0001, h4100) + B(conj ∘ (h1001), h4000) + 4B(h1001, h3100)
                                     + 4B(h1100, h3001) + 4B(h1101, h3000) + 6B(h2000, h2101)
                                     + 6B(h2001, h2100) + 4B₁(ϕ, h3100, K01) + B₁(conj ∘ (ϕ), h4000, K01)
                                     + 4B₁(h1100, h3000, K01) + 6B₁(h2000, h2100, K01) + 6C(ϕ, ϕ, h2101)
                                     + 4C(ϕ, conj ∘ (ϕ), h3001) + 4C(ϕ, h0001, h3100) + 4C(ϕ, conj ∘ (h1001), h3000)
                                     + 12C(ϕ, h1001, h2100) + 12C(ϕ, h1100, h2001) + 12C(ϕ, h1101, h2000)
                                     + C(conj ∘ (ϕ), h0001, h4000) + 4C(conj ∘ (ϕ), h1001, h3000) + 6C(conj ∘ (ϕ), h2000, h2001)
                                     + 4C(h0001, h1100, h3000) + 6C(h0001, h2000, h2100)
                                     + 3C(conj ∘ (h1001), h2000, h2000) + 12C(h1001, h1100, h2000)
                                     + 6C₁(ϕ, ϕ, h2100, K01) + 4C₁(ϕ, conj ∘ (ϕ), h3000, K01)
                                     + 12C₁(ϕ, h1100, h2000, K01) + 3C₁(conj ∘ (ϕ), h2000, h2000, K01)
                                     + 4D(ϕ, ϕ, ϕ, h1101) + 6D(ϕ, ϕ, conj ∘ (ϕ), h2001) + 6D(ϕ, ϕ, h0001, h2100)
                                     + 6D(ϕ, ϕ, conj ∘ (h1001), h2000) + 12D(ϕ, ϕ, h1001, h1100)
                                     + 4D(ϕ, conj ∘ (ϕ), h0001, h3000) + 12D(ϕ, conj ∘ (ϕ), h1001, h2000)
                                     + 12D(ϕ, h0001, h1100, h2000) + 3D(conj ∘ (ϕ), h0001, h2000, h2000)
                                     + 4D₁(ϕ, ϕ, ϕ, h1100, K01) + 6D₁(ϕ, ϕ, conj ∘ (ϕ), h2000, K01)
                                     + E(ϕ, ϕ, ϕ, ϕ, conj ∘ (h1001)) + 4E(ϕ, ϕ, ϕ, conj ∘ (ϕ), h1001) + 4E(ϕ, ϕ, ϕ, h0001, h1100)
                                     + 6E(ϕ, ϕ, conj ∘ (ϕ), h0001, h2000) + E₁(ϕ, ϕ, ϕ, ϕ, conj ∘ (ϕ), K01) + K(ϕ, ϕ, ϕ, ϕ, conj ∘ (ϕ), h0001))
              -
              12 * (1 + im * b201) * (Δ(3λ) \ (Δ′(3λ) - θ * Δ(3λ))) * h3000(θ)
              -
              12 * c₁ * exp(3λ * θ) * (Δ(3λ) \ ((Δ′(3λ) - θ * Δ(3λ)) * h3001(0) + (3 / 2) * im * b101 * (Δ′′(3λ) - θ^2 * Δ(3λ)) * h3000(0)))
              -
              3 * im * b101 * exp(3λ * θ) * (Δ(3λ) \ ((Δ′(3λ) - θ * Δ(3λ)) * h4100(0) + 6c₁ * (Δ′′(3λ) - θ^2 * Δ(3λ)) * h3000(0))))



  #K02
  M0002 = (2A₁(h0001, K01) + B(h0001, h0001) + J₂(K01, K01))
  M1002 = (2A₁(h1001, K01) + 2B(h0001, h1001) + A₂(ϕ, K01, K01) + 2B₁(ϕ, h0001, K01) + C(ϕ, h0001, h0001))
  ℳ1002 = (B(ϕ, (θ -> Δ(0) \ (M0002))) + M1002)
  B1002(θ) = Aᴵᴺⱽ(λ, ℳ1002, 0)(θ) - 2 * im * b101 * A2ᴵᴺⱽ(λ, h1001(0), -im * b101)(θ)
  M2002 = (2A₁(h2001, K01) + 2B(h0001, h2001) + 2B(h1001, h1001) + A₂(h2000, K01, K01) + 4B₁(ϕ, h1001, K01)
           + 2B₁(h0001, h2000, K01) + 4C(ϕ, h0001, h1001) + C(h0001, h0001, h2000) + B₂(ϕ, ϕ, K01, K01) + 2C₁(ϕ, ϕ, h0001, K01)
           + D(ϕ, ϕ, h0001, h0001))
  M1102 = (2A₁(h1101, K01) + 2B(h0001, h1101) + 2B(conj ∘ h1001, h1001) + A₂(h1100, K01, K01) + 4real(B₁(conj ∘ ϕ, h1001, K01))
           + 2B₁(h0001, h1100, K01) + 4real(C(conj ∘ ϕ, h0001, h1001)) + C(h0001, h0001, h1100) + B₂(ϕ, conj ∘ ϕ, K01, K01)
           + 2C₁(ϕ, conj ∘ ϕ, h0001, K01) + D(ϕ, conj ∘ ϕ, h0001, h0001))

  ℳ2002(θ) = (M2002 + 2B(ϕ, B1002) + B(h2000, (θ -> Δ(0) \ M0002)) + C(ϕ, ϕ, (θ -> Δ(0) \ M0002))
              -
              4 * im * b101 * ((Δ′(2λ) - θ * Δ(2λ)) * h2001(0) + im * b101 * (Δ′′(2λ) - θ^2 * Δ(2λ)) * h2000(0)))
  ℳ1102(_) = (M1102 + 2real(B(conj ∘ ϕ, B1002)) + B(h1100, (θ -> Δ(0) \ M0002)) + C(ϕ, conj ∘ ϕ, (θ -> Δ(0) \ M0002)))

  M2102 = (2A₁(h2101, K01) + 2B(h0001, h2101) + 2B(conj ∘ h1001, h2001) + 4B(h1001, h1101) + A₂(h2100, K01, K01)
           + 4B₁(ϕ, h1101, K01) + 2B₁(conj ∘ ϕ, h2001, K01) + 2B₁(h0001, h2100, K01) + 2B₁(conj ∘ h1001, h2000, K01)
           + 4B₁(h1001, h1100, K01) + 4C(ϕ, h0001, h1101) + 4C(ϕ, conj ∘ h1001, h1001) + 2C(conj ∘ ϕ, h0001, h2001)
           + 2C(conj ∘ ϕ, h1001, h1001) + C(h0001, h0001, h2100) + 2C(h0001, conj ∘ h1001, h2000) + 4C(h0001, h1001, h1100)
           + 2B₂(ϕ, h1100, K01, K01) + B₂(conj ∘ ϕ, h2000, K01, K01) + 2C₁(ϕ, ϕ, conj ∘ h1001, K01) + 4C₁(ϕ, conj ∘ ϕ, h1001, K01)
           + 4C₁(ϕ, h0001, h1100, K01) + 2C₁(conj ∘ ϕ, h0001, h2000, K01) + 2D(ϕ, ϕ, h0001, conj ∘ h1001) + 4D(ϕ, conj ∘ ϕ, h0001, h1001)
           + 2D(ϕ, h0001, h0001, h1100) + D(conj ∘ ϕ, h2000, h0001, h0001) + C₂(ϕ, ϕ, conj ∘ ϕ, K01, K01) + 2D₁(ϕ, ϕ, conj ∘ ϕ, h0001, K01)
           + E(ϕ, ϕ, conj ∘ ϕ, h0001, h0001))


  Q102 = -real(first(transpose(p) * ℳ1002))
  Q202 = -0.5 * (real(first(transpose(p) * (2B(ϕ, (θ -> Δ(0) \ ℳ1102(θ))) + B(conj ∘ ϕ, (θ -> exp(2 * λ * θ) * (Δ(2λ) \ ℳ2002(θ))))
                                            + B(h2100, (θ -> Δ(0) \ M0002)) + B(h2000, conj ∘ B1002) + 2B(h1100, B1002) + C(ϕ, ϕ, conj ∘ B1002)
                                            + 2C(ϕ, conj ∘ ϕ, B1002) + 2C(ϕ, h1100, (θ -> Δ(0) \ M0002)) + C(conj ∘ ϕ, h2000, (θ -> Δ(0) \ M0002)) + D(ϕ, ϕ, conj ∘ ϕ, (θ -> Δ(0) \ M0002))
                                            + M2102)))
                 +
                 imag(first(transpose(p) * ℳ1002)) * rP2)

  K02 = P \ [Q102; Q202]

  h0002(_) = Δ(0) \ (J₁ * K02 + M0002)
  b102 = imag(first(transpose(p) * (A₁(ϕ, K02) + B(ϕ, h0002) + M1002)))


  h1002(θ) = Aᴵᴺⱽ(λ, A₁(ϕ, K02) + B(ϕ, h0002) + M1002, -im * b102)(θ) - 2 * im * b101 * A2ᴵᴺⱽ(λ, h1001(0), -im * b101)(θ)

  h2002(θ) = (exp(2λ * θ) * Δ(2λ) \ (A₁(h2000, K02) + 2B(ϕ, h1002) + B(h0002, h2000) + B₁(ϕ, ϕ, K02) + C(ϕ, ϕ, h0002) + M2002) - 2 * im * b102 * (Δ(2λ) \ (Δ′(2λ) - θ * Δ(2λ))) * h2001(θ)
              -
              4 * im * b101 * exp(2λ * θ) * (Δ(2λ) \ ((Δ′(2λ) - θ * Δ(2λ)) * h2001(0) + im * b101 * (Δ′′(2λ) - θ^2 * Δ(2λ)) * h2000(0))))
  h1102(θ) = Δ(0) \ (A₁(h1100, K02) + 2real(B(conj ∘ ϕ, h1002)) + B(h0002, h1100) + B₁(ϕ, conj ∘ ϕ, K02) + C(ϕ, conj ∘ ϕ, h0002) + M1102)


  R2102 = (A₁(h2100, K02) + 2B(ϕ, h1102) + B(conj ∘ ϕ, h2002) + B(h0002, h2100) + B(conj ∘ h1002, h2000) + 2B(h1002, h1100)
           + 2B₁(ϕ, h1100, K02) + B₁(conj ∘ ϕ, h2000, K02) + C(ϕ, ϕ, conj ∘ h1002) + 2C(ϕ, conj ∘ ϕ, h1002) + 2C(ϕ, h0002, h1100)
           + C(conj ∘ ϕ, h0002, h2000) + C₁(ϕ, ϕ, conj ∘ ϕ, K02) + D(ϕ, ϕ, conj ∘ ϕ, h0002) + M2102)

  b202 = 0.5 * imag(first(transpose(p) * (R2102)))

  h2102(θ) = (Aᴵᴺⱽ(λ, R2102, -2 * im * b202)(θ) - im * b102 * A2ᴵᴺⱽ(λ, h2100(0), -2c₁)(θ)
              -
              2c₁ * A3ᴵᴺⱽ(λ, h1002(0), -(im * b102 * q + 2 * im * b101 * h1001(0)), -2b101^2)(θ) - 4(1 + im * b201) * A2ᴵᴺⱽ(λ, h1001(0), -im * b101)(θ)
              -
              2 * im * b101 * A3ᴵᴺⱽ(λ, h2101(0), -(2 * (1 + im * b201) * q + im * b101 * h2100(0) + 2c₁ * h1001(0)), 4 * im * c₁ * b101)(θ))

  h3002(θ) = (exp(3λ * θ) * Δ(3λ) \ (2A₁(h3001, K01) + A₁(h3000, K02) + 3B(ϕ, h2002)
                                     + 2B(h0001, h3001) + B(h0002, h3000) + 6B(h1001, h2001)
                                     + 3B(h1002, h2000) + A₂(h3000, K01, K01) + 6B₁(ϕ, h2001, K01)
                                     + 3B₁(ϕ, h2000, K02) + 2B₁(h0001, h3000, K01) + 6B₁(h1001, h2000, K01)
                                     + 3C(ϕ, ϕ, h1002) + 6C(ϕ, h0001, h2001) + 3C(ϕ, h0002, h2000)
                                     + 6C(ϕ, h1001, h1001) + C(h0001, h0001, h3000)
                                     + 6C(h0001, h1001, h2000) + 3B₂(ϕ, h2000, K01, K01) + C₁(ϕ, ϕ, ϕ, K02)
                                     + 6C₁(ϕ, ϕ, h1001, K01) + 6C₁(ϕ, h0001, h2000, K01) + D(ϕ, ϕ, ϕ, h0002)
                                     + 6D(ϕ, ϕ, h0001, h1001) + 3D(ϕ, h0001, h0001, h2000)
                                     + C₂(ϕ, ϕ, ϕ, K01, K01) + 2D₁(ϕ, ϕ, ϕ, h0001, K01)
                                     + E(ϕ, ϕ, ϕ, h0001, h0001))
              -
              3 * im * b102 * (Δ(3λ) \ (Δ′(3λ) - θ * Δ(3λ))) * h3000(θ)
              -
              6 * im * b101 * exp(3λ * θ) * (Δ(3λ) \ ((Δ′(3λ) - θ * Δ(3λ)) * h3001(0) + (3 / 2) * im * b101 * (Δ′′(3λ) - θ^2 * Δ(3λ)) * h3000(0))))


  #K11
  M0011 = (A₁(h0010, K01) + A₁(h0001, K10) + B(h0001, h0010) + J₂(K01, K10))
  M1011 = (A₁(h1010, K01) + A₁(h1001, K10) + B(h0001, h1010)
           + B(h0010, h1001) + A₂(ϕ, K01, K10) + B₁(ϕ, h0010, K01)
           + B₁(ϕ, h0001, K10) + C(ϕ, h0001, h0010))
  ℳ1011 = (B(ϕ, (θ -> Δ(0) \ (M0011))) + M1011)
  B1011(θ) = Aᴵᴺⱽ(λ, ℳ1011, 0)(θ) - (1 + im * b110) * A2ᴵᴺⱽ(λ, h1001(0), -im * b101)(θ) - im * b101 * A2ᴵᴺⱽ(λ, h1010(0), -(1 + im * b110))(θ)
  M2011 = (A₁(h2010, K01) + A₁(h2001, K10) + B(h0001, h2010)
           + B(h0010, h2001) + 2B(h1001, h1010) + A₂(h2000, K01, K10)
           + 2B₁(ϕ, h1010, K01) + 2B₁(ϕ, h1001, K10) + B₁(h0010, h2000, K01)
           + B₁(h0001, h2000, K10) + 2C(ϕ, h0001, h1010) + 2C(ϕ, h0010, h1001)
           + C(h0001, h0010, h2000) + B₂(ϕ, ϕ, K10, K01) + C₁(ϕ, ϕ, h0010, K01)
           + C₁(ϕ, ϕ, h0001, K10) + D(ϕ, ϕ, h0001, h0010))
  M1111 = (A₁(h1110, K01) + A₁(h1101, K10) + B(h0001, h1110)
           + B(h0010, h1101) + 2 * real(B(conj ∘ (h1001), h1010)) + A₂(h1100, K01, K10)
           + 2 * real(B₁(conj ∘ (ϕ), h1010, K01)) + 2 * real(B₁(conj ∘ (ϕ), h1001, K10)) + B₁(h0010, h1100, K01)
           + B₁(h0001, h1100, K10) + 2 * real(C(conj ∘ (ϕ), h0001, h1010)) + 2 * real(C(conj ∘ (ϕ), h0010, h1001))
           + C(h0001, h0010, h1100) + B₂(ϕ, conj ∘ (ϕ), K01, K10) + C₁(ϕ, conj ∘ (ϕ), h0010, K01)
           + C₁(ϕ, conj ∘ (ϕ), h0001, K10) + D(ϕ, conj ∘ (ϕ), h0001, h0010))

  ℳ2011(θ) = (M2011 + 2B(ϕ, B1011) + B(h2000, (θ -> Δ(0) \ M0011)) + C(ϕ, ϕ, (θ -> Δ(0) \ M0011))
              -
              2(1 + im * b110) * ((Δ′(2λ) - θ * Δ(2λ)) * h2001(0) + im * b101 * (Δ′′(2λ) - θ^2 * Δ(2λ)) * h2000(0))
              -
              2 * im * b101 * ((Δ′(2λ) - θ * Δ(2λ)) * h2010(0) + (1 + im * b110) * (Δ′′(2λ) - θ^2 * Δ(2λ)) * h2000(0)))
  ℳ1111(θ) = (M1111 + 2real(B(conj ∘ ϕ, B1011)) + B(h1100, (θ -> Δ(0) \ M0011)) + C(ϕ, conj ∘ ϕ, (θ -> Δ(0) \ M0011))
              -
              (Δ′(0) - θ * Δ(0)) * h1101(θ))

  M2111 = (A₁(h2110, K01) + A₁(h2101, K10) + B(h0001, h2110) + B(h0010, h2101)
           + B(conj ∘ (h1001), h2010) + B(conj ∘ (h1010), h2001) + 2B(h1001, h1110) + 2B(h1010, h1101)
           + A₂(h2100, K01, K10) + 2B₁(ϕ, h1110, K01) + 2B₁(ϕ, h1101, K10) + B₁(conj ∘ (ϕ), h2010, K01)
           + B₁(conj ∘ (ϕ), h2001, K10) + B₁(h0010, h2100, K01) + B₁(conj ∘ (h1010), h2000, K01) + 2B₁(h1010, h1100, K01)
           + B₁(h0001, h2100, K10) + B₁(conj ∘ (h1001), h2000, K10) + 2B₁(h1001, h1100, K10) + 2C(ϕ, h0001, h1110)
           + 2C(ϕ, h0010, h1101) + 2C(ϕ, conj ∘ (h1001), h1010) + 2C(ϕ, conj ∘ (h1010), h1001) + C(conj ∘ (ϕ), h0001, h2010)
           + C(conj ∘ (ϕ), h0010, h2001) + 2C(conj ∘ (ϕ), h1001, h1010) + C(h0001, h0010, h2100) + C(h0001, conj ∘ (h1010), h2000)
           + 2C(h0001, h1010, h1100) + C(h0010, conj ∘ (h1001), h2000) + 2C(h0010, h1001, h1100)
           + 2B₂(ϕ, h1100, K01, K10) + B₂(conj ∘ (ϕ), h2000, K01, K10) + C₁(ϕ, ϕ, conj ∘ (h1010), K01) + C₁(ϕ, ϕ, conj ∘ (h1001), K10)
           + 2C₁(ϕ, conj ∘ (ϕ), h1010, K01) + 2C₁(ϕ, conj ∘ (ϕ), h1001, K10) + 2C₁(ϕ, h0010, h1100, K01)
           + 2C₁(ϕ, h0001, h1100, K10) + C₁(conj ∘ (ϕ), h0010, h2000, K01) + C₁(conj ∘ (ϕ), h0001, h2000, K10)
           + D(ϕ, ϕ, h0001, conj ∘ (h1010)) + D(ϕ, ϕ, h0010, conj ∘ (h1001)) + 2D(ϕ, conj ∘ (ϕ), h0001, h1010)
           + 2D(ϕ, conj ∘ (ϕ), h0010, h1001) + 2D(ϕ, h0001, h0010, h1100) + D(conj ∘ (ϕ), h0001, h0010, h2000)
           + C₂(ϕ, ϕ, conj ∘ (ϕ), K01, K10) + D₁(ϕ, ϕ, conj ∘ (ϕ), h0010, K01) + D₁(ϕ, ϕ, conj ∘ (ϕ), h0001, K10) + E(ϕ, ϕ, conj ∘ (ϕ), h0001, h0010))


  Q111 = -real(first(transpose(p) * ℳ1011))
  Q211 = -0.5 * (real(first(transpose(p) * (2B(ϕ, (θ -> Δ(0) \ ℳ1111(θ))) + B(conj ∘ ϕ, (θ -> exp(2 * λ * θ) * (Δ(2λ) \ ℳ2011(θ))))
                                            + B(h2100, (θ -> Δ(0) \ M0011)) + B(h2000, conj ∘ B1011) + 2B(h1100, B1011) + C(ϕ, ϕ, conj ∘ B1011)
                                            + 2C(ϕ, conj ∘ ϕ, B1011) + 2C(ϕ, h1100, (θ -> Δ(0) \ M0011)) + C(conj ∘ ϕ, h2000, (θ -> Δ(0) \ M0011)) + D(ϕ, ϕ, conj ∘ ϕ, (θ -> Δ(0) \ M0011))
                                            + M2111)))
                 +
                 imag(first(transpose(p) * ℳ1011)) * rP2)

  K11 = P \ [Q111; Q211]

  h0011(θ) = Δ(0) \ (J₁ * K11 + M0011)
  b111 = imag(first(transpose(p) * (A₁(ϕ, K11) + B(ϕ, h0011) + M1011)))
  h1011(θ) = Aᴵᴺⱽ(λ, A₁(ϕ, K11) + B(ϕ, h0011) + M1011, -im * b111)(θ) - (1 + im * b110) * A2ᴵᴺⱽ(λ, h1001(0), -im * b101)(θ) - im * b101 * A2ᴵᴺⱽ(λ, h1010(0), -(1 + im * b110))(θ)

  #K03
  M0003 = (3A₁(h0002, K01) + 3A₁(h0001, K02) + 3B(h0001, h0002)
           + 3J₂(K01, K02) + 3A₂(h0001, K01, K01) + 3B₁(h0001, h0001, K01)
           + J₃(K01, K01, K01) + C(h0001, h0001, h0001))
  M1003 = (3A₁(h1002, K01) + 3A₁(h1001, K02) + 3B(h0001, h1002)
           + 3B(h0002, h1001) + 3A₂(ϕ, K01, K02) + 3A₂(h1001, K01, K01)
           + 3B₁(ϕ, h0002, K01) + 3B₁(ϕ, h0001, K02) + 6B₁(h0001, h1001, K01)
           + 3C(ϕ, h0001, h0002) + 3C(h0001, h0001, h1001) + A₃(ϕ, K01, K01, K01)
           + 3B₂(ϕ, h0001, K01, K01) + 3C₁(ϕ, h0001, h0001, K01)
           + D(ϕ, h0001, h0001, h0001))
  ℳ1003 = (B(ϕ, (θ -> Δ(0) \ (M0003))) + M1003)
  B1003(θ) = Aᴵᴺⱽ(λ, ℳ1003, 0)(θ) - 3 * im * b102 * A2ᴵᴺⱽ(λ, h1001(0), -im * b101)(θ) - 3 * im * b101 * A3ᴵᴺⱽ(λ, h1002(0), -(im * b102 * q + 2 * im * b101 * h1001(0)), -2b101^2)(θ)

  M2003 = (3 * A₁(h2002, K01) + 3 * A₁(h2001, K02) + 3 * B(h0001, h2002) + 3 * B(h0002, h2001)
           + 6 * B(h1001, h1002) + 3 * A₂(h2001, K01, K01) + 3 * A₂(h2000, K01, K02)
           + 6 * B₁(ϕ, h1002, K01) + 6 * B₁(ϕ, h1001, K02) + 6 * B₁(h0001, h2001, K01)
           + 3 * B₁(h0002, h2000, K01) + 6 * B₁(h1001, h1001, K01) + 3 * B₁(h0001, h2000, K02)
           + 6 * C(ϕ, h0001, h1002) + 6 * C(ϕ, h0002, h1001) + 3 * C(h0001, h0001, h2001)
           + 3 * C(h0001, h0002, h2000) + 6 * C(h0001, h1001, h1001) + A₃(h2000, K01, K01, K01)
           + 3 * B₂(ϕ, ϕ, K01, K02) + 6 * B₂(ϕ, h1001, K01, K01) + 3 * B₂(h0001, h2000, K01, K01)
           + 3 * C₁(ϕ, ϕ, h0002, K01) + 3 * C₁(ϕ, ϕ, h0001, K02) + 12 * C₁(ϕ, h0001, h1001, K01)
           + 3 * C₁(h0001, h0001, h2000, K01) + 3 * D(ϕ, ϕ, h0001, h0002) + 6 * D(ϕ, h0001, h0001, h1001)
           + D(h0001, h0001, h0001, h2000) + B₃(ϕ, ϕ, K01, K01, K01) + 3 * C₂(ϕ, ϕ, h0001, K01, K01)
           + 3 * D₁(ϕ, ϕ, h0001, h0001, K01) + E(ϕ, ϕ, h0001, h0001, h0001))

  M1103 = (3 * A₁(h1102, K01) + 3 * A₁(h1101, K02) + 3 * B(h0001, h1102) + 3 * B(h0002, h1101)
           + 3 * B(conj ∘ (h1001), h1002) + 3 * B(conj ∘ (h1002), h1001) + 3 * A₂(h1101, K01, K01)
           + 3 * A₂(h1100, K01, K02) + 6 * real(B₁(conj ∘ (ϕ), h1002, K01)) + 6 * real(B₁(conj ∘ (ϕ), h1001, K02))
           + 6 * B₁(h0001, h1101, K01) + 3 * B₁(h0002, h1100, K01) + 6 * B₁(conj ∘ (h1001), h1001, K01)
           + 3 * B₁(h0001, h1100, K02) + 6 * real(C(conj ∘ (ϕ), h0001, h1002)) + 6 * real(C(conj ∘ (ϕ), h0002, h1001))
           + 3 * C(h0001, h0001, h1101) + 3 * C(h0001, h0002, h1100) + 6 * C(h0001, conj ∘ (h1001), h1001)
           + A₃(h1100, K01, K01, K01) + 3 * B₂(ϕ, conj ∘ (ϕ), K01, K02) + 6 * real(B₂(conj ∘ (ϕ), h1001, K01, K01))
           + 3 * B₂(h0001, h1100, K01, K01) + 3 * C₁(ϕ, conj ∘ (ϕ), h0002, K01) + 3 * C₁(ϕ, conj ∘ (ϕ), h0001, K02)
           + 12 * real(C₁(conj ∘ (ϕ), h0001, h1001, K01)) + 3 * C₁(h0001, h0001, h1100, K01)
           + 3 * D(ϕ, conj ∘ (ϕ), h0001, h0002) + 6 * real(D(conj ∘ (ϕ), h0001, h0001, h1001)) + D(h0001, h0001, h0001, h1100)
           + B₃(ϕ, conj ∘ (ϕ), K01, K01, K01) + 3 * C₂(ϕ, conj ∘ (ϕ), h0001, K01, K01) + 3 * D₁(ϕ, conj ∘ (ϕ), h0001, h0001, K01)
           + E(ϕ, conj ∘ (ϕ), h0001, h0001, h0001))

  ℳ2003(θ) = (M2003 + 2B(ϕ, B1003) + B(h2000, (θ -> Δ(0) \ M0003)) + C(ϕ, ϕ, (θ -> Δ(0) \ M0003))
              -
              6 * im * b102 * ((Δ′(2λ) - θ * Δ(2λ)) * h2001(0) + im * b101 * (Δ′′(2λ) - θ^2 * Δ(2λ)) * h2000(0))
              -
              6 * im * b101 * ((Δ′(2λ) - θ * Δ(2λ)) * h2002(0) + (Δ′′(2λ) - θ^2 * Δ(2λ)) * (im * b102 * h2000(0) + 2 * im * b101 * h2001(0)) - (4 / 3) * b101^2 * (Δ′′′(2λ) - θ^3 * Δ(2λ)) * h2000(0)))
  ℳ1103(θ) = (M1103 + 2real(B(conj ∘ ϕ, B1003)) + B(h1100, (θ -> Δ(0) \ M0003)) + C(ϕ, conj ∘ ϕ, (θ -> Δ(0) \ M0003)))

  M2103 = (3A₁(h2102, K01) + 3A₁(h2101, K02) + 3B(h0001, h2102) + 3B(h0002, h2101) + 3B(conj ∘ (h1001), h2002)
           + 3B(conj ∘ (h1002), h2001) + 6B(h1001, h1102) + 6B(h1002, h1101) + 3A₂(h2101, K01, K01)
           + 3A₂(h2100, K01, K02) + 6 * B₁(ϕ, h1102, K01) + 6 * B₁(ϕ, h1101, K02) + 3 * B₁(conj ∘ (ϕ), h2002, K01)
           + 3 * B₁(conj ∘ (ϕ), h2001, K02) + 6 * B₁(h0001, h2101, K01) + 3 * B₁(h0002, h2100, K01) + 6 * B₁(conj ∘ (h1001), h2001, K01)
           + 3 * B₁(conj ∘ (h1002), h2000, K01) + 12 * B₁(h1001, h1101, K01) + 6 * B₁(h1002, h1100, K01) + 3 * B₁(h0001, h2100, K02)
           + 3 * B₁(conj ∘ (h1001), h2000, K02) + 6 * B₁(h1001, h1100, K02) + 6 * C(ϕ, h0001, h1102) + 6 * C(ϕ, h0002, h1101)
           + 6 * C(ϕ, conj ∘ (h1001), h1002) + 6 * C(ϕ, conj ∘ (h1002), h1001) + 3 * C(conj ∘ (ϕ), h0001, h2002) + 3 * C(conj ∘ (ϕ), h0002, h2001)
           + 6 * C(conj ∘ (ϕ), h1001, h1002) + 3 * C(h0001, h0001, h2101) + 3 * C(h0001, h0002, h2100) + 6 * C(h0001, conj ∘ (h1001), h2001)
           + 3 * C(h0001, conj ∘ (h1002), h2000) + 12 * C(h0001, h1001, h1101) + 6 * C(h0001, h1002, h1100) + 3 * C(h0002, conj ∘ (h1001), h2000)
           + 6 * C(h0002, h1001, h1100) + 6 * C(conj ∘ (h1001), h1001, h1001) + A₃(h2100, K01, K01, K01) + 6 * B₂(ϕ, h1101, K01, K01)
           + 6 * B₂(ϕ, h1100, K01, K02) + 3 * B₂(conj ∘ (ϕ), h2001, K01, K01) + 3 * B₂(conj ∘ (ϕ), h2000, K01, K02)
           + 3 * B₂(h0001, h2100, K01, K01) + 3 * B₂(conj ∘ (h1001), h2000, K01, K01) + 6 * B₂(h1001, h1100, K01, K01)
           + 3 * C₁(ϕ, ϕ, conj ∘ (h1002), K01) + 3 * C₁(ϕ, ϕ, conj ∘ (h1001), K02) + 6 * C₁(ϕ, conj ∘ (ϕ), h1002, K01)
           + 6 * C₁(ϕ, conj ∘ (ϕ), h1001, K02) + 12 * C₁(ϕ, h0001, h1101, K01) + 6 * C₁(ϕ, h0002, h1100, K01)
           + 12 * C₁(ϕ, conj ∘ (h1001), h1001, K01) + 6 * C₁(ϕ, h0001, h1100, K02) + 6 * C₁(conj ∘ (ϕ), h0001, h2001, K01)
           + 3 * C₁(conj ∘ (ϕ), h0002, h2000, K01) + 6 * C₁(conj ∘ (ϕ), h1001, h1001, K01) + 3 * C₁(conj ∘ (ϕ), h0001, h2000, K02)
           + 3 * C₁(h0001, h0001, h2100, K01) + 6 * C₁(h0001, conj ∘ (h1001), h2000, K01) + 12 * C₁(h0001, h1001, h1100, K01)
           + 3 * D(ϕ, ϕ, h0001, conj ∘ (h1002)) + 3 * D(ϕ, ϕ, h0002, conj ∘ (h1001)) + 6 * D(ϕ, conj ∘ (ϕ), h0001, h1002)
           + 6 * D(ϕ, conj ∘ (ϕ), h0002, h1001) + 6 * D(ϕ, h0001, h0001, h1101) + 6 * D(ϕ, h0001, h0002, h1100)
           + 12 * D(ϕ, h0001, conj ∘ (h1001), h1001) + 3 * D(conj ∘ (ϕ), h0001, h0001, h2001) + 3 * D(conj ∘ (ϕ), h0001, h0002, h2000)
           + 6 * D(conj ∘ (ϕ), h0001, h1001, h1001) + D(h0001, h0001, h0001, h2100) + 3 * D(h0001, h0001, conj ∘ (h1001), h2000)
           + 6 * D(h0001, h0001, h1001, h1100) + 2 * B₃(ϕ, h1100, K01, K01, K01) + B₃(conj ∘ (ϕ), h2000, K01, K01, K01)
           + 3 * C₂(ϕ, ϕ, conj ∘ (ϕ), K01, K02) + 3 * C₂(ϕ, ϕ, conj ∘ (h1001), K01, K01) + 6 * C₂(ϕ, conj ∘ (ϕ), h1001, K01, K01)
           + 6 * C₂(ϕ, h0001, h1100, K01, K01) + 3 * C₂(conj ∘ (ϕ), h0001, h2000, K01, K01) + 3 * D₁(ϕ, ϕ, conj ∘ (ϕ), h0002, K01)
           + 3 * D₁(ϕ, ϕ, conj ∘ (ϕ), h0001, K02) + 6 * D₁(ϕ, ϕ, h0001, conj ∘ (h1001), K01) + 12 * D₁(ϕ, conj ∘ (ϕ), h0001, h1001, K01)
           + 6 * D₁(ϕ, h0001, h0001, h1100, K01) + 3 * D₁(conj ∘ (ϕ), h0001, h0001, h2000, K01) + 3 * E(ϕ, ϕ, conj ∘ (ϕ), h0001, h0002)
           + 3 * E(ϕ, ϕ, h0001, h0001, conj ∘ (h1001)) + 6 * E(ϕ, conj ∘ (ϕ), h0001, h0001, h1001) + 2 * E(ϕ, h0001, h0001, h0001, h1100)
           + E(conj ∘ (ϕ), h0001, h0001, h0001, h2000) + C₃(ϕ, ϕ, conj ∘ (ϕ), K01, K01, K01)
           + 3 * E₁(ϕ, ϕ, conj ∘ (ϕ), h0001, h0001, K01) + K(ϕ, ϕ, conj ∘ (ϕ), h0001, h0001, h0001)
  )


  Q103 = -real(first(transpose(p) * ℳ1003))
  Q203 = -0.5 * (real(first(transpose(p) * (2B(ϕ, (θ -> Δ(0) \ ℳ1103(θ))) + B(conj ∘ ϕ, (θ -> exp(2 * λ * θ) * (Δ(2λ) \ ℳ2003(θ))))
                                            + B(h2100, (θ -> Δ(0) \ M0003)) + B(h2000, conj ∘ B1003) + 2B(h1100, B1003) + C(ϕ, ϕ, conj ∘ B1003)
                                            + 2C(ϕ, conj ∘ ϕ, B1003) + 2C(ϕ, h1100, (θ -> Δ(0) \ M0003)) + C(conj ∘ ϕ, h2000, (θ -> Δ(0) \ M0003)) + D(ϕ, ϕ, conj ∘ ϕ, (θ -> Δ(0) \ M0003))
                                            + M2103)))
                 +
                 imag(first(transpose(p) * ℳ1003)) * rP2)

  K03 = P \ [Q103; Q203]

  h0003(θ) = Δ(0) \ (J₁ * K03 + M0003)
  b103 = imag(first(transpose(p) * (A₁(ϕ, K03) + B(ϕ, h0003) + M1003)))
  h1003(θ) = Aᴵᴺⱽ(λ, A₁(ϕ, K03) + B(ϕ, h0003) + M1003, -im * b103)(θ) - 3 * im * b102 * A2ᴵᴺⱽ(λ, h1001(0), -im * b101)(θ) - 3 * im * b101 * A3ᴵᴺⱽ(λ, h1002(0), -(im * b102 * q + 2 * im * b101 * h1001(0)), -2b101^2)(θ)

  nmfm = GenHopfNormalformHigherOrder(
    K10,
    K01,
    K02,
    K11,
    K03,
    c₁,
    c₂,
    c₃,
    ℓ₁,
    ℓ₂,
    ℓ₃,
    a3201,
    q,
    h2000,
    h1100,
    h3000,
    h2100,
    h4000,
    h3100,
    h2200,
    h5000,
    h4100,
    h3200,
    h6000,
    h5100,
    h4200,
    h3300,
    h4300,
    h7000,
    h6100,
    h5200,
    h0010,
    h0001,
    h0002,
    h0011,
    h0003,
    h1002,
    h1011,
    h1003,
    h1001,
    h1010,
    h2001,
    h2010,
    h2002,
    h1101,
    h1110,
    h1102,
    h3001,
    h3010,
    h2101,
    h2110,
    h4001,
    h2201,
    h3101,
    h4101,
    h5001,
    h3201,
    h2102,
    h3002,
    b101,
    b110,
    b201,
    b102)

  hopf = @set hopf.nmfm = nmfm
end

# function to create initial guess for LPC branch emanating from Generalized Hopf
function generalizedHopfToPsol(jet, genh, ϵ, ntst, ncol, τs)
  # extract normal form coefficients
  q = genh.nmfm.q
  h0010 = genh.nmfm.h0010(0)
  h0001 = genh.nmfm.h0001(0)
  h1100 = real(genh.nmfm.h1100(0))
  h2000 = genh.nmfm.h2000(0)
  γ110 = genh.nmfm.γ110
  γ101 = genh.nmfm.γ101
  γ210 = genh.nmfm.γ210
  γ201 = genh.nmfm.γ201
  c₁ = genh.nmfm.c₁
  c₂ = genh.nmfm.c₂
  ω₂ = genh.nmfm.ω₂

  α = real([γ110 γ101; γ210 γ201]) \ [0.0; -2real(c₂)]

  # construct initial guess
  t = range(0.0, 1.0, ntst * ncol + 1)
  profile = [2 * real(exp(2 * pi * t * im) * q) * ϵ + (real(h1100) + h0010 * α[1] + h0001 * α[2] + real(exp(4 * pi * t * im) * h2000)) * ϵ^2 for t ∈ t]
  pars = genh.parameters + α * ϵ^2
  T = 2pi / (abs(genh.ω) + (imag(c₁) - 2real(c₂) * ω₂) * ϵ^2)
  psol(profile, pars, collect(t), T, ncol, nothing, nothing)
end

# TODO: implement this
# function to create initial guess for LPC branch emanating from Generalized Hopf using higher order predictor
function generalizedHopfToPsolHigherOrder(jet, genh, ϵ, ntst, ncol, τs)
  # extract normal form coefficients
  q = genh.nmfm.q
  h0010 = genh.nmfm.h0010(0)
  h0001 = genh.nmfm.h0001(0)
  h0002 = genh.nmfm.h0002(0)
  h0011 = genh.nmfm.h0011(0)
  h0003 = genh.nmfm.h0003(0)
  h1100 = real(genh.nmfm.h1100(0))
  h2000 = genh.nmfm.h2000(0)
  h1010 = genh.nmfm.h1010(0)
  h1001 = genh.nmfm.h1001(0)
  h1002 = genh.nmfm.h1002(0)
  h1011 = genh.nmfm.h1011(0)
  h1003 = genh.nmfm.h1003(0)
  h2001 = genh.nmfm.h2001(0)
  h2010 = genh.nmfm.h2010(0)
  h2002 = genh.nmfm.h2002(0)
  h1101 = genh.nmfm.h1101(0)
  h1110 = genh.nmfm.h1110(0)
  h1102 = genh.nmfm.h1102(0)
  h3001 = genh.nmfm.h3001(0)
  h2101 = genh.nmfm.h2101(0)
  h4001 = genh.nmfm.h4001(0)
  h3101 = genh.nmfm.h3101(0)
  h2201 = genh.nmfm.h2201(0)
  h5001 = genh.nmfm.h5001(0)
  h4101 = genh.nmfm.h4101(0)
  h3201 = genh.nmfm.h3201(0)
  h3000 = genh.nmfm.h3000(0)
  h2100 = genh.nmfm.h2100(0)
  h4000 = genh.nmfm.h4000(0)
  h3100 = genh.nmfm.h3100(0)
  h2200 = genh.nmfm.h2200(0)
  h5000 = genh.nmfm.h5000(0)
  h4100 = genh.nmfm.h4100(0)
  h3200 = genh.nmfm.h3200(0)
  h6000 = genh.nmfm.h6000(0)
  h5100 = genh.nmfm.h5100(0)
  h4200 = genh.nmfm.h4200(0)
  h3300 = genh.nmfm.h3300(0)
  h7000 = genh.nmfm.h7000(0)
  h6100 = genh.nmfm.h6100(0)
  h5200 = genh.nmfm.h5200(0)
  h4300 = genh.nmfm.h4300(0)
  K10 = genh.nmfm.K10
  K01 = genh.nmfm.K01
  K02 = genh.nmfm.K02
  K11 = genh.nmfm.K11
  K03 = genh.nmfm.K03
  c₁ = genh.nmfm.c₁
  c₂ = genh.nmfm.c₂
  c₃ = genh.nmfm.c₃
  b101 = genh.nmfm.b101
  b110 = genh.nmfm.b110
  b102 = genh.nmfm.b102
  b201 = genh.nmfm.b201
  a3201 = genh.nmfm.a3201

  #approximmations of the parameters
  β₁ = real(c₂) * ϵ^4 + 2(real(c₃) - a3201 * real(c₂)) * ϵ^6
  β₂ = -2real(c₂) * ϵ^2 + (4a3201 * real(c₂) - 3 * real(c₃)) * ϵ^4
  β₂₂ = 4 * (real(c₂))^2 * ϵ^4 - 4real(c₂) * (4 * a3201 * real(c₂) - 3 * real(c₃)) * ϵ^6  #\beta_2^2
  β₁₂ = -2 * (real(c₂))^2 * ϵ^6       #\beta_1 * \beta_2
  β₂₂₂ = -8 * (real(c₂))^3 * ϵ^6 #\beta_2^3

  # construct initial guess
  t = range(0.0, 1.0, ntst * ncol + 1)
  profile = [(2 * real(exp(2 * pi * t * im) * q) * ϵ + h0010 * β₁ + h0001 * β₂ + (1 / 2) * β₂₂ * h0002 + β₁₂ * h0011 + (1 / 6) * h0003 * β₂₂₂
              + 2real(exp(2 * pi * t * im) * h1010) * ϵ * β₁ + real(exp(4 * pi * t * im) * h2010) * ϵ^2 * β₁ + real(h1110) * ϵ^2 * β₁
              + 2real(exp(2 * pi * t * im) * h1002) * ϵ * β₂₂ * (1 / 2) + real(exp(4 * pi * t * im) * h2002) * β₂₂ * (1 / 2) + real(h1102) * β₂₂ * (1 / 2)
              + 2real(exp(2 * pi * t * im) * h1011) * ϵ * β₁₂ + 2real(exp(2 * pi * t * im) * h1003) * ϵ * β₂₂₂ * (1 / 6)
              + 2real(exp(2 * pi * t * im) * h1001) * ϵ * β₂ + real(exp(4 * pi * t * im) * h2001) * ϵ^2 * (1 / 2) * β₂₂
              + real(h1101) * ϵ^2 * β₂ + (2 / 6)real(exp(6 * pi * t * im) * h3001) * ϵ^3 * β₂ + real(exp(2 * pi * t * im) * h2101) * ϵ^3 * (1 / 2) * β₂
              + (2 / 24)real(exp(8 * pi * t * im) * h4001) * ϵ^4 * β₂ + (2 / 6)real(exp(4 * pi * t * im) * h3101) * ϵ^4 * β₂ + (1 / 4) * real(h2201) * ϵ^4 * β₂
              + (2 / 120)real(exp(10 * pi * t * im) * h5001) * ϵ^5 * β₂ + (2 / 24)real(exp(6 * pi * t * im) * h4101) * ϵ^5 * β₂ + (2 / 12)real(exp(2 * pi * t * im) * h3201) * ϵ^5 * β₂
              + real(exp(4 * pi * t * im) * h2000) * ϵ^2 + real(h1100) * ϵ^2
              + (2 / 6)real(exp(6 * pi * t * im) * h3000) * ϵ^3 + real(exp(2 * pi * t * im) * h2100) * ϵ^3
              + (2 / 24)real(exp(8 * pi * t * im) * h4000) * ϵ^4 + (2 / 6)real(exp(4 * pi * t * im) * h3100) * ϵ^4 + (1 / 4) * real(h2200) * ϵ^4
              + (2 / 120)real(exp(10 * pi * t * im) * h5000) * ϵ^5 + (2 / 24)real(exp(6 * pi * t * im) * h4100) * ϵ^5 + (2 / 12)real(exp(4 * pi * t * im) * h3200) * ϵ^5
              + (2 / 720)real(exp(12 * pi * t * im) * h6000) * ϵ^6 + (2 / 120)real(exp(8 * pi * t * im) * h5100) * ϵ^6 + (2 / 48)real(exp(4 * pi * t * im) * h4200) * ϵ^6 + (1 / 36)real(h3300) * ϵ^6
              + (2 / 5040)real(exp(14 * pi * t * im) * h7000) * ϵ^7 + (2 / 720)real(exp(10 * pi * t * im) * h6100) * ϵ^7 + (2 / 240)real(exp(6 * pi * t * im) * h5200) * ϵ^7 + (2 / 144)real(exp(2 * pi * t * im) * h4300) * ϵ^7) for t ∈ t]
  pars = genh.parameters - 2real(c₂) * K01 * ϵ^2 + (real(c₂) * K10 + (4a3201 * real(c₂) - 3 * real(c₃)) * K01 + 2 * (real(c₂))^2 * K02) * ϵ^4 + (2(real(c₃) - a3201 * real(c₂)) * K10 - 2 * (real(c₂))^2 * K11 - 2real(c₂) * (4 * a3201 * real(c₂) - 3 * real(c₃)) * K02 - (8 / 6) * (real(c₂))^3 * K03) * ϵ^6
  T = 2pi / (abs(genh.ω) + (imag(c₁) - 2real(c₂) * b101) * ϵ^2 + (real(c₂) * b110 + (4a3201 * real(c₂) - 3 * real(c₃)) * b101 + 2 * (real(c₂))^2 * b102 - 2real(c₂) * b201 + imag(c₂)) * ϵ^4)
  psol(profile, pars, collect(t), T, ncol, nothing, nothing)
end
