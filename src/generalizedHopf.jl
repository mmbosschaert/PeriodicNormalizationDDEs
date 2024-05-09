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

# TODO: implement this
function normalform_beta(jet, hopf::GenHopf, τs)
  m = length(τs)
  φ = repeat(hopf.coords, 1, m)
  α = hopf.parameters

  λ = hopf.ω * im
  Δ(λ) = jet.Δ(φ, α, λ)
  Δ′(λ) = jet.Δ′(φ, α, λ)
  Δ′′(λ) = jet.Δ′′(φ, α, λ)

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

  h2010(θ) = exp(2λ * θ) * Δ(2λ) \ (A₁(h2000, v10) + 2B(ϕ, h1010) + B(h2000, h0010) + B₁(ϕ, ϕ, v10) + C(ϕ, ϕ, h0010)) - 2γ110 * Δ(2λ) \ (Δ′(2λ) - θ * Δ(2λ)) * h2000(θ)
  h2001(θ) = exp(2λ * θ) * Δ(2λ) \ (A₁(h2000, v01) + 2B(ϕ, h1001) + B(h2000, h0001) + B₁(ϕ, ϕ, v01) + C(ϕ, ϕ, h0001)) - 2γ101 * Δ(2λ) \ (Δ′(2λ) - θ * Δ(2λ)) * h2000(θ)

  h1110(θ) = Δ(0.0) \ (A₁(h1100, v10) + 2real(B(conj ∘ ϕ, h1010)) + B(h1100, h0010) + B₁(ϕ, conj ∘ ϕ, v10) + C(ϕ, conj ∘ ϕ, h0010)) - 2 * real(γ110) * Δ(0.0) \ (Δ′(0.0) - θ * Δ(0.0)) * h1100(θ)
  h1101(θ) = Δ(0.0) \ (A₁(h1100, v01) + 2real(B(conj ∘ ϕ, h1001)) + B(h1100, h0001) + B₁(ϕ, conj ∘ ϕ, v01) + C(ϕ, conj ∘ ϕ, h0001)) - 2 * real(γ101) * Δ(0.0) \ (Δ′(0.0) - θ * Δ(0.0)) * h1100(θ)

  γ210 = 0.5 * transpose(p) * (A₁(h2100, v10) + B(conj ∘ ϕ, h2010) + 2B(ϕ, h1110) + B(h2100, h0010) + B(h2000, conj ∘ h1010) + 2B(h1100, h1010) + B₁(h2000, conj ∘ ϕ, v10) + 2B₁(ϕ, h1100, v10) + 2C(ϕ, conj ∘ ϕ, h1010) + C(h2000, conj ∘ ϕ, h0010) + C(ϕ, ϕ, conj ∘ h1010) + 2C(ϕ, h1100, h0010) + C₁(ϕ, ϕ, conj ∘ ϕ, v10) + D(ϕ, ϕ, conj ∘ ϕ, h0010)) |> first
  γ201 = 0.5 * transpose(p) * (A₁(h2100, v01) + B(conj ∘ ϕ, h2001) + 2B(ϕ, h1101) + B(h2100, h0001) + B(h2000, conj ∘ h1001) + 2B(h1100, h1001) + B₁(h2000, conj ∘ ϕ, v01) + 2B₁(ϕ, h1100, v01) + 2C(ϕ, conj ∘ ϕ, h1001) + C(h2000, conj ∘ ϕ, h0001) + C(ϕ, ϕ, conj ∘ h1001) + 2C(ϕ, h1100, h0001) + C₁(ϕ, ϕ, conj ∘ ϕ, v01) + D(ϕ, ϕ, conj ∘ ϕ, h0001)) |> first

  # (;γ110,γ101,γ210,γ201,c₁,c₂,ℓ₁,ℓ₂,q,h2000,h1100,h2010,h2001,h1101,h1110,h2100,h2200,h3100,h1010,h1001,h0010,h0001)

  ω₂ = imag([γ110 γ101] * (real([γ110 γ101; γ210 γ201]) \ [0.0; 1.0])) |> first

  nmfm = GenHopfNormalformHigherOrder(
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
