mutable struct Hopf
    coords::Vector{Float64}
    parameters::Vector{Float64}
    v::Vector{ComplexF64}
    ω::Float64
    stability::Union{Vector{ComplexF64},Nothing}
    nmfm::Union{ComplexF64,Nothing}
end

function vec(hopf::Hopf)::Vector{Float64}
    vcat(hopf.coords, hopf.parameters,real(hopf.v),imag(hopf.v),hopf.ω)
end

function vec_to_point(::Type{Hopf},v::Vector{Float64},dims)
    Hopf(v[1:dims], v[dims+1:2dims], v[2dims+1:3dims] + v[3dims+1:4dims]*im, v[4dims+1], nothing, nothing)
end

function Hopf_res!(res, model, τs, Δre, Δim, xx, xx_prev, n)

    x, α, vre, vim , ω = xx[1:n], xx[n+1:n+2], xx[n+3:2n+2], xx[2n+3:3n+2], xx[3n+3]
    vre_prev, vim_prev = xx_prev[n+3:2n+2], xx_prev[2n+3:3n+2]

    res[1:n] = model(repeat(x,1,length(τs)),α)
    res[n+1:2n]  = Δre(repeat(x,1,length(τs)),α,ω)*vre - Δim(repeat(x,1,length(τs)),α,ω)*vim
    res[2n+1:3n] = Δre(repeat(x,1,length(τs)),α,ω)*vim + Δim(repeat(x,1,length(τs)),α,ω)*vre
    res[3n+1] = dot(vre,vre_prev) - dot(vim,vim_prev) - 1
    res[3n+2] = dot(vim,vre_prev) + dot(vre,vim_prev)
end

function Hopf_res!(res, model, τs, Δ, xx, xx_prev, n)

    x, α, v, ω = xx[1:n], xx[n+1:n+2], xx[n+3:2n+2] + xx[2n+3:3n+2]*im, xx[3n+3]
    v_prev = xx_prev[n+3:2n+2] + xx_prev[2n+3:3n+2]*im

    res[1:n] = model(repeat(x,1,length(τs)),α)
    res[n+1:2n]  = real(Δ(repeat(x,1,length(τs)),α,ω*im)*v)
    res[2n+1:3n] = imag(Δ(repeat(x,1,length(τs)),α,ω*im)*v)
    res[3n+1] = real(dot(v,v_prev)) - 1
    res[3n+2] = imag(dot(v,v_prev))
end

function normalform(jet, hopf::Hopf, τs)

    m = length(τs)
    φ = repeat(hopf.coords,1,m)
    α = hopf.parameters

    λ = hopf.ω*im
    Δ(λ) = jet.Δ(φ,α,λ)
    Δ′(λ)= jet.Δ′(φ,α,λ)
    # q = ([D rand(2,1); rand(1,3)]\[0;0;1])[1:2]
    q = nullspace(Δ(λ); atol=1e-14)
    p = nullspace(Δ(λ)';atol=1e-14)

    # normalize
    p /=  transpose(p)*(Δ′(λ)*q)
    # @show first(transpose(p)*(Δ′(λ)*q)) ≈ 1.0

    # multi-linear forms at bt point
    Ξ(h) = vcat([h(-τ) for τ ∈ τs]...)
    B(v₁,v₂)  = jet.D2(φ,α,Ξ(v₁),Ξ(v₂))
    C(v₁,v₂,v₃)  = jet.D3(φ,α,Ξ(v₁),Ξ(v₂),Ξ(v₃))

    ϕ(θ) = exp(λ*θ)*q
    h2000(θ) = exp(2*λ*θ)*(Δ(2*λ)\B(ϕ,ϕ))
    h1100(θ) = Δ(0)\B(ϕ,conj∘ϕ)

    c₁ = first((1/2)*transpose(p)*(B(conj∘ϕ,h2000) + 2*B(ϕ,h1100) + C(ϕ,ϕ,conj∘ϕ)))
    ℓ₁ = real(c₁)/hopf.ω
    ℓ₁
end

function locate_genh(jet, branch, τs)
    for hopf in branch
        hopf.nmfm = normalform(jet, hopf, τs)
    end
    nmfm = [hopf.nmfm for hopf in branch]
    findall(!iszero,sign.(nmfm)[2:end] - sign.(nmfm)[1:end-1])
end
