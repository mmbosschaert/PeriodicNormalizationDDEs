function Dγ(φ₁::Function,φ₂::Function,jet,τ,periodicsolution,τs,ncol,ap)
    T = periodicsolution["period"]
    par = periodicsolution["parameter"][ap]
    ts = T*vec(periodicsolution["mesh"])
    profile = periodicsolution["profile"][1:2,:]
    # nodes = first(legendre(ncol))
    # colpoints = hcat([ts[i*ncol+1] .+ (ts[(i+1)*ncol+1] - ts[i*ncol+1])/2 .* (nodes .+ 1) for i in 0:ntst-1]...)
    testintervals = ts[1:ncol:end]
    # γ = [profile[:,i] for i ∈ 1:size(profile,2)]
    γ = Vector{eltype(profile)}[eachcol(profile)...] 
    γζ = γζs(τ,τs,ncol,T,testintervals,ts,γ)
    Φ₁ = hcat(φ₁.(τ,-τs)...)
    Φ₂ = hcat(φ₂.(τ,-τs)...)
    jet.D2(γζ,par,Φ₁,Φ₂)
end

function Dγ(γ::Function,jet,τ,τs,par)
    γ = hcat(γ.(τ,-τs)...)
    jet.D01(γ,par)
end

# for branch continued with Julia
function Dγ(γ::Function,φ₁::Function,φ₂::Function,jet,τ,τs,par)
    γ = hcat(γ.(τ,-τs)...)
    Φ₁ = hcat(φ₁.(τ,-τs)...)
    Φ₂ = hcat(φ₂.(τ,-τs)...)
    jet.D2(γ,par,Φ₁,Φ₂)
end

# for branch continued with Julia
function Dγ(γ::Function,φ₁::Function,φ₂::Function,φ₃::Function,jet,τ,τs,par)
    γ = hcat(γ.(τ,-τs)...)
    Φ₁ = hcat(φ₁.(τ,-τs)...)
    Φ₂ = hcat(φ₂.(τ,-τs)...)
    Φ₃ = hcat(φ₃.(τ,-τs)...)
    jet.D3(γ,par,Φ₁,Φ₂,Φ₃)
end

function Dγ(φ₁::Function,φ₂::Function,φ₃::Function,jet,τ,periodicsolution,τs,ncol,ap)
    T = periodicsolution["period"]
    par = periodicsolution["parameter"][ap]
    ts = T*vec(periodicsolution["mesh"])
    profile = periodicsolution["profile"][1:2,:]
    # nodes = first(legendre(ncol))
    # colpoints = hcat([ts[i*ncol+1] .+ (ts[(i+1)*ncol+1] - ts[i*ncol+1])/2 .* (nodes .+ 1) for i in 0:ntst-1]...)
    testintervals = ts[1:ncol:end]
    # γ = [profile[:,i] for i ∈ 1:size(profile,2)]
    γ = Vector{eltype(profile)}[eachcol(profile)...] 
    γζ = γζs(τ,τs,ncol,T,testintervals,ts,γ)
    Φ₁ = hcat(φ₁.(τ,-τs)...)
    Φ₂ = hcat(φ₂.(τ,-τs)...)
    Φ₃ = hcat(φ₃.(τ,-τs)...)
    jet.D3(γζ,par,Φ₁,Φ₂,Φ₃)
end
