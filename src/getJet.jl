function getJet(model, n, τs)

# calculate derivatives
jac = Symbolics.jacobian
m = length(τs)

@variables φ[1:n,1:m] α[1:2] λ
φ = reshape(Symbolics.get_variables(φ), n, m) # cast into a Matrix of symbols
α = vec(reshape(Symbolics.get_variables(α), 2, 1))

@variables v₁[1:n*length(τs)] v₂[1:n*length(τs)] v₃[1:n*length(τs)] v₄[1:n*length(τs)] v₅[1:n*length(τs)] p₁[1:2] p₂[1:2] p₃[1:2]

# TODO: replace vcat with Symbolics.scalarize
v₁ = vcat(v₁)
v₂ = vcat(v₂)
v₃ = vcat(v₃)
v₄ = vcat(v₄)
v₅ = vcat(v₅)
p₁ = vcat(p₁)
p₂ = vcat(p₂)
p₃ = vcat(p₃)

D1 = build_function(jac(model(φ, α), φ[:]), φ, α, expression=Val{false})[1]
D2 = build_function(jac(jac(model(φ, α), φ[:])*v₁, φ[:])*v₂, φ, α, v₁, v₂, expression=Val{false})[1]
D2v1 = build_function(jac(jac(model(φ, α), φ[:])*v₁, φ[:]), φ, α, v₁, expression=Val{false})[1]
D11 = build_function(jac(jac(model(φ, α), φ[:])*v₁, α)*p₁, φ, α, v₁, p₁, expression=Val{false})[1]
D11v1 = build_function(jac(jac(model(φ, α), φ[:])*v₁, α), φ, α, v₁, expression=Val{false})[1]
D12 = build_function(jac(jac(jac(model(φ, α), φ[:])*v₁, α)*p₁, α)*p₂, φ, α, v₁, p₁, p₂, expression=Val{false})[1]
D21 = build_function(jac(jac(jac(model(φ, α), φ[:])*v₁, φ[:])*v₂, α)*p₁, φ, α, v₁, v₂, p₁, expression=Val{false})[1]
D31 = build_function(jac(jac(jac(jac(model(φ, α), φ[:])*v₁, φ[:])*v₂, φ[:])*v₃, α)*p₁, φ, α, v₁, v₂, v₃, p₁, expression=Val{false})[1]
D3 =  build_function(jac(jac(jac(model(φ, α), φ[:])*v₁, φ[:])*v₂, φ[:])*v₃, φ, α, v₁, v₂, v₃, expression=Val{false})[1]
D4 =  build_function(jac(jac(jac(jac(model(φ, α), φ[:])*v₁, φ[:])*v₂, φ[:])*v₃, φ[:])*v₄, φ, α, v₁, v₂, v₃, v₄, expression=Val{false})[1]
D5 = build_function(jac(jac(jac(jac(jac(model(φ, α), φ[:])*v₁, φ[:])*v₂, φ[:])*v₃, φ[:])*v₄, φ[:])*v₅, φ, α, v₁, v₂, v₃, v₄, v₅, expression=Val{false})[1]
D01 = build_function(jac(model(φ, α), α), φ, α, expression = Val{false})[1]
D02 = build_function(jac(jac(model(φ, α), α)*p₁, α)*p₂, φ, α, p₁, p₂, expression=Val{false})[1]
D03 = build_function(jac(jac(jac(model(φ, α), α)*p₁, α)*p₂, α)*p₃, φ, α, p₁, p₂, p₃, expression=Val{false})[1]

M = [D1(φ, α)[:,j*n+1:(j+1)*n] for j in 0:length(τs)-1]
Δ = build_function(Matrix(I*λ,n,n) - sum(M.*vec([exp(-τ*λ) for τ ∈ τs])), φ, α, λ, expression=Val{false})[1]
Δ′   = build_function(Matrix(I,n,n) + sum(M.*vec([τ*exp(-τ*λ) for τ ∈ τs])), φ, α, λ, expression=Val{false})[1]
Δ′′  = build_function(sum(M.*vec([-τ^2*exp(-τ*λ) for τ ∈ τs])), φ, α, λ, expression=Val{false})[1]
Δ′′′ = build_function(sum(M.*vec([ τ^3*exp(-τ*λ) for τ ∈ τs])), φ, α, λ, expression=Val{false})[1]
Δ⁴   = build_function(sum(M.*vec([-τ^4*exp(-τ*λ) for τ ∈ τs])), φ, α, λ, expression=Val{false})[1]
Δ⁵   = build_function(sum(M.*vec([ τ^5*exp(-τ*λ) for τ ∈ τs])), φ, α, λ, expression=Val{false})[1]

# matrices for stability calculation
Ms = build_function(M, φ, α, expression=Val{false})[1]

# better for memory (used in own BVP solvers)
M = []
for j=0:length(τs) - 1
    push!(M,build_function(D1(φ, α)[:,j*n+1:(j+1)*n], φ, α, expression=Val{false})[1])
end

jet = (
    system=model,
    M=M,
    Ms=Ms,
    Δ=Δ, 
    Δ′ = Δ′,
    Δ′′ = Δ′′,
    Δ′′′ = Δ′′′,
    Δ⁴ = Δ⁴, 
    Δ⁵ = Δ⁵,
    D1=D1,
    D01=D01,
    D11=D11,
    D11v1=D11v1,
    D2=D2,
    D2v1=D2v1,
    D02=D02,
    D3=D3,
    D21=D21,
    D12=D12,
    D03=D03,
    D4=D4,
    D31=D31,
    D5=D5)

jet

end
